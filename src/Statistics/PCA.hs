{-# LANGUAGE BangPatterns #-}

module Statistics.PCA where

import           Data.Array            as Arr
import           Data.Array.Repa       as R
import           Data.Binary
import           Data.List             as L
import           Data.Vector           as V
import           Data.Vector.Unboxed   as VU
import           Foreign.Storable
import           Numeric.LinearAlgebra as LA
import           Utils.Parallel

data PCAMatrix a = PCAMatrix
  { pcaMean   :: !(VU.Vector a)
  , pcaMatrix :: !(V.Vector (VU.Vector a))
  } deriving (Show)

instance (Binary a, Unbox a) => Binary (PCAMatrix a) where
  put (PCAMatrix m mat) = do
    put . VU.toList $ m
    put . V.toList . V.map VU.toList $ mat
  get = do
    m <- get
    mat <- get
    return $! PCAMatrix (VU.fromList m) (V.fromList . L.map VU.fromList $ mat)

{-# INLINE computeRemoveMean #-}
computeRemoveMean ::
     (Num e, Unbox e, NFData e, Fractional e)
  => ParallelParams
  -> [VU.Vector e]
  -> (VU.Vector e, [VU.Vector e])
computeRemoveMean parallelParams xs =
  (mean, parMapChunk parallelParams rdeepseq (removeMean mean) xs)
  where
    !s = L.foldl1' (VU.zipWith (+)) xs
    !mean = VU.map (/ fromIntegral (L.length xs)) s

{-# INLINE computeRemoveMeanS #-}
computeRemoveMeanS ::
     (Num e, Unbox e, NFData e, Fractional e)
  => [VU.Vector e]
  -> (VU.Vector e, [VU.Vector e])
computeRemoveMeanS xs = (mean, L.map (removeMean mean) xs)
  where
    !s = L.foldl1' (VU.zipWith (+)) xs
    !mean = VU.map (/ fromIntegral (L.length xs)) s

{-# INLINE removeMean #-}
removeMean :: (Num e, Unbox e) => VU.Vector e -> VU.Vector e -> VU.Vector e
removeMean mean xs = VU.zipWith (-) xs mean

-- xs are mean-removed
{-# INLINE covarianceMatrix #-}
covarianceMatrix ::
     (Unbox e, NFData e, Fractional e, Storable e)
  => ParallelParams
  -> [VU.Vector e]
  -> Matrix e
covarianceMatrix parallelParams xs = (len >< len) . R.toList $ rArr
  where
    ys =
      parMapChunk
        parallelParams
        rdeepseq
        (\(i, j) ->
           ( (i, j)
           , (L.sum . L.map (\vec -> (vec VU.! i) * (vec VU.! j)) $ xs) /
             fromIntegral (L.length xs) -- .
            ))
      -- L.filter (\(i, j) -> j >= i) $
        [(i, j) | i <- [0 .. len - 1], j <- [0 .. len - 1]]
    len = VU.length . L.head $ xs
    arr = array ((0, 0), (len - 1, len - 1)) ys
    rArr =
      fromFunction
        (Z :. len :. len)
        (\(Z :. i :. j) ->
           if i <= j
             then arr Arr.! (i, j)
             else arr Arr.! (j, i))

-- xs are mean-removed
{-# INLINE covarianceMatrixP #-}
covarianceMatrixP ::
     (Unbox e, NFData e, Fractional e, Storable e)
  => ParallelParams
  -> [VU.Vector e]
  -> IO (Matrix e)
covarianceMatrixP parallelParams xs = do
  zs <- computeUnboxedP rArr
  return $! (len >< len) . R.toList $ zs
  where
    !ys =
      L.map
        (\(i, j) ->
           ( (i, j)
           , (L.sum . L.map (\vec -> (vec VU.! i) * (vec VU.! j)) $ xs) /
             (fromIntegral (L.length xs)))) .
      L.filter (\(i, j) -> j >= i) $
      [(i, j) | i <- [0 .. len - 1], j <- [0 .. len - 1]]
    len = VU.length . L.head $ xs
    !arr = array ((0, 0), (len - 1, len - 1)) ys -- :: AU.Array (Int, Int) (Complex Double)
    rArr =
      fromFunction
        (Z :. len :. len)
        (\(Z :. i :. j) ->
           if i <= j
             then arr Arr.! (i, j)
             else arr Arr.! (j, i))

-- compute PCAMatrix and using it to reduce the dimensions of input data
{-# INLINE pcaCovariance #-}
pcaCovariance ::
     (Unbox e, NFData e, Fractional e, Storable e, Field e)
  => ParallelParams
  -> Int
  -> [VU.Vector e]
  -> IO (PCAMatrix e, VU.Vector Double, [VU.Vector e])
pcaCovariance parallelParams n xs = do
  let (mean, meanRemovedVecs) = computeRemoveMean parallelParams xs
  covMat <- covarianceMatrixP parallelParams meanRemovedVecs
  let (val', vec') = eigSH $ trustSym covMat
      val = LA.toList val'
      vec = L.map (VU.fromList . LA.toList) . toColumns $ vec'
      pairs = L.unzip . L.take n . L.reverse . L.sortOn fst $ L.zip val vec
      mat = V.fromList . snd $ pairs
      eigenValVec = VU.fromList . fst $ pairs
      pcaMat = PCAMatrix mean mat
      reducedVecs =
        parMapChunk
          parallelParams
          rdeepseq
          (pcaReduction pcaMat)
          meanRemovedVecs
  return (pcaMat, eigenValVec, reducedVecs)

{-# INLINE pcaSVD #-}
pcaSVD ::
     (Num e, Unbox e, NFData e, Fractional e, Field e)
  => ParallelParams
  -> Int
  -> [VU.Vector e]
  -> (PCAMatrix e, VU.Vector Double, [VU.Vector e])
pcaSVD parallelParams n xs = (pcaMat, eigenValVec, reducedVecs)
  where
    (mean, meanRemovedVecs) = computeRemoveMean parallelParams xs
    d'' = fromRows . L.map (LA.fromList . VU.toList) $ meanRemovedVecs
    (_, vec', uni') = thinSVD d''
    vec = L.map abs . LA.toList $ vec'
    uni = L.map (VU.fromList . LA.toList) . toColumns $ uni'
    pairs = L.unzip . L.take n . L.reverse $ sortOn fst $ L.zip vec uni
    mat = V.fromList . snd $ pairs
    eigenValVec = VU.fromList . fst $ pairs
    pcaMat = PCAMatrix mean mat
    reducedVecs =
      parMapChunk parallelParams rdeepseq (pcaReduction pcaMat) meanRemovedVecs

{-# INLINE pcaSVDS #-}
pcaSVDS ::
     (Num e, Unbox e, NFData e, Fractional e, Field e)
  => Int
  -> [VU.Vector e]
  -> (PCAMatrix e, VU.Vector Double, [VU.Vector e])
pcaSVDS n xs = (pcaMat, eigenValVec, reducedVecs)
  where
    (mean, meanRemovedVecs) = computeRemoveMeanS xs
    d'' = fromRows . L.map (LA.fromList . VU.toList) $ meanRemovedVecs
    (_, vec', uni') = thinSVD d''
    vec = L.map abs . LA.toList $ vec'
    uni = L.map (VU.fromList . LA.toList) . toColumns $ uni'
    pairs = L.unzip . L.take n . L.reverse $ sortOn fst $ L.zip vec uni
    mat = V.fromList . snd $ pairs
    eigenValVec = VU.fromList . fst $ pairs
    pcaMat = PCAMatrix mean mat
    reducedVecs = L.map (pcaReduction pcaMat) meanRemovedVecs

-- perform pcaReduction on uncentered data
{-# INLINE pcaReduction #-}
pcaReduction :: (Num e, Unbox e) => PCAMatrix e -> VU.Vector e -> VU.Vector e
pcaReduction (PCAMatrix mean mat) x
  -- VU.zipWith (+) mean .
 = V.convert . V.map (VU.sum . VU.zipWith (*) y) $ mat
  where
    y = removeMean mean x

{-# INLINE pcaReductionP #-}
pcaReductionP ::
     (NFData e, Num e, Unbox e)
  => ParallelParams
  -> PCAMatrix e
  -> [VU.Vector e]
  -> [VU.Vector e]
pcaReductionP parallelParams pcaMat =
  parMapChunk parallelParams rdeepseq (pcaReduction pcaMat)


{-# INLINE zcaWhiten #-}
zcaWhiten ::
     (Num e, Unbox e, Floating e)
  => PCAMatrix e
  -> VU.Vector e
  -> VU.Vector e
  -> VU.Vector e
zcaWhiten (PCAMatrix mean mat) eigVec x =
  VU.fromList $ L.map (VU.sum . VU.zipWith (*) z) mat'
  where
    y = removeMean mean x
    z =
      VU.zipWith (\y z -> z / (sqrt y)) eigVec .
      V.convert . V.map (VU.sum . VU.zipWith (*) y) $
      mat
    mat' = L.map VU.fromList . L.transpose . L.map VU.toList . V.toList $ mat
