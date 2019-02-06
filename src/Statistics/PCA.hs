{-# LANGUAGE DeriveFunctor #-}
module Statistics.PCA where

import           Control.DeepSeq
import           Data.Array                 as Arr
import           Data.Array.Repa            as R
import           Data.Binary
import           Data.List                  as L
import           Data.Vector                as V
import           Data.Vector.Unboxed        as VU
import           Foreign.Storable
import           Numeric.LinearAlgebra      as LA
import           Numeric.LinearAlgebra.Data as LA
import           Utils.Parallel

data PCAData a
  = PCARowData { getPCARowData :: [[a]] } -- [a] is in a row and one feature per row
  | PCAColData { getPCAColData :: [[a]] } -- [a] is in a column and N feature per col
  deriving (Functor)

instance (NFData a) => NFData (PCAData a) where
  rnf (PCARowData x) = rnf x
  rnf (PCAColData x) = rnf x


data PCAMatrix a = PCAMatrix
  { pcaMean   :: [a]         -- length is N
  , pcaMatrix :: LA.Matrix a -- MxN, where N is the length of input features)
  , pcaEigVal :: [Double]
  } deriving (Show)

instance (Binary a, Unbox a, Element a) => Binary (PCAMatrix a) where
  put (PCAMatrix m mat eig) = do
    put m
    put . LA.toLists $ mat
    put eig
  get = do
    m <- get
    mat <- get
    e <- get
    return $! PCAMatrix m (LA.fromLists mat) e

{-# INLINE removeMean #-}
removeMean :: (Num e, Unbox e) => [e] -> PCAData e -> PCAData e
removeMean mean (PCARowData xss) =
  PCARowData $ L.zipWith (\xs m -> L.map (\x -> x - m) xs) xss mean
removeMean mean (PCAColData xss) =
  PCAColData $ L.map (\xs -> L.zipWith (-) xs mean) xss

{-# INLINE removeMeanP #-}
removeMeanP ::
     (Num e, Unbox e, NFData e)
  => ParallelParams
  -> [e]
  -> PCAData e
  -> PCAData e
removeMeanP parallelParams mean (PCARowData xss) =
  PCARowData $
  parZipWithChunk
    parallelParams
    rdeepseq
    (\xs m -> L.map (\x -> x - m) xs)
    xss
    mean
removeMeanP parallelParams mean (PCAColData xss) =
  PCAColData $
  parMapChunk parallelParams rdeepseq (\xs -> L.zipWith (-) xs mean) xss


{-# INLINE computeRemoveMean #-}
computeRemoveMean ::
     (Num e, Unbox e, NFData e, Fractional e)
  => ParallelParams
  -> PCAData e
  -> ([e],PCAData e)
computeRemoveMean parallelParams ys@(PCARowData xs) =
  let mean =
        parMapChunk
          parallelParams
          rdeepseq
          ((/ fromIntegral (L.length . L.head $ xs)) . L.sum)
          xs
   in (mean, removeMeanP parallelParams mean ys)
computeRemoveMean parallelParams ys@(PCAColData xs) =
  (mean, removeMeanP parallelParams mean ys)
  where
    s = L.foldl1' (L.zipWith (+)) xs
    mean = L.map (/ fromIntegral (L.length xs)) s

{-# INLINE computeRemoveMeanS #-}
computeRemoveMeanS ::
     (Num e, Unbox e, Fractional e) => PCAData e -> ([e], PCAData e)
computeRemoveMeanS ys@(PCARowData xs) =
  let mean = L.map ((/ fromIntegral (L.length . L.head $ xs)) . L.sum) xs
   in (mean, removeMean mean ys)
computeRemoveMeanS ys@(PCAColData xs) = (mean, removeMean mean ys)
  where
    s = L.foldl1' (L.zipWith (+)) xs
    mean = L.map (/ fromIntegral (L.length xs)) s

{-# INLINE pcaReduction #-}
pcaReduction ::
     (Unbox e, Element e, Numeric e) => PCAMatrix e -> PCAData e -> PCAData e
pcaReduction (PCAMatrix mean mat eig) xs@(PCARowData _) =
  let ys = LA.fromLists . getPCARowData . removeMean mean $ xs
   in PCARowData . LA.toLists $ mat LA.<> ys
pcaReduction (PCAMatrix mean mat eig) xs@(PCAColData _) =
  let ys =
        fromColumns . L.map LA.fromList . getPCAColData . removeMean mean $ xs
   in PCAColData . L.map LA.toList . toColumns $ mat LA.<> ys

{-# INLINE pcaReductionP #-}
pcaReductionP ::
     (Unbox e, Element e, Numeric e, NFData e)
  => ParallelParams
  -> PCAMatrix e
  -> PCAData e
  -> PCAData e
pcaReductionP parallelParams (PCAMatrix mean mat _) xs@(PCARowData _) =
  let ys =
        LA.fromLists . getPCAColData . removeMeanP parallelParams mean $
        xs
   in PCARowData . LA.toLists $ mat LA.<> ys
pcaReductionP parallelParams (PCAMatrix mean mat _) xs@(PCAColData _) =
  let ys =
        fromColumns .
        L.map LA.fromList . getPCAColData . removeMeanP parallelParams mean $
        xs
   in PCAColData . L.map LA.toList . toColumns $ mat LA.<> ys

-- Use HMatrix to do the matrix operation
{-# INLINE pcaSVDS #-}
pcaSVDS ::
     (Unbox e, Element e, Numeric e, Num e, Fractional e, Field e)
  => Int
  -> PCAData e
  -> (PCAMatrix e, [Double], PCAData e)
pcaSVDS n xs = (pcaMat, (fst . L.unzip $ sortedPairs), reducedVecs)
  where
    (mean, meanRemovedVecs) = computeRemoveMeanS xs
    d'' =
      case meanRemovedVecs of
        PCARowData ys -> LA.fromLists . L.transpose $ ys
        PCAColData ys -> LA.fromLists $ ys
    (_, vec', uni') = thinSVD d''
    vec = L.map abs . LA.toList $ vec'
    uni = toColumns $ uni'
    sortedPairs = L.reverse $ sortOn fst $ L.zip vec uni
    pairs = L.unzip . L.take n $ sortedPairs
    mat = fromRows . snd $ pairs
    eigenValVec = fst $ pairs
    pcaMat = PCAMatrix mean mat eigenValVec
    reducedVecs = pcaReduction pcaMat meanRemovedVecs
    
{-# INLINE pcaSVDS' #-}
pcaSVDS' ::
     (Unbox e, Element e, Numeric e, Num e, Fractional e, Field e)
  => Int
  -> PCAData e
  -> (PCAMatrix e, [Double], PCAData e)
pcaSVDS' _ xs = (pcaMat, eigenValVec, reducedVecs)
  where
    (mean, meanRemovedVecs) = computeRemoveMeanS xs
    d'' =
      case meanRemovedVecs of
        PCARowData ys -> LA.fromLists . L.transpose $ ys
        PCAColData ys -> LA.fromLists $ ys
    (_, vec', uni') = thinSVD d''
    vec = L.map abs . LA.toList $ vec'
    uni = toColumns $ uni'
    sortedPairs = L.reverse $ sortOn fst $ L.zip vec uni
    eigList = fst . L.unzip $ sortedPairs
    n = div (L.length . L.filter (> 1) $ eigList) 4
    pairs = L.unzip . L.take n $ sortedPairs
    mat = fromRows . snd $ pairs
    eigenValVec = fst $ pairs
    pcaMat = PCAMatrix mean mat eigenValVec
    reducedVecs = pcaReduction pcaMat meanRemovedVecs

{-# INLINE pcaSVDP #-}
pcaSVDP ::
     (Num e, Unbox e, NFData e, Fractional e, Field e)
  => ParallelParams
  -> Int
  -> PCAData e
  -> (PCAMatrix e, PCAData e)
pcaSVDP parallelParams n xs = (pcaMat, reducedVecs)
  where
    (mean, meanRemovedVecs) = computeRemoveMean parallelParams xs
    d'' =
      case meanRemovedVecs of
        PCARowData ys -> LA.fromLists . L.transpose $ ys
        PCAColData ys -> LA.fromLists $ ys
    (_, vec', uni') = thinSVD d''
    vec = L.map abs . LA.toList $ vec'
    uni = toColumns $ uni'
    pairs = L.unzip . L.take n . L.reverse $ sortOn fst $ L.zip vec uni
    mat = fromRows . snd $ pairs
    eigenValVec = fst $ pairs
    pcaMat = PCAMatrix mean mat eigenValVec
    reducedVecs = pcaReductionP parallelParams pcaMat meanRemovedVecs

{-# INLINE zcaWhiten #-}
zcaWhiten :: PCAMatrix Double -> PCAData Double -> PCAData Double
zcaWhiten (PCAMatrix mean mat eig) xs@(PCARowData _) =
  let ys = LA.fromLists . getPCARowData . removeMean mean $ xs
      diagMat = diagl . L.map (\x -> 1 / sqrt x) $ eig
   in PCARowData . LA.toLists $ (tr' mat) LA.<> diagMat LA.<> mat LA.<> ys
zcaWhiten (PCAMatrix mean mat eig) xs@(PCAColData _) =
  let ys =
        fromColumns . L.map LA.fromList . getPCAColData . removeMean mean $ xs
      diagMat = diagl . L.map (\x -> 1 / sqrt x) $ eig
   in PCAColData . L.map LA.toList . toColumns $
      (tr' mat) LA.<> diagMat LA.<> mat LA.<> ys

-- -- perform pcaReduction on uncentered data
-- {-# INLINE pcaReduction #-}
-- pcaReduction :: (Num e, Unbox e) => PCAMatrix e -> VU.Vector e -> VU.Vector e
-- pcaReduction (PCAMatrix mean mat _) x
--   -- VU.zipWith (+) mean .
--  = V.convert . V.map (VU.sum . VU.zipWith (*) y) $ mat
--   where
--     y = removeMean mean x

-- {-# INLINE pcaReductionP #-}
-- pcaReductionP ::
--      (NFData e, Num e, Unbox e)
--   => ParallelParams
--   -> PCAMatrix e
--   -> [VU.Vector e]
--   -> [VU.Vector e]
-- pcaReductionP parallelParams pcaMat =
--   parMapChunk parallelParams rdeepseq (pcaReduction pcaMat)




-- {-# INLINE pcaReduction' #-}
-- pcaReduction' ::
--      (Num e, Unbox e, Storable e, Element e, Numeric e)
--   => PCAMatrix' e
--   -> [VU.Vector e]
--   -> [VU.Vector e]
-- pcaReduction' (PCAMatrix' mean mat) xs =
--   L.map (VU.fromList . LA.toList) . toColumns $ mat LA.<> ys
--   where
--     ys = fromColumns . L.map (LA.fromList . VU.toList . removeMean mean) $ xs

-- {-# INLINE pcaReductionP' #-}
-- pcaReductionP' ::
--      (Num e, Unbox e, Storable e, Element e, Numeric e)
--   => ParallelParams
--   -> PCAMatrix' e
--   -> [[VU.Vector e]]
--   -> [[VU.Vector e]]
-- pcaReductionP' parallelParams pcaMat =
--   parMapChunk parallelParams rdeepseq (pcaReduction' pcaMat)

-- {-# INLINE pcaReduction2Rows' #-}
-- pcaReduction2Rows' ::
--      (Num e, Unbox e, Storable e, Element e, Numeric e)
--   => PCAMatrix' e
--   -> [VU.Vector e]
--   -> [VU.Vector e]
-- pcaReduction2Rows' (PCAMatrix' mean mat) xs =
--   L.map (VU.fromList . LA.toList) . toRows $ mat LA.<> ys
--   where
--     means = VU.toList mean
--     ys =
--       fromRows $
--       L.zipWith (\m -> LA.fromList . VU.toList . VU.map (\y -> y - m)) means xs

-- {-# INLINE pcaReduction2RowsP' #-}
-- pcaReduction2RowsP' ::
--      (Num e, Unbox e, Storable e, Element e, Numeric e)
--   => ParallelParams
--   -> PCAMatrix' e
--   -> [[VU.Vector e]]
--   -> [[VU.Vector e]]
-- pcaReduction2RowsP' parallelParams pcaMat =
--   parMapChunk parallelParams rdeepseq (pcaReduction2Rows' pcaMat)


-- {-# INLINE zcaWhiten #-}
-- zcaWhiten ::
--      (Num e, Unbox e, Floating e)
--   => PCAMatrix e
--   -> VU.Vector e
--   -> VU.Vector e
--   -> VU.Vector e
-- zcaWhiten (PCAMatrix mean mat _) eigVec x =
--   VU.fromList $ L.map (VU.sum . VU.zipWith (*) z) mat'
--   where
--     y = removeMean mean x
--     z =
--       VU.zipWith (\y z -> z / (sqrt y)) eigVec .
--       V.convert . V.map (VU.sum . VU.zipWith (*) y) $
--       mat
--     mat' = L.map VU.fromList . L.transpose . L.map VU.toList . V.toList $ mat
