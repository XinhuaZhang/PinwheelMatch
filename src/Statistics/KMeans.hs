{-# LANGUAGE BangPatterns #-}
module Statistics.KMeans
  ( ClusterCenter
  , KMeansModel(..)
  , kmeans
  , computeSoftAssignment
  , computeSoftAssignmentP
  , vlad
  ) where
 
import           Control.Monad       as M
import           Utils.Parallel
import           Data.Binary
import           Data.List           as L
import           Data.Vector         as V
import           Data.Vector.Unboxed as VU
import           System.Random
import           Text.Printf

type ClusterCenter = V.Vector (VU.Vector Double)

data KMeansModel = KMeansModel
  { clusterSize :: V.Vector Int
  , center      :: ClusterCenter
  }

instance Binary KMeansModel where
  put (KMeansModel cs c) = do
    put . V.toList $ cs
    put . L.map VU.toList . V.toList $ c
  get = do
    cs <- get
    xs <- get
    return $! KMeansModel (V.fromList cs) . V.fromList . L.map VU.fromList $ xs

type Assignment = V.Vector Int

randomClusterCenterPP
  :: [VU.Vector Double]
  -> Int
  -> [V.Vector (VU.Vector Double)]
  -> IO [VU.Vector Double]
randomClusterCenterPP !centers 0 _ = return centers
randomClusterCenterPP [] !n xs = do
  randomGen <- M.replicateM (L.length xs) newStdGen
  let ys =
        V.concat $
        parZipWith
          rdeepseq
          (\x' g ->
             V.map fst .
             V.filter (\(_, p') -> p' > 0.5) .
             V.zip x' . getRandomVector g . V.length . L.head $
             xs)
          xs
          randomGen
  randomClusterCenterPP
    [VU.map (/ fromIntegral (V.length ys)) . V.foldl1' (VU.zipWith (+)) $ ys]
    (n - 1)
    xs
randomClusterCenterPP !centers !n xs = do
  randomGen <- M.replicateM (L.length xs) newStdGen
  let ys =
        parMap
          rdeepseq
          (V.map (\x' -> L.minimum $ L.map (distFunc x') centers))
          xs
      zs = parMap rdeepseq V.sum ys
      s = L.sum zs
      as =
        V.concat $
        parZipWith3
          rdeepseq
          (\xs' ys' g ->
             V.map fst . V.filter snd $
             V.zipWith3
               (\x' y' r' ->
                  if y' / s > r'
                    then (x', True)
                    else (x', False))
               xs'
               ys' .
             getRandomVector g . V.length . L.head $
             xs)
          xs
          ys
          randomGen
      center' =
        VU.map (/ fromIntegral (V.length as)) . V.foldl1' (VU.zipWith (+)) $ as
  if V.length as == 0
    then randomClusterCenterPP centers n xs
    else randomClusterCenterPP (center' : centers) (n - 1) xs
    
{-# INLINE getRandomVector #-}

getRandomVector :: StdGen -> Int -> V.Vector Double
getRandomVector gen n = V.unfoldrN n (\g -> Just . randomR (0, 1) $ g) gen

{-# INLINE computeAssignmentP #-}
computeAssignmentP :: ClusterCenter
                   -> [V.Vector (VU.Vector Double)]
                   -> [Assignment]
computeAssignmentP = parMap rdeepseq . computeAssignment

{-# INLINE computeMeanP #-}
computeMeanP
  :: Int
  -> Int
  -> [Assignment]
  -> [V.Vector (VU.Vector Double)]
  -> [V.Vector (VU.Vector Double)]
computeMeanP k nf = parZipWith rdeepseq (computeMean k nf)

{-# INLINE computeDistortionP #-}
computeDistortionP :: ClusterCenter
                   -> [Assignment]
                   -> [V.Vector (VU.Vector Double)]
                   -> Double
computeDistortionP clusterCenter assignments =
  L.sum . parZipWith rdeepseq (computeDistortion clusterCenter) assignments

{-# INLINE meanList2ClusterCenter #-}
meanList2ClusterCenter
  :: [V.Vector (VU.Vector Double)]
  -> [V.Vector (VU.Vector Double)]
  -> IO ClusterCenter
meanList2ClusterCenter vecs =
  fmap V.fromList . randomClusterCenterPP zs (k - L.length zs)
  where
    xs =
      L.map
        (V.map
           (\x ->
               if VU.any isNaN x
                 then (undefined, 0)
                 else (x, 1)))
        vecs
    ys =
      L.foldl'
        (V.zipWith
           (\(s, count) (vec, n) ->
               if n == 0
                 then (s, count)
                 else (VU.zipWith (+) s vec, count + 1)))
        (V.replicate k (VU.replicate nf 0, 0 :: Double))
        xs
    zs =
      V.toList .
      V.map (\(s, count) -> VU.map (/ count) s) .
      V.filter (\(_, count) -> count /= 0) $
      ys
    nf = VU.length . V.head . L.head $ vecs
    k = V.length . L.head $ vecs

computeClusterSize :: Int -> [Assignment] -> V.Vector Int
computeClusterSize k xs =
  V.accumulate (+) (V.replicate k 0) $ V.zip vec (V.replicate (V.length vec) 1)
  where
    vec = V.concat xs

kmeans :: ParallelParams -> Int -> FilePath -> Double -> [VU.Vector Double] -> IO KMeansModel
kmeans parallelParams k filePath threshold xs = do
  when
    (k > L.length xs)
    (error $
     "kmeans: number of center is greater than number of points.\n" L.++ show k L.++
     " vs " L.++
     (show . L.length $ xs))
  printf "Looking for %d centers in %d samples.\n" k (L.length xs)
  print (VU.length . L.head $ xs)
  randomCenter <- randomClusterCenterPP [] k ys
  putStrLn "Done."
  go
    (V.fromListN k randomCenter)
    (V.fromListN k randomCenter)
    (fromIntegral (maxBound :: Int))
    ys
  where
    nf = VU.length . L.head $ xs
    ys -- Because the data will be partitioned again and again,
     =
      L.map V.fromList . -- it is better partitioning them first, then using parMap
      splitList (div (L.length xs) (numThread parallelParams)) $
      xs
    go oldCenter !center' lastDistortion zs = do
      let assignment = computeAssignmentP center' zs
          distortion = computeDistortionP center' assignment zs
          ratio = (lastDistortion - distortion) / distortion
      printf "%e %.6f\n" distortion ratio
         -- distortion >= lastDistortion ||
      if abs ratio < threshold
        then do
          putStrLn ""
          return $! KMeansModel (computeClusterSize k assignment) oldCenter
        else do
          newCenter <-
            meanList2ClusterCenter (computeMeanP k nf assignment zs) zs
          encodeFile filePath (KMeansModel (computeClusterSize k assignment) newCenter)
          go center' newCenter distortion zs

{-# INLINE computeSoftAssignment #-}
computeSoftAssignment :: ClusterCenter -> VU.Vector Double -> VU.Vector Double
computeSoftAssignment center' x =
  normalizeVec . V.convert . V.map (\y -> max 0 (mean - y)) $ dist
  where
    dist = V.map (distFunc x) center'
    mean = V.sum dist / fromIntegral (V.length center')

computeSoftAssignmentP :: ParallelParams
                       -> ClusterCenter
                       -> [VU.Vector Double]
                       -> [VU.Vector Double]
computeSoftAssignmentP parallelParams center' =
  parMapChunk parallelParams rdeepseq (computeSoftAssignment center')

{-# INLINE computeAssignment #-}
computeAssignment :: ClusterCenter -> V.Vector (VU.Vector Double) -> Assignment
computeAssignment cluster =
  V.map (\x -> V.minIndex . V.map (distFunc x) $ cluster)

{-# INLINE computeMean #-}
computeMean
  :: Int
  -> Int
  -> Assignment
  -> V.Vector (VU.Vector Double)
  -> V.Vector (VU.Vector Double)
computeMean k nf assignmet =
  V.map (\(s, count) -> VU.map (/ count) s) .
  V.accumulate
    (\(s, count) vec -> (VU.zipWith (+) s vec, count + 1))
    (V.replicate k (VU.replicate nf 0, 0)) .
  V.zip assignmet

{-# INLINE computeDistortion #-}
computeDistortion :: ClusterCenter
                  -> Assignment
                  -> V.Vector (VU.Vector Double)
                  -> Double
computeDistortion clusterCenter assignments =
  V.sum .
  V.zipWith
    (\assignment vec -> distFunc vec (clusterCenter V.! assignment))
    assignments

{-# INLINE distFunc #-}
distFunc :: VU.Vector Double -> VU.Vector Double -> Double
distFunc vec1 vec2 =
  VU.sum $ VU.zipWith (\a b -> (a - b) ^ (2 :: Int)) vec1 vec2

{-# INLINE splitList #-}
splitList :: Int -> [a] -> [[a]]
splitList _ [] = []
splitList n ys = as : splitList n bs
  where
    (as, bs) = L.splitAt n ys

{-# INLINE normalizeVec #-}

normalizeVec :: VU.Vector Double -> VU.Vector Double
normalizeVec vec
  | s == 0 = VU.replicate (VU.length vec) 0
  | otherwise = VU.map (/ s) vec
  where
    s = sqrt . VU.sum . VU.map (^ (2 :: Int)) $ vec

{-# INLINE vlad #-}

vlad :: ClusterCenter -> [VU.Vector Double] -> VU.Vector Double
vlad center' =
  normalizeVec .
  VU.map (\x -> (signum x) * sqrt (abs x)) .
  VU.concat .
  V.toList .
  V.accum (VU.zipWith (+)) (V.replicate k (VU.replicate vecLen 0)) .
  L.map
    (\vec ->
       let idx = V.minIndex . V.map (distFunc vec) $ center'
       in (idx, VU.zipWith (-) vec (center' V.! idx)))
  where
    k = V.length center'
    vecLen = VU.length . V.head $ center'
    

-- {-# INLINE vlad #-}

-- vlad :: ClusterCenter -> [VU.Vector Double] -> VU.Vector Double
-- vlad center' =
--   normalizeVec .
--   VU.map (\x -> (signum x) * sqrt (abs x)) .
--   VU.concat .
--   V.toList .
--   V.accum (VU.zipWith (+)) (V.replicate k (VU.replicate vecLen 0)) .
--   L.concatMap
--     (\vec ->
--        let xs =
--              L.take 4 .
--              fst .
--              L.unzip .
--              L.sortOn snd . V.toList . V.imap (\i c -> (i, distFunc vec $ c)) $
--              center'
--        in L.map (\i -> (i, VU.zipWith (-) vec (center' V.! i))) xs)
--   where
--     k = V.length center'
--     vecLen = VU.length . V.head $ center'
