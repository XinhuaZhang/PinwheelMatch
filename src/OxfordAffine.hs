module OxfordAffine
  ( module OxfordAffine
  , module Statistics.KMeans
  , module Statistics.PCA
  , module Types
  , module Utils
  ) where

import           Control.Monad       as M
import           Data.Array.Repa     as R
import           Data.Conduit        as C
import           Data.List           as L
import           Data.Set            as S
import           Data.Vector.Unboxed as VU
import           OxfordAffine.Parser
import           OxfordAffine.Region
import           Statistics.KMeans
import           Statistics.PCA
import           System.FilePath
import           System.IO           as IO
import           Text.Printf
import           Types
import           Utils
import           Utils.Parallel


{-# INLINE computeAffineRegionParameterFromFile #-}
computeAffineRegionParameterFromFile ::
     FilePath -> FilePath -> Double -> Double -> IO AffineRegionParameter
computeAffineRegionParameterFromFile imagePath detectedRegionFilePath scale delta = do
  xs <- parseDetectedRegionFile detectedRegionFilePath
  return . AffineData imagePath $ xs

{-# INLINE computeAffineRegionFromFile #-}
computeAffineRegionFromFile ::
     FilePath -> FilePath -> Double -> Double -> IO AffineRegion
computeAffineRegionFromFile imagePath detectedRegionFilePath scale delta = do
  xs <- parseDetectedRegionFile detectedRegionFilePath
  let ys = parMap rdeepseq (detectedRegionSet scale delta) xs
  return . AffineData imagePath $ ys

{-# INLINE writeAffineDescriptor #-}
writeAffineDescriptor :: FilePath -> AffineDescriptor -> IO ()
writeAffineDescriptor folderPath (AffineData imagePath xs) = do
  withFile (folderPath </> (takeBaseName imagePath) L.++ ".pinwheel") WriteMode $ \h -> do
    IO.hPutStrLn h (show . VU.length . detectedRegionFeature . L.head $ xs)
    IO.hPutStrLn h (show . L.length $ xs)
    M.mapM_
      (\(DetectedRegion x y a b c vec) ->
         let str =
               printf
                 "%.8f %.8f %.8f %.8f %.8f%s\n"
                 x
                 y
                 a
                 b
                 c
                 (L.concatMap (\v -> printf " %.6f" v :: String) . VU.toList $
                  vec)
          in IO.hPutStr h str)
      xs

{-# INLINE getFeature #-}
getFeature :: (Unbox e) => R.Array U DIM3 e -> [(Int, Int)] -> [VU.Vector e]
getFeature arr =
  L.map (\(x, y) -> toUnboxed . computeS . R.slice arr $ (Z :. R.All :. x :. y))

{-# INLINE getFeatureFromAffineRegion #-}
getFeatureFromAffineRegion ::
     (Int,Int) ->  Int -> R.Array U DIM3 Double -> AffineRegion -> AffineRegionFeature
getFeatureFromAffineRegion boundary stride arr =
  fmap
    (L.map l2norm .
     getFeature arr . downsampleIndex boundary stride . S.toList)

{-# INLINE getFeatureFromAffineRegionP #-}
getFeatureFromAffineRegionP ::
     ParallelParams
  -> (Int, Int)
  -> Int
  -> R.Array U DIM3 Double
  -> AffineRegion
  -> AffineRegionFeature
getFeatureFromAffineRegionP parallelParams boundary stride arr (AffineData path xs) =
  AffineData path .
  parMapChunk
    parallelParams
    rdeepseq
    (fmap
       (L.map l2norm .
        getFeature arr . downsampleIndex boundary stride . S.toList)) $
  xs

{-# INLINE computeVLAD #-}
computeVLAD :: KMeansModel -> AffineRegionFeature -> AffineDescriptor
computeVLAD (KMeansModel _ c) = fmap (vlad c)

{-# INLINE computeVLADP #-}
computeVLADP ::
     ParallelParams -> KMeansModel -> AffineRegionFeature -> AffineDescriptor
computeVLADP parallelParams (KMeansModel _ c) (AffineData path xs) =
  AffineData path . parMapChunk parallelParams rdeepseq (fmap (vlad c)) $ xs

{-# INLINE descriptorPCA #-}
descriptorPCA :: PCAMatrix Double -> AffineDescriptor -> AffineDescriptor
descriptorPCA pcaMat = fmap (l2norm . VU.map (\x -> (signum x) * sqrt (abs x)) . pcaReduction pcaMat)

{-# INLINE descriptorPCAP #-}
descriptorPCAP ::
     ParallelParams -> PCAMatrix Double -> AffineDescriptor -> AffineDescriptor
descriptorPCAP parallelParams pcaMat (AffineData path xs) =
  AffineData path .
  parMapChunk parallelParams rdeepseq (fmap (l2norm . VU.map (\x -> (signum x) * sqrt (abs x)) . pcaReduction pcaMat)) $
  xs

{-# INLINE descriptorPCAP' #-}
descriptorPCAP' ::
     ParallelParams -> PCAMatrix Double -> VU.Vector Double  -> AffineDescriptor -> AffineDescriptor
descriptorPCAP' parallelParams pcaMat eigVec (AffineData path xs) =
  AffineData path .
  parMapChunk
    parallelParams
    rdeepseq
    (fmap
       (l2norm .
        VU.map (\x -> (signum x) * sqrt (abs x)) . zcaWhiten pcaMat eigVec)) $
  xs

{-# INLINE downsampleIndex #-}
downsampleIndex :: (Int, Int) -> Int -> [(Int, Int)] -> [(Int, Int)]
downsampleIndex (maxX, maxY) stride xs =
  L.filter
    (\(x, y) ->
       if (mod x stride == 0) &&
          (mod y stride == 0) &&
          (x >= 0) && (x < maxX) && (y >= 0) && (y < maxY)
         then True
         else False)
    xs

{-# INLINE norm8bit #-}
norm8bit :: VU.Vector Double -> VU.Vector Double
norm8bit vec =
  let minV = VU.minimum vec
      maxV = VU.maximum vec
   in VU.map (\x -> (x - minV) / (maxV - minV) * 255) vec

{-# INLINE getFeaturePatch #-}
getFeaturePatch :: Int -> AffineData () -> R.Array U DIM3 Double -> AffineRegionFeature
getFeaturePatch stride region arr =
  fmap (\_ -> sliceArrayPositionStride stride arr) region
  
{-# INLINE getFeaturePatchP #-}
getFeaturePatchP ::
     ParallelParams
  -> Int
  -> AffineRegionParameter
  -> [R.Array U DIM3 Double]
  -> AffineRegionFeature
getFeaturePatchP parallelParams stride (AffineData path xs) arrs =
  AffineData path $
  parZipWithChunk
    parallelParams
    rdeepseq
    (\x arr -> fmap (\_ -> L.map l2norm $ sliceArrayPositionStride stride arr) x)
    xs
    arrs
