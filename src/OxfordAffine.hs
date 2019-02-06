module OxfordAffine
  ( module OxfordAffine
  , module Statistics.KMeans
  , module Statistics.PCA
  , module Types
  , module Utils
  ) where

import           Control.Monad                as M
import           Control.Monad.IO.Class       (liftIO)
import           Control.Monad.Parallel       as MP
import           Control.Monad.Trans.Resource
import           Data.Array.Repa              as R
import           Data.Conduit                 as C
import           Data.Conduit.List            as CL
import           Data.List                    as L
import           Data.Set                     as S
import           Data.Vector.Generic          as VG
import           Data.Vector.Storable         as VS
import           Data.Vector.Unboxed          as VU
import           Descriptor.Pinwheel
import           DFT.Plan
import           Image.IO
import           Image.Transform
import           OxfordAffine.Parser
import           OxfordAffine.Region
import           Statistics.KMeans
import           Statistics.PCA
import           System.FilePath
import           System.IO                    as IO
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

-- {-# INLINE descriptorPCA #-}
-- descriptorPCA :: PCAMatrix Double -> AffineDescriptor -> AffineDescriptor
-- descriptorPCA pcaMat = fmap (l2norm . VU.map (\x -> (signum x) * sqrt (abs x)) . pcaReduction pcaMat)

-- {-# INLINE descriptorPCAP #-}
-- descriptorPCAP ::
--      ParallelParams -> PCAMatrix Double -> AffineDescriptor -> AffineDescriptor
-- descriptorPCAP parallelParams pcaMat (AffineData path xs) =
--   AffineData path .
--   parMapChunk parallelParams rdeepseq (fmap (l2norm . VU.map (\x -> (signum x) * sqrt (abs x)) . pcaReduction pcaMat)) $
--   xs

-- {-# INLINE descriptorPCAP' #-}
-- descriptorPCAP' ::
--      ParallelParams -> PCAMatrix Double -> VU.Vector Double  -> AffineDescriptor -> AffineDescriptor
-- descriptorPCAP' parallelParams pcaMat eigVec (AffineData path xs) =
--   AffineData path .
--   parMapChunk
--     parallelParams
--     rdeepseq
--     (fmap
--        (l2norm .
--         VU.map (\x -> (signum x) * sqrt (abs x)) . zcaWhiten pcaMat eigVec)) $
--   xs

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
getFeaturePatch ::
     Int
  -> AffineRegionParameter
  -> R.Array U DIM3 Double
  -> AffineRegionFeature
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
    (\x arr -> fmap (\_ -> sliceArrayPositionStride stride arr) x)
    xs
    arrs

{-# INLINE getFeaturePatchP' #-}
getFeaturePatchP' ::
     ParallelParams
  -> Int
  -> AffineRegionParameter
  -> [R.Array U DIM3 Double]
  -> AffineRegionFeature
getFeaturePatchP' parallelParams stride (AffineData path xs) arrs =
  AffineData path $
  parZipWithChunk
    parallelParams
    rdeepseq
    (\x arr -> fmap (\_ ->  sliceArrayPositionStride' stride arr) x)
    xs
    arrs

{-# INLINE getFeatureFromArray #-}
getFeatureFromArray :: Int -> R.Array U DIM3 Double -> [VU.Vector Double]
getFeatureFromArray stride arr =
  let (Z :. _ :. cols :. rows) = extent arr
      idx =
        downsampleIndex
          (cols, rows)
          stride
          [(i, j) | i <- [0 .. cols - 1], j <- [0 .. rows - 1]]
   in L.map l2norm . getFeature arr $ idx


{-# INLINE getFeatureFromArrayP #-}
getFeatureFromArrayP ::
     ParallelParams -> Int -> R.Array U DIM3 Double -> [VU.Vector Double]
getFeatureFromArrayP parallelParams stride arr =
  let (Z :. _ :. cols :. rows) = extent arr
      idx =
        downsampleIndex
          (cols, rows)
          stride
          [(i, j) | i <- [0 .. cols - 1], j <- [0 .. rows - 1]]
   in parMapChunk parallelParams rdeepseq l2norm . getFeature arr $ idx

-- {-# INLINE applyPinwheelConvolutionCascade #-}
-- applyPinwheelConvolutionCascade ::
--      ParallelParams
--   -> Int
--   -> Int
--   -> Int
--   -> Int
--   -> DFTPlan
--   -> Int
--   -> Int
--   -> [[VS.Vector (Complex Double)]]
--   -> [[VS.Vector Double]]
--   -> [[VS.Vector Double]]
--   -> IO ([[VS.Vector Double]], [[VS.Vector Double]])
-- applyPinwheelConvolutionCascade parallelParams 1 stride _ descriptorLength plan rows cols filters xs ys = do
--   filteredPatchVec1 <-
--     MP.mapM
--       (applyPinwheelConvolution plan rows cols filters .
--        L.map (VS.map (\x -> x :+ 0)))
--       xs
--   filteredPatchVec2 <-
--     MP.mapM
--       (applyPinwheelConvolution plan rows cols filters .
--        L.map (VS.map (\x -> x :+ 0)))
--       ys
--   let pinwheelFeature1 =
--         parMapChunk
--           parallelParams
--           rdeepseq
--           (L.map l2norm .
--            sliceArrayPositionStride 1 .
--            fromUnboxed
--              (Z :. (L.length . L.head . L.head $ filteredPatchVec1) :. cols :.
--               rows) .
--            VS.convert . VS.concat . L.head) $
--         filteredPatchVec1
--       pinwheelFeature2 =
--         parMapChunk
--           parallelParams
--           rdeepseq
--           (L.map l2norm .
--            sliceArrayPositionStride 1 .
--            fromUnboxed
--              (Z :. (L.length . L.head . L.head $ filteredPatchVec2) :. cols :.
--               rows) .
--            VS.convert . VS.concat . L.head) $
--         filteredPatchVec2
--       pinwheelFeatureStride1 =
--         parMapChunk
--           parallelParams
--           rdeepseq
--           (L.map l2norm .
--            sliceArrayPositionStride stride .
--            fromUnboxed
--              (Z :. (L.length . L.head . L.head $ filteredPatchVec1) :. cols :.
--               rows) .
--            VS.convert . VS.concat . L.head) $
--         filteredPatchVec1
--       pinwheelFeatureStride2 =
--         parMapChunk
--           parallelParams
--           rdeepseq
--           (L.map l2norm .
--            sliceArrayPositionStride stride .
--            fromUnboxed
--              (Z :. (L.length . L.head . L.head $ filteredPatchVec2) :. cols :.
--               rows) .
--            VS.convert . VS.concat . L.head) $
--         filteredPatchVec2
--       (pcaMat, _, _) =
--         pcaSVDS'
--           descriptorLength
--           (L.concat pinwheelFeatureStride1 L.++ L.concat pinwheelFeatureStride2)
--       dimReducedPinwheelFeature1 =
--         pcaReduction2RowsP' parallelParams pcaMat pinwheelFeature1
--       dimReducedPinwheelFeature2 =
--         pcaReduction2RowsP' parallelParams pcaMat pinwheelFeature2
--   return
--     ( L.map (L.map VU.convert) dimReducedPinwheelFeature1
--     , L.map (L.map VU.convert) dimReducedPinwheelFeature2)
-- applyPinwheelConvolutionCascade parallelParams n stride pcaLength descriptorLength plan rows cols filters xs ys = do
--   filteredPatchVec1 <-
--     MP.mapM
--       (applyPinwheelConvolution plan rows cols filters .
--        L.map (VS.map (\x -> x :+ 0)))
--       xs
--   filteredPatchVec2 <-
--     MP.mapM
--       (applyPinwheelConvolution plan rows cols filters .
--        L.map (VS.map (\x -> x :+ 0)))
--       ys
--   let pinwheelFeature1 =
--         parMapChunk
--           parallelParams
--           rdeepseq
--           (L.map l2norm .
--            sliceArrayPositionStride 1 .
--            fromUnboxed
--              (Z :. (L.length . L.head . L.head $ filteredPatchVec1) :. cols :.
--               rows) .
--            VS.convert . VS.concat . L.head) $
--         filteredPatchVec1
--       pinwheelFeature2 =
--         parMapChunk
--           parallelParams
--           rdeepseq
--           (L.map l2norm .
--            sliceArrayPositionStride 1 .
--            fromUnboxed
--              (Z :. (L.length . L.head . L.head $ filteredPatchVec2) :. cols :.
--               rows) .
--            VS.convert . VS.concat . L.head) $
--         filteredPatchVec2
--       pinwheelFeatureStride1 =
--         parMapChunk
--           parallelParams
--           rdeepseq
--           (L.map l2norm .
--            sliceArrayPositionStride stride .
--            fromUnboxed
--              (Z :. (L.length . L.head . L.head $ filteredPatchVec1) :. cols :.
--               rows) .
--            VS.convert . VS.concat . L.head) $
--         filteredPatchVec1
--       pinwheelFeatureStride2 =
--         parMapChunk
--           parallelParams
--           rdeepseq
--           (L.map l2norm .
--            sliceArrayPositionStride stride .
--            fromUnboxed
--              (Z :. (L.length . L.head . L.head $ filteredPatchVec2) :. cols :.
--               rows) .
--            VS.convert . VS.concat . L.head) $
--         filteredPatchVec2
--       (pcaMat, _, _) =
--         pcaSVDS'
--           pcaLength
--           (L.concat pinwheelFeatureStride1 L.++ L.concat pinwheelFeatureStride2)
--       dimReducedPinwheelFeature1 =
--         pcaReduction2RowsP' parallelParams pcaMat pinwheelFeature1
--       dimReducedPinwheelFeature2 =
--         pcaReduction2RowsP' parallelParams pcaMat pinwheelFeature2
--   applyPinwheelConvolutionCascade
--     parallelParams
--     (n - 1)
--     stride
--     pcaLength
--     descriptorLength
--     plan
--     rows
--     cols
--     filters
--     (L.map (L.map VU.convert) dimReducedPinwheelFeature1)
--     (L.map (L.map VU.convert) dimReducedPinwheelFeature2)


{-# INLINE firstLayerConduit #-}
firstLayerConduit ::
     ParallelParams
  -> DFTPlan
  -> [[VS.Vector (Complex Double)]]
  -> Int
  -> Int
  -> ConduitT (ImageRepa Double) [VS.Vector Double] (ResourceT IO) ()
firstLayerConduit parallelParams plan filters rows cols = do
  xs <- CL.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do ys <-
          liftIO $
          MP.mapM
            (fmap (L.head) .
             applyPinwheelConvolution plan rows cols filters .
             L.map (VU.convert . VU.map (\x -> x :+ 0)) .
             sliceArrayChannel . imageContent)
            xs
        -- zs <-
        --   liftIO $
        --   MP.mapM
        --     (fmap (l2norm' . L.head) .
        --      applyPinwheelConvolution plan rows cols filters .
        --      L.map (VS.map (\x -> x :+ 0)))
        --     ys
        sourceList ys
        firstLayerConduit parallelParams plan filters rows cols)

{-# INLINE concatConduit #-}
concatConduit :: ConduitT (a, a) a (ResourceT IO) ()
concatConduit = mapFoldable (\(x,y) -> [x,y])

{-# INLINE getPCAColDataConduit #-}
getPCAColDataConduit ::
     ParallelParams
  -> Int
  -> Int
  -> Int
  -> ConduitT [VS.Vector Double] (PCAData Double) (ResourceT IO) ()
getPCAColDataConduit parallelParams stride rows cols = do
  xs <- CL.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do let ys =
              parMapChunk
                parallelParams
                rdeepseq
                (PCAColData .
                 L.map VU.toList .
                 sliceArrayPositionStride stride .
                 fromUnboxed (Z :. (L.length . L.head $ xs) :. cols :. rows) .
                 VS.convert . VS.concat) $
              xs
        sourceList ys
        getPCAColDataConduit parallelParams stride rows cols)


{-# INLINE getPCARowDataConduit #-}
getPCARowDataConduit ::
     ParallelParams
  -> Int
  -> Int
  -> ConduitT [VS.Vector Double] (PCAData Double) (ResourceT IO) ()
getPCARowDataConduit parallelParams rows cols = do
  xs <- CL.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do let ys =
              parMapChunk parallelParams rdeepseq (PCARowData . L.map VS.toList) $
              xs
        sourceList ys
        getPCARowDataConduit parallelParams rows cols)

{-# INLINE pcaSink #-}
pcaSink ::
     ParallelParams
  -> Int
  -> ConduitT (PCAData Double) Void (ResourceT IO) (PCAMatrix Double)
pcaSink parallelarams k = do
  xs <- CL.consume
  if L.null xs
    then error "pcaSink: empty stream."
    else case L.head xs of
           PCARowData _ ->
             error $ "pcaSink: row data is inefficent to concatenate."
           PCAColData _ -> do
             let (pcaMat, eigList, _) =
                   pcaSVDS k . (PCAColData . L.concatMap getPCAColData) $ xs
             liftIO . print . L.length $ eigList
             liftIO . print $ eigList
             return pcaMat

-- return position-wise features
{-# INLINE pcaConduit #-}
pcaConduit ::
     ParallelParams
  -> PCAMatrix Double
  -> ConduitT (PCAData Double) (PCAData Double) (ResourceT IO) ()
pcaConduit parallelParams pcaMat = do
  xs <- CL.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do let ys = parMapChunk parallelParams rdeepseq (pcaReduction pcaMat) xs
        sourceList ys
        pcaConduit parallelParams pcaMat)

{-# INLINE zcaConduit #-}
zcaConduit ::
     ParallelParams
  -> PCAMatrix Double
  -> ConduitT (PCAData Double) (PCAData Double) (ResourceT IO) ()
zcaConduit parallelParams pcaMat = do
  xs <- CL.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do let ys = parMapChunk parallelParams rdeepseq (zcaWhiten pcaMat) xs
        sourceList ys
        zcaConduit parallelParams pcaMat)

-- {-# INLINE pcaRowConduit #-}
-- pcaRowConduit ::
--      ParallelParams
--   -> PCAMatrix' Double
--   -> ConduitT [VS.Vector Double] [VS.Vector Double] (ResourceT IO) ()
-- pcaRowConduit parallelParams pcaMat = do
--   xs <- CL.take (batchSize parallelParams)
--   unless
--     (L.null xs)
--     (do let ys =
--               L.map (L.map VU.convert) .
--               pcaReduction2RowsP' parallelParams pcaMat .
--               L.map (L.map VS.convert) $
--               xs
--         sourceList ys
--         pcaRowConduit parallelParams pcaMat)


cascadeConvolutionConduit ::
     ParallelParams
  -> DFTPlan
  -> [[VS.Vector (Complex Double)]]
  -> Int
  -> Int
  -> ConduitT [VS.Vector Double] [VS.Vector Double] (ResourceT IO) ()
cascadeConvolutionConduit parallelParams plan filters rows cols = do
  xs <- CL.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do filteredVec <-
          liftIO $
          MP.mapM
            (fmap (L.head) .
             applyPinwheelConvolution plan rows cols filters .
             L.map (VS.map (\x -> x :+ 0)))
            xs
        sourceList filteredVec
        cascadeConvolutionConduit parallelParams plan filters rows cols)


cascadeConvolutionConduit' ::
     ParallelParams
  -> DFTPlan
  -> [[VS.Vector (Complex Double)]]
  -> Int
  -> Int
  -> ConduitT [VS.Vector Double] [VS.Vector Double] (ResourceT IO) ()
cascadeConvolutionConduit' parallelParams plan filters rows cols = do
  xs <- CL.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do filteredVec <-
          liftIO $
          MP.mapM
            (fmap (l2norm' . L.head) .
             applyPinwheelConvolution plan rows cols filters .
             L.map (VS.map (\x -> x :+ 0)))
            xs
        sourceList filteredVec
        cascadeConvolutionConduit' parallelParams plan filters rows cols)

{-# INLINE zca #-}
zca ::
     ParallelParams
  -> Int
  -> AffineDescriptor
  -> AffineDescriptor
  -> (AffineDescriptor, AffineDescriptor)
zca parallelParams n xs'@(AffineData xPath as) ys'@(AffineData yPath bs) =
  let xs = L.map (VU.toList . detectedRegionFeature) . affineDataFeature $ xs'
      ys = L.map (VU.toList . detectedRegionFeature) . affineDataFeature $ ys'
      (pcaMat, _, _) = pcaSVDS n . PCAColData $ xs L.++ ys
      xsWhiten =
        L.map VU.fromList . getPCAColData . zcaWhiten pcaMat . PCAColData $ xs
      ysWhiten =
        L.map VU.fromList . getPCAColData . zcaWhiten pcaMat . PCAColData $ ys
   in ( AffineData xPath $
        parZipWithChunk
          parallelParams
          rdeepseq
          (\region vec ->
             fmap
               (\_ -> l2norm . VU.map (\x -> (signum x) * sqrt (abs x)) $ vec)
               region)
          as
          xsWhiten
      , AffineData yPath $
        parZipWithChunk
          parallelParams
          rdeepseq
          (\region vec ->
             fmap
               (\_ -> l2norm . VU.map (\x -> (signum x) * sqrt (abs x)) $ vec)
               region)
          bs
          ysWhiten)
