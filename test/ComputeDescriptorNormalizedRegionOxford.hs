module ComputeDescriptorNormalizedRegionOxford where

import           Control.Monad                 as M
import           Control.Monad.Parallel        as MP
import           Data.Array.Repa               as R
import           Data.Complex
import           Data.List                     as L
import           Data.Set                      as S
import           Data.Vector.Storable          as VS
import           Data.Vector.Unboxed           as VU
import           Descriptor.EllipticalPinwheel
import           Descriptor.Pinwheel
import           Image.IO
import           OxfordAffine
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils
import           Utils.Parallel

main = do
  args@(imagePath1:patchFolderPath1:detectedRegionFilePath1:imagePath2:patchFolderPath2:detectedRegionFilePath2:scaleStr:deltaStr:radialFreqStr:angularFreqStr:alphaStr:numKmeansCenterStr:kmeansThresholdStr:descriptorLengthStr:strideStr:deltaThetaStr:deltaAStr:maxAStr:numThreadStr:batchSizeStr:_) <-
    getArgs
  let folderPath = "output/test/ComputeDescriptorNormalizedRegionOxford"
      scale = read scaleStr :: Double
      delta = read deltaStr :: Double
      radialFreq = read radialFreqStr :: Double
      angularFreq = read angularFreqStr :: Double
      alpha = read alphaStr :: Double
      numKmeansCenter = read numKmeansCenterStr :: Int
      kmeansThreshold = read kmeansThresholdStr :: Double
      descriptorLength = read descriptorLengthStr :: Int
      stride = read strideStr :: Int
      deltaTheta = (read deltaThetaStr :: Double) / 180 * pi
      deltaA = read deltaAStr :: Double
      maxA = read maxAStr :: Double
      numThread = read numThreadStr :: Int
      batchSize = read batchSizeStr :: Int
      parallelParams = ParallelParams numThread batchSize
  createDirectoryIfMissing True folderPath
  -- read image
  (ImageRepa 8 img1) <- readImageRepa imagePath1 False
  (ImageRepa 8 img2) <- readImageRepa imagePath2 False
  -- read image list
  patchPathList1' <- listDirectory patchFolderPath1
  patchPathList2' <- listDirectory patchFolderPath2
  let patchPathList1 =
        L.map
          (printf "%s/%05d.png" patchFolderPath1)
          [1 .. L.length patchPathList1']
      patchPathList2 =
        L.map
          (printf "%s/%05d.png" patchFolderPath2)
          [1 .. L.length patchPathList2']
  -- read affine regions
  affineRegion1 <-
    computeAffineRegionParameterFromFile
      imagePath1
      detectedRegionFilePath1
      scale
      delta
  affineRegion2 <-
    computeAffineRegionParameterFromFile
      imagePath2
      detectedRegionFilePath2
      scale
      delta
  when
    ((L.length patchPathList1 /= (L.length . affineDataFeature $ affineRegion1)) ||
     (L.length patchPathList2 /= (L.length . affineDataFeature $ affineRegion2)))
    (error $
     printf
       "%d patches in %s\n%d regions in %s\n%d patches in %s\n%d regions in %s"
       (L.length patchPathList1)
       (L.length . affineDataFeature $ affineRegion1)
       (L.length patchPathList2)
       (L.length . affineDataFeature $ affineRegion2))
  -- read patches
  patches1 <- MP.mapM (\path -> readImageRepa path False) patchPathList1
  patches2 <- MP.mapM (\path -> readImageRepa path False) patchPathList2
  let (Z :. nf :. cols :. rows) = extent . imageContent . L.head $ patches1
      -- mean1 = meanRepa img1
      -- mean2 = meanRepa img2
      -- patches1 =
      --   parMapChunk
      --     parallelParams
      --     rseq
      --     (\(ImageRepa d arr) -> let m = meanRepa arr
      --                            in ImageRepa d . removeMeanRepa m $ arr)
      --     patches1'
      -- patches2 =
      --   parMapChunk
      --     parallelParams
      --     rseq
      --     (\(ImageRepa d arr) -> let m = meanRepa arr
      --                            in ImageRepa d . removeMeanRepa m $ arr)
      --     patches2'
      pinwheelParams =
        PinwheelParams
          rows
          cols
          [0 .. radialFreq]
          [-angularFreq .. angularFreq]
          alpha
      -- canonicalEllipseList =
      --   [ CanonicalEllipse a theta
      --   | a <- [1 + deltaA,1 + 2 * deltaA .. maxA]
      --   , theta <- [0,deltaTheta .. pi - deltaTheta]
      --   ]
      -- ellipticalPinwheelParamsList =
      --   L.map
      --     (EllipticalPinwheelParams
      --        rows
      --        cols
      --        [0 .. radialFreq]
      --        [-angularFreq .. angularFreq]
      --        alpha)
      --     canonicalEllipseList
  print pinwheelParams
  -- create fftw plan and filters
  (plan, filters) <-
    makePinwheelConvolution emptyPlan pinwheelParams Convolution
  -- ellipticalPinwheelFilter <-
  --   MP.mapM
  --     (\params -> makeEllipticalPinwheelConvolution' plan params Convolution)
  --     ellipticalPinwheelParamsList
  -- apply filters to patches
  filteredPatchVec1 <-
    MP.mapM
      (applyPinwheelConvolution plan rows cols filters .
       L.map (VU.convert . VU.map (\x -> x :+ 0)) .
       sliceArrayChannel . imageContent)
      patches1
  filteredPatchVec2 <-
    MP.mapM
      (applyPinwheelConvolution plan rows cols filters .
       L.map (VU.convert . VU.map (\x -> x :+ 0)) .
       sliceArrayChannel . imageContent)
      patches2
  -- ellipticalPinwheelFilteredVec1 <-
  --   MP.mapM
  --     (\patch ->
  --        M.mapM
  --          (\f ->
  --             applyPinwheelConvolution plan rows cols f .
  --             L.map (VU.convert . VU.map (\x -> x :+ 0)) .
  --             sliceArrayChannel . imageContent $
  --             patch)
  --          ellipticalPinwheelFilter)
  --     patches1
  -- ellipticalPinwheelFilteredVec2 <-
  --   MP.mapM
  --     (\patch ->
  --        M.mapM
  --          (\f ->
  --             applyPinwheelConvolution plan rows cols f .
  --             L.map (VU.convert . VU.map (\x -> x :+ 0)) .
  --             sliceArrayChannel . imageContent $
  --             patch)
  --          ellipticalPinwheelFilter)
  --     patches2
  let -- pooledVec1 =
      --   parZipWithChunk
      --     parallelParams
      --     rdeepseq
      --     (\a b -> L.foldl1' (L.zipWith (L.zipWith (VS.zipWith (max)))) (a : b))
      --     filteredPatchVec1
      --     ellipticalPinwheelFilteredVec1
      -- pooledVec2 =
      --   parZipWithChunk
      --     parallelParams
      --     rdeepseq
      --     (\a b -> L.foldl1' (L.zipWith (L.zipWith (VS.zipWith (max)))) (a : b))
      --     filteredPatchVec2
      --     ellipticalPinwheelFilteredVec2
      pooledVec1 = filteredPatchVec1
      pooledVec2 = filteredPatchVec2
      pinwheelFeature1 =
        getFeaturePatchP parallelParams 1 affineRegion1 .
        L.map
          (fromUnboxed
             (Z :. (L.length . L.head . L.head $ filteredPatchVec1) :. cols :.
              rows) .
           VS.convert . VS.concat . L.head) $
        pooledVec1
        -- filteredPatchVec1
      pinwheelFeature2 =
        getFeaturePatchP parallelParams 1 affineRegion2 .
        L.map
          (fromUnboxed
             (Z :. (L.length . L.head . L.head $ filteredPatchVec2) :. cols :.
              rows) .
           VS.convert . VS.concat . L.head) $
        pooledVec2
        -- filteredPatchVec2
      pinwheelFeatureStride1 =
        getFeaturePatchP parallelParams stride affineRegion1 .
        L.map
          (fromUnboxed
             (Z :. (L.length . L.head . L.head $ filteredPatchVec1) :. cols :.
              rows) .
           VS.convert . VS.concat . L.head) $
        pooledVec1
        -- filteredPatchVec1
      pinwheelFeatureStride2 =
        getFeaturePatchP parallelParams stride affineRegion2 .
        L.map
          (fromUnboxed
             (Z :. (L.length . L.head . L.head $ filteredPatchVec2) :. cols :.
              rows) .
           VS.convert . VS.concat . L.head) $
        pooledVec2
        -- filteredPatchVec2
  -- compute K-means model
  kmeansModel <-
    kmeans
      parallelParams
      numKmeansCenter
      (folderPath </> "test.kmeans")
      kmeansThreshold .
    L.concatMap detectedRegionFeature $
    affineDataFeature pinwheelFeatureStride1 L.++
    affineDataFeature pinwheelFeatureStride2
  let vladFeatures1 = computeVLADP parallelParams kmeansModel pinwheelFeature1
      vladFeatures2 = computeVLADP parallelParams kmeansModel pinwheelFeature2
      -- (pcaMat, eigValVec, _) =
      --   pcaSVDS
      --     descriptorLength
      --     ((L.map (detectedRegionFeature) . affineDataFeature $ vladFeatures1) L.++
      --      (L.map (detectedRegionFeature) . affineDataFeature $ vladFeatures2))
      -- vladFeaturesPCA1 = descriptorPCAP parallelParams pcaMat vladFeatures1
      -- vladFeaturesPCA2 = descriptorPCAP parallelParams pcaMat vladFeatures2
      -- vladFeaturesPCA1 = descriptorPCAP' parallelParams pcaMat eigValVec vladFeatures1
      -- vladFeaturesPCA2 = descriptorPCAP' parallelParams pcaMat eigValVec vladFeatures2
  writeAffineDescriptor folderPath vladFeatures1
  writeAffineDescriptor folderPath vladFeatures2
  -- writeAffineDescriptor folderPath $ vladFeaturesPCA1
  -- writeAffineDescriptor folderPath $ vladFeaturesPCA2
