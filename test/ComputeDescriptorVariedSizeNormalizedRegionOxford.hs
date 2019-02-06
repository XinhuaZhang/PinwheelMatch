module ComputeDescriptorVariedSizeNormalizedRegionOxford where

import           Control.Monad          as M
import           Control.Monad.Parallel as MP
import           Data.Array.Repa        as R
import           Data.Complex
import           Data.List              as L
import           Data.Set               as S
import           Data.Vector.Storable   as VS
import           Data.Vector.Unboxed    as VU
import           Descriptor.Pinwheel
import           Image.IO
import           OxfordAffine
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Parallel
import           Utils

main = do
  args@(imagePath1:patchFolderPath1:detectedRegionFilePath1:imagePath2:patchFolderPath2:detectedRegionFilePath2:scaleStr:deltaStr:radialFreqStr:angularFreqStr:alphaStr:numKmeansCenterStr:kmeansThresholdStr:descriptorLengthStr:strideStr:numThreadStr:batchSizeStr:_) <-
    getArgs
  let folderPath =
        "output/test/ComputeDescriptorVariedSizeNormalizedRegionOxford"
      scale = read scaleStr :: Double
      delta = read deltaStr :: Double
      radialFreq = read radialFreqStr :: Double
      angularFreq = read angularFreqStr :: Double
      alpha = read alphaStr :: Double
      numKmeansCenter = read numKmeansCenterStr :: Int
      kmeansThreshold = read kmeansThresholdStr :: Double
      descriptorLength = read descriptorLengthStr :: Int
      stride = read strideStr :: Int
      numThread = read numThreadStr :: Int
      batchSize = read batchSizeStr :: Int
      parallelParams = ParallelParams numThread batchSize
  createDirectoryIfMissing True folderPath
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
  let pinwheelParamsList1 =
        L.map
          (\(ImageRepa _ img) ->
             let (Z :. nf :. cols :. rows) = extent img
              in PinwheelParams
                   rows
                   cols
                   [0 .. radialFreq]
                   [-angularFreq .. angularFreq]
                   alpha) $
        patches1
      pinwheelParamsList2 =
        L.map
          (\(ImageRepa _ img) ->
             let (Z :. nf :. cols :. rows) = extent img
              in PinwheelParams
                   rows
                   cols
                   [0 .. radialFreq]
                   [-angularFreq .. angularFreq]
                   alpha) $
        patches2
  -- create fftw plan and filters
  plan <-
    M.foldM (makePinwheelConvolutionPlan Convolution) emptyPlan .
    S.toList . S.fromList $
    pinwheelParamsList1 L.++ pinwheelParamsList2
  -- apply filters to patches
  filteredPatchVec1 <-
    MP.mapM
      (\(patch, params@(PinwheelParams rows cols _ _ _)) -> do
         filters <- makePinwheelConvolution' plan params Convolution
         applyPinwheelConvolution plan rows cols filters .
           L.map (VU.convert . VU.map (\x -> x :+ 0)) .
           sliceArrayChannel . imageContent $
           patch) $
    L.zip patches1 pinwheelParamsList1 :: IO [[[VS.Vector Double]]]
  filteredPatchVec2 <-
    MP.mapM
      (\(patch, params@(PinwheelParams rows cols _ _ _)) -> do
         filters <- makePinwheelConvolution' plan params Convolution
         applyPinwheelConvolution plan rows cols filters .
           L.map (VU.convert . VU.map (\x -> x :+ 0)) .
           sliceArrayChannel . imageContent $
           patch) $
    L.zip patches2 pinwheelParamsList2
  let pinwheelFeature1 =
        getFeaturePatchP parallelParams 1 affineRegion1 $
        L.zipWith
          (\filteredPatchVec (PinwheelParams rows cols _ _ _) ->
             fromUnboxed
               (Z :. (L.length . L.head $ filteredPatchVec) :. cols :. rows) .
             VS.convert . VS.concat . L.head $
             filteredPatchVec)
          filteredPatchVec1
          pinwheelParamsList1
      pinwheelFeature2 =
        getFeaturePatchP parallelParams 1 affineRegion2 $
        L.zipWith
          (\filteredPatchVec (PinwheelParams rows cols _ _ _) ->
             fromUnboxed
               (Z :. (L.length . L.head $ filteredPatchVec) :. cols :. rows) .
             VS.convert . VS.concat . L.head $
             filteredPatchVec)
          filteredPatchVec2
          pinwheelParamsList2
      pinwheelFeatureStride1 =
        getFeaturePatchP parallelParams stride affineRegion1 $
        L.zipWith
          (\filteredPatchVec (PinwheelParams rows cols _ _ _) ->
             fromUnboxed
               (Z :. (L.length . L.head $ filteredPatchVec) :. cols :. rows) .
             VS.convert . VS.concat . L.head $
             filteredPatchVec)
          filteredPatchVec1
          pinwheelParamsList1
      pinwheelFeatureStride2 =
        getFeaturePatchP parallelParams stride affineRegion2 $
        L.zipWith
          (\filteredPatchVec (PinwheelParams rows cols _ _ _) ->
             fromUnboxed
               (Z :. (L.length . L.head $ filteredPatchVec) :. cols :. rows) .
             VS.convert . VS.concat . L.head $
             filteredPatchVec)
          filteredPatchVec2
          pinwheelParamsList2
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
