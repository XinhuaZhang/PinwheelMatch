module ComputeDescriptor where

import           Control.Monad        as M
import           Data.Array.Repa      as R
import           Data.Complex
import           Data.List            as L
import           Data.Set             as S
import           Data.Vector.Storable as VS
import           Data.Vector.Unboxed  as VU
import           Descriptor.Pinwheel
import           Image.IO
import           OxfordAffine
import           System.Directory
import           System.Environment
import           System.FilePath
import           Utils.Parallel

main = do
  args@(imagePath1:detectedRegionFilePath1:imagePath2:detectedRegionFilePath2:scaleStr:deltaStr:radialFreqStr:angularFreqStr:alphaStr:numKmeansCenterStr:kmeansThresholdStr:descriptorLengthStr:strideStr:numThreadStr:batchSizeStr:_) <-
    getArgs
  (Image 8 img1) <- readImageRepa imagePath1 False
  (Image 8 img2) <- readImageRepa imagePath2 False
  let folderPath = "output/test/ComputeDescriptor"
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
      (Z :. nf :. cols :. rows) = extent img1
      pinwheelParams =
        PinwheelParams
          rows
          cols
          [0 .. radialFreq]
          [-angularFreq .. angularFreq]
          alpha
      imgVec1 =
        L.map
          (\i ->
             VU.convert .
             VU.map (\x -> x :+ 0) . toUnboxed . computeS . R.slice img1 $
             (Z :. i :. R.All :. R.All))
          [0 .. nf - 1]
      imgVec2 =
        L.map
          (\i ->
             VU.convert .
             VU.map (\x -> x :+ 0) . toUnboxed . computeS . R.slice img2 $
             (Z :. i :. R.All :. R.All))
          [0 .. nf - 1]
  createDirectoryIfMissing True folderPath
  print pinwheelParams
  -- read affine regions
  affineRegion1 <-
    computeAffineRegionFromFile imagePath1 detectedRegionFilePath1 scale delta
  affineRegion2 <-
    computeAffineRegionFromFile imagePath2 detectedRegionFilePath2 scale delta
  -- create fftw plan and filters
  (plan, filters) <-
    makePinwheelConvolution emptyPlan pinwheelParams Convolution
  -- apply filters to image
  filteredVec1 <- applyPinwheelConvolutionP plan rows cols filters imgVec1
  filteredVec2 <- applyPinwheelConvolutionP plan rows cols filters imgVec2
  let filteredImage1 =
        fromUnboxed (Z :. (L.length . L.head $ filteredVec1) :. cols :. rows) .
        VS.convert . VS.concat . L.head $
        filteredVec1
      filteredImage2 =
        fromUnboxed (Z :. (L.length . L.head $ filteredVec2) :. cols :. rows) .
        VS.convert . VS.concat . L.head $
        filteredVec2
      pinwheelFeature1 =
        getFeatureFromAffineRegionP
          parallelParams
          (cols - 1, rows - 1)
          1
          filteredImage1
          affineRegion1
      pinwheelFeature2 =
        getFeatureFromAffineRegionP
          parallelParams
          (cols - 1, rows - 1)
          1
          filteredImage2
          affineRegion2
      pinwheelFeature1' =
        getFeatureFromAffineRegionP
          parallelParams
          (cols - 1, rows - 1)
          stride
          filteredImage1
          affineRegion1
      pinwheelFeature2' =
        getFeatureFromAffineRegionP
          parallelParams
          (cols - 1, rows - 1)
          stride
          filteredImage2
          affineRegion2
  -- compute K-means model
  kmeansModel <-
    kmeans
      parallelParams
      numKmeansCenter
      (folderPath </> "test.kmeans")
      kmeansThreshold .
    L.concatMap detectedRegionFeature $
    affineDataFeature pinwheelFeature1' L.++ affineDataFeature pinwheelFeature2'
  -- compute VLAD features
  let vladFeatures1 = computeVLADP parallelParams kmeansModel pinwheelFeature1
      vladFeatures2 = computeVLADP parallelParams kmeansModel pinwheelFeature2
      (pcaMat, eigValVec, _) =
        pcaSVDS
          descriptorLength
          ((L.map (detectedRegionFeature) . affineDataFeature $ vladFeatures1) L.++
           (L.map (detectedRegionFeature) . affineDataFeature $ vladFeatures2))
      vladFeaturesPCA1 = descriptorPCAP parallelParams pcaMat vladFeatures1
      vladFeaturesPCA2 = descriptorPCAP parallelParams pcaMat vladFeatures2
      -- vladFeaturesPCA1 = descriptorPCAP' parallelParams pcaMat eigValVec vladFeatures1
      -- vladFeaturesPCA2 = descriptorPCAP' parallelParams pcaMat eigValVec vladFeatures2
  writeAffineDescriptor folderPath vladFeatures1
  writeAffineDescriptor folderPath vladFeatures2
  -- writeAffineDescriptor folderPath $ vladFeaturesPCA1
  -- writeAffineDescriptor folderPath $ vladFeaturesPCA2
