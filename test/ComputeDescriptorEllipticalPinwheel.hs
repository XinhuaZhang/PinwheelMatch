module ComputeDescriptorEllipticalPinwheel where

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
import           Utils.Parallel

main = do
  args@(imagePath1:detectedRegionFilePath1:imagePath2:detectedRegionFilePath2:scaleStr:deltaStr:radialFreqStr:angularFreqStr:alphaStr:numKmeansCenterStr:kmeansThresholdStr:descriptorLengthStr:strideStr:deltaThetaStr:deltaAStr:maxAStr:numThreadStr:batchSizeStr:_) <-
    getArgs
  (ImageRepa 8 img1) <- readImageRepa imagePath1 False
  (ImageRepa 8 img2) <- readImageRepa imagePath2 False
  let folderPath = "output/test/ComputeDescriptorEllipticalPinwheel"
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
      (Z :. nf :. cols :. rows) = extent img1
      pinwheelParams =
        PinwheelParams
          rows
          cols
          [0 .. radialFreq]
          [-angularFreq .. angularFreq]
          alpha
      canonicalEllipseList =
        [ CanonicalEllipse a theta
        | a <- [1 + deltaA,1 + 2 * deltaA .. maxA]
        , theta <- [0,deltaTheta .. pi - deltaTheta]
        ]
      ellipticalPinwheelParamsList =
        L.map
          (EllipticalPinwheelParams
             rows
             cols
             [0 .. radialFreq]
             [-angularFreq .. angularFreq]
             alpha)
          canonicalEllipseList
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
  print pinwheelParams
  print
    [ (a, theta)
    | a <- [1 + deltaA,1 + 2 * deltaA .. maxA]
    , theta <- [0,deltaTheta .. pi - deltaTheta]
    ]
  createDirectoryIfMissing True folderPath
  -- read affine regions
  affineRegion1 <-
    computeAffineRegionFromFile imagePath1 detectedRegionFilePath1 scale delta
  affineRegion2 <-
    computeAffineRegionFromFile imagePath2 detectedRegionFilePath2 scale delta
  -- create fftw plan and filters
  plan <- makePinwheelConvolutionPlan Convolution emptyPlan pinwheelParams
  pinwheelFilter <- makePinwheelConvolution' plan pinwheelParams Convolution
  ellipticalPinwheelFilter <-
    MP.mapM
      (\params -> makeEllipticalPinwheelConvolution' plan params Convolution)
      ellipticalPinwheelParamsList
  -- apply filters to image
  pinwheelFilteredVec1 <-
    applyPinwheelConvolutionP plan rows cols pinwheelFilter imgVec1
  pinwheelFilteredVec2 <-
    applyPinwheelConvolutionP plan rows cols pinwheelFilter imgVec2
  ellipticalPinwheelFilteredVec1 <-
    M.mapM
      (\filter -> applyPinwheelConvolutionP plan rows cols filter imgVec1)
      ellipticalPinwheelFilter
  ellipticalPinwheelFilteredVec2 <-
    M.mapM
      (\filter -> applyPinwheelConvolutionP plan rows cols filter imgVec2)
      ellipticalPinwheelFilter
  -- pooling
  let pooledVec1 =
        L.foldl1'
          (L.zipWith (parZipWith rdeepseq (VS.zipWith (+))))
          (pinwheelFilteredVec1 : ellipticalPinwheelFilteredVec1)
      pooledVec2 =
        L.foldl1'
          (L.zipWith (parZipWith rdeepseq (VS.zipWith (+))))
          (pinwheelFilteredVec2 : ellipticalPinwheelFilteredVec2)
      filteredImage1 =
        fromUnboxed
          (Z :. (L.length . L.head . L.head $ ellipticalPinwheelFilteredVec1) :.
           cols :.
           rows) .
        VS.convert . VS.concat . L.head $
        pooledVec1
      filteredImage2 =
        fromUnboxed
          (Z :. (L.length . L.head . L.head $ ellipticalPinwheelFilteredVec2) :.
           cols :.
           rows) .
        VS.convert . VS.concat . L.head $
        pooledVec2
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
        getFeatureFromArrayP parallelParams stride filteredImage1
      pinwheelFeature2' =
        getFeatureFromArrayP parallelParams stride filteredImage2
  -- compute K-means model
  kmeansModel <-
    kmeans
      parallelParams
      numKmeansCenter
      (folderPath </> "test.kmeans")
      kmeansThreshold $
    pinwheelFeature1' L.++ pinwheelFeature2'
  -- compute VLAD features
  let vladFeatures1 = computeVLADP parallelParams kmeansModel pinwheelFeature1
      vladFeatures2 = computeVLADP parallelParams kmeansModel pinwheelFeature2
      -- (pcaMat, eigValVec, _) =
      --   pcaSVDS
      --     descriptorLength
      --     ((L.map detectedRegionFeature . affineDataFeature $ vladFeatures1) L.++
      --      (L.map detectedRegionFeature . affineDataFeature $ vladFeatures2))
      -- vladFeaturesPCA1 = descriptorPCAP parallelParams pcaMat vladFeatures1
      -- vladFeaturesPCA2 = descriptorPCAP parallelParams pcaMat vladFeatures2
      -- vladFeaturesPCA1 = descriptorPCAP' parallelParams pcaMat eigValVec vladFeatures1
      -- vladFeaturesPCA2 = descriptorPCAP' parallelParams pcaMat eigValVec vladFeatures2
  writeAffineDescriptor folderPath vladFeatures1
  writeAffineDescriptor folderPath vladFeatures2
  -- writeAffineDescriptor folderPath $ vladFeaturesPCA1
  -- writeAffineDescriptor folderPath $ vladFeaturesPCA2
