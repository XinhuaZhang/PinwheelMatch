module ComputeDescriptorNormalizedRegionOxfordCascade where

import           Control.Monad                 as M
import           Numeric.LinearAlgebra.Data as LA
import           Control.Monad.Parallel        as MP
import           Data.Array.Repa               as R
import           Data.Complex
import           Data.Conduit                  as C
import           Data.Conduit.List             as CL
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
import           Utils.Parallel                hiding ((.|))

main = do
  args@(imagePath1:patchFolderPath1:detectedRegionFilePath1:imagePath2:patchFolderPath2:detectedRegionFilePath2:scaleStr:deltaStr:radialFreqStr:angularFreqStr:alphaStr:numKmeansCenterStr:kmeansThresholdStr:descriptorLengthStr:pcaLengthStr:strideStr:nStr:numThreadStr:batchSizeStr:_) <-
    getArgs
  let folderPath = "output/test/ComputeDescriptorNormalizedRegionOxfordCascade"
      scale = read scaleStr :: Double
      delta = read deltaStr :: Double
      radialFreq = read radialFreqStr :: Double
      angularFreq = read angularFreqStr :: Double
      alpha = read alphaStr :: Double
      numKmeansCenter = read numKmeansCenterStr :: Int
      kmeansThreshold = read kmeansThresholdStr :: Double
      descriptorLength = read descriptorLengthStr :: Int
      pcaLength = read pcaLengthStr :: Int
      stride = read strideStr :: Int
      n = read nStr :: Int
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
       patchFolderPath1
       (L.length . affineDataFeature $ affineRegion1)
       detectedRegionFilePath1
       (L.length patchPathList2)
       patchFolderPath2
       (L.length . affineDataFeature $ affineRegion2)
       detectedRegionFilePath2)
  img <- readImageRepa (L.head patchPathList1) False
  let (Z :. nf :. cols :. rows) = extent . imageContent $ img
      pinwheelParams =
        PinwheelParams
          rows
          cols
          [0 .. radialFreq]
          [-angularFreq .. angularFreq]
          alpha
  print pinwheelParams
  -- create fftw plan and filters
  (plan, filters) <-
    makePinwheelConvolution emptyPlan pinwheelParams Convolution
  xs1 <-
    runConduitRes $
    CL.sourceList (patchPathList1 L.++ patchPathList2) .| readImageConduit False .|
    firstLayerConduit parallelParams plan filters rows cols .|
    CL.consume
  pcaMat1 <-
    runConduitRes $
    CL.sourceList xs1 .| getPCAColDataConduit parallelParams stride rows cols .|
    pcaSink parallelParams pcaLength
  xs2 <-
    runConduitRes $
    CL.sourceList xs1 .| getPCARowDataConduit parallelParams rows cols .|
    pcaConduit parallelParams pcaMat1 .|
    CL.map (L.map VS.fromList . getPCARowData) .|
    cascadeConvolutionConduit' parallelParams plan filters rows cols .|
    CL.consume
  pcaMat2 <-
    runConduitRes $
    CL.sourceList xs2 .| getPCAColDataConduit parallelParams stride rows cols .|
    pcaSink parallelParams descriptorLength
  -- xs3 <-
  --   runConduitRes $
  --   CL.sourceList xs2 .| getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat2 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit' parallelParams plan filters rows cols .|
  --   CL.consume
  -- pcaMat3 <-
  --   runConduitRes $
  --   CL.sourceList xs3 .| getPCAColDataConduit parallelParams stride rows cols .|
  --   pcaSink parallelParams descriptorLength
  let (as, bs) = L.splitAt (L.length patchPathList1) xs2
  filteredPatchVec1 <-
    runConduitRes $
    CL.sourceList as .| getPCARowDataConduit parallelParams rows cols .|
    pcaConduit parallelParams pcaMat2 .|
    CL.consume
  filteredPatchVec2 <-
    runConduitRes $
    CL.sourceList bs .| getPCARowDataConduit parallelParams rows cols .|
    pcaConduit parallelParams pcaMat2 .|
    CL.consume
  -- pcaMat1 <-
  --   runConduitRes $
  --   CL.sourceList (patchPathList1 L.++ patchPathList2) .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCAColDataConduit parallelParams stride rows cols .|
  --   pcaSink parallelParams pcaLength
  -- print . pcaEigVal $ pcaMat1
  -- pcaMat2 <-
  --   runConduitRes $
  --   CL.sourceList (patchPathList1 L.++ patchPathList2) .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat1 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit parallelParams plan filters rows cols .|
  --   getPCAColDataConduit parallelParams stride rows cols .|
  --   pcaSink parallelParams (pcaLength ^ (2 :: Int))
  -- print . pcaEigVal $ pcaMat2
  -- pcaMat3 <-
  --   runConduitRes $
  --   CL.sourceList (patchPathList1 L.++ patchPathList2) .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat1 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat2 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit' parallelParams plan filters rows cols .|
  --   getPCAColDataConduit parallelParams stride rows cols .|
  --   pcaSink parallelParams descriptorLength
  -- print . pcaEigVal $ pcaMat3
  -- filteredPatchVec1 <-
  --   runConduitRes $
  --   CL.sourceList patchPathList1 .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat1 .|
  --   CL.consume
  -- filteredPatchVec2 <-
  --   runConduitRes $
  --   CL.sourceList patchPathList2 .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat1 .|
  --   CL.consume
  -- filteredPatchVec1 <-
  --   runConduitRes $
  --   CL.sourceList patchPathList1 .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat1 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat2 .|
  --   CL.consume
  -- filteredPatchVec2 <-
  --   runConduitRes $
  --   CL.sourceList patchPathList2 .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat1 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat2 .|
  --   CL.consume
  -- filteredPatchVec1 <-
  --   runConduitRes $
  --   CL.sourceList patchPathList1 .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat1 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat2 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit' parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat3 .|
  --   CL.consume
  -- filteredPatchVec2 <-
  --   runConduitRes $
  --   CL.sourceList patchPathList2 .| readImageConduit False .|
  --   firstLayerConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat1 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat2 .|
  --   CL.map (L.map VS.fromList . getPCARowData) .|
  --   cascadeConvolutionConduit' parallelParams plan filters rows cols .|
  --   getPCARowDataConduit parallelParams rows cols .|
  --   pcaConduit parallelParams pcaMat3 .|
  --   CL.consume
  let pinwheelFeature1 =
        getFeaturePatchP parallelParams 1 affineRegion1 .
        L.map
          (fromListUnboxed
             (Z :. (L.length . getPCARowData . L.head $ filteredPatchVec1) :.
              cols :.
              rows) .
           L.concat . getPCARowData) $
        filteredPatchVec1
      pinwheelFeature2 =
        getFeaturePatchP parallelParams 1 affineRegion2 .
        L.map
          (fromListUnboxed
             (Z :. (L.length . getPCARowData . L.head $ filteredPatchVec1) :.
              cols :.
              rows) .
           L.concat . getPCARowData) $
        filteredPatchVec2
      pinwheelFeatureStride1 =
        getFeaturePatchP parallelParams stride affineRegion1 .
        L.map
          (fromListUnboxed
             (Z :. (L.length . getPCARowData . L.head $ filteredPatchVec1) :.
              cols :.
              rows) .
           L.concat . getPCARowData) $
        filteredPatchVec1
      pinwheelFeatureStride2 =
        getFeaturePatchP parallelParams stride affineRegion2 .
        L.map
          (fromListUnboxed
             (Z :. (L.length . getPCARowData . L.head $ filteredPatchVec1) :.
              cols :.
              rows) .
           L.concat . getPCARowData) $
        filteredPatchVec2
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
  let vladFeatures1' = computeVLADP parallelParams kmeansModel pinwheelFeature1
      vladFeatures2' = computeVLADP parallelParams kmeansModel pinwheelFeature2
      (vladFeatures1, vladFeatures2) =
        zca
          parallelParams
          ((div descriptorLength 5) * numKmeansCenter)
          vladFeatures1'
          vladFeatures2'
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
