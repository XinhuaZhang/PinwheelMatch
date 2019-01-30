module PlotEllipsePinwheel where

import           Control.Monad              as M
import           Data.Array.Repa            as R
import           Data.Array.ST
import           Data.Array.Unboxed         as Arr
import           Data.List                  as L
import           Data.Set                   as S
import           Descriptor.EllipsePinwheel
import           Image.IO
import           OxfordAffine.Parser
import           OxfordAffine.Region
import           System.Directory
import           System.Environment
import           System.FilePath

main = do
  args@(detectedRegionFilePath:scaleStr:deltaStr:radialFreqStr:angularFreqStr:alphaStr:idxStr:_) <-
    getArgs
  xs <- parseDetectedRegionFile detectedRegionFilePath
  let folderPath = "output/test/PlotEllipsePinwheel/"
      scale = read scaleStr :: Double
      delta = read deltaStr :: Double
      radialFreq = read radialFreqStr :: Double
      angularFreq = read angularFreqStr :: Double
      alpha = read alphaStr :: Double
      idx = read idxStr :: Int
      z@(DetectedRegion x y _ _ _ _) = (L.reverse xs) !! idx
      ellipse = computeCanonicalEllipse z
      ellipseSet = detectedRegionBoundary scale delta z
      ((xMin, xMax), (yMin, yMax)) =
        boundary . detectedRegionFeature $ ellipseSet
      cols = xMax - xMin + 1
      rows = yMax - yMin + 1
      arr =
        (runSTUArray $ do
           arr <-
             newArray
               ((0 :: Int, 0 :: Int, 0 :: Int), (0 :: Int, cols - 1, rows - 1))
               0
           M.mapM_ (\(i, j) -> writeArray arr (0, i, j) 255) .
             L.filter (inRange ((0, 0), (cols - 1, rows - 1))) .
             L.map
               (\(i, j) ->
                  (i + (div cols 2) - round x, j + (div rows 2) - round y)) $
             S.toList . detectedRegionFeature $ ellipseSet
           return arr) :: UArray (Int, Int, Int) Double
  createDirectoryIfMissing True folderPath
  plotImageRepa (folderPath </> "ellipse.png") .
    ImageRepa 8 . fromListUnboxed (Z :. (1 :: Int) :. cols :. rows) . Arr.elems $
    arr
  let params =
        EllpisePinwheelParams 64 64 [radialFreq] [angularFreq] alpha ellipse
      filter = L.head . L.head $ makeEllpisePinwheelExpansion params 32 32
  print ellipse
  plotImageRepaComplex (folderPath </> "ellipsePinwheel.png") .
    ImageRepa 8 . fromUnboxed (Z :. (1 :: Int) :. 64 :. 64) $
    filter
