module PlotDetectedRegion where
import           Control.Monad               as M
import           Control.Parallel.Strategies
import           Data.Array.Repa             as R
import           Data.Array.ST
import           Data.Array.Unboxed          as Arr
import           Data.List                   as L
import           Data.Set                    as S
import           Image.IO
import           OxfordAffine.Parser
import           OxfordAffine.Region
import           System.Directory
import           System.Environment
import           System.FilePath

main = do
  args@(imagePath:detectedRegionFilePath:scaleStr:deltaStr:_) <- getArgs
  (ImageRepa 8 img) <- readImageRepa imagePath True
  xs <- parseDetectedRegionFile detectedRegionFilePath
  let folderPath = "output/test/PlotDetectedRegion/"
      scale = read scaleStr :: Double
      delta = read deltaStr :: Double
      (Z :. nf :. cols :. rows) = extent img
      ys =
        L.concat $
        parMap
          rdeepseq
          (S.toList . detectedRegionFeature . detectedRegionBoundary scale delta)
          xs
      zs = L.filter (\idx -> inRange ((0, 0), (cols - 1, rows - 1)) idx) ys
      newArray =
        (runSTUArray $ do
           arr <-
             newListArray ((0, 0, 0), (nf - 1, cols - 1, rows - 1)) . R.toList $
             img
           M.mapM_
             (\(x, y) -> do
                writeArray arr (0, x, y) 255
                writeArray arr (1, x, y) 255
                writeArray arr (2, x, y) 255)
             zs
           return arr) :: UArray (Int, Int, Int) Double
      ys1 =
        L.concat $
        parMap
          rdeepseq
          (S.toList . detectedRegionFeature . detectedRegionSet scale delta)
          xs
      zs1 = L.filter (\idx -> inRange ((0, 0), (cols - 1, rows - 1)) idx) ys1
      newArray1 =
        (runSTUArray $ do
           arr <-
             newListArray ((0, 0, 0), (nf - 1, cols - 1, rows - 1)) . R.toList $
             img
           M.mapM_
             (\(x, y) -> do
                writeArray arr (0, x, y) 255
                writeArray arr (1, x, y) 255
                writeArray arr (2, x, y) 255)
             zs1
           return arr) :: UArray (Int, Int, Int) Double
  createDirectoryIfMissing True folderPath
  plotImageRepa (folderPath L.++ takeBaseName imagePath L.++ ".png") .
    ImageRepa 8 . fromListUnboxed (extent img) . Arr.elems $
    newArray
  plotImageRepa (folderPath L.++ takeBaseName imagePath L.++ "_fill.png") .
    ImageRepa 8 . fromListUnboxed (extent img) . Arr.elems $
    newArray1
