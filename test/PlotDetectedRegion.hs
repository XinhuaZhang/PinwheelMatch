import           Control.Monad               as M
import           Control.Parallel.Strategies
import           Data.Array.Repa             as R
import           Data.Array.ST
import           Data.Array.Unboxed          as Arr
import           Data.List                   as L
import           Data.Set                    as S
import           Image.IO
import           Preprocess.Parser
import           Preprocess.Region
import           System.Environment
import           System.FilePath

main = do
  args@(imagePath:detectedRegionFilePath:scaleStr:deltaStr:_) <- getArgs
  (Image 8 img) <- readImageRepa imagePath True
  xs <- parseDetectedRegionFile detectedRegionFilePath
  let scale = read scaleStr :: Double
      delta = read deltaStr :: Double
      (Z :. nf :. rows :. cols) = extent img
      ys =
        L.concat $
        parMap rdeepseq (S.toList . detectedRegionBoundary scale delta) xs
      zs = L.filter (\idx -> inRange ((0, 0), (rows - 1, cols - 1)) idx) ys
      newArray =
        (runSTUArray $ do
           arr <-
             newListArray ((0, 0, 0), (nf - 1, rows - 1, cols - 1)) . R.toList $
             img
           M.mapM_
             (\(x, y) -> do
                writeArray arr (0, x, y) 255
                writeArray arr (1, x, y) 255
                writeArray arr (2, x, y) 255)
             zs
           return arr) :: UArray (Int, Int, Int) Double
      ys1 =
        L.concat $ parMap rdeepseq (S.toList . detectedRegionSet scale delta) xs
      zs1 = L.filter (\idx -> inRange ((0, 0), (rows - 1, cols - 1)) idx) ys1
      newArray1 =
        (runSTUArray $ do
           arr <-
             newListArray ((0, 0, 0), (nf - 1, rows - 1, cols - 1)) . R.toList $
             img
           M.mapM_
             (\(x, y) -> do
                writeArray arr (0, x, y) 255
                writeArray arr (1, x, y) 255
                writeArray arr (2, x, y) 255)
             zs1
           return arr) :: UArray (Int, Int, Int) Double
  plotImageRepa ("output/" L.++ takeBaseName imagePath L.++ ".png") .
    Image 8 . fromListUnboxed (extent img) . Arr.elems $
    newArray
  plotImageRepa ("output/" L.++ takeBaseName imagePath L.++ "_fill.png") .
    Image 8 . fromListUnboxed (extent img) . Arr.elems $
    newArray1
