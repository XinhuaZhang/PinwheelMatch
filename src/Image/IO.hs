module Image.IO
  ( ImageRepa(..)
  , readImagePathList
  , imagePathSource
  , readImageRepa
  , readImageConduit
  , plotImageRepa
  ) where

import           Codec.Picture                hiding (Image)
import           Control.Monad.IO.Class       (liftIO)
import           Control.Monad.Trans.Resource
import           Data.Array.Repa              as R
import           Data.Conduit                 as C
import           Data.Conduit.List            as CL
import           Data.List                    as L
import           Data.Vector.Unboxed          as VU
import           Data.Word
import           GHC.Float

data ImageRepa = Image
  { imageDepth :: Int
  , imageContent :: (Array U DIM3 Double)
  } 

readImagePathList :: FilePath -> IO [String]
readImagePathList = fmap lines . readFile

imagePathSource :: FilePath -> ConduitT () FilePath (ResourceT IO) ()
imagePathSource filePath = do
  pathList <- liftIO $ readImagePathList filePath
  sourceList pathList

readImageRepa :: FilePath -> Bool -> IO ImageRepa
readImageRepa filePath isColor = do
  buffer <- liftIO $ readImage filePath
  case buffer of
    Left msg -> error msg
    Right dImg ->
      return $
      if isColor
        then case dImg of
               ImageY8 img ->
                 Image 8 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) -> fromIntegral $ pixelAt img i j :: Double)
               ImageY16 img ->
                 Image 16 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) -> fromIntegral $ pixelAt img i j :: Double)
               ImageYF img ->
                 Image 64 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) -> float2Double $ pixelAt img i j :: Double)
               ImageRGB8 img ->
                 Image 8 . computeS $
                 fromFunction
                   (Z :. (3 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. k :. j :. i) ->
                       let (PixelRGB8 r g b) = pixelAt img i j
                       in case k of
                            0 -> fromIntegral r
                            1 -> fromIntegral g
                            2 -> fromIntegral b
                            _ -> error "readImageConduit: dimension error.")
               ImageRGB16 img ->
                 Image 16 . computeS $
                 fromFunction
                   (Z :. (3 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. k :. j :. i) ->
                       let (PixelRGB16 r g b) = pixelAt img i j
                       in case k of
                            0 -> fromIntegral r
                            1 -> fromIntegral g
                            2 -> fromIntegral b
                            _ -> error "readImageConduit: dimension error.")
               ImageRGBF img ->
                 Image 64 . computeS $
                 fromFunction
                   (Z :. (3 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. k :. j :. i) ->
                       let (PixelRGBF r g b) = pixelAt img i j
                       in case k of
                            0 -> float2Double r
                            1 -> float2Double g
                            2 -> float2Double b
                            _ -> error "readImageConduit: dimension error.")
               img ->
                 let rgbImg = convertRGB8 img
                 in Image 8 . computeS $
                    fromFunction
                      (Z :. (3 :: Int) :. imageHeight rgbImg :. imageWidth rgbImg)
                      (\(Z :. k :. j :. i) ->
                          let (PixelRGB8 r g b) = pixelAt rgbImg i j
                          in case k of
                               0 -> fromIntegral r
                               1 -> fromIntegral g
                               2 -> fromIntegral b
                               _ -> error "readImageConduit: dimension error.")
        else case dImg of
               ImageY8 img ->
                 Image 8 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) -> fromIntegral $ pixelAt img i j :: Double)
               ImageY16 img ->
                 Image 16 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) -> fromIntegral $ pixelAt img i j :: Double)
               ImageYF img ->
                 Image 64 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) -> float2Double $ pixelAt img i j :: Double)
               ImageRGB8 img ->
                 Image 8 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) ->
                       let (PixelRGB8 r g b) = pixelAt img i j
                       in rgb2Gray
                            (fromIntegral (maxBound :: Word8))
                            (fromIntegral r)
                            (fromIntegral g)
                            (fromIntegral b))
               ImageRGB16 img ->
                 Image 16 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) ->
                       let (PixelRGB16 r g b) = pixelAt img i j
                       in rgb2Gray
                            (fromIntegral (maxBound :: Word16))
                            (fromIntegral r)
                            (fromIntegral g)
                            (fromIntegral b))
               ImageRGBF img ->
                 Image 64 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageHeight img :. imageWidth img)
                   (\(Z :. _ :. j :. i) ->
                       let (PixelRGBF r g b) = pixelAt img i j
                       in rgb2Gray
                            1
                            (float2Double r)
                            (float2Double g)
                            (float2Double b))
               img ->
                 let rgbImg = convertRGB8 img
                 in Image 8 . computeS $
                    fromFunction
                      (Z :. (1 :: Int) :. imageHeight rgbImg :. imageWidth rgbImg)
                      (\(Z :. _ :. j :. i) ->
                          let (PixelRGB8 r g b) = pixelAt rgbImg i j
                          in rgb2Gray
                               (fromIntegral (maxBound :: Word8))
                               (fromIntegral r)
                               (fromIntegral g)
                               (fromIntegral b))


readImageConduit :: Bool -> ConduitT FilePath ImageRepa (ResourceT IO) ()
readImageConduit isColor =
  awaitForever
    (\filePath -> do
       img <- liftIO $ readImageRepa filePath isColor
       yield img)

{-# INLINE rgb2Gray #-}
rgb2Gray :: Double -> Double -> Double -> Double -> Double
rgb2Gray bound r g b
  | yLinear <= 0.0031308 = 12.92 * yLinear * bound
  | otherwise = (1.055 * (yLinear ** (1 / 2.4)) - 0.055) * bound
  where
    yLinear =
      0.2126 * func bound r + 0.7152 * func bound g + 0.0722 * func bound b

{-# INLINE func #-}
func :: Double -> Double -> Double
func bound x
  | y < 0.04045 = y / 12.92
  | otherwise = ((y + 0.055) / 1.055) ** 2.4
  where y = x / bound


{-# INLINE normalizeImageRepa #-}

normalizeImageRepa
  :: ImageRepa -> ImageRepa
normalizeImageRepa i@(Image depth img)
  | maxV == minV =
    Image depth . computeS . R.map (\x -> (fromIntegral $ 2 ^ depth - 1)) $ img
  | otherwise =
    Image depth .
    computeS .
    R.map (\x -> (x - minV) / (maxV - minV) * (fromIntegral $ 2 ^ depth - 1)) $
    img
  where
    maxV = VU.maximum . toUnboxed $ img
    minV = VU.minimum . toUnboxed $ img

plotImageRepa :: FilePath -> ImageRepa -> IO ()
plotImageRepa filePath img@(Image depth x) = do
  let Z :. nfp' :. nyp' :. nxp' = extent x
      normalizedImg = imageContent . normalizeImageRepa $ img
      w =
        case nfp' of
          1 ->
            ImageY8 $
            generateImage
              (\i j ->
                  let v =
                        fromIntegral . round $
                        normalizedImg R.! (Z :. 0 :. j :. i)
                  in v)
              nxp'
              nyp'
          3 ->
            ImageRGB8 $
            generateImage
              (\i j ->
                  let r =
                        fromIntegral . round $
                        normalizedImg R.! (Z :. 0 :. j :. i)
                      g =
                        fromIntegral . round $
                        normalizedImg R.! (Z :. 1 :. j :. i)
                      b =
                        fromIntegral . round $
                        normalizedImg R.! (Z :. 2 :. j :. i)
                  in PixelRGB8 r g b)
              nxp'
              nyp'
          _ ->
            error $
            "Image is neither a gray image nor a color image. There are " L.++
            show nfp' L.++
            " channels."
  savePngImage filePath w
