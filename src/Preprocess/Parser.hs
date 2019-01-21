module Preprocess.Parser
  ( parseAffineMatrixFile
  , parseDetectedRegionFile
  ) where

import           Data.Char
import           Data.List      as L
import           Data.Text      as T
import           Data.Text.IO   as T
import           Data.Text.Read as T
import           System.IO
import           Types

{-# INLINE  parseDouble #-}
parseDouble :: Text -> [Double]
parseDouble txt
  | T.null xs = []
  | otherwise =
    case (T.double xs) of
      Left msg -> error msg
      Right (y, ys) -> y : parseDouble ys
  where
    xs =
      T.dropWhile
        (\c ->
           if (not . isDigit $ c) && (c /= '-')
             then True
             else False)
        txt

{-# INLINE parseAffineMatrix #-}
parseAffineMatrix :: Text -> AffineMatrix
parseAffineMatrix txt =
  let xs = parseDouble $ txt
      size = L.length xs
   in if size == 9
        then AffineMatrix . (3 >< 3) $ xs
        else error $
             "parseAffineMatrix: Matrix size " L.++ show size L.++ " /= 9."

{-# INLINE parseAffineMatrixFile #-}
parseAffineMatrixFile :: FilePath -> IO AffineMatrix
parseAffineMatrixFile filePath =
  withFile filePath ReadMode $ \h -> do
    txt <- T.hGetContents h
    return . parseAffineMatrix $ txt

{-# INLINE parseDetectedRegion #-}
parseDetectedRegion :: Text -> [DetectedRegion]
parseDetectedRegion =
  L.map ((\(x:y:a:b:c:[]) -> DetectedRegion x y a b c) . parseDouble) . T.lines

{-# INLINE parseDetectedRegionFile #-}
parseDetectedRegionFile :: FilePath -> IO [DetectedRegion]
parseDetectedRegionFile filePath =
  withFile filePath ReadMode $ \h -> do
    txt <- T.hGetContents h
    return . parseDetectedRegion $ txt
