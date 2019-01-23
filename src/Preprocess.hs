module Preprocess where

import           Control.Monad               as M
import           Control.Parallel.Strategies
import           Data.List                   as L
import           Data.Vector.Unboxed         as VU
import           Preprocess.Parser
import           Preprocess.Region
import           System.FilePath
import           System.IO                   as IO
import           Text.Printf

{-# INLINE computeAffineRegionFromFile #-}
computeAffineRegionFromFile ::
     FilePath -> FilePath -> Double -> Double -> IO AffineRegion
computeAffineRegionFromFile imagePath detectedRegionFilePath scale delta = do
  xs <- parseDetectedRegionFile detectedRegionFilePath
  let ys = parMap rdeepseq (detectedRegionSet scale delta) xs
  return . AffineRegion imagePath $ L.zip xs ys

{-# INLINE writeAffineDescriptor #-}
writeAffineDescriptor :: FilePath -> AffineDescriptor -> IO ()
writeAffineDescriptor folderPath (AffineDescriptor imagePath xs) = do
  withFile (folderPath </> takeBaseName imagePath L.++ ".pinwheel") WriteMode $ \h -> do
    IO.hPutStrLn h (show . VU.length . detectedRegionDescriptor . L.head $ xs)
    IO.hPutStrLn h (show . L.length $ xs)
    M.mapM_
      (\(DetectedRegion x y a b c vec) ->
         let str =
               printf
                 "%f %f %f %f %f%s\n"
                 x
                 y
                 a
                 b
                 c
                 (L.concatMap (\v -> " " L.++ show v) . VU.toList $ vec)
          in IO.hPutStr h str)
      xs
