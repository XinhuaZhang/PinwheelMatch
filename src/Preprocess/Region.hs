module Preprocess.Region
  ( module Types
  , detectedRegionSet
  , detectedRegionBoundary
  ) where

import           Data.List as L
import           Data.Set  as S
import           Types

{-# INLINE radius #-}
radius :: DetectedRegion -> (Int, Int) -> Double
radius (DetectedRegion x y _ _ _) (a, b) =
  (fromIntegral a - x) ^ (2 :: Int) + (fromIntegral b - y) ^ (2 :: Int)

{-# INLINE angleFunctionRad #-}
angleFunctionRad
  :: (Floating a, Ord a)
  => a -> a -> a
angleFunctionRad 0 0 = 0
angleFunctionRad i 0
  | i > 0 = 0.0
  | otherwise = 1.0 * pi
angleFunctionRad 0 j
  | j > 0 = pi / 2.0
  | otherwise = 3.0 * pi / 2.0
angleFunctionRad i j
  | i > 0
  , ratio > 0 = ratio
  | i > 0
  , ratio < 0 = 2.0 * pi + ratio
  | otherwise = pi + ratio
  where
    ratio = atan $ j / i


{-# INLINE boundaryRadius #-}
boundaryRadius :: DetectedRegion -> Double -> Double -> Double
boundaryRadius (DetectedRegion _ _ a b c) scale theta =
  let deltaX = ((a * cos theta) + b * sin theta) * scale
      deltaY = ((b * cos theta) + c * sin theta) * scale
   in deltaX ^ (2 :: Int) + deltaY ^ (2 :: Int)

-- ((minX,maxX),(minY,maxY))
{-# INLINE boundary #-}
boundary :: DetectedRegion -> Double -> ((Double, Double), (Double, Double))
boundary (DetectedRegion x y a b c) scale =
  let thetaX = atan $ b / a
      thetaY = atan $ c / b
   in ( sort
          (x + (a * cos thetaX + b * sin thetaX) * scale)
          (x + (a * cos (thetaX + pi) + b * sin (thetaX + pi)) * scale)
      , sort
          (y + (b * cos thetaY + c * sin thetaY) * scale)
          (y + (b * cos (thetaY + pi) + c * sin (thetaY + pi)) * scale))
  where
    sort m n =
      if m > n
        then (n, m)
        else (m, n)

{-# INLINE inRange #-}
inRange :: Double -> DetectedRegion -> (Int,Int) -> Bool
inRange scale (DetectedRegion x y a b c) (i, j) =
  let x' = fromIntegral i - x
      y' = fromIntegral j - y
      theta = angleFunctionRad x' y'
      t' = atan ((a * tan theta - b) / (c - b * tan theta))
      t =
        if x' >= 0
          then if y' < 0
                  then t' + pi
                  else t'
          else t' + pi
      normalizedR = ((x' - b * sin t) / a) ^ 2 + ((y' - b * cos t) / c) ^ 2
   in if normalizedR <= 1
        then True
        else False
        
{-# INLINE detectedRegionSet' #-}
detectedRegionSet' :: Double -> Double -> DetectedRegion -> [(Int, Int)]
detectedRegionSet' scale delta (DetectedRegion x y a b c) =
  L.map
    (\theta ->
       ( round $ y + (b * (cos theta) + c * sin theta) * scale
       , round $ x + (a * (cos theta) + b * sin theta) * scale)) $
  [0,delta .. 2 * pi]

{-# INLINE detectedRegionSet #-}
detectedRegionSet :: Double -> Double -> DetectedRegion -> Set (Int, Int)
detectedRegionSet scale delta detectedRegion@(DetectedRegion x y _ _ _) =
  let ((minX, maxX), (minY, maxY)) = boundary detectedRegion scale
   in S.fromList . L.concatMap (\s -> detectedRegionSet' s delta detectedRegion) $
      [0,delta .. 1]

{-# INLINE detectedRegionBoundary #-}
detectedRegionBoundary :: Double -> Double -> DetectedRegion -> Set (Int, Int)
detectedRegionBoundary scale delta (DetectedRegion x y a b c) =
  S.fromList .
  L.map
    (\theta ->
       ( round (y + (b * (cos theta) + c * sin theta) * scale)
       , round (x + (a * (cos theta) + b * sin theta) * scale) )) $
  [0,delta .. 2 * pi]
