module Preprocess.Region
  ( module Types
  , detectedRegionSet
  , detectedRegionBoundary
  ) where

import           Data.List             as L
import           Data.Set              as S
import           Numeric.LinearAlgebra as NL
import           Types
import           Utils

-- ((minX,maxX),(minY,maxY))
{-# INLINE boundary #-}
boundary :: Set (Int,Int) -> ((Int, Int), (Int, Int))
boundary set =
  let (xs, ys) = L.unzip . S.toList $ set
   in ((L.minimum xs, L.maximum xs), (L.minimum ys, L.maximum ys))

{-# INLINE inRange #-}
inRange :: Double -> DetectedRegion -> (Int,Int) -> Bool
inRange scale (DetectedRegion x y a' b' c' _) (i, j) =
  let x' = fromIntegral i - x
      y' = fromIntegral j - y
      a = a' / (scale^2)
      b = b' / (scale^2)
      c = c' / (scale^2)
      normalizedR = a * x' ^ 2 + 2 * b * x' * y' + c * y' ^ 2
   in if normalizedR <= 1
        then True
        else False

{-# INLINE detectedRegionSet #-}
detectedRegionSet :: Double -> Double -> DetectedRegion -> Set (Int, Int)
detectedRegionSet scale delta detectedRegion@(DetectedRegion x y _ _ _ _) =
  let ((minY, maxY), (minX, maxX)) =
        boundary $ detectedRegionBoundary scale delta detectedRegion
   in S.fromList .
      L.map (\(a, b) -> (b, a)) . L.filter (inRange scale detectedRegion) $
      [(i, j) | i <- [minX .. maxX], j <- [minY .. maxY]]

{-# INLINE detectedRegionBoundary #-}
detectedRegionBoundary :: Double -> Double -> DetectedRegion -> Set (Int, Int)
detectedRegionBoundary scale delta (DetectedRegion x y a b c _) =
  let (eigVal, eigVec) = eigSH . trustSym . (2 >< 2) $ [a, b, b, c]
      (l2:l1:[]) = L.map (\x -> 1 / sqrt x) $ NL.toList eigVal
      alpha = atan2 (eigVec `atIndex` (1, 1)) (eigVec `atIndex` (1, 0))
      xs =
        L.map
          (\t -> (scale * l1 * cos t, scale * l2 * sin t))
          [0,delta .. 2 * pi]
      ys =
        L.map
          (\(i, j) ->
             ( x + i * cos alpha + j * sin alpha
             , y + j * cos alpha - i * sin alpha))
          xs
   in S.fromList . L.map (\(i, j) -> (round j, round i)) $ ys
