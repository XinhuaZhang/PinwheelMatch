module Utils where

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
