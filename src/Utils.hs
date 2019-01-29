module Utils where

import           Data.Array.Repa     as R
import           Data.List           as L
import           Data.Vector.Unboxed as VU

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

{-# INLINE l2norm #-}
l2norm :: VU.Vector Double -> VU.Vector Double
l2norm vec
  | s == 0 = vec
  | otherwise = VU.map (/ s) vec
  where
    s = sqrt . VU.sum . VU.map (^ (2 :: Int)) $ vec

{-# INLINE sliceArrayChannel #-}
sliceArrayChannel :: (Unbox e) => R.Array U DIM3 e -> [VU.Vector e]
sliceArrayChannel arr =
  let (Z :. nf :. _ :. _) = extent arr
   in L.map
        (\i -> toUnboxed . computeS . R.slice arr $ (Z :. i :. R.All :. R.All))
        [0 .. nf - 1]
        
{-# INLINE sliceArrayPosition #-}
sliceArrayPosition :: (Unbox e) => R.Array U DIM3 e -> [VU.Vector e]
sliceArrayPosition arr =
  let (Z :. _ :. cols :. rows) = extent arr
   in L.map
        (\(i, j) -> toUnboxed . computeS . R.slice arr $ (Z :. R.All :. i :. j))
        [(i, j) | i <- [0 .. cols - 1], j <- [0 .. rows - 1]]
        
{-# INLINE sliceArrayPositionStride #-}
sliceArrayPositionStride ::
     (Unbox e) => Int -> R.Array U DIM3 e -> [VU.Vector e]
sliceArrayPositionStride stride arr =
  let (Z :. _ :. cols :. rows) = extent arr
   in L.map
        (\(i, j) -> toUnboxed . computeS . R.slice arr $ (Z :. R.All :. i :. j))
        [(i, j) | i <- [0,stride .. cols - 1], j <- [0,stride .. rows - 1]]
