module Types
  ( module Types
  , module Numeric.LinearAlgebra
  ) where

import Numeric.LinearAlgebra

newtype AffineMatrix =
  AffineMatrix (Matrix Double)
  deriving (Show)

data DetectedRegion = DetectedRegion
  { detectedRegionX :: Double
  , detectedRegionY :: Double
  , detectedRegionA :: Double
  , detectedRegionB :: Double
  , detectedRegionC :: Double
  } deriving (Show)
