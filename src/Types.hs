module Types
  ( module Types
  , module Numeric.LinearAlgebra
  ) where

import           Data.Set
import           Data.Vector.Unboxed   as VU
import           Numeric.LinearAlgebra

newtype AffineMatrix =
  AffineMatrix (Matrix Double)
  deriving (Show)

data DetectedRegion = DetectedRegion
  { detectedRegionX          :: Double
  , detectedRegionY          :: Double
  , detectedRegionA          :: Double
  , detectedRegionB          :: Double
  , detectedRegionC          :: Double
  , detectedRegionDescriptor :: VU.Vector Double
  } deriving (Show)

data AffineRegion = AffineRegion
  { affineRegionFilePath :: FilePath
  , affineRegionSets     :: [(DetectedRegion, Set (Int, Int))]
  }

data AffineDescriptor = AffineDescriptor
  { affineDescriptorImagePath :: FilePath
  , affineDescriptorRegion    :: [DetectedRegion]
  }
