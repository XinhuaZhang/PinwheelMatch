{-# LANGUAGE DeriveFunctor #-}
module Types
  ( module Types
  , module Numeric.LinearAlgebra
  ) where

import           Control.DeepSeq
import           Data.Set
import           Data.Vector.Unboxed   as VU
import           Numeric.LinearAlgebra

newtype AffineMatrix =
  AffineMatrix (Matrix Double)
  deriving (Show)

data DetectedRegion a = DetectedRegion
  { detectedRegionX       :: Double
  , detectedRegionY       :: Double
  , detectedRegionA       :: Double
  , detectedRegionB       :: Double
  , detectedRegionC       :: Double
  , detectedRegionFeature :: a
  } deriving (Show,Functor)

instance (NFData a) => NFData (DetectedRegion a) where
  rnf (DetectedRegion x y a b c z) =
    rnf x `seq` rnf y `seq` rnf a `seq` rnf b `seq` rnf c `seq` rnf z `seq` ()

data AffineData a = AffineData
  { affineDataImagePath :: FilePath
  , affineDataFeature   :: [DetectedRegion a]
  } deriving (Functor)

instance (NFData a) => NFData (AffineData a) where
  rnf (AffineData x xs) = rnf x `seq` rnf xs `seq` ()

type AffineRegionParameter = AffineData ()                   -- x,y,a,b,c
type AffineRegion = AffineData (Set (Int, Int))              -- location
type AffineRegionFeature = AffineData [VU.Vector Double]     -- pinwheel feature
type AffineDescriptor = AffineData (VU.Vector Double)        -- VLAD feature

data CanonicalEllipse = CanonicalEllipse
  { canonicalEllipseA     :: Double
  , canonicalEllipseAlpha :: Double
  } deriving (Read, Show)

data PolarEllipseIndex = PolarEllipseIndex
  { polarEllipseIndexR     :: Double
  , polarEllipseIndexTheta :: Double
  } deriving (Read, Show)
