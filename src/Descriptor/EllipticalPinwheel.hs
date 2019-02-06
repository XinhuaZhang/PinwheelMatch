module Descriptor.EllipticalPinwheel
  ( ConvolutionalType(..)
  , EllipticalPinwheelParams(..)
  , makeEllipticalPinwheelExpansion
  , makeEllipticalPinwheelConvolution
  , makeEllipticalPinwheelConvolution'
  , module DFT.Plan
  ) where

import           Control.Arrow
import           Control.Monad               as M
import           Control.Monad.Parallel      as MP
import           Control.Parallel.Strategies
import           Data.Complex                as C
import           Data.List                   as L
import           Data.Vector.Storable        as VS
import           Data.Vector.Unboxed         as VU
import           Descriptor.Pinwheel         (ConvolutionalType (..))
import           DFT.Plan
import           OxfordAffine.Region
import           Types
import           Utils

data EllipticalPinwheelParams = EllipticalPinwheelParams
  { ellipticalPinwheelParamsRows             :: Int
  , ellipticalPinwheelParamsCols             :: Int
  , ellipticalPinwheelParamsRadialFreq       :: [Double]
  , ellipticalPinwheelParamsAngularFreq      :: [Double]
  , ellipticalPinwheelParamsAlpha            :: Double
  , ellipticalPinwheelParamsCanonicalEllipse :: CanonicalEllipse
  } deriving (Read, Show)

{-# INLINE ellipticalPinwheel #-}
ellipticalPinwheel ::
     CanonicalEllipse
  -> Double
  -> Double
  -> Double
  -> Int
  -> Int
  -> Complex Double
ellipticalPinwheel canonicalEllipse rf af alpha x y
  | r <= (0.75 / pi * (abs af)) || r <= (0.75 / pi * (abs rf)) = 0
  | otherwise = (((r) :+ 0) ** (alpha :+ (-rf))) * exp (0 :+ ((-af) * theta))
  where
    PolarEllipseIndex r theta = computePolarEllipse (x, y) canonicalEllipse

-- For display purpose
{-# INLINE makeFilter #-}
makeFilter :: Int -> Int -> Int -> Int -> (Int -> Int -> a) -> [a]
makeFilter rows cols rCenter cCenter f =
  [f (r - rCenter) (c - cCenter) | c <- [0 .. cols - 1], r <- [0 .. rows - 1]]

{-# INLINE makeEllipticalPinwheelExpansion #-}
makeEllipticalPinwheelExpansion ::
     EllipticalPinwheelParams -> Int -> Int -> [[VU.Vector (Complex Double)]]
makeEllipticalPinwheelExpansion (EllipticalPinwheelParams rows cols rfs afs alpha canonicalEllipse) rCenter cCenter =
  [ [ VU.fromList $
  makeFilter
    rows
    cols
    rCenter
    cCenter
    (ellipticalPinwheel canonicalEllipse rf af alpha)
  | af <- afs
  ]
  | rf <- rfs
  ]

-- For convolution using convolution theorem
{-# INLINE makeFilterConvolution #-}
makeFilterConvolution :: Int -> Int -> (Int -> Int -> a) -> [a]
makeFilterConvolution rows cols f =
  [ let x =
          if r < (rows `div` 2)
            then r
            else r - rows
        y =
          if c < (cols `div` 2)
            then c
            else c - cols
     in f x y
  | c <- [0 .. cols - 1]
  , r <- [0 .. rows - 1]
  ]

{-# INLINE conjugateFunc #-}
conjugateFunc ::
     ConvolutionalType -> ([Complex Double] -> [Complex Double])
conjugateFunc x =
  case x of
    Convolution      -> id
    Crosscorrelation -> Prelude.map conjugate

-- make both fftw plan and EllipticalPinwheel filters
{-# INLINE makeEllipticalPinwheelConvolution #-}
makeEllipticalPinwheelConvolution ::
     DFTPlan
  -> EllipticalPinwheelParams
  -> ConvolutionalType
  -> IO (DFTPlan, [[VS.Vector (Complex Double)]])
makeEllipticalPinwheelConvolution plan (EllipticalPinwheelParams rows cols rfs afs alpha canonicalEllipse) filterType = do
  let filterTemp =
        VS.fromList . conjugateFunc filterType $
        makeFilterConvolution
          rows
          cols
          (ellipticalPinwheel canonicalEllipse (L.last rfs) (L.last afs) alpha)
      filterList =
        L.map
          (\rf ->
             L.map
               (\af ->
                  VS.fromList . conjugateFunc filterType $
                  makeFilterConvolution
                    rows
                    cols
                    (ellipticalPinwheel canonicalEllipse rf af alpha))
               afs)
          rfs
  lock <- getFFTWLock
  (p1, vec) <- dft1dGPlan lock plan [cols, rows] [0, 1] filterTemp
  (p2, _) <- idft1dGPlan lock p1 [cols, rows] [0, 1] vec
  filters <-
    MP.mapM
      (dftExecuteBatch p2 (DFTPlanID DFT1DG [cols, rows] [0, 1]))
      filterList
  return (p2, filters)

{-# INLINE makeEllipticalPinwheelConvolution' #-}
makeEllipticalPinwheelConvolution' ::
     DFTPlan
  -> EllipticalPinwheelParams
  -> ConvolutionalType
  -> IO [[VS.Vector (Complex Double)]]
makeEllipticalPinwheelConvolution' plan (EllipticalPinwheelParams rows cols rfs afs alpha canonicalEllipse) filterType = do
  let filterTemp =
        VS.fromList . conjugateFunc filterType $
        makeFilterConvolution
          rows
          cols
          (ellipticalPinwheel canonicalEllipse (L.last rfs) (L.last afs) alpha)
      filterList =
        L.map
          (\rf ->
             L.map
               (\af ->
                  VS.fromList . conjugateFunc filterType $
                  makeFilterConvolution
                    rows
                    cols
                    (ellipticalPinwheel canonicalEllipse rf af alpha))
               afs)
          rfs
  filters <-
    M.mapM
      (dftExecuteBatch plan (DFTPlanID DFT1DG [cols, rows] [0, 1]))
      filterList
  return filters
