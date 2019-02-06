module Descriptor.Pinwheel
  ( ConvolutionalType(..)
  , PinwheelParams(..)
  , makePinwheelExpansion
  , makePinwheelConvolutionPlan
  , makePinwheelConvolution
  , makePinwheelConvolution'
  , applyPinwheelConvolution
  , applyPinwheelConvolutionP
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
import           DFT.Plan
import           Types
import           Utils

data ConvolutionalType
  = Convolution
  | Crosscorrelation
  deriving (Show, Read)

data PinwheelParams = PinwheelParams
  { pinwheelParamsRows        :: Int
  , pinwheelParamsCols        :: Int
  , pinwheelParamsRadialFreq  :: [Double]
  , pinwheelParamsAngularFreq :: [Double]
  , pinwheelParamsAlpha       :: Double
  } deriving (Read, Show, Eq, Ord)

{-# INLINE pinwheel #-}
pinwheel :: Double -> Double -> Double -> Int -> Int -> Complex Double
pinwheel rf af alpha x y
  -- | r <= (2 / pi * (abs af)) = 0
  | r <= (0.75 / pi * (abs af)) || r <= (0.75 / pi * (abs rf)) = 0
  | otherwise = (((r) :+ 0) ** (alpha :+ (-rf))) * exp (0 :+ ((-af) * theta))
  where
    r = sqrt . fromIntegral $ x ^ (2 :: Int) + y ^ (2 :: Int)
    theta = angleFunctionRad (fromIntegral x) (fromIntegral y)

-- For display purpose
{-# INLINE makeFilter #-}
makeFilter :: Int -> Int -> Int -> Int -> (Int -> Int -> a) -> [a]
makeFilter rows cols rCenter cCenter f =
  [f (r - rCenter) (c - cCenter) | c <- [0 .. cols - 1], r <- [0 .. rows - 1]]

{-# INLINE makePinwheelExpansion #-}
makePinwheelExpansion ::
     PinwheelParams -> Int -> Int -> [[VU.Vector (Complex Double)]]
makePinwheelExpansion (PinwheelParams rows cols rfs afs alpha) rCenter cCenter =
  [ [ VU.fromList $ makeFilter rows cols rCenter cCenter (pinwheel rf af alpha)
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

{-# INLINE makePinwheelConvolutionPlan #-}
makePinwheelConvolutionPlan ::
     ConvolutionalType -> DFTPlan -> PinwheelParams -> IO DFTPlan
makePinwheelConvolutionPlan filterType plan (PinwheelParams rows cols rfs afs alpha) = do
  let filterTemp =
        VS.fromList . conjugateFunc filterType $
        makeFilterConvolution
          rows
          cols
          (pinwheel (L.last rfs) (L.last afs) alpha)
  lock <- getFFTWLock
  (p1, vec) <- dft1dGPlan lock plan [cols, rows] [0, 1] filterTemp
  (p2, _) <- idft1dGPlan lock p1 [cols, rows] [0, 1] vec
  return p2

-- make both fftw plan and pinwheel filters
{-# INLINE makePinwheelConvolution #-}
makePinwheelConvolution ::
     DFTPlan
  -> PinwheelParams
  -> ConvolutionalType
  -> IO (DFTPlan, [[VS.Vector (Complex Double)]])
makePinwheelConvolution plan (PinwheelParams rows cols rfs afs alpha) filterType = do
  let filterTemp =
        VS.fromList . conjugateFunc filterType $
        makeFilterConvolution
          rows
          cols
          (pinwheel (L.last rfs) (L.last afs) alpha)
      filterList =
        L.map
          (\rf ->
             L.map
               (\af ->
                  VS.fromList . conjugateFunc filterType $
                  makeFilterConvolution rows cols (pinwheel rf af alpha))
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

{-# INLINE makePinwheelConvolution' #-}
makePinwheelConvolution' ::
     DFTPlan
  -> PinwheelParams
  -> ConvolutionalType
  -> IO [[VS.Vector (Complex Double)]]
makePinwheelConvolution' plan (PinwheelParams rows cols rfs afs alpha) filterType = do
  let filterTemp =
        VS.fromList . conjugateFunc filterType $
        makeFilterConvolution
          rows
          cols
          (pinwheel (L.last rfs) (L.last afs) alpha)
      filterList =
        L.map
          (\rf ->
             L.map
               (\af ->
                  VS.fromList . conjugateFunc filterType $
                  makeFilterConvolution rows cols (pinwheel rf af alpha))
               afs)
          rfs
  filters <-
    M.mapM
      (dftExecuteBatch plan (DFTPlanID DFT1DG [cols, rows] [0, 1]))
      filterList
  return filters

{-# INLINE applyPinwheelConvolution #-}
applyPinwheelConvolution ::
     DFTPlan
  -> Int
  -> Int
  -> [[VS.Vector (Complex Double)]]
  -> [VS.Vector (Complex Double)]
  -> IO [[VS.Vector Double]]
applyPinwheelConvolution plan rows cols filters xs = do
  ys <- dftExecuteBatch plan (DFTPlanID DFT1DG [cols,rows] [0,1]) xs
  fmap (\x -> [L.map (VS.map magnitude) x]) .
    dftExecuteBatch plan (DFTPlanID IDFT1DG [cols,rows] [0,1]) .
    L.concatMap (\x -> L.concatMap (L.map (VS.zipWith (*) x)) filters) $
    ys

{-# INLINE applyPinwheelConvolutionP #-}
applyPinwheelConvolutionP ::
     DFTPlan
  -> Int
  -> Int
  -> [[VS.Vector (Complex Double)]]
  -> [VS.Vector (Complex Double)]
  -> IO [[VS.Vector Double]]
applyPinwheelConvolutionP plan rows cols filters xs = do
  ys <- dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0,1]) xs
  fmap (\x -> [parMap rdeepseq (VS.map magnitude) x]) .
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0,1]) .
    L.concatMap (\x -> L.concatMap (L.map (VS.zipWith (*) x)) filters) $
    ys
