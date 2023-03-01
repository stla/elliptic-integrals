module Math.EllipticIntegrals.Elliptic
  where
import Math.EllipticIntegrals.Carlson  ( carlsonRF', carlsonRD', carlsonRJ' )
import Data.Complex                    ( realPart, Complex )
import Math.EllipticIntegrals.Internal ( toCplx, getPhiK )

-- | Elliptic integral of the first kind.
ellipticF' :: 
     Double -- ^ bound on the relative error passed to `carlsonRF'`
  -> Complex Double -- ^ amplitude
  -> Complex Double -- ^ parameter
  -> Complex Double
ellipticF' err phi m
  | phi == 0 =
    toCplx 0
  | m == 1 && abs(realPart phi) == pi/2 =
    toCplx (0/0)
  | m == 1 && abs(realPart phi) < pi/2 =
    atanh(sin phi)
  | abs(realPart phi) <= pi/2 =
    if m == 0
      then
        phi
      else
        let sine = sin phi in
        let sine2 = sine*sine in
        let (cosine2, oneminusmsine2) = (1 - sine2, 1 - m*sine2) in
        sine * carlsonRF' err cosine2 oneminusmsine2 1
  | otherwise =
    let (phi', k) = getPhiK phi in
    2 * fromIntegral k * ellipticF' err (pi/2) m + ellipticF' err phi' m

-- | Elliptic integral of the first kind.
ellipticF :: 
     Complex Double -- ^ amplitude
  -> Complex Double -- ^ parameter
  -> Complex Double
ellipticF = ellipticF' 1e-15

-- | Elliptic integral of the second kind.
ellipticE' :: 
     Double -- ^ bound on the relative error passed to `carlsonRF'` and `carlsonRD'` 
  -> Complex Double -- ^ amplitude
  -> Complex Double -- ^ parameter
  -> Complex Double
ellipticE' err phi m
  | phi == 0 =
    toCplx 0
  | abs(realPart phi) <= pi/2 =
    case m of
      0 -> phi
      1 -> sin phi
      _ ->
        let sine = sin phi in
        let sine2 = sine*sine in
        let (cosine2, oneminusmsine2) = (1 - sine2, 1 - m*sine2) in
        sine * (carlsonRF' err cosine2 oneminusmsine2 1 -
          m * sine2 / 3 * carlsonRD' err cosine2 oneminusmsine2 1)
  | otherwise =
    let (phi', k) = getPhiK phi in
    2 * fromIntegral k * ellipticE' err (pi/2) m + ellipticE' err phi' m

-- | Elliptic integral of the second kind.
ellipticE :: 
     Complex Double -- ^ amplitude
  -> Complex Double -- ^ parameter
  -> Complex Double
ellipticE = ellipticE' 1e-15

-- | Elliptic integral of the third kind.
ellipticPI' :: 
     Double -- ^ bound on the relative error passed to `carlsonRF'` and `carlsonRJ'` 
  -> Complex Double -- ^ amplitude
  -> Complex Double -- ^ characteristic
  -> Complex Double -- ^ parameter
  -> Complex Double
ellipticPI' err phi n m
  | phi == 0 =
    toCplx 0
  | phi == pi/2 && n == 1 =
    0/0
  | phi == pi/2 && m == 0 =
    pi/2/sqrt(1-n)
  | phi == pi/2 && m == n =
    ellipticE' err (pi/2) m / (1-m)
  | phi == pi/2 && n == 0 =
    ellipticF' err (pi/2) m
  | abs(realPart phi) <= pi/2 =
    let sine = sin phi in
    let sine2 = sine*sine in
    let (cosine2, oneminusmsine2) = (1 - sine2, 1 - m*sine2) in
    sine * (carlsonRF' err cosine2 oneminusmsine2 1 +
      n * sine2 / 3 * carlsonRJ' err cosine2 oneminusmsine2 1 (1-n*sine2))
  | otherwise =
    let (phi', k) = getPhiK phi in
    2 * fromIntegral k * ellipticPI' err (pi/2) n m + ellipticPI' err phi' n m

-- | Elliptic integral of the third kind.
ellipticPI ::
     Complex Double -- ^ amplitude
  -> Complex Double -- ^ characteristic
  -> Complex Double -- ^ parameter
  -> Complex Double
ellipticPI = ellipticPI' 1e-15

-- | Jacobi zeta function.
jacobiZeta' ::
     Double -- ^ bound on the relative error passed to `ellipticF'` and `ellipticE'` 
  -> Complex Double -- ^ amplitude
  -> Complex Double -- ^ parameter
  -> Complex Double
jacobiZeta' err phi m =
  if m == 1
    then
      if abs(realPart phi) <= pi/2
        then sin phi
        else let (phi',_) = getPhiK phi in sin phi'
    else
      ellipticE' err phi m -
        ellipticE' err (pi/2) m / ellipticF' err (pi/2) m *
        ellipticF' err phi m

-- | Jacobi zeta function.
jacobiZeta ::
     Complex Double -- ^ amplitude
  -> Complex Double -- ^ parameter
  -> Complex Double
jacobiZeta = jacobiZeta' 1e-15

