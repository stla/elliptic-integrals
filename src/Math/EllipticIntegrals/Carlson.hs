module Math.EllipticIntegrals.Carlson
  (carlsonRF, carlsonRF',
  carlsonRD, carlsonRD',
  carlsonRJ, carlsonRJ',
  carlsonRC, carlsonRC',
  carlsonRG, carlsonRG')
  where
import Data.Complex                    ( magnitude, Complex )
import Math.EllipticIntegrals.Internal ( Cplx, toCplx, atanC )

rf_ :: Cplx -> Cplx -> Cplx -> Double -> ((Double,Double,Double), Cplx)
rf_ x y z err =
  let a = (x+y+z)/3 in
  let delta = map (\u -> magnitude(1-u/a)) [x,y,z] in
  if maximum delta < err
    then ((delta !! 0, delta !! 1, delta !! 2), a)
    else
      let (sqrtx, sqrty, sqrtz) = (sqrt x, sqrt y, sqrt z) in
      let lambda = sqrtx*sqrty + sqrty*sqrtz + sqrtz*sqrtx in
      rf_ ((x+lambda)/4) ((y+lambda)/4) ((z+lambda)/4) err

-- | Carlson integral RF.
carlsonRF' :: 
     Double -- ^ bound on the relative error
  -> Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double -- ^ third variable 
  -> Complex Double
carlsonRF' err x y z =
  if zeros > 1
    then error "At most one of x, y, z can be 0"
    else
      let ((dx,dy,dz), a) = rf_ x y z err in
      let (e2,e3) = (dx*dy + dy*dz + dz*dx, dx*dy*dz) in
      toCplx(1 - e2/10 + e3/14 + e2*e2/24 - 3*e2*e3/44 - 5*e2*e2*e2/208 +
          3*e3*e3/104 + e2*e2*e3/16) / sqrt a
  where
    zeros = sum (map (\u -> fromEnum (u == 0)) [x,y,z])

-- | Carlson integral RF.
carlsonRF :: 
     Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double -- ^ third variable 
  -> Complex Double
carlsonRF = carlsonRF' 1e-15

rd_ :: Cplx -> Cplx -> Cplx -> Cplx -> Cplx -> Double ->
       ((Double,Double,Double), Cplx, Cplx, Cplx)
rd_ x y z s fac err =
  let a = (x+y+z+z+z)/5 in
  let delta = map (\u -> magnitude(1-u/a)) [x,y,z] in
  if maximum delta < err
    then ((delta !! 0, delta !! 1, delta !! 2), a, s, fac)
    else
      let (sqrtx, sqrty, sqrtz) = (sqrt x, sqrt y, sqrt z) in
      let lambda = sqrtx*sqrty + sqrty*sqrtz + sqrtz*sqrtx in
      let s' = s + fac / (sqrt z * (z + lambda)) in
      rd_ ((x+lambda)/4) ((y+lambda)/4) ((z+lambda)/4) s' (fac/4) err

-- | Carlson integral RD.
carlsonRD' ::
     Double -- ^ bound on the relative error
  -> Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double -- ^ third variable 
  -> Complex Double
carlsonRD' err x y z =
  if zeros > 1
    then error "At most one of x, y, z can be 0"
    else
      let ((dx,dy,dz), a, s, fac) = rd_ x y z 0 1 err in
      let
        (e2,e3,e4,e5) = (dx*dy + dy*dz + 3*dz*dz + 2*dz*dx + dx*dz + 2*dy*dz,
          dz*dz*dz + dx*dz*dz + 3*dx*dy*dz + 2*dy*dz*dz + dy*dz*dz + 2*dx*dz*dz,
          dy*dz*dz*dz + dx*dz*dz*dz + dx*dy*dz*dz + 2*dx*dy*dz*dz,
          dx*dy*dz*dz*dz) in
      3*s + fac * toCplx(1 - 3*e2/14 + e3/6 + 9*e2*e2/88 - 3*e4/22 - 9*e2*e3/52 +
        3*e5/26 - e2*e2*e2/16 + 3*e3*e3/40 + 3*e2*e4/20 + 45*e2*e2*e3/272 -
        9*(e3*e4 + e2*e5)/68) / a / sqrt a
  where
    zeros = sum (map (\u -> fromEnum (u == 0)) [x,y,z])

-- | Carlson integral RD.
carlsonRD ::
     Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double -- ^ third variable 
  -> Complex Double
carlsonRD = carlsonRD' 1e-15

rj_ :: Cplx -> Cplx -> Cplx -> Cplx -> Cplx -> Double -> Cplx -> Int ->
       Double -> [Cplx] -> [Cplx] -> Double -> (Cplx, Int, [Cplx], [Cplx])
rj_ x y z p a maxmagns delta f fac d e err =
  let q = (4/err)**(1/6) * maxmagns / fromIntegral f in
  if magnitude a > q
    then (a, f, d, e)
    else
      let dnew = (sqrt p + sqrt x)*(sqrt p + sqrt y)*(sqrt p + sqrt z)
          d' = (fromIntegral f * dnew) : d
          e' = (toCplx fac * delta / dnew / dnew) : e
          lambda = sqrt x * sqrt y + sqrt y * sqrt z + sqrt z * sqrt x
          x' = (x + lambda) / 4
          y' = (y + lambda) / 4
          z' = (z + lambda) / 4
          p' = (p + lambda) / 4
          a' = (a + lambda) / 4
      in
      rj_ x' y' z' p' a' maxmagns delta (4*f) (fac/64) d' e' err

-- | Carlson integral RJ.
carlsonRJ' ::
     Double -- ^ bound on the relative error
  -> Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double -- ^ third variable 
  -> Complex Double -- ^ fourth variable 
  -> Complex Double
carlsonRJ' err x y z p =
  if zeros > 1
    then error "At most one of x, y, z, p can be 0"
    else
      let a0 = (x + y + z + p + p) / 5
          maxmagns = maximum $ map (\u -> magnitude(a0-u)) [x, y, z, p]
          delta = (p-x)*(p-y)*(p-z)
      in
      let (a, f, d, e) = rj_ x y z p a0 maxmagns delta 1 1 [] [] err
          f' = fromIntegral f
      in
      let x' = (a0 - x) / f' / a
          y' = (a0 - y) / f' / a
          z' = (a0 - z) / f' / a
          p' = -(x'+y'+z') / 2
          e2 = x'*y' + x'*z' + y'*z' - 3*p'*p'
          e3 = x'*y'*z' + 2*e2*p' + 4*p'*p'*p'
          e4 = p'*(2*x'*y'*z' + e2*p' + 3*p'*p'*p')
          e5 = x'*y'*z'*p'*p'
          h = zipWith (\u v -> atanx_over_x(sqrt u) / v) e d
      in
      (1 - 3*e2/14 + e3/6 + 9*e2*e2/88 - 3*e4/22 - 9*e2*e3/52 + 3*e5/26) /
        f' / a / sqrt a + 6 * sum h
  where
    zeros = sum (map (\u -> fromEnum (u == 0)) [x,y,z,p])
    atanx_over_x w = if w == 0 then 1 else atanC w / w

-- | Carlson integral RJ.
carlsonRJ ::
     Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double -- ^ third variable 
  -> Complex Double -- ^ fourth variable 
  -> Complex Double
carlsonRJ = carlsonRJ' 1e-15


rc_ :: Cplx -> Cplx -> Cplx -> Double -> Int -> Double -> (Cplx, Int)
rc_ x y a magn f err =
  let q = (1/3/err)**(1/8) * magn / fromIntegral f in
  if magnitude a > q
    then (a, f)
    else
      let lambda = 2 * sqrt x * sqrt y + y
          a' = (a + lambda) / 4
          x' = (x + lambda) / 4
          y' = (y + lambda) / 4
      in
      rc_ x' y' a' magn (4*f) err

-- | Carlson integral RC.
carlsonRC' ::
     Double -- ^ bound on the relative error
  -> Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double
carlsonRC' err x y =
  if y == 0
    then error "y cannot be 0"
    else
      let a0 = (x + y + y) / 3
          magn = magnitude(a0-x)
      in
      let (a, f) = rc_ x y a0 magn 1 err
          f' = fromIntegral f
      in
      let s = (y - a0) / f' / a in
      (1 + 3*s*s/10 + s*s*s/7 + 3*s*s*s*s/8 + 9*s*s*s*s*s/22 +
        159*s*s*s*s*s*s/208 + 9*s*s*s*s*s*s*s/8) / sqrt a

-- | Carlson integral RC.
carlsonRC ::
     Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double
carlsonRC = carlsonRC' 1e-15


-- | Carlson integral RG.
carlsonRG' ::
     Double -- ^ bound on the relative error passed to `CarlsonRD'`
  -> Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double -- ^ third variable 
  -> Complex Double
carlsonRG' err x y z =
  if zeros > 1
    then sqrt(x+y+z) / 2
    else
      if z == 0
        then carlsonRG' err z x y
        else
          (z * carlsonRF' err x y z -
            (x-z) * (y-z) * carlsonRD' err x y z / 3 +
            sqrt x * sqrt y / sqrt z) / 2
  where
    zeros = sum (map (\u -> fromEnum (u == 0)) [x,y,z])

-- | Carlson integral RG.
carlsonRG ::
     Complex Double -- ^ first variable
  -> Complex Double -- ^ second variable 
  -> Complex Double -- ^ third variable 
  -> Complex Double
carlsonRG = carlsonRG' 1e-15

