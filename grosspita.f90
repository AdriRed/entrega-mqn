PROGRAM Main
   IMPLICIT NONE
   INTEGER :: n1, iter, it, i, itw
   REAL(8) :: a0, step, num_atoms, time, alpha, cequ, as3n
   REAL(8) :: pi, piin, pi2in, alpha2, cvar
   REAL(8) :: xr2, xnorm, ene0, radious, xkin, poth0, potself, chem, xaver, xnormden, radious2, pot
   REAL(8), DIMENSION(1000) :: xr, frev, freo, fred, xmu, fren, den, u

   ! Compute constants
   pi = 4.0d0 * DATAN(1.0d0)
   piin = 1.0d0 / (4.0d0 * pi)
   pi2in = DSQRT(piin)

   ! Read input values
   READ(*,*) a0, n1, step, num_atoms, time, alpha, iter

   alpha2 = alpha * alpha
   cvar = 2.0d0 * DSQRT(alpha)**3 / DSQRT(DSQRT(pi))

   ! Initialize wave function
   DO i = 1, n1
      xr(i) = step * DBLE(i - 1)
      xr2 = xr(i) * xr(i)
      frev(i) = cvar * xr(i) * DEXP(-0.5d0 * alpha2 * xr2)
      freo(i) = frev(i)
   END DO

   ! Define potential parameters
   cequ = a0 * num_atoms
   as3n = num_atoms * a0 * a0 * a0
   itw = 0

   ! Iterative process
   DO it = 1, iter
      itw = itw + 1
      xnorm = 0.0d0
      ene0 = 0.0d0

      DO i = 2, n1 - 1
         fred(i) = (freo(i - 1) + freo(i + 1) - 2.0d0 * freo(i)) / (step * step)
      END DO
      fred(n1) = (freo(n1 - 1) - 2.0d0 * freo(n1)) / (step * step)

      DO i = 1, n1
         xr2 = xr(i) * xr(i)
         IF (i == 1) THEN
            xmu(i) = 0.0d0
         ELSE
            ene0 = ene0 - freo(i) * fred(i) * 0.5d0 + 0.5d0 * xr2 * freo(i) * freo(i) + 0.5d0 * cequ * xr2 * (freo(i) / xr(i))**4
            xmu(i) = -0.5d0 * fred(i) / freo(i) + 0.5d0 * xr2 + cequ * (freo(i) / xr(i))**2
         END IF
         fren(i) = freo(i) - time * xmu(i) * freo(i)
         xnorm = xnorm + fren(i) * fren(i)
      END DO

      xnorm = DSQRT(xnorm * step)
      ene0 = ene0 * step

      IF (itw == 200) THEN
         WRITE(*,*) 'ene0 =', ene0
         itw = 0
      END IF

      DO i = 1, n1
         freo(i) = fren(i) / xnorm
      END DO

      IF (it == iter) THEN
         WRITE(10, '(2E15.5)') (xr(i), xmu(i), i = 2, n1)
      END IF
   END DO

   ! Calculation of radius, potential and kinetic energy, density, and single particle potential
   DO i = 2, n1 - 1
      fred(i) = (freo(i - 1) + freo(i + 1) - 2.0d0 * freo(i)) / (step * step)
   END DO
   fred(n1) = (freo(n1 - 1) - 2.0d0 * freo(n1)) / (step * step)

   radious = 0.0d0
   xkin = 0.0d0
   poth0 = 0.0d0
   potself = 0.0d0
   chem = 0.0d0
   xaver = 0.0d0
   xnormden = 0.0d0

   DO i = 2, n1
      xr2 = xr(i) * xr(i)
      radious = radious + xr2 * freo(i) * freo(i)
      xkin = xkin + freo(i) * fred(i)
      poth0 = poth0 + xr2 * freo(i) * freo(i)
      potself = potself + xr2 * (freo(i) / xr(i))**4
      chem = chem + xmu(i) * freo(i) * freo(i)
      u(i) = 0.5d0 * xr2 + cequ * (freo(i) / xr(i))**2
      den(i) = (freo(i) / xr(i))**2
      xnormden = xnormden + den(i) * xr2
      xaver = xaver + freo(i) * freo(i) * as3n * den(i)
   END DO

   radious2 = radious * step
   radious = DSQRT(radious * step)
   xaver = xaver * step
   chem = chem * step
   xkin = -xkin * step * 0.5D0
   poth0 = 0.5d0 * poth0 * step
   potself = potself * step * cequ * 0.5d0
   pot = potself + poth0
   xnormden = xnormden * step

   WRITE(*,*) 'xnormden =', xnormden
   WRITE(*, '(A, E12.5, A, E12.5)') ' * ener', ene0, '   average chemical=', chem
   WRITE(*, '(A, E12.5, A, E12.5)') ' * kin-ener=', xkin, '     total-pot= ', pot
   WRITE(*, '(A, E12.4, A, E12.4)') ' * poth0=', poth0, '     potint =', potself
   WRITE(*, '(A, E12.4, A, E15.7)') ' * radious =', radious, ' radious2 =', radious2

   DO i = 2, n1
      WRITE(9, '(2E15.5)') xr(i), den(i)
   END DO

   STOP
END PROGRAM Main

