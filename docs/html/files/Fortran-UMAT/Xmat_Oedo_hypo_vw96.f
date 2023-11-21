program Oedometric_test
   implicit none

   integer, parameter :: dp = selected_real_kind(15)

   integer, parameter :: ndi = 3
   integer, parameter :: nshr = 3
   integer, parameter :: ntens = 6
   integer, parameter :: nstatv = 20
   integer, parameter :: nprops = 16
   real(dp), dimension(ntens) :: ddsddt, drplde, stran
   real(dp), dimension(ntens, ntens) :: ddsdde
   real(dp), dimension(3, 3) :: drot, dfgrd0, dfgrd1
   real(dp), dimension(3) :: coords
   real(dp), dimension(2) :: time
   real(dp), dimension(1) :: predef, dpred
   real(dp) :: sse, spd, scd, rpl, drpldt, temp, dtemp, celent, pnewdt
   integer, dimension(4) :: jstep
   integer :: noel, npt, layer, kspt, kinc

   character(len=80) :: materialname
   real(dp), dimension(nprops) :: materialparameters
   real(dp), dimension(nstatv) :: statevariables, inoutstate
   real(dp), dimension(ntens) :: stress, inoutstress
   real(dp), dimension(ntens) :: strain, inpstrain
   real(dp), dimension(6) :: intergranular_strain
   real(dp), dimension(4) :: oedo_pressure
   real(dp) :: dt, sigma1, voidratio, nu, numbersign
   integer :: istep, idx, filenr, stat, ixx, jxx
   integer, parameter :: maxiter = 50000
   logical :: breakall
   real(dp) :: starttime, endtime

   ddsddt = 0.0_dp
   drplde = 0.0_dp
   stran = 0.0_dp
   ddsdde = 0.0_dp
   drot = reshape([(1.0_dp, (0.0_dp, ixx = 1, 3), jxx = 1, 2), 1.0_dp], [3, 3])
   dfgrd0 = drot
   dfgrd1 = drot
   coords = 0.0_dp

   predef = 0.0_dp
   dpred = 0.0_dp
   celent = 0.0_dp
   pnewdt = 0.0_dp
   sse = 0.0_dp
   spd = 0.0_dp
   scd = 0.0_dp
   rpl = 0.0_dp
   drpldt = 0.0_dp
   temp = 0.0_dp
   dtemp = 0.0_dp

   jstep = 0
   noel = 0
   npt = 0
   layer = 0
   kspt = 0
   kinc = 0


   ! Decide which constitutive model to use by materialname
   materialname = 'Hypo-VW96_Hostun Sand 01'
   nu = 0.25_dp
   voidratio = 0.7_dp
   materialparameters = [0.576_dp, nu, 1.0e6_dp, 0.29_dp, 0.63_dp, 1.0_dp, 1.15_dp, 0.13_dp, 2.0_dp, &
      2.0_dp, 5.0_dp, 1.0e-4_dp, 0.5_dp, 6.0_dp, 0.0_dp, 0.0_dp]
   intergranular_strain = [0.0_dp, -0.0001_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

   statevariables = 0.0_dp
   statevariables(1) = voidratio
   statevariables(2:4) = intergranular_strain(1:3)
   statevariables(5:6) = intergranular_strain(4)
   statevariables(7:8) = intergranular_strain(5)
   statevariables(9:10) = intergranular_strain(6)

   ! Time increment to be used
   time = [0.0001_dp, 0.0001_dp]
   dt = 0.0001_dp

   ! State loads
   oedo_pressure = [ &
   -10.0_dp, -100.0_dp, -50.0_dp, -200.0_dp]

   ! Give initial stress and strain matrices (in vector form)
   sigma1 = oedo_pressure(1)
   stress = [sigma1/2.0_dp, sigma1, sigma1/2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   strain = [0.0_dp, -0.05_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]*dt

   ! Adjust input stress and input strain to UMAT defaults
   stress = [stress(1:4), stress(6), stress(5)]
   strain = [strain(1:3), 2.0_dp*strain(4), 2.0_dp*strain(6), 2.0_dp*strain(5)]

   inoutstress = stress
   inoutstate = statevariables


   filenr = 20
   open(filenr, file="Xmat_Oedo_Fortran-UMAT.csv", iostat=stat)
   if (stat /= 0) then
      write(*, *) "Access problem during write attempt"
      close(filenr)
      stop
   end if

   ! Write output for start point as well
   write(filenr, '(f13.6, a, f13.6)') -inoutstress(2), '   ', voidratio

   idx = 1
   breakall = .False.

   call cpu_time(starttime)
   loading_cycle: &
   do istep = 2, size(oedo_pressure)
      numbersign = (-1.0_dp)**istep
      inpstrain = numbersign*strain
      loading_loop: &
      do
         if (numbersign*inoutstress(2) < numbersign*oedo_pressure(istep)) then
            exit loading_loop
         end if
         call UMAT(inoutstress, inoutstate, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
                   stran, inpstrain, time, dt, temp, dtemp, predef, dpred, materialname, ndi, nshr, &
                   ntens, nstatv, materialparameters, nprops, coords, drot, pnewdt, celent, dfgrd0, &
                   dfgrd1, noel, npt, layer, kspt, jstep, kinc)

         voidratio = inoutstate(1)
         write(filenr, '(f13.6, a, f13.6)') -inoutstress(2), '   ', voidratio

         idx = idx + 1
         if (idx > maxiter) then
            write(*, *) 'Maximum number of specified iterations reached'
            breakall = .True.
            exit loading_loop
         end if
      end do loading_loop
      if (breakall) then
         exit loading_cycle
      end if
   end do loading_cycle
   call cpu_time(endtime)
   print '("Xmat_Oedo_Fortran: ",f6.3,"s")', endtime - starttime
   close(filenr)
end program Oedometric_test
