program Oedometric_test
   implicit none

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: nstatv = 20
   integer, parameter :: nprops = 4
   character(len=80) :: materialname
   real(dp), dimension(nprops) :: materialparameters
   real(dp), dimension(nstatv) :: statevariables, oldstate, newstate
   real(dp), dimension(6) :: stress, oldstress, newstress, strain, inpstrain
   real(dp), dimension(6, 6) :: jacobian
   real(dp), dimension(4) :: oedo_pressure
   real(dp) :: timeincrement, totaltime, voidratio, sigma1, numbersign
   real(dp) :: starttime, endtime
   integer, parameter :: maxiter = 50000
   integer :: idx, istep, filenr, stat
   logical :: breakall

   ! Decide which constitutive model to use by materialname
   materialname = 'Hypo-Wu92_Karlsruher Sand dicht 01'
   voidratio = 0.65_dp
   materialparameters = [-106.5_dp, -801.5_dp, -797.1_dp, 1077.7_dp]

   statevariables = 0.0_dp
   statevariables(1) = voidratio

   ! Time increment to be used
   totaltime = 0.0001_dp
   timeincrement = 0.0001_dp

   ! State loads
   oedo_pressure = [ &
   -10.0_dp, -100.0_dp, -50.0_dp, -200.0_dp]

   ! Give initial stress and strain matrices (in vector form)
   sigma1 = oedo_pressure(1)
   stress = [sigma1/2.0_dp, sigma1, sigma1/2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   strain = [0.0_dp, -0.05_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]*timeincrement

   oldstress = stress
   oldstate = statevariables


   filenr = 20
   open(filenr, file="Xmat_Oedo_Fortran.csv", iostat=stat)
   if (stat /= 0) then
      write(*, *) "Access problem during write attempt"
      close(filenr)
      stop
   end if

   ! Write output for start point as well
   write(filenr, '(f13.6, a, f13.6)') -oldstress(2), '   ', voidratio

   idx = 1
   breakall = .False.

   call cpu_time(starttime)
   loading_cycle: &
   do istep = 2, size(oedo_pressure)
      numbersign = (-1.0_dp)**istep
      inpstrain = numbersign*strain
      loading_loop: &
      do
         if (numbersign*oldstress(2) < numbersign*oedo_pressure(istep)) then
            exit loading_loop
         end if
         call xmat_console(materialname, size(materialparameters), materialparameters, &
            size(oldstate), oldstate, size(inpstrain), oldstress, inpstrain, timeincrement, totaltime, &
            newstress, newstate, jacobian)

         voidratio = newstate(1)
         oldstress = newstress
         oldstate = newstate
         write(filenr, '(f13.6, a, f13.6)') -oldstress(2), '   ', voidratio

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
