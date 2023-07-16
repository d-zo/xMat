! ------------------------------------------------------------------ ! -----------------------------
program xmat_tester                                                  ! Test program for xmat_console
! ------------------------------------------------------------------ !
   implicit none
   
   integer, parameter :: dp = selected_real_kind(15)
   character(len=80) :: ela_name, hyw_name, wo2_name, bak_name, bad_name, ba2_name, ne2_name
   real(dp), dimension(2) :: ela_params
   real(dp), dimension(4) :: hyw_params
   real(dp), dimension(16) :: wo2_params
   real(dp), dimension(7) :: bak_params
   real(dp), dimension(9) :: bad_params
   real(dp), dimension(10) :: ba2_params
   real(dp), dimension(15) :: ne2_params
   real(dp), dimension(1) :: ela_oldstate, ela_newstate
   real(dp), dimension(1) :: hyw_oldstate, hyw_newstate
   real(dp), dimension(11) :: wo2_oldstate, wo2_newstate
   real(dp), dimension(1) :: bak_oldstate, bak_newstate
   real(dp), dimension(1) :: bad_oldstate, bad_newstate
   real(dp), dimension(12) :: ba2_oldstate, ba2_newstate
   real(dp), dimension(11) :: ne2_oldstate, ne2_newstate
   real(dp), dimension(6) :: initial_stress, stress, dot_strain, newstress
   real(dp), dimension(6, 6) :: jacobian
   real(dp) :: timeincrement, totaltime, voidratio
   real(dp) :: starttime, endtime

   interface
      subroutine xmat_console(materialname, nparams, materialparameters, nstatevar, statevariables, &
         ncomponents, oldstress, oldstrain, timeincrement, totaltime, newstress, newstate, jacobian)
         integer, parameter :: dp = selected_real_kind(15)
         character(len=80), intent(in) :: materialname
         integer, intent(in) :: nparams, nstatevar
         real(dp), dimension(nparams), intent(in) :: materialparameters
         real(dp), dimension(nstatevar), intent(in) :: statevariables
         integer, intent(in) :: ncomponents
         real(dp), dimension(ncomponents), intent(in) :: oldstress, oldstrain
         real(dp), intent(in) :: timeincrement, totaltime
         real(dp), dimension(ncomponents), intent(out) :: newstress
         real(dp), dimension(nstatevar), intent(out) :: newstate
         real(dp), dimension(ncomponents, ncomponents), intent(out) :: jacobian
      end subroutine xmat_console
   end interface

   initial_stress = [-250.0_dp, -100.0_dp, -100.0_dp, 10.0_dp, 5.0_dp, 0.0_dp]
   dot_strain = [-0.001_dp, 0.0003_dp, 0.0003_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   voidratio = 0.699_dp
   timeincrement = 0.01_dp
   totaltime = 0.0_dp

   call printvector(name='(INI) stress', vector=initial_stress)
   call printvector(name='(INI)dstrain', vector=dot_strain)
   write(*, *) '---'

   ! ---

   ela_name = 'ELAS-0000_Test'
   ela_params = [20000.0_dp, 0.25_dp]
   call cpu_time(starttime)
   call xmat_console(ela_name, size(ela_params), ela_params, size(ela_oldstate), ela_oldstate, &
                     size(dot_strain), initial_stress, dot_strain, timeincrement, totaltime, &
                     newstress, ela_newstate, jacobian)
   call cpu_time(endtime)
   call printvector(vector=newstress, name='(ELAS-0000) stress')
   call printvector(vector=ela_newstate, name='(ELAS-0000) state ')
   write(*, '(a, f16.7, a)') '(ELAS-0000): ', endtime - starttime, 's'
   write(*, *) '---'
   
   ! ---

   hyw_name = 'HYPO-WU92_Test'
   hyw_params = [-106.5_dp, -801.5_dp, -797.1_dp, 1077.7_dp]
   hyw_oldstate = 0.0_dp
   hyw_oldstate(1) = voidratio
   call cpu_time(starttime)
   call xmat_console(hyw_name, size(hyw_params), hyw_params, size(hyw_oldstate), hyw_oldstate, &
                     size(dot_strain), initial_stress, dot_strain, timeincrement, totaltime, &
                     newstress, hyw_newstate, jacobian)
   call cpu_time(endtime)
   call printvector(vector=newstress, name='(HYPO-WU92) stress')
   call printvector(vector=hyw_newstate, name='(HYPO-WU92) state ')
   write(*, '(a, f16.7, a)') '(HYPO-WU92): ', endtime - starttime, 's'
   write(*, *) '---'
   
   ! ---

   wo2_name = 'HYPO-VW96_Test'
   wo2_params = [0.524_dp, 0.0_dp, 2600000.0_dp, 0.27_dp, &
                 0.61_dp, 0.98_dp, 1.10_dp, 0.18_dp, 1.0_dp, &
                 2.0_dp, 5.0_dp, 0.0001_dp, 0.5_dp, 6.0_dp, 0.0_dp, voidratio]
   wo2_oldstate = 0.0_dp
   wo2_oldstate(1) = voidratio
   call cpu_time(starttime)
   call xmat_console(wo2_name, size(wo2_params), wo2_params, size(wo2_oldstate), wo2_oldstate, &
                     size(dot_strain), initial_stress, dot_strain, timeincrement, totaltime, &
                     newstress, wo2_newstate, jacobian)
   call cpu_time(endtime)
   call printvector(vector=newstress, name='(HYPO-VW96) stress')
   call printvector(vector=wo2_newstate, name='(HYPO-VW96) state ')
   write(*, '(a, f16.7, a)') '(HYPO-VW96): ', endtime - starttime, 's'
   write(*, *) '---'

   ! ---

   bak_name = 'BARO-KO15_Test'
   bak_params = [0.590_dp, 1.0_dp, 2.5076_dp, 1000.0_dp, &
                 40.0_dp, 0.91_dp, 0.5_dp]
   bak_oldstate = 0.0_dp
   bak_oldstate(1) = voidratio
   call cpu_time(starttime)
   call xmat_console(bak_name, size(bak_params), bak_params, size(bak_oldstate), bak_oldstate, &
                     size(dot_strain), initial_stress, dot_strain, timeincrement, totaltime, &
                     newstress, bak_newstate, jacobian)
   call cpu_time(endtime)
   call printvector(vector=newstress, name='(BARO-KO15) stress')
   call printvector(vector=bak_newstate, name='(BARO-KO15) state ')
   write(*, '(a, f16.7, a)') '(BARO-KO15): ', endtime - starttime, 's'
   write(*, *) '---'

   ! ---

   bad_name = 'BARO-SC18_Test'
   bad_params = [0.590_dp, 430.0_dp, 0.74_dp, 3.0_dp, &
                 0.904_dp, 1.093_dp, 3.0_dp, 5.0_dp, 10.0_dp]
   bad_oldstate = 0.0_dp
   bad_oldstate(1) = voidratio
   call cpu_time(starttime)
   call xmat_console(bad_name, size(bad_params), bad_params, size(bad_oldstate), bad_oldstate, &
                     size(dot_strain), initial_stress, dot_strain, timeincrement, totaltime, &
                     newstress, bad_newstate, jacobian)
   call cpu_time(endtime)
   call printvector(vector=newstress, name='(BARO-SC18) stress')
   call printvector(vector=bad_newstate, name='(BARO-SC18) state ')
   write(*, '(a, f16.7, a)') '(BARO-SC18): ', endtime - starttime, 's'
   write(*, *) '---'

   ! ---

   ba2_name = 'BARO-KO21_Test'
   ba2_params = [0.590_dp, 0.50_dp, -0.20_dp, 5000.0_dp, 0.50_dp, &
      30.0_dp, 0.0_dp, 1.0_dp, 0.02_dp, 0.9_dp]
   ba2_oldstate = 0.0_dp
   ba2_oldstate(1) = voidratio
   call cpu_time(starttime)
   call xmat_console(ba2_name, size(ba2_params), ba2_params, size(ba2_oldstate), ba2_oldstate, &
                     size(dot_strain), initial_stress, dot_strain, timeincrement, totaltime, &
                     newstress, ba2_newstate, jacobian)
   call cpu_time(endtime)
   call printvector(vector=newstress, name='(BARO-KO21) stress')
   call printvector(vector=ba2_newstate, name='(BARO-KO21) state ')
   write(*, '(a, f16.7, a)') '(BARO-KO21): ', endtime - starttime, 's'
   write(*, *) '---'

   ! ---

   ! FIXME: The valid stress bounds seem to be different than the other ones.
   !        Clarify and adjust so that one initial_stress can be used for all cases
   stress = [-100.0_dp, -100.0_dp, -100.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   
   ne2_name = 'VIHY-NI03_Test'
   ne2_params = [0.6545_dp, 0.3_dp, 0.0656_dp, 0.0062_dp, 0.5_dp, 0.016_dp, 0.0000008329_dp, 0.436_dp, &
                 0.0_dp, 2.0_dp, 5.0_dp, 0.0001_dp, 0.5_dp, 6.0_dp, 1.0_dp]
   ne2_oldstate = 0.0_dp
   ne2_oldstate(1) = voidratio
   call cpu_time(starttime)
   call xmat_console(ne2_name, size(ne2_params), ne2_params, size(ne2_oldstate), ne2_oldstate, &
                     size(dot_strain), stress, dot_strain, timeincrement, totaltime, &
                     newstress, ne2_newstate, jacobian)
   call cpu_time(endtime)
   call printvector(vector=newstress, name='(VIHY-NI03) stress')
   call printvector(vector=ne2_newstate, name='(VIHY-NI03) state ')
   write(*, '(a, f16.7, a)') '(VIHY-NI03): ', endtime - starttime, 's'
   write(*, *) '---'


   contains


   ! -----------------------------------------------------------------------------------------------
   subroutine printvector(vector, name)
   ! -----------------------------------------------------------------------------------------------
      implicit none

      real(dp), dimension(:), intent(in) :: vector
      character(len=*), intent(in) :: name
      character(len=16) :: formatstring

      formatstring = '(a,a,  f16.7)'
      write(formatstring(6:7), fmt='(i2)') size(vector)
      write(*, formatstring) name, ' = ', vector
   end subroutine printvector
end program xmat_tester
