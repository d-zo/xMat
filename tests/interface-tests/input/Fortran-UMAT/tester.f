! ------------------------------------------------------------------ ! -----------------------------
program xmat_tester                                                  ! Test program for xmat_console
! ------------------------------------------------------------------ !
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
   !
   character(len=80) :: materialname
   real(dp), dimension(16) :: materialparameters
   real(dp), dimension(20) :: inoutstate
   real(dp), dimension(6) :: inoutstress, inpstrain
   real(dp) :: timeincrement, totaltime, voidratio
   integer :: ixx, jxx

   ddsddt = 0.0_dp
   drplde = 0.0_dp
   stran = 0.0_dp
   ddsdde = 0.0_dp
   drot = reshape([(1.0_dp, (0.0_dp, ixx = 1, 3), jxx = 1, 2), 1.0_dp], [3, 3])
   dfgrd0 = drot
   dfgrd1 = drot
   coords = 0.0_dp
   !
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
   !
   jstep = 0
   noel = 0
   npt = 0
   layer = 0
   kspt = 0
   kinc = 0

   inoutstress = [-100.0_dp, -250.0_dp, -250.0_dp, 50.0_dp, 10.0_dp, 0.0_dp]
   inpstrain = [-0.00005_dp, 0.000015_dp, 0.000015_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   voidratio = 0.8_dp
   timeincrement = 0.001_dp
   totaltime = 0.0_dp
   time = [totaltime, totaltime]

   materialname = 'HYPO-VW96_Test'
   materialparameters = [0.5777_dp, 0.0_dp, 4000000.0_dp, 0.27_dp, &
                         0.677_dp, 1.054_dp, 1.212_dp, 0.14_dp, 2.5_dp, &
                         1.1_dp, 2.2_dp, 0.0001_dp, 0.1_dp, 5.5_dp, 0.0_dp, voidratio]
   inoutstate = 0.0_dp
   inoutstate(1) = voidratio

   ! Adjust input stress and input strain to UMAT defaults
   inoutstress = [inoutstress(1:4), inoutstress(6), inoutstress(5)]
   inpstrain = [inpstrain(1:3), 2.0_dp*inpstrain(4), 2.0_dp*inpstrain(6), 2.0_dp*inpstrain(5)]

   call UMAT(inoutstress, inoutstate, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
      stran, inpstrain, time, timeincrement, temp, dtemp, predef, dpred, &
      materialname, ndi, nshr, ntens, nstatv, materialparameters, nprops, &
      coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, jstep, kinc)

   ! Adjust stress and jacobian received from UMAT
   inoutstress = [inoutstress(1:4), inoutstress(6), inoutstress(5)]
   ddsdde(:, 5:6) = reshape([ddsdde(:, 6), ddsdde(:, 5)], [6, 2])
   ddsdde(5:6, :) = transpose(reshape([ddsdde(6, :), ddsdde(5, :)], [6, 2]))

   call printvector(vector=inoutstress, name='newstress')
   call printvector(vector=inoutstate, name='newstate')
   call printmatrix(matrix=ddsdde, name='jacobian')


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


   ! -----------------------------------------------------------------------------------------------
   subroutine printmatrix(matrix, name)
   ! -----------------------------------------------------------------------------------------------
      ! Die Ausgabe weicht von der Darstellung im Speicher ab. Fortran speichert Matrizen
      ! spaltenorientiert, die Ausgabe erfolgt zeilenorientiert.
      implicit none

      real(dp), dimension(:, :), intent(in) :: matrix
      character(len=*), intent(in) :: name
      character(len=16) :: formatstring
      integer :: idx

      write(*, '(a,a)') name, ' = '
      formatstring = '(a,i2,a,  f16.7)'
      write(formatstring(9:10), fmt='(i2)') size(matrix, 2)

      do idx = 1, size(matrix, 1)
         write(*, formatstring) '(', idx, ', :) ', matrix(idx, :)
      end do
      write(*, '(a)') ''
   end subroutine printmatrix
end program xmat_tester
