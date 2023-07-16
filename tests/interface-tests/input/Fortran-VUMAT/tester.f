! ------------------------------------------------------------------ ! -----------------------------
program xmat_tester                                                  ! Test program for xmat_console
! ------------------------------------------------------------------ !
   implicit none
   
   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: nblock = 1
   integer, parameter :: ndir = 3
   integer, parameter :: nshr = 3
   integer, parameter :: nstatev = 20
   integer, parameter :: nfieldv = 1
   integer, parameter :: nprops = 16
   integer :: lanneal
   real(dp) :: stepTime, totalTime
   real(dp), dimension(1) :: charLength, density, tempOld, enerInternOld, enerInelasOld, tempNew, &
                             enerInternNew, enerInelasNew
   real(dp), dimension(1, 6) :: stretchOld, stretchNew
   real(dp), dimension(1, 3) :: relSpinInc
   real(dp), dimension(1, 9) :: defgradOld, defgradNew
   real(dp), dimension(1, 1) :: fieldOld, fieldNew
   real(dp), dimension(1, 1) :: coordMp
   !
   character(len=80) :: materialname
   real(dp), dimension(nprops) :: materialparameters
   real(dp), dimension(nstatev) :: statevariables
   real(dp), dimension(ndir + nshr) :: stress, strain
   real(dp), dimension(1, nstatev) :: inpstate, outstate
   real(dp), dimension(1, ndir + nshr) :: inpstress, outstress
   real(dp), dimension(1, ndir + nshr) :: inpstrain
   real(dp) :: timeincrement, voidratio
   integer :: ixx, jxx

   lanneal = 0
   stepTime = 0.0_dp
   totalTime = 0.0_dp

   coordMp = 0.0_dp
   charLength = 0.0_dp
   density = 0.0_dp
   relSpinInc = 0.0_dp
   tempOld = 0.0_dp
   stretchOld = 0.0_dp
   defgradOld = 0.0_dp
   fieldOld = 0.0_dp
   enerInternOld = 0.0_dp
   enerInelasOld = 0.0_dp
   tempNew = 0.0_dp
   stretchNew = 0.0_dp
   defgradNew = 0.0_dp
   fieldNew = 0.0_dp
   enerInternNew = 0.0_dp
   enerInelasNew = 0.0_dp

   inpstress(1, :) = [-100.0_dp, -250.0_dp, -250.0_dp, 50.0_dp, 10.0_dp, 0.0_dp]
   inpstrain(1, :) = [-0.00005_dp, 0.000015_dp, 0.000015_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   voidratio = 0.8_dp
   timeincrement = 0.001_dp

   materialname = 'HYPO-VW96_Test'
   materialparameters = [0.5777_dp, 0.0_dp, 4000000.0_dp, 0.27_dp, &
                         0.677_dp, 1.054_dp, 1.212_dp, 0.14_dp, 2.5_dp, &
                         1.1_dp, 2.2_dp, 0.0001_dp, 0.1_dp, 5.5_dp, 0.0_dp, voidratio]
   inpstate = 0.0_dp
   inpstate(1, 1) = voidratio

   call VUMAT(nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, stepTime, totalTime, &
      timeincrement, materialname, coordMp, charLength, materialparameters, density, inpstrain, &
      relSpinInc, tempOld, stretchOld, defgradOld, fieldOld, inpstress, inpstate, &
      enerInternOld, enerInelasOld, tempNew, stretchNew, defgradNew, fieldNew, &
      outstress, outstate, enerInternNew, enerInelasNew)

   call printvector(vector=outstress(1, :), name='newstress')
   call printvector(vector=outstate(1, :), name='newstate')


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
