! ------------------------------------------------------------------ ! -----------------------------
program xmat_tester                                                  ! Test program for xmat_console
! ------------------------------------------------------------------ !
   implicit none
   
   integer, parameter :: dp = selected_real_kind(15)
   character(len=80) :: materialname
   real(dp), dimension(16) :: materialparameters
   real(dp), dimension(20) :: oldstate, newstate, dotstate
   real(dp), dimension(6) :: initial_stress, stress, strain, dotstress, newstress
   real(dp), dimension(6, 6) :: jacobian
   real(dp) :: timeincrement, totaltime, voidratio

   initial_stress = [-100.0_dp, -250.0_dp, -250.0_dp, 50.0_dp, 10.0_dp, 0.0_dp]
   strain = [-0.00005_dp, 0.000015_dp, 0.000015_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   voidratio = 0.8_dp
   timeincrement = 0.001_dp
   totaltime = 0.0_dp

   materialname = 'HYPO-VW96_Test'
   materialparameters = [0.5777_dp, 0.0_dp, 4000000.0_dp, 0.27_dp, &
                         0.677_dp, 1.054_dp, 1.212_dp, 0.14_dp, 2.5_dp, &
                         1.1_dp, 2.2_dp, 0.0001_dp, 0.1_dp, 5.5_dp, 0.0_dp, voidratio]
   oldstate = 0.0_dp
   oldstate(1) = voidratio

   call Dot_Values(materialname, size(materialparameters), materialparameters, &
      size(oldstate), oldstate, size(strain), initial_stress, strain, timeincrement, totaltime, &
      dotstress, dotstate, jacobian)

   call printvector(vector=dotstress, name='dotstress')
   call printvector(vector=dotstate, name='dotstate')
   call printmatrix(matrix=jacobian, name='jacobian')

   call xmat_console(materialname, size(materialparameters), materialparameters, &
      size(oldstate), oldstate, size(strain), initial_stress, strain, timeincrement, totaltime, &
      newstress, newstate, jacobian)

   call printvector(vector=newstress, name='newstress')
   call printvector(vector=newstate, name='newstate')
   call printmatrix(matrix=jacobian, name='jacobian')


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
