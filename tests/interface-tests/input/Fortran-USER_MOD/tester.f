! ------------------------------------------------------------------ ! -----------------------------
program xmat_tester                                                  ! Test program for xmat_console
! ------------------------------------------------------------------ !
   implicit none
   
   integer, parameter :: dp = selected_real_kind(15)
   real(dp), dimension(6) :: initial_stress, dot_strain, newstress
   real(dp), dimension(20) :: sigvec
   real(dp), dimension(12) :: strainvec
   real(dp), dimension(6, 6) :: jacobian
   real(dp) :: timeincrement, totaltime, voidratio
   integer, parameter :: nStat = 20
   real(dp), dimension(nStat) :: inpstate, newstate
   integer :: IDTask, iMod, IsUndr, iStep, iTer, Iel, Jnt, ipl
   integer :: NonSym, iStrsDep, iTimeDep, iTang, iPrjDir, iPrjLen, iAbort
   real(dp) :: X, Y, Z, Swp0, Bulk_W, Swp
   real(dp), dimension(50) :: Props

   initial_stress = [-100.0_dp, -250.0_dp, -250.0_dp, 50.0_dp, 10.0_dp, 0.0_dp]
   dot_strain = [-0.00005_dp, 0.000015_dp, 0.000015_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   voidratio = 0.8_dp
   timeincrement = 0.001_dp
   totaltime = 0.0_dp

   iMod = 3 ! Defined in PlaxisInformationPool as 'HYPO-VW96'
   IsUndr = 0
   iStep = 0
   iTer = 0
   Iel = 1
   Jnt = 1
   X = 0.0_dp
   Y = 0.0_dp
   Z = 0.0_dp
   Props = 0.0_dp
   Props(1:16) = [0.5777_dp, 0.0_dp, 4000000.0_dp, 0.27_dp, &
                  0.677_dp, 1.054_dp, 1.212_dp, 0.14_dp, 2.5_dp, &
                  1.1_dp, 2.2_dp, 0.0001_dp, 0.1_dp, 5.5_dp, 0.0_dp, voidratio]
   strainvec = 0.0_dp
   strainvec(1:6) = dot_strain
   sigvec = 0.0_dp
   sigvec(1:6) = initial_stress
   inpstate = 0.0_dp
   inpstate(1) = voidratio

   IDTask = 2
   call USER_MOD(IDTask, iMod, IsUndr, iStep, iTer, Iel, Jnt, X, Y, Z, &
      totaltime, timeincrement, Props, sigvec, Swp0, inpstate, &
      strainvec, jacobian, Bulk_W, newstress, Swp, newstate, &
      ipl, nStat, NonSym, iStrsDep, iTimeDep, iTang, iPrjDir, iPrjLen, iAbort)

   call printvector(vector=newstress, name='newstress')
   call printvector(vector=newstate, name='newstate')

   IDTask = 3
   call USER_MOD(IDTask, iMod, IsUndr, iStep, iTer, Iel, Jnt, X, Y, Z, &
      totaltime, timeincrement, Props, sigvec, Swp0, inpstate, &
      strainvec, jacobian, Bulk_W, newstress, Swp, newstate, &
      ipl, nStat, NonSym, iStrsDep, iTimeDep, iTang, iPrjDir, iPrjLen, iAbort)

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
