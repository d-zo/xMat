! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                x m a t _ c o n s o l e _ C I n t e r               ! xmat_console wrapper function for interoperability with C
! ------------------------------------------------------------------ !
subroutine xmat_console_CInter( &
   ! ============= variables passed in for information ============= !
      c_materialname, c_nparams, c_materialparameters, c_nstatevar, c_statevariables, &
      c_ncomponents, c_oldstress, c_oldstrain, c_timeincrement, c_totaltime, &
   ! ==================== user coding to define ==================== !
      c_newstress, c_newstate, c_jacobian) bind(C, name='xmat_console_CInter')
! ------------------------------------------------------------------ !
   use iso_c_binding, only: c_int, c_double, c_char, c_null_char
   use General_Settings, only: dp, C_To_Fortran_String
   implicit none

   interface
      subroutine xmat_console(materialname, nparams, materialparameters, nstatevar, statevariables, &
         ncomponents, oldstress, oldstrain, timeincrement, totaltime, newstress, newstate, jacobian)
         use General_Settings, only: dp
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

   character(kind=c_char, len=1), dimension(80), intent(in) :: c_materialname
   integer(c_int), intent(in) :: c_nparams, c_nstatevar, c_ncomponents
   real(c_double), dimension(c_nparams), intent(in) :: c_materialparameters
   real(c_double), dimension(c_nstatevar), intent(in) ::c_statevariables
   real(c_double), dimension(c_ncomponents), intent(in) :: c_oldstress, c_oldstrain
   real(c_double), intent(in) :: c_timeincrement, c_totaltime
   real(c_double), dimension(c_ncomponents), intent(out) :: c_newstress
   real(c_double), dimension(c_nstatevar), intent(out) :: c_newstate
   real(c_double), dimension(c_ncomponents*c_ncomponents), intent(out) :: c_jacobian
   ! --------------------------------------------------------------- !
   character(len=80) :: materialname
   real(dp), dimension(c_nparams) :: materialparameters
   real(dp), dimension(c_nstatevar) :: statevariables, newstate
   real(dp), dimension(c_ncomponents) :: oldstress, oldstrain, newstress
   real(dp), dimension(c_ncomponents, c_ncomponents) :: jacobian
   real(dp) :: timeincrement, totaltime
   integer :: nparams, nstatevar, ncomponents, idx, jdx

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall :: xmat_console_CInter            ! Export function name in dll
#endif

   materialname       = C_To_Fortran_String(length=80, c_string=c_materialname)
   nparams            = c_nparams
   materialparameters = real(c_materialparameters, dp)
   nstatevar          = c_nstatevar
   statevariables     = real(c_statevariables, dp)
   ncomponents        = c_ncomponents
   oldstress          = real(c_oldstress, dp)
   oldstrain          = real(c_oldstrain, dp)
   timeincrement      = real(c_timeincrement, dp)
   totaltime          = real(c_totaltime, dp)

   call xmat_console(materialname=materialname, nparams=nparams, materialparameters=materialparameters, &
      nstatevar=nstatevar, statevariables=statevariables, ncomponents=ncomponents, oldstress=oldstress, &
      oldstrain=oldstrain, timeincrement=timeincrement, totaltime=totaltime, newstress=newstress, &
      newstate=newstate, jacobian=jacobian)

   c_newstress = real(newstress, c_double)
   c_newstate = real(newstate, c_double)
   do idx = 1, ncomponents
      do jdx = 1, ncomponents
         c_jacobian((jdx-1)*ncomponents+idx) = real(jacobian(jdx, idx), c_double)
      end do
   end do
end subroutine xmat_console_CInter
