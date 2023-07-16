   ! --------------------------------------------------------------- !
   subroutine Abort_If_Not_In_Interval(name, number, limits)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Write_Error_And_Exit
      use Debug, only: Formatval
      !
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: number
      real(dp), dimension(2), intent(in) :: limits
      ! ------------------------------------------------------------ !
      if (.not. Value_In_Interval(number, limits)) then
         call Write_Error_And_Exit('Abort_If_Not_In_Interval: ' // &
            Formatval(name, number) // ' not in valid interval')
      end if
   end subroutine Abort_If_Not_In_Interval
