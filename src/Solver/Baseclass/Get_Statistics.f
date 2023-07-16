   ! --------------------------------------------------------------- !
   function Get_Statistics(this) result(statistics)
   ! --------------------------------------------------------------- !
      class(Solver), intent(in) :: this
      integer, dimension(3) :: statistics
      ! ------------------------------------------------------------ !
      integer :: integ_success

      integ_success = 0
      if (this%last_integ_success) then
         integ_success = 1
      end if
      statistics = [integ_success, this%last_num_steps_accepted, this%last_num_steps_rejected]
   end function Get_Statistics
