   ! --------------------------------------------------------------- !
   subroutine Set_Values(this, overall_dt, dot_strain)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(inout) :: this
      real(dp), intent(in) :: overall_dt
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      ! ------------------------------------------------------------ !
      this%overall_dt = overall_dt
      this%dot_strain = dot_strain
   end subroutine Set_Values
