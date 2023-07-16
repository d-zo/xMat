   ! --------------------------------------------------------------- !
   pure function Dimensionless(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(6) :: Dimensionless
      ! ------------------------------------------------------------ !
      Dimensionless = Nonzero_Division(val=vec6, fac=Trace(vec6))    ! Dimensionless `\hat{\mathbf{M}} = \frac{\mathbf{M}}{\tr{(\mathbf{M})}}`
   end function Dimensionless
