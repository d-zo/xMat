   ! --------------------------------------------------------------- !
   pure function Dimensionless(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(9) :: Dimensionless
      ! ------------------------------------------------------------ !
      Dimensionless = Nonzero_Division(val=vec9, fac=Trace(vec9))    ! Dimensionless `\hat{\mathbf{M}} = \frac{\mathbf{M}}{\tr{(\mathbf{M})}}`
   end function Dimensionless
