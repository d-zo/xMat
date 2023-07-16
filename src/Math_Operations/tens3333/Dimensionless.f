   ! --------------------------------------------------------------- !
   pure function Dimensionless(mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33
      real(dp), dimension(3, 3) :: Dimensionless
      ! ------------------------------------------------------------ !
      Dimensionless = Nonzero_Division(val=mat33, fac=Trace(mat33))  ! Dimensionless `\hat{\mathbf{M}} = \frac{\mathbf{M}}{\tr{(\mathbf{M})}}`
   end function Dimensionless
