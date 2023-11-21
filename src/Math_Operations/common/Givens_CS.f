   ! --------------------------------------------------------------- !
   pure function Givens_CS(refval, zeroval)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: refval, zeroval
      real(dp), dimension(2, 2) :: Givens_CS
      ! ------------------------------------------------------------ !
      real(dp) :: tau, cos_theta, sin_theta

      if (abs(zeroval) < setting_epsilon_extra) then
         cos_theta = 1.0_dp
         sin_theta = 0.0_dp
      else if (abs(zeroval) > abs(refval)) then
         tau = -refval/zeroval;
         sin_theta = 1.0_dp/sqrt(1.0_dp + tau**2)
         cos_theta = tau*sin_theta
      else
         tau = -zeroval/refval;
         cos_theta = 1.0_dp/sqrt(1.0_dp + tau**2)
         sin_theta = tau*cos_theta
      end if

      Givens_CS = reshape([cos_theta, sin_theta, -sin_theta, cos_theta], [2, 2])
   end function Givens_CS
