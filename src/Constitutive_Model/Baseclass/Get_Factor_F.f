   ! --------------------------------------------------------------- !
   pure function Get_Factor_F(stress_dless_dev)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: const_root2, const_root3, const_root6, Norm, Determinant, Double_Contraction22
      !
      real(dp), dimension(__matrix__), intent(in) :: stress_dless_dev
      real(dp) :: Get_Factor_F
      ! ------------------------------------------------------------ !
      real(dp) :: tan_psi, stress_contraction, cos_3theta

      tan_psi = const_root3*Norm(stress_dless_dev)                   ! (2.67) (1) of Niemunis (2003): `\tan{(\psi)} = \sqrt{3}||\hat{\mathbf{T}}^\mathrm{dev}||`
      stress_contraction = Double_Contraction22(stress_dless_dev, stress_dless_dev)

      cos_3theta = 1.0_dp                                            ! Assume term equals one in case denominator `\hat{\mathbf{T}}^\mathrm{dev}:\hat{\mathbf{T}}^\mathrm{dev}` is zero
      if (stress_contraction > setting_epsilon) then
         ! For the given definition of `\hat{\mathbf{T}}^\mathrm{dev}`, the terms `\tr{(\hat{\mathbf{T}}^\mathrm{dev}\cdot\hat{\mathbf{T}}^\mathrm{dev}\cdot\hat{\mathbf{T}}^\mathrm{dev})}` and `3\det{(\hat{\mathbf{T}}^\mathrm{dev})}` are equal
         cos_3theta = -const_root6 * 3.0_dp &                        ! (2.67) (2) of Niemunis (2003): `\cos{(3\theta)} = -\sqrt{6}\frac{\tr{(\hat{\mathbf{T}}^\mathrm{dev}\cdot\hat{\mathbf{T}}^\mathrm{dev}\cdot\hat{\mathbf{T}}^\mathrm{dev})}}{\left[\hat{\mathbf{T}}^\mathrm{dev}:\hat{\mathbf{T}}^\mathrm{dev}\right]^{1.5}}`
                    * Determinant(stress_dless_dev) / (stress_contraction**1.5_dp)
      end if
      cos_3theta = max(-1.0_dp, min(1.0_dp, cos_3theta))             ! Restrict `\cos{(3\theta)}` in interval `[-1, 1]`

      Get_Factor_F = sqrt(abs(tan_psi*tan_psi/8.0_dp &               ! Modified (2.66) of Niemunis (2003):
                   + (2.0_dp - tan_psi*tan_psi) &                    ! `F = \sqrt{|\frac{1}{8}\tan^2{(\psi)} + \frac{2-\tan^2{(\psi)}}{2+\sqrt{2}\tan{(\psi)}\cos{(3\theta)}}|} - \frac{1}{2\sqrt{2}} \tan{(\psi)}`
                   /(2.0_dp + const_root2*tan_psi*cos_3theta))) - tan_psi/(2.0_dp*const_root2)
   end function Get_Factor_F
