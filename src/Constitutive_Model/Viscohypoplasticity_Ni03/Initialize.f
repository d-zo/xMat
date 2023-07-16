   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculateJacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Write_Error_And_Exit
      use Debug, only: Formatval
      use Math_Operations, only: Abort_If_Not_In_Interval, const_root2, const_root3
      !
      class(Viscohypoplasticity_Ni03), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculateJacobian, firstcall
      ! ------------------------------------------------------------ !
      real(dp) :: K_0

      ! --- Initialize base class variables
      call this%Base_Initialization(calculateJacobian=calculateJacobian, provideJacobian=.True.)
      this%direct_variables_mask(11) = .True.

      ! --- Standard viscohypoplastic parameters
      this%param_e_100  = params(1)                                  ! Void ratio `e_{100}` for `p_{e_{100}}=100 kPa` and chosen `D_r`
      this%param_nu     = params(2)                                  ! Poisson ratio `\nu` for increased shear stiffness
      this%param_lambda = params(3)                                  ! Slope of normal consolidated line in void ratio-pressure diagram `\lambda`
      this%param_kappa  = params(4)                                  ! Slope of reloading branch in void ratio-pressure diagram `\kappa`
      this%param_beta_b = params(5)                                  ! Factor `\beta_b` to describe loading shape
      this%param_I_v    = params(6)                                  ! Viscosity index `I_v`
      this%param_D_r    = params(7)                                  ! Reference creep rate `D_r`
      this%param_phi_c  = params(8)                                  ! Friction angle `\varphi_c` (in radians)

      ! --- Extended parameters for intergranular strain
      this%param_m_T    = params(10)                                 ! Orthogonal loading history multiplier `m_T`
      this%param_m_R    = params(11)                                 ! Reverse loading history multiplier `m_R`
      this%param_R_max  = params(12)                                 ! Elastic strain range `R_\mathrm{max}`
      this%param_beta_R = params(13)                                 ! Parameter to control evolution of intergranular strain `\beta_R`
      this%param_chi    = params(14)                                 ! Parameter to control degradation of stiffness `\chi`

      ! --- Additional viscohypoplastic parameters
      this%param_OCR    = params(15)                                 ! Overconsolidation ratio `OCR` (currently not used)

      ! --- Check valid parameter range
      if (firstcall) then
         ! Limits taken/adjusted from user routine of Niemunis
         call Abort_If_Not_In_Interval('e_100', this%param_e_100, [setting_epsilon, 20.0_dp])
         call Abort_If_Not_In_Interval('nu', this%param_nu, [0.0_dp, 0.5_dp])
         call Abort_If_Not_In_Interval('lambda', this%param_lambda, [0.0_dp, 2.0_dp])
         call Abort_If_Not_In_Interval('kappa', this%param_kappa, [0.0_dp, this%param_lambda])
         call Abort_If_Not_In_Interval('beta_b', this%param_beta_b, [0.0_dp, 1.0_dp])
         call Abort_If_Not_In_Interval('I_v', this%param_I_v, [0.0_dp, 1.0_dp])
         call Abort_If_Not_In_Interval('D_r', this%param_D_r, [setting_epsilon, 1.0e9_dp])
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
         call Abort_If_Not_In_Interval('m_T', this%param_m_T, [0.4_dp, 10.0_dp])
         call Abort_If_Not_In_Interval('m_R', this%param_m_R, [0.4_dp, 20.0_dp])
         call Abort_If_Not_In_Interval('R_max', this%param_R_max, [setting_epsilon, 0.1_dp])
         call Abort_If_Not_In_Interval('beta_R', this%param_beta_R, [setting_epsilon, 100.0_dp])
         call Abort_If_Not_In_Interval('chi', this%param_chi, [0.0_dp, 100.0_dp])
      end if

      ! --- Derived parameters
      this%derived_a = const_root3*(3.0_dp - sin(this%param_phi_c)) &! Precalculate (2.65) of Niemunis (2003): `a = \frac{\sqrt{3}\left(3-\sin{(\varphi_c)}\right)}{2\sqrt{2}\sin{(\varphi_c)}}`
                     / (2.0_dp*const_root2*sin(this%param_phi_c))
      this%derived_a2 = this%derived_a**2                            ! `a^2`

      this%derived_pq_incl = 6.0_dp*sin(this%param_phi_c) &          ! inclination in p-q diagram `\frac{6\sin{(\varphi_c)}}{3-\sin{(\varphi_c)}}`
                           / (3.0_dp - sin(this%param_phi_c))

      K_0 = (-2.0_dp - this%derived_a2 + sqrt(36.0_dp &              ! (4.94) of Niemunis (2003): `K_0^{up} = K_0 = \frac{-2-a^2+\sqrt{36+36a^2+a^4}}{16}`
          + 36.0_dp*this%derived_a2 + this%derived_a2**2))/16.0_dp   ! which resembles the upper limit of stress ratio `K_0`
      this%derived_beta_b = 1.0_dp/(this%param_kappa * (1.0_dp &     ! (4.72) of Niemunis (2003): `\beta_b = \frac{1}{\kappa_0\left(1 + a^2/(1 + 2K_0)\right)}`
                          + this%derived_a2 / (1.0_dp + 2.0_dp*K_0)))

      this%derived_b2 = (1.0_dp + this%derived_a2/3.0_dp) &          ! Modified (4.183) of Niemunis (2003): `b^2 = \left(1+\frac{a^2}{3}\right)\frac{1-2\nu}{1+\nu} - 1`
                      * (1.0_dp - 2.0_dp*this%param_nu) &            ! (needed for increased shear stiffness)
                      / (1.0_dp + this%param_nu) - 1.0_dp
   end subroutine Initialize
