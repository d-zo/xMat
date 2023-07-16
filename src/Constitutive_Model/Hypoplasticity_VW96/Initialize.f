   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculateJacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Write_Error_And_Exit
      use Debug, only: Formatval
      use Math_Operations, only: Abort_If_Not_In_Interval, const_root2, const_root3
      !
      class(Hypoplasticity_VW96), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculateJacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculateJacobian=calculateJacobian, provideJacobian=.True.)
      this%direct_variables_mask(11) = .True.

      ! --- Standard hypoplastic parameters
      this%param_phi_c   = params(1)                                 ! Friction angle `\varphi_c` (in radians)
      this%param_nu      = params(2)                                 ! Poisson ratio `\nu` for increased shear stiffness
      this%param_h_s     = params(3)                                 ! Granulate hardness `h_s` (in kPa)
      this%param_n_H     = params(4)                                 ! Exponent `n_H`
      this%param_e_d0    = params(5)                                 ! Densest void ratio for zero stress `e_{d0}`
      this%param_e_c0    = params(6)                                 ! Critical void ratio for zero stress `e_{c0}`
      this%param_e_i0    = params(7)                                 ! Loosest void ratio for zero stress `e_{i0}`
      this%param_alpha_H = params(8)                                 ! Pycnotropic exponent `\alpha_H`
      this%param_beta_H  = params(9)                                 ! Barotropic exponent `\beta_H`

      ! --- Extended parameters for intergranular strain
      this%param_m_T    = params(10)                                 ! Orthogonal loading history multiplier `m_T`
      this%param_m_R    = params(11)                                 ! Reverse loading history multiplier `m_R`
      this%param_R_max  = params(12)                                 ! Elastic strain range `R_\mathrm{max}`
      this%param_beta_R = params(13)                                 ! Parameter to control evolution of intergranular strain `\beta_R`
      this%param_Chi    = params(14)                                 ! Parameter to control degradation of stiffness `\chi`

      ! --- Additional parameters
      ! param(15) empty due to params structure compatibility with some other routines
      this%param_e_ini  = params(16)                                 ! Void ratio at K0-state `e_\mathrm{ini}`


      ! --- Check valid parameter range
      if (firstcall) then
         ! Limits taken/adjusted from user routine of Niemunis
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
         call Abort_If_Not_In_Interval('nu', this%param_nu, [0.0_dp, 0.5_dp])
         call Abort_If_Not_In_Interval('h_s', this%param_h_s, [100.0_dp, 1.0e9_dp])
         call Abort_If_Not_In_Interval('n_H', this%param_n_H, [0.0_dp, 10.0_dp])
         call Abort_If_Not_In_Interval('e_d0', this%param_e_d0, [setting_epsilon, 20.0_dp])
         call Abort_If_Not_In_Interval('e_c0', this%param_e_c0, [this%param_e_d0, 20.0_dp])
         call Abort_If_Not_In_Interval('e_i0', this%param_e_i0, [this%param_e_c0, 20.0_dp])
         call Abort_If_Not_In_Interval('m_T', this%param_m_T, [0.4_dp, 10.0_dp])
         call Abort_If_Not_In_Interval('m_R', this%param_m_R, [0.4_dp, 20.0_dp])
         call Abort_If_Not_In_Interval('R_max', this%param_R_max, [setting_epsilon, 0.1_dp])
         call Abort_If_Not_In_Interval('beta_R', this%param_beta_R, [setting_epsilon, 100.0_dp])
         call Abort_If_Not_In_Interval('Chi', this%param_Chi, [0.0_dp, 100.0_dp])
      end if

      ! --- Derived parameters
      this%derived_a = const_root3*(3.0_dp - sin(this%param_phi_c)) &! Precalculate (2.65) of Niemunis (2003): `a = \frac{\sqrt{3}\left(3-\sin{(\varphi_c)}\right)}{2\sqrt{2}\sin{(\varphi_c)}}`
                     / (2.0_dp*const_root2*sin(this%param_phi_c))
      this%derived_a2 = this%derived_a**2                            ! `a^2`

      this%derived_f_b_lastterm = 3.0_dp + this%derived_a2 &         ! Precalculate `f_{b,\mathrm{lastterm}} = 3+a^2-a\sqrt{3}\left(\frac{e_{i0}-e_{d0}}{e_{c0}-e_{d0}}\right)^{\alpha}`
                                - this%derived_a*const_root3 &       ! `f_{b,\mathrm{lastterm}}` is the last term of (2.72) in Niemunis (2003)
                                * (((this%param_e_i0 - this%param_e_d0) &
                                / (this%param_e_c0 - this%param_e_d0))**this%param_alpha_H)
      if (this%derived_f_b_lastterm <= 0.0_dp) then
         call Write_Error_And_Exit('Param_Hypoplasticity: Last term of f_b is not greater than zero ' // &
            '(sign of stresses would change)' // char(10) // Formatval('f_b_lastterm', this%derived_f_b_lastterm))
      end if

      this%derived_b2 = (1.0_dp + this%derived_a2/3.0_dp &           ! From (4.183) of Niemunis (2003): `b^2 = \left(1+\frac{a^2}{3} + \frac{a}{\sqrt{3}}\right)\frac{1-2\nu}{1+\nu} - 1`
                      + this%derived_a/const_root3) &                ! (needed for increased shear stiffness)
                      * (1.0_dp - 2.0_dp*this%param_nu)/(1.0_dp + this%param_nu) - 1.0_dp
   end subroutine Initialize
