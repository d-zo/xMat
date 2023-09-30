   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: Abort_If_Not_In_Interval, const_root2, const_root3
      !
      class(Barodesy_Sc18), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      ! Auxiliary terms `f_1`, `f_2`, `f_3` and `f_4`
      real(dp), parameter :: fak1 = log((1.0_dp + const_root2)/const_root2)
      real(dp), parameter :: fak2 = log((3.0_dp*const_root2 - const_root3)/(3.0_dp*const_root2))
      real(dp), parameter :: fak3 = log((1.0_dp + const_root2)/(2.0_dp + const_root2))
      real(dp), parameter :: fak4 = fak1*log((3.0_dp + 3.0_dp*const_root2) &
                                  / (3.0_dp + 3.0_dp*const_root2 - const_root3)) - fak3*fak2
      ! NOTE: Consider making cmin a class variable since it is also used in Calculate_Dot_State()
      real(dp), parameter :: cmin = -30.0_dp
      real(dp), parameter :: p_r = 1.0_dp                            ! Reference pressure `p_r` for parameters of barodey is 1 kPa
      real(dp) :: K_c, alpha_p, alpha_c, alpha_0

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.False.)

      ! --- Parameters
      this%param_phi_c = params(1)                                   ! Friction angle `\varphi_c` (in radians)
      this%param_K_r   = params(2)                                   ! Reference stiffness `K_r` at 1 kPa
      this%param_xi    = params(3)                                   ! Exponent `\xi` in stiffness term
      this%param_kappa = params(4)                                   ! Stiffness releation `\kappa` between loading and unloading
      this%param_e_c0  = params(5)                                   ! Critical void ratio for zero stress `e_{c0}`
      this%param_ce    = params(6)                                   ! Ratio `c_e` of crit. and comparable void ratio at normal compression
      this%param_c5    = params(7)                                   ! Weighting factor `c_5` for pycnotropy
      this%param_c6    = params(8)                                   ! Control factor `c_6` for splitting `f` and `g`
      this%param_c7    = params(9)                                   ! Interpolation factor `c_7` of `b` between compression and extension

      ! --- Check valid parameter range
      if (firstcall) then
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
         call Abort_If_Not_In_Interval('K_r', this%param_K_r, [0.1_dp, 1.0e9_dp])
         call Abort_If_Not_In_Interval('xi', this%param_xi, [0.2_dp, 2.0_dp])
         call Abort_If_Not_In_Interval('kappa', this%param_kappa, [setting_epsilon, 1.0e9_dp])
         call Abort_If_Not_In_Interval('e_c0', this%param_e_c0, [setting_epsilon, 20.0_dp])
         call Abort_If_Not_In_Interval('c_e', this%param_ce, [1.0_dp, 5.0_dp])
         call Abort_If_Not_In_Interval('c_5', this%param_c5, [0.1_dp, 100.0_dp])
         call Abort_If_Not_In_Interval('c_6', this%param_c6, [0.1_dp, 100.0_dp])
         call Abort_If_Not_In_Interval('c_7', this%param_c7, [0.1_dp, 100.0_dp])
      end if

      ! --- Derived parameters (taken/adapted from F. Schranz (2018), p.166 and p.175)
      K_c = (1.0_dp - sin(this%param_phi_c)) &                       ! `K_c = \frac{1-\sin{(\varphi_c)}}{1+\sin{(\varphi_c)}}`
          / (1.0_dp + sin(this%param_phi_c))

      alpha_p = log((9.0_dp - 11.0_dp*sin(this%param_phi_c)) &       ! `\alpha_p = \log{\left(\frac{9-11\sin{(\varphi_c)}}{9+13\sin{(\varphi_c)}}\right)}`
              / (9.0_dp + 13.0_dp*sin(this%param_phi_c)))*const_root3/2.0_dp
      alpha_c = const_root2/const_root3*log(K_c)                     ! `\alpha_c = \sqrt{\frac{2}{3}}\log{\left(\frac{1-\sin{(\varphi_c)}}{1+\sin{(\varphi_c)}}\right)}`
      alpha_0 = log(1.0_dp - sin(this%param_phi_c))                  ! `\alpha_0 = \log{\left(1-\sin{(\varphi_c)}\right)}`

      this%derived_c1 = (log((alpha_p - cmin) / (alpha_c - cmin)) &  ! `c_1 = \frac{f_1}{f_4}\log{\left(\frac{\alpha_p+30}{\alpha_c+30}\right)} - \frac{f_2}{f_4}\log{\left(\frac{\alpha_0+30}{\alpha_c+30}\right)}`
                      * fak1 - log((alpha_0 - cmin)/(alpha_c - cmin)) * fak2)/fak4
      this%derived_c2 = (log((alpha_0 - cmin)/(alpha_c - cmin)) &    ! `c_2 = \frac{\log{\left(\frac{\alpha_0+30}{\alpha_c+30}\right)} - f_3 c_1}{f_1}`
                      - fak3*this%derived_c1)/fak1
      this%derived_c3 = (alpha_c - cmin) &                           ! `c_3 = \left(\alpha_c + 30\right)\frac{\left(1+\sqrt{2}\right)^{c_1}}{\sqrt{2}^{c_2}}`
                      * (1.0_dp + const_root2)**this%derived_c1 / const_root2**this%derived_c2
      this%derived_c4 = this%param_K_r &                             ! `c_4 = K_r\left(\sqrt{3} p_r\right)^{-\xi}`
                      * (const_root3*p_r)**(-this%param_xi)
      this%derived_zeta = -1.0_dp/K_c                                ! `\zeta = -\frac{1}{K_c} = -\frac{1+\sin{(\varphi_c)}}{1-\sin{(\varphi_c)}}`

      this%derived_b_comp = (this%param_c5 &                         ! (4.50) of Schranz (2018): `b_\mathrm{comp} = \frac{c_5\left(c_e^\zeta - 1\right) -3}{\sqrt{3}}`
                          * (this%param_ce**this%derived_zeta - 1.0_dp) - 3.0_dp)/const_root3
      this%derived_b_ext = (this%param_c5 &                          ! (4.56) of Schranz (2018): `b_\mathrm{ext} = \frac{c_5\left(1 - c_e^\zeta\right) -3\kappa}{\sqrt{3}}`
                         * (1.0_dp - this%param_ce**this%derived_zeta) - 3.0_dp*this%param_kappa)/const_root3

      this%derived_c8 = -this%param_c6 &                             ! (4.65) of Schranz (2018): `c_8 = \frac{-c_6}{\sqrt{3} b_\mathrm{ext}}`
                      / (const_root3*this%derived_b_ext)
   end subroutine Initialize
