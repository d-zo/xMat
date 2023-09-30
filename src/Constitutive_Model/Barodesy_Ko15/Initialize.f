   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: Abort_If_Not_In_Interval
      !
      class(Barodesy_Ko15), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      real(dp) :: K_c

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.False.)

      ! --- Parameters
      this%param_phi_c = params(1)                                   ! Friction angle `\varphi_c` (in radians)
      this%param_c2    = params(2)                                   ! Parameter `c_2`
      this%param_c3    = params(3)                                   ! Parameter `c_3`
      this%param_c4    = params(4)                                   ! Parameter `c_4`
      this%param_c5    = params(5)                                   ! Parameter `c_5`
      this%param_e_c0  = params(6)                                   ! Critical void ratio `e_{c0}`
      this%param_e_min = params(7)                                   ! Minimal void ratio `e_\text{min}`

      ! --- Check valid parameter range
      if (firstcall) then
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
         call Abort_If_Not_In_Interval('e_c0', this%param_e_c0, [setting_epsilon, 20.0_dp])
         call Abort_If_Not_In_Interval('e_min', this%param_e_min, [setting_epsilon, 20.0_dp])
      end if

      K_c = (1.0_dp - sin(this%param_phi_c)) &                       ! `K_c = \frac{1-\sin{(\varphi_c)}}{1+\sin{(\varphi_c)}}`
          / (1.0_dp + sin(this%param_phi_c))

      this%derived_c1 = sqrt(2.0_dp/3.0_dp)*log(K_c)                 ! (32) of Kolymbas (2015): `c_1 = \sqrt{\frac{2}{3}}\log{(K_c)}`
   end subroutine Initialize
