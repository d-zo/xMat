   ! --------------------------------------------------------------- !
   pure subroutine Elasticity(this, youngs_modulus, nu, dot_strain, dot_stress, jacobian)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_numerical_jacobian
      use Math_Operations, only: Double_Contraction42, const_identity4d_sym, const_identity4d_tr
      !
      class(Constitutive_Model), intent(in) :: this
      real(dp), intent(in) :: youngs_modulus, nu
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      real(dp), dimension(__matrix__), intent(out) :: dot_stress
      real(dp), dimension(__tensor__), intent(out) :: jacobian
      ! ------------------------------------------------------------ !
      associate(Youngsmat => jacobian)
         ! Using the formula is not the most efficient way, but it is independent of the used representation
         Youngsmat = youngs_modulus/(1.0_dp + nu) &                  ! `\mathcal{E} = \frac{E}{1+\nu}\left(\mathcal{I}^\mathrm{sym} + \frac{\nu}{1-2\nu}\mathcal{I}^\mathrm{tr}\right)`
                   * const_identity4d_sym + youngs_modulus*nu/(1.0_dp - nu - 2.0_dp*nu**2)*const_identity4d_tr
         dot_stress = Double_Contraction42(Youngsmat, dot_strain)    ! `\mathbf{\dot{T}} = \mathcal{E}:\mathbf{D}`
      end associate

      if ((.not. this%calculateJacobian) .or. (setting_numerical_jacobian)) then
         jacobian = 0.0_dp
      end if
   end subroutine Elasticity
