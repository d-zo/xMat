   ! --------------------------------------------------------------- !
   pure subroutine Intergranular_Strain(this, L_mat, fdN_mat, D_vis, hypoplastic, igran_strain, &
      R_max, m_T, m_R, beta_R, chi, dt, dot_strain, dot_stress, dot_igran_strain, jacobian)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, setting_numerical_jacobian
      use Math_Operations, only: Nonzero_Division, Norm, Double_Contraction22, Double_Contraction42, &
                                 Double_Contraction44, Dyadic_Product22, const_identity4d_sym
      !
      class(Constitutive_Model), intent(in) :: this
      real(dp), dimension(__tensor__), intent(in) :: L_mat
      real(dp), dimension(__matrix__), intent(in) :: fdN_mat, D_vis, igran_strain
      logical, intent(in) :: hypoplastic
      real(dp), intent(in) :: R_max, m_T, m_R, beta_R, chi, dt
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      real(dp), dimension(__matrix__), intent(out) :: dot_stress, dot_igran_strain
      real(dp), dimension(__tensor__), intent(out) :: jacobian
      ! ------------------------------------------------------------ !
      real(dp) :: norm_h, norm_D, rho, load_dir, r_c, ref_dt, sub_dt, sub_dtmin, done_dt
      real(dp), dimension(__matrix__) :: h_dir, new_igran_strain
      real(dp), dimension(__tensor__) :: L_mod
      logical :: precheck, breakloop

      ref_dt = dt
      ! In case dt is zero, use ref_dt = 1 to prevent division by zero
      if (ref_dt < setting_epsilon) then
         ref_dt = 1.0_dp
      end if

      dot_stress = 0.0_dp
      jacobian = 0.0_dp

      done_dt = 0.0_dp
      new_igran_strain = igran_strain

      norm_D = Norm(dot_strain)
      ! Use given dt (which might be zero) only here to calculate precheck
      precheck = (norm_D*dt <= 0.1_dp*R_max)
      breakloop = .False.

      associate(M_mat => L_mod)
         igranloop: &
         do
            norm_h = Norm(new_igran_strain)
            rho = norm_h/R_max                                       ! (4.5) of Niemunis (2003): `\rho = \frac{||\mathbf{h}||}{R_\mathrm{max}}` with `0 \leq \rho \leq 1`

            load_dir = Double_Contraction22(new_igran_strain, &      ! Proportional to loading direction `\vec{\mathbf{h}}:\mathbf{D}`
                       dot_strain)

            if (precheck) then
               sub_dt = ref_dt
               breakloop = .True.
            else
               ! This branch can't be reached if either dt or norm_D is zero
               sub_dtmin = 0.1_dp*R_max/norm_D

               if (load_dir > 0.0_dp) then
                  if (rho < 0.99_dp) then
                     sub_dt = min(1.01_dp*ref_dt, abs(0.1_dp*R_max/((1.0_dp - rho**beta_R)*norm_D)))
                  else
                     sub_dt = 1.01_dp*ref_dt
                  end if
               else
                  sub_dt = sub_dtmin
               end if

               if ((rho > 0.01_dp) .and. (load_dir >= 0.0_dp)) then
                  ! if norm_h is zero, so is rho
                  sub_dt = sub_dtmin + (sub_dt - sub_dtmin)*load_dir/(norm_h*norm_D)
               end if

               if (done_dt + sub_dt >= ref_dt) then
                  sub_dt = ref_dt - done_dt
                  breakloop = .True.
               end if
            end if

            h_dir = Nonzero_Division(val=new_igran_strain, &         ! Strain direction `\vec{\mathbf{h}} = \frac{\mathbf{h}}{||\mathbf{h}||}`
               fac=norm_h)

            r_c = rho**Chi
            L_mod = Dyadic_Product22(Double_Contraction42(L_mat, &   ! `\mathcal{L}_\mathrm{mod} = \mathcal{L}:\vec{\mathbf{h}}\otimes\vec{\mathbf{h}}`
                    h_dir), h_dir)

            if (load_dir > 0.0_dp) then                              ! if `\vec{\mathbf{h}}:\mathbf{D} > 0` (loading branch) then use
               M_mat = (r_c*m_T + (1.0_dp - r_c)*m_R)*L_mat &        ! upper branch of (4.11) of Niemunis (2003):
                     + r_c*(1.0_dp - m_T)*L_mod                      ! `\mathcal{M} = \left[\rho^\chi m_T + (1 - \rho^\chi)m_R\right]\mathcal{L} + \rho^\chi(1 - m_T)\mathcal{L}_\mathrm{mod}`
               if (hypoplastic) then
                  M_mat = M_mat &                                    ! with hypoplasticity adjust `\mathcal{M} = \mathcal{M} + \rho^\chi \mathbf{N}\otimes \vec{\mathbf{h}}`
                        + r_c*Dyadic_Product22(fdN_mat, h_dir)
               end if

               dot_igran_strain = Double_Contraction42( &            ! upper branch of (4.12) of Niemunis (2003): `\overset{\circ}{\mathbf{h}} = \left(\mathcal{I} - \vec{\mathbf{h}}\otimes\vec{\mathbf{h}} \rho^{\beta_R}\right):\mathbf{D}`
                                  const_identity4d_sym - rho**beta_R*Dyadic_Product22(h_dir, h_dir), dot_strain)
            else
               M_mat = (r_c*m_T + (1.0_dp - r_c)*m_R)*L_mat &        ! lower branch of (4.11) of Niemunis (2003):
                     + r_c*(m_R - m_T)*L_mod                         ! `\mathcal{M} = \left[\rho^\chi m_T + (1 - \rho^\chi)m_R\right]\mathcal{L} + \rho^\chi(m_R - m_T)\mathcal{L}_\mathrm{mod}`
               dot_igran_strain = dot_strain                         ! lower branch of (4.12) of Niemunis (2003): `\overset{\circ}{\mathbf{h}} = \mathbf{D}`
            end if

            ! Estimate the integrated intergranular strain and determine its norm. Adjust it, if it is greater than R_max
            new_igran_strain = new_igran_strain + dot_igran_strain*sub_dt
            norm_h = Norm(new_igran_strain)

            if (norm_h > R_max) then
               new_igran_strain = new_igran_strain * R_max/norm_h
            end if

            dot_stress = dot_stress &                                ! Incremental variant of (4.6) of Niemunis (2003): `\overset{\circ}{\mathbf{T}} = \mathcal{M}:\mathbf{D}`
                       + Double_Contraction42(M_mat, dot_strain)*sub_dt/ref_dt
            if (.not. hypoplastic) then
               dot_stress = dot_stress &                             ! Incremental variant of (4.121) of Niemunis (2003): `\overset{\circ}{\mathbf{T}} = \mathcal{M}:\mathbf{D} - \mathcal{L}:\mathbf{D}^\mathrm{vis}`
                          - Double_Contraction42(L_mat, D_vis)*sub_dt/ref_dt
            end if

            done_dt = done_dt + sub_dt

            if ((this%calculate_jacobian) .and. (.not. setting_numerical_jacobian)) then
               jacobian = jacobian + M_mat*sub_dt/ref_dt
            end if

            if (breakloop) then
               exit igranloop
            end if
         end do igranloop
         dot_igran_strain = (new_igran_strain - igran_strain)/ref_dt
      end associate
   end subroutine Intergranular_Strain
