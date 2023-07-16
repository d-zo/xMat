   ! --------------------------------------------------------------- !
   pure function Deviatoric_Part(vec6) result(dev6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(6) :: dev6
      ! ------------------------------------------------------------ !
      dev6 = vec6                                                    ! Deviatoric part (3D) `\mathbf{M}^\mathrm{dev} = \mathbf{M} - \frac{\tr{(\mathbf{M})}}{n}\mathbf{I}`
      associate(ndi => global_num_direct_components)
         dev6(1:ndi) = dev6(1:ndi) - Trace(vec6)/ndi
      end associate
   end function Deviatoric_Part
