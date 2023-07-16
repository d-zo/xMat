   ! --------------------------------------------------------------- !
   pure function Deviatoric_Part(vec9) result(dev9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(9) :: dev9
      ! ------------------------------------------------------------ !
      dev9 = vec9                                                    ! Deviatoric part (3D) `\mathbf{M}^\mathrm{dev} = \mathbf{M} - \frac{\tr{(\mathbf{M})}}{n}\mathbf{I}`
      associate(ndi => global_num_direct_components)
         dev9(1:ndi) = dev9(1:ndi) - Trace(vec9)/ndi
      end associate
   end function Deviatoric_Part
