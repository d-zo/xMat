   ! --------------------------------------------------------------- !
   pure function Deviatoric_Part(mat33) result(dev33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33
      real(dp), dimension(3, 3) :: dev33
      ! ------------------------------------------------------------ !
      real(dp) :: modtr
      integer :: idx

      dev33 = mat33                                                  ! Deviatoric part (3D) `\mathbf{M}^\mathrm{dev} = \mathbf{M} - \frac{\tr{(\mathbf{M})}}{n}\mathbf{I}`
      modtr = Trace(mat33)/global_num_direct_components
      do idx = 1, global_num_direct_components
         dev33(idx, idx) = dev33(idx, idx) - modtr
      end do
   end function Deviatoric_Part
