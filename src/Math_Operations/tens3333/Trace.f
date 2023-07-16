   ! --------------------------------------------------------------- !
   pure function Trace(mat33) result(tr33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33
      real(dp) :: tr33
      ! ------------------------------------------------------------ !
      integer :: idx

      tr33 = 0.0_dp
      do idx = 1, global_num_direct_components
         tr33 = tr33 + mat33(idx, idx)                               ! Trace `\tr{(\mathbf{M})} = \sum_{i=1}^n M_{ii}`
      end do
   end function Trace
