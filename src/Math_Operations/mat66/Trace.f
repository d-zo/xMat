   ! --------------------------------------------------------------- !
   pure function Trace(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp) :: Trace
      ! ------------------------------------------------------------ !
      Trace = sum(vec6(1:global_num_direct_components))              ! Trace `\tr{(\mathbf{M})} = \sum_{i=1}^n M_{ii}`
   end function Trace
