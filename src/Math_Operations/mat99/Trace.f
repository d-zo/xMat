   ! --------------------------------------------------------------- !
   pure function Trace(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp) :: Trace
      ! ------------------------------------------------------------ !
      Trace = sum(vec9(1:global_num_direct_components))              ! Trace `\tr{(\mathbf{M})} = \sum_{i=1}^n M_{ii}`
   end function Trace
