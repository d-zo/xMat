   ! --------------------------------------------------------------- !
   pure function Squared(vec6) result(sq_vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(6) :: sq_vec6
      ! ------------------------------------------------------------ !
      sq_vec6(1) = vec6(1)**2 + vec6(4)**2 + vec6(6)**2              ! Square `\mathbf{M}^2 = \mathbf{M}\cdot \mathbf{M}`
      sq_vec6(2) = vec6(4)**2 + vec6(2)**2 + vec6(5)**2              ! (square of a symmetric matrix is also symmetric)
      sq_vec6(3) = vec6(6)**2 + vec6(5)**2 + vec6(3)**2
      sq_vec6(4) = vec6(4)*(vec6(1) + vec6(2)) + vec6(5)*vec6(6)
      sq_vec6(5) = vec6(5)*(vec6(2) + vec6(3)) + vec6(4)*vec6(6)
      sq_vec6(6) = vec6(6)*(vec6(1) + vec6(3)) + vec6(4)*vec6(5)
   end function Squared
