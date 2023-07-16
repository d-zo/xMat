   ! --------------------------------------------------------------- !
   pure function Determinant(mat33) result(det33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33
      real(dp) :: det33
      ! ------------------------------------------------------------ !
      if (global_num_direct_components == 3) then
         det33 = mat33(1, 1)*(mat33(2, 2)*mat33(3, 3) - mat33(2, 3)*mat33(3, 2)) &
               + mat33(1, 2)*(mat33(2, 3)*mat33(3, 1) - mat33(2, 1)*mat33(3, 3)) &
               + mat33(1, 3)*(mat33(2, 1)*mat33(3, 2) - mat33(2, 2)*mat33(3, 1))
      else
         det33 = mat33(1, 1)*mat33(2, 2) - mat33(1, 2)*mat33(2, 1)
      end if
   end function Determinant
