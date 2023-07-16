   ! --------------------------------------------------------------- !
   pure function Format_dp_dim4(name, tensor)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      real(dp), dimension(:, :, :, :), intent(in) :: tensor
      character(len=((len(submatrixname)+2+mat_identifierlength+2 + (onenumberlength+2)*size(tensor, 4)) &
         *size(tensor, 3) + len(name)+tens_identifierlength+1)*size(tensor, 2)*size(tensor, 1)) :: Format_dp_dim4
      ! ------------------------------------------------------------ !
      character(len=tens_identifierlength) :: identifier
      integer :: matrixlaenge, idx1, idx2, position

      matrixlaenge = (len(submatrixname)+2+mat_identifierlength+2 + (onenumberlength+2)*size(tensor, 4)) &
                   * size(tensor, 3) + len(name)+tens_identifierlength+1
      identifier = '(:, :, :, :)'
      do idx1 = 1, size(tensor, 1)
         write(identifier(2:2), '(i1)') idx1
         do idx2 = 1, size(tensor, 2)
            position = ((idx1-1)*size(tensor, 2) + (idx2-1))*matrixlaenge
            write(identifier(5:5), '(i1)') idx2

            write(Format_dp_dim4(position+1:position+matrixlaenge), '(a, a, a, a)') &
               name, identifier, char(10), Format_dp_dim2(name=submatrixname, matrix=tensor(idx1, idx2, :, :))
         end do
      end do
   end function Format_dp_dim4
