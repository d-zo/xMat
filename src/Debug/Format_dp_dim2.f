   ! --------------------------------------------------------------- !
   pure function Format_dp_dim2(name, matrix)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      real(dp), dimension(:, :), intent(in) :: matrix
      character(len=(len(name)+2 + mat_identifierlength+2 + &
         (onenumberlength+2)*size(matrix, 1))*size(matrix, 2)) :: Format_dp_dim2
      ! ------------------------------------------------------------ !
      character(len=mat_identifierlength) :: identifier
      integer :: zeilenlaenge, idx

      zeilenlaenge = len(name)+2 + mat_identifierlength+2 + (onenumberlength+2)*size(matrix, 1)
      identifier = '(:, :)'
      do idx = 1, size(matrix, 2)
         write(identifier(2:2), '(i1)') idx
         write(Format_dp_dim2((idx-1)*zeilenlaenge+1:idx*zeilenlaenge), '(a, a, a, a)') &
            name, identifier, '  ', Format_dp_dim1(name='', vector=matrix(idx, :))
      end do
   end function Format_dp_dim2
