   ! --------------------------------------------------------------- !
   pure function Double_Contraction44(tens3333a, tens3333b) result(tens_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3, 3, 3), intent(in) :: tens3333a, tens3333b
      real(dp), dimension(3, 3, 3, 3) :: tens_out
      ! ------------------------------------------------------------ !
      real(dp), dimension(9) :: matveca, matvecb
      integer :: idx, jdx, kdx, ldx

      do idx = 1, 3                                                  ! Double contraction `\mathcal{T} : \mathcal{F} = \sum_{kl} T_{ijkl}\cdot F_{klmn}`
         do jdx = 1, 3
            matveca = reshape(tens3333a(idx, jdx, :, :), [9])

            do kdx = 1, 3
               do ldx = 1, 3
                  matvecb = reshape(tens3333b(:, :, kdx, ldx), [9])
                  tens_out(idx, jdx, kdx, ldx) = dot_product(matveca, matvecb)
               end do
            end do
         end do
      end do
   end function Double_Contraction44
