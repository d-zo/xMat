   ! --------------------------------------------------------------- !
   pure function Mat99_To_Tens(mat99) result(tens3333)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9), intent(in) :: mat99
      real(dp), dimension(3, 3, 3, 3) :: tens3333
      ! ------------------------------------------------------------ !
      integer :: m_idx, m_jdx, t_idx, t_jdx, t_kdx, t_ldx

      do m_jdx = 1, 9
         call Ref_Index(ref_idx=m_jdx, idx=t_kdx, jdx=t_ldx)
         do m_idx = 1, 9
            call Ref_Index(ref_idx=m_idx, idx=t_idx, jdx=t_jdx)
            tens3333(t_idx, t_jdx, t_kdx, t_ldx) = mat99(m_idx, m_jdx)
         end do
      end do
   end function Mat99_To_Tens
