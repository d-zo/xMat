   ! --------------------------------------------------------------- !
   pure function Postprocess_State_Matrices(identifier, intername, nstates, states) result(mod_states)
   ! --------------------------------------------------------------- !
      character(len=setting_len_id), intent(in) :: identifier
      character(len=13), intent(in) :: intername
      integer, intent(in) :: nstates
      real(dp), dimension(nstates), intent(in) :: states
      real(dp), dimension(nstates) :: mod_states
      ! ------------------------------------------------------------ !
      real(dp), dimension(9) :: ig_vec

      mod_states = states

      if ((identifier == setting_id_hypo_vw96) .or. (identifier == setting_id_viscohypo_ni03)) then
         ! Adjust intergranular strain
         ig_vec = [mod_states(2:4), 2.0_dp*mod_states(5:10)]

         if (intername == 'umat_abq_std ') then
            ! Switch order corresponding to Abaqus/Standard based on internal representation
            ! (1, 1), (2, 2), (3, 3), (1, 2), (2, 1), (2, 3), (3, 2), (1, 3), (3, 1)
            ig_vec(6:9) = [ig_vec(8), ig_vec(9), ig_vec(6), ig_vec(7)]
         end if
         mod_states(2:10) = ig_vec
      end if
   end function Postprocess_State_Matrices
