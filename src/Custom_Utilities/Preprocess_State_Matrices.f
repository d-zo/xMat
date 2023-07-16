   ! --------------------------------------------------------------- !
   function Preprocess_State_Matrices(identifier, intername, nstates, states, rotation) result(mod_states)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Write_Error_And_Exit
      use Math_Operations, only: Vec9_To_Mat, Matrix_Rotation, Mat_To_Vec9
      !
      character(len=setting_len_id), intent(in) :: identifier
      character(len=13), intent(in) :: intername
      integer, intent(in) :: nstates
      real(dp), dimension(nstates), intent(in) :: states
      real(dp), dimension(3, 3), intent(in) :: rotation
      real(dp), dimension(nstates) :: mod_states
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: save_mat
      real(dp), dimension(9) :: save_vec

      mod_states = states

      if ((identifier == setting_id_hypo_vw96) .or. (identifier == setting_id_viscohypo_ni03)) then
         ! Adjust intergranular strain
         if (nstates < 10) then
            call Write_Error_And_Exit('Preprocess_State_Matrices: Expecting at least ' // &
               '10 (internal) states for this constitutive model')
         end if

         save_vec = [mod_states(2:4), 0.5_dp*mod_states(5:10)]

         ! Rotate previous intergranular strain to new coordinate orientation
         save_mat = Vec9_To_Mat(vec9=save_vec)
         save_mat = Matrix_Rotation(save_mat, rotation)
         save_vec = Mat_To_Vec9(mat=save_mat)

         if (intername == 'umat_abq_std ') then
            ! Switch order corresponding to Abaqus/Standard after transformation based on internal representation
            ! (1, 1), (2, 2), (3, 3), (1, 2), (2, 1), (2, 3), (3, 2), (1, 3), (3, 1)
            save_vec(6:9) = [save_vec(8), save_vec(9), save_vec(6), save_vec(7)]
         end if

         mod_states(2:10) = save_vec
      else if (identifier == setting_id_barodesy_ko21) then
         ! Adjust last saved strain
         if (nstates < 12) then
            call Write_Error_And_Exit('Preprocess_State_Matrices: Expecting at least ' // &
               '12 (internal) states for this constitutive model')
         end if

         ! Rotate last saved strain to new coordinate orientation
         save_mat = Vec9_To_Mat(vec9=mod_states(4:12))
         save_mat = Matrix_Rotation(save_mat, rotation)
         mod_states(4:12) = Mat_To_Vec9(mat=save_mat)
      end if
   end function Preprocess_State_Matrices
