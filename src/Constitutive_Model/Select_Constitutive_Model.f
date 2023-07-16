   ! --------------------------------------------------------------- !
   subroutine Select_Constitutive_Model(new_constitutive_model, identifier, params, calculateJacobian, firstcall)
   ! --------------------------------------------------------------- !
      use Elasticity_Class
      use Hypoplasticity_Wu92_Class
      use Hypoplasticity_VW96_Class
      use Viscohypoplasticity_Ni03_Class
      use Barodesy_Ko15_Class
      use Barodesy_Sc18_Class
      use Barodesy_Ko21_Class
      use Test_DGL_Class
      !
      class(Constitutive_Model), allocatable, intent(out) :: new_constitutive_model
      character(len=setting_len_id), intent(in) :: identifier
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculateJacobian, firstcall
      ! ------------------------------------------------------------ !
      if (identifier == setting_id_elasticity) then
         allocate(Elasticity::new_constitutive_model)
      else if (identifier == setting_id_hypo_wu92) then
         allocate(Hypoplasticity_Wu92::new_constitutive_model)
      else if (identifier == setting_id_hypo_vw96) then
         allocate(Hypoplasticity_VW96::new_constitutive_model)
      else if (identifier == setting_id_viscohypo_ni03) then
         allocate(Viscohypoplasticity_Ni03::new_constitutive_model)
      else if (identifier == setting_id_barodesy_ko15) then
         allocate(Barodesy_Ko15::new_constitutive_model)
      else if (identifier == setting_id_barodesy_sc18) then
         allocate(Barodesy_Sc18::new_constitutive_model)
      else if (identifier == setting_id_barodesy_ko21) then
         allocate(Barodesy_Ko21::new_constitutive_model)
      else if (identifier == setting_id_test_dgl) then
         allocate(Test_DGL::new_constitutive_model)
      else
         call Write_Error_And_Exit('Select_Constitutive_Model: Identifier >' // identifier // '< unknown')
      end if
      call new_constitutive_model%Initialize(params=params, calculateJacobian=calculateJacobian, &
         firstcall=firstcall)
   end subroutine Select_Constitutive_Model
