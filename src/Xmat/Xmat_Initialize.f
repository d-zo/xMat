   ! --------------------------------------------------------------- !
   subroutine Xmat_Initialize(xmat_obj, solver_name, constitutive_model_name, material_parameters, &
      calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      class(Xmat), allocatable, intent(out) :: xmat_obj
      character(len=*), intent(in) :: solver_name
      character(len=setting_len_id), intent(in) :: constitutive_model_name
      real(dp), dimension(:), intent(in) :: material_parameters
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      if (.not. allocated(xmat_obj)) then
         allocate(Xmat::xmat_obj)
      end if

      call Select_Solver(new_solver=xmat_obj%xmat_solver, name=solver_name)
      call Select_Constitutive_Model(new_constitutive_model=xmat_obj%xmat_constitutive_model, &
         identifier=constitutive_model_name, params=material_parameters, calculate_jacobian=calculate_jacobian, &
         firstcall=firstcall)
      call xmat_obj%xmat_solver%Set_State_Mask(statemask=xmat_obj%xmat_constitutive_model%Get_Direct_Variables_Mask())
   end subroutine Xmat_Initialize
