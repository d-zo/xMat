! ==================================================================================================================== !
module Math_Operations
   use General_Settings, only: dp, setting_epsilon, setting_epsilon_extra, global_num_direct_components, &
                               global_num_shear_components
   implicit none

   real(dp), parameter :: const_root2 = sqrt(2.0_dp)
   real(dp), parameter :: const_root3 = sqrt(3.0_dp)
   real(dp), parameter :: const_root6 = sqrt(6.0_dp)
   !
   real(dp), parameter, dimension(3, 3) :: const_identity2d = &      ! 2dim identity matrix `\mathbf{I}` where `I_{ij} = \delta_{ij}`
      reshape([ &                                                    ! with Kronecker-symbol `\delta_{ij} = 1` for `i=j` and 0 else
      1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
   ! 4dim symmetric identity tensor `\mathcal{I}^\mathrm{sym}` with `I^\mathrm{sym}_{ijkl} = \frac{1}{2}\left(\delta_{ik}\otimes\delta_{jl} + \delta_{il}\otimes\delta_{jk}\right)` which returns the symmetric part of the matrix: `\mathcal{I}^\mathrm{sym}:\mathbf{A} = \mathrm{sym}\left(\mathbf{A}\right)`
   real(dp), parameter, dimension(3, 3, 3, 3) :: const_identity4d_sym = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3, 3, 3])
   ! 4dim identity tensor representing `\mathcal{I}^\mathrm{tr}` with `I^\mathrm{tr}_{ijkl} = \delta_{ij}\otimes\delta_{kl}` which maps the trace of a matrix to each element of a diagonal matrix:
   ! `\mathcal{I}^\mathrm{tr}:\mathbf{A} = \tr{(\mathbf{A})}\cdot\delta_{ii}`
   real(dp), parameter, dimension(3, 3, 3, 3) :: const_identity4d_tr = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3, 3, 3])
   ! This index ordering is used to switch between tensor and matrix notation
   ! and is not the order of the tensor elements
   integer, parameter, dimension(2, 9) :: ref_indices = reshape([ &
      1, 1,   2, 2,   3, 3,   1, 2,   2, 1,   2, 3,   3, 2,   1, 3,   3, 1], [2, 9])

   private
   public const_root2, const_root3, const_root6, &
          const_identity2d, const_identity4d_sym, const_identity4d_tr, &
          Nonzero_Division, Trace, Dimensionless, Deviatoric_Part, Norm, Determinant, Squared, &
          Inverse_Tensor, Double_Contraction22, Double_Contraction42, Double_Contraction44, Dyadic_Product22, &
          Tensor_Partialtrace, Matrix_Rotation, Matrix_Exponential, Sqrt_Of_Sum_Of_Squares, &
          Is_Diagonal, Value_In_Interval, Abort_If_Not_In_Interval, Set_Element_In_Tensor, &
          Get_Matrix_From_Tensorpart, Set_Matrix_In_Tensorpart, Get_Elements_From_Matrixlist, &
          Set_Elements_In_Matrixlist, Ref_Index, Perturbate, Vec9_To_Mat, Mat_To_Vec9, &
          Pack_States, Unpack_States, Import_Matrix, Export_Matrix, Export_Tensor


   contains


#addfile function Nonzero_Division


#addfile function Trace


#addfile function Dimensionless


#addfile function Deviatoric_Part


#addfile function Norm


#addfile function Determinant


#addfile function Squared


#addfile function LU_Decomposition


#addfile function LU_Solve


#addfile function Inverse_Internal


#addfile function Inverse_Tensor


#addfile function Double_Contraction22


#addfile function Double_Contraction42


#addfile function Double_Contraction44


#addfile function Dyadic_Product22


#addfile function Tensor_Partialtrace


#addfile function Matrix_Rotation


#addfile function Sqrt_Of_Sum_Of_Squares


#addfile subroutine Tridiagonal_Matrix


#addfile function Givens_CS


#addfile subroutine QL_Decomposition


#addfile subroutine Eigendecomposition


#addfile function Matrix_Exponential


#addfile function Is_Diagonal


#addfile subroutine Value_In_Interval


#addfile subroutine Abort_If_Not_In_Interval


#addfile subroutine Set_Element_In_Tensor


#addfile function Get_Matrix_From_Tensorpart


#addfile subroutine Set_Matrix_In_Tensorpart


#addfile function Get_Elements_From_Matrixlist


#addfile subroutine Set_Elements_In_Matrixlist


#addfile subroutine Ref_Index


#addfile subroutine Perturbate


#addfile function Mat99_To_Tens


#addfile function Tens_To_Mat99


#addfile function Vec9_To_Mat


#addfile function Mat_To_Vec9


#addfile function Pack_States


#addfile subroutine Unpack_States


#addfile function Import_Matrix


#addfile function Export_Matrix


#addfile function Export_Tensor
end module Math_Operations
