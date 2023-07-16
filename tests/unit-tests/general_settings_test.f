! ------------------------------------------------------------------ ! --------------------------- !
module General_Settings_Test
   use General_Settings
   use fruit

   implicit none

   integer, parameter, dimension(3, 4) :: nel_direct_shear = reshape([ &
      3, 2, 1, &
      4, 3, 1, &
      5, 3, 2, &
      6, 3, 3], [3, 4])

   private
   public :: General_Settings_Testing


   contains


   ! --------------------------------------------------------------- !
   subroutine General_Settings_Testing()
   ! --------------------------------------------------------------- !
      call Estimate_Components_Test()
      call Check_Input_Dimensions_Test()
      call Is_Nan_Test()
      ! Excluding the following messaging functions
      ! - Write_Message()
      ! - Write_Warning()
      ! - Write_Error_And_Exit()
   end subroutine General_Settings_Testing


   ! --------------------------------------------------------------- !
   subroutine Estimate_Components_Test()
   ! --------------------------------------------------------------- !
      integer :: idx, nel, num_direct, num_shear
      integer, dimension(2) :: res, ref_comp
      character(len=6) :: id_num

      id_num = 'nel=  '
      do nel = -2, 9
         write(id_num(5:6), '(i2)') nel
         call Estimate_Components(nel=nel, num_dimensions=num_direct, num_shear=num_shear)
         res = [num_direct, num_shear]

         ref_comp = [0, 0]
         do idx = 1, size(nel_direct_shear, 2)
            if (nel == nel_direct_shear(1, idx)) then
               ref_comp = nel_direct_shear(2:3, idx)
               exit
            end if
         end do

         call assert_equals(var1=ref_comp, var2=res, n=2, message='Estimate_Components(' // id_num // ')')
      end do
   end subroutine Estimate_Components_Test


   ! --------------------------------------------------------------- !
   subroutine Check_Input_Dimensions_Test()
   ! --------------------------------------------------------------- !
      integer :: num_direct, num_shear, idx
      logical :: check_dim, check_valid
      character(len=31) :: id_num

      id_num = 'num_dimensions=  , num_shear=  '
      do num_direct = -2, 6
         write(id_num(16:17), '(i2)') num_direct
         do num_shear = -2, 6
            write(id_num(30:31), '(i2)') num_shear
            check_dim = Check_Input_Dimensions(num_dimensions=num_direct, num_shear=num_shear)

            check_valid = .False.
            do idx = 1, size(nel_direct_shear, 2)
               if ((num_direct == nel_direct_shear(2, idx)) .and. (num_shear == nel_direct_shear(3, idx))) then
                  check_valid = .True.
                  exit
               end if
            end do
            call assert_equals(var1=check_dim, var2=check_valid, message='Check_Input_Dimensions(' // id_num // ')')
         end do
      end do

      ! Reset to a valid state
      check_dim = Check_Input_Dimensions(num_dimensions=3, num_shear=3)
   end subroutine Check_Input_Dimensions_Test

   ! --------------------------------------------------------------- !
   subroutine Is_Nan_Test()
   ! --------------------------------------------------------------- !
      real(dp) :: x_null, x_nan, x_inf

      x_null = 0.0_dp
      x_nan = 0.0_dp/x_null
      x_inf = 1.0_dp/x_null

      call assert_equals(var1=.False., var2=Is_Nan(1.0_dp), message='Is_Nan(1.0_dp)')
      call assert_equals(var1=.False., var2=Is_Nan(-0.01_dp), message='Is_Nan(-0.01_dp)')
      call assert_equals(var1=.False., var2=Is_Nan(0.0_dp), message='Is_Nan(0.0_dp)')
      call assert_equals(var1=.False., var2=Is_Nan(x_inf), message='Is_Nan(1.0_dp/0.0_dp)')
      call assert_equals(var1=.False., var2=Is_Nan(0.0_dp/1.0_dp), message='Is_Nan(0.0_dp/1.0_dp)')
      call assert_equals(var1=.True., var2=Is_Nan(x_nan), message='Is_Nan(0.0_dp/0.0_dp)')
   end subroutine Is_Nan_Test
end module General_Settings_Test
