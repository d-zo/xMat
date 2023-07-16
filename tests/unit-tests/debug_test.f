! ------------------------------------------------------------------ ! --------------------------- !
module Debug_Test
   use General_Settings, only: dp
   use Debug
   use fruit

   implicit none

   private
   public :: Debug_Testing


   contains


   ! --------------------------------------------------------------- !
   subroutine Debug_Testing()
   ! --------------------------------------------------------------- !
      call Number_To_String_Test()
      call Format_logical_Test()
      call Format_dp_dim0_Test()
      ! Excluding the following debug functions
      ! - Format_dp_dim1()
      ! - Format_dp_dim2()
      ! - Format_dp_dim4()
   end subroutine Debug_Testing


   ! --------------------------------------------------------------- !
   subroutine Number_To_String_Test()
   ! --------------------------------------------------------------- !
      ! FIXME: Assuming onenumberlength=12. If this or any output format is changed, the code here has to change as well
      call assert_equals(var1='         0.0, ', var2=Number_To_String(number=0.0_dp), &
         message='Number_To_String(number=0.0_dp)')

      call assert_equals(var1='    12.00000, ', var2=Number_To_String(number=12.0_dp), &
         message='Number_To_String(number=12.0_dp)')

      call assert_equals(var1='    -0.03456, ', var2=Number_To_String(number=-0.03456_dp), &
         message='Number_To_String(number=-0.03456_dp)')

      call assert_equals(var1='   1.200E+03, ', var2=Number_To_String(number=1200.0_dp), &
         message='Number_To_String(number=1200.0_dp)')

      call assert_equals(var1='  -3.254E-06, ', var2=Number_To_String(number=-0.0000032544_dp), &
         message='Number_To_String(number=-0.0000032544_dp)')
   end subroutine Number_To_String_Test


   ! --------------------------------------------------------------- !
   subroutine Format_logical_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1='test:  Yes', var2=Format_logical(name='test', bool=.True.), &
         message='Format_logical(name=''test'', bool=.True.)')

      call assert_equals(var1=':  No ', var2=Format_logical(name='', bool=.False.), &
         message='Format_logical(name='''', bool=.False.)')

      call assert_equals(var1='Some very long name:  Yes', &
         var2=Format_logical(name='Some very long name', bool=.True.), &
         message='Format_logical(name=''Some very long name'', bool=.True.)')
   end subroutine Format_logical_Test


   ! --------------------------------------------------------------- !
   subroutine Format_dp_dim0_Test()
   ! --------------------------------------------------------------- !
      ! FIXME: Assuming onenumberlength=12. If this or any output format is changed, the code here has to change as well
      call assert_equals(var1='test           0.0', var2=Format_dp_dim0(name='test', number=0.0_dp), &
         message='Format_dp_dim0(name=''test'', number=0.0_dp)')

      call assert_equals(var1='      12.00000', var2=Format_dp_dim0(name='', number=12.0_dp), &
         message='Format_dp_dim0(name='''', number=12.0_dp)')

      call assert_equals(var1='Some very long name      -0.03456', &
         var2=Format_dp_dim0(name='Some very long name', number=-0.03456_dp), &
         message='Format_dp_dim0(name=''Some very long name'', number=-0.03456_dp)')

      call assert_equals(var1='short     1.200E+03', var2=Format_dp_dim0(name='short', number=1200.0_dp), &
         message='Format_dp_dim0(name='''', number=1200.0_dp)')

      call assert_equals(var1='special char' // char(10) // '    -3.254E-06', &
         var2=Format_dp_dim0(name='special char' // char(10), number=-0.0000032544_dp), &
         message='Format_dp_dim0(name=''special char'' // char(10), number=-0.0000032544_dp)')
   end subroutine Format_dp_dim0_Test
end module Debug_Test
