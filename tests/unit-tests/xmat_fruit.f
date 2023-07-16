program xmat_fruit
   use fruit
   use general_settings_test
   use debug_test
   use math_operations_test

   call init_fruit()

   call run_test_case(tc=General_Settings_Testing, tc_name='Testing "General Settings"')
   write(*, '(a)') ' <- General Settings'

   call run_test_case(tc=Debug_Testing, tc_name='Testing "Debug"')
   write(*, '(a)') ' <- Debug'

   call run_test_case(tc=Math_Operations_Testing, tc_name='Testing "Math Operations"')
   write(*, '(a)') ' <- Math Operations'

   call fruit_summary()
end program xmat_fruit 
