! ------------------------------------------------------------------ ! --------------------------- !
module Math_Operations_Test
   use General_Settings, only: dp, setting_epsilon
   use Math_Operations
   use fruit

   implicit none

   real(dp), parameter, dimension(3, 3, 3, 3) :: t_tens_id = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3, 3, 3])
   real(dp), parameter, dimension(3, 3, 3, 3) :: t_tens_tr = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3, 3, 3])

   real(dp), parameter, dimension(3, 3) :: t_mat33_rot01 = reshape([ &
      1.0_dp, 0.0_dp,        0.0_dp, &
      0.0_dp, sqrt(0.75_dp), 0.5_dp, &
      0.0_dp, -0.5_dp,       sqrt(0.75_dp)], [3, 3])
   real(dp), parameter, dimension(3, 3) :: t_mat33_rot02 = reshape([ &
      -0.5_dp,         0.0_dp,       sqrt(0.75_dp), &
      sqrt(0.375_dp),  sqrt(0.5_dp), sqrt(0.125_dp), &
      -sqrt(0.375_dp), sqrt(0.5_dp), -sqrt(0.125_dp)], [3, 3])

   real(dp), parameter, dimension(9) :: t_vec9_bsp01 = [ &
      7.0_dp, 9.0_dp, 5.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 1.0_dp, 1.0_dp]
   real(dp), parameter, dimension(9) :: t_vec9_bsp02 = [ &
      3.0_dp, -2.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, 5.0_dp, 5.0_dp, -6.0_dp, -6.0_dp]

   real(dp), parameter, dimension(9) :: t_vec9_bsp03 = [ &
      5.0_dp, 3.0_dp, 7.0_dp, 1.0_dp, 2.0_dp, 0.5_dp, 1.5_dp, -0.5_dp, -1.0_dp]

   real(dp), parameter, dimension(9, 9) :: t_mat99_bsp01 = reshape([ &
      100.0_dp,  500.0_dp,  900.0_dp,  2.0_dp,  4.0_dp,  6.0_dp,  8.0_dp,  3.0_dp,  7.0_dp,  &
      3700.0_dp, 4100.0_dp, 4500.0_dp, 38.0_dp, 40.0_dp, 42.0_dp, 44.0_dp, 39.0_dp, 43.0_dp, &
      7300.0_dp, 7700.0_dp, 8100.0_dp, 74.0_dp, 76.0_dp, 78.0_dp, 80.0_dp, 75.0_dp, 79.0_dp, &
      10.0_dp,   14.0_dp,   18.0_dp,   11.0_dp, 13.0_dp, 15.0_dp, 17.0_dp, 12.0_dp, 16.0_dp, &
      28.0_dp,   32.0_dp,   36.0_dp,   29.0_dp, 31.0_dp, 33.0_dp, 35.0_dp, 30.0_dp, 34.0_dp, &
      46.0_dp,   50.0_dp,   54.0_dp,   47.0_dp, 49.0_dp, 51.0_dp, 53.0_dp, 48.0_dp, 52.0_dp, &
      64.0_dp,   68.0_dp,   72.0_dp,   65.0_dp, 67.0_dp, 69.0_dp, 71.0_dp, 66.0_dp, 70.0_dp, &
      19.0_dp,   23.0_dp,   27.0_dp,   20.0_dp, 22.0_dp, 24.0_dp, 26.0_dp, 21.0_dp, 25.0_dp, &
      55.0_dp,   59.0_dp,   63.0_dp,   56.0_dp, 58.0_dp, 60.0_dp, 62.0_dp, 57.0_dp, 61.0_dp], [9, 9])

   real(dp), parameter, dimension(9, 9) :: t_mat99_bsp02 = reshape([ &
      10.0_dp, 40.0_dp, 60.0_dp, 0.1_dp, 0.1_dp, 0.2_dp, 0.2_dp, 0.3_dp, 0.3_dp, &
      40.0_dp, 20.0_dp, 50.0_dp, 0.4_dp, 0.4_dp, 0.5_dp, 0.5_dp, 0.6_dp, 0.6_dp, &
      60.0_dp, 50.0_dp, 30.0_dp, 0.7_dp, 0.7_dp, 0.8_dp, 0.8_dp, 0.9_dp, 0.9_dp, &
      0.1_dp,  0.4_dp,  0.7_dp,  1.0_dp, 1.0_dp, 4.0_dp, 4.0_dp, 6.0_dp, 6.0_dp, &
      0.1_dp,  0.4_dp,  0.7_dp,  1.0_dp, 1.0_dp, 4.0_dp, 4.0_dp, 6.0_dp, 6.0_dp, &
      0.2_dp,  0.5_dp,  0.8_dp,  4.0_dp, 4.0_dp, 2.0_dp, 2.0_dp, 5.0_dp, 5.0_dp, &
      0.2_dp,  0.5_dp,  0.8_dp,  4.0_dp, 4.0_dp, 2.0_dp, 2.0_dp, 5.0_dp, 5.0_dp, &
      0.3_dp,  0.6_dp,  0.9_dp,  6.0_dp, 6.0_dp, 5.0_dp, 5.0_dp, 3.0_dp, 3.0_dp, &
      0.3_dp,  0.6_dp,  0.9_dp,  6.0_dp, 6.0_dp, 5.0_dp, 5.0_dp, 3.0_dp, 3.0_dp], [9, 9])

   real(dp), parameter, dimension(3, 3) :: t_mat33_id = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_id = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_id = [1.0_dp, 1.0_dp, 1.0_dp, &
                                                     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_id = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_zero = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_zero = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_zero = [0.0_dp, 0.0_dp, 0.0_dp, &
                                                       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_zero = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_bsp01 = reshape([ &
      5.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 50.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_bsp01 = [5.0_dp, 0.5_dp, 50.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_bsp01 = [5.0_dp, 0.5_dp, 50.0_dp, &
                                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_bsp01 = reshape([ &
      5.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 50.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_bsp02 = reshape([ &
      -0.3_dp, 0.0_dp,  0.0_dp, &
      0.0_dp, -30.0_dp, 0.0_dp, &
      0.0_dp,  0.0_dp, -3.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_bsp02 = [-0.3_dp, -30.0_dp, -3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_bsp02 = [-0.3_dp, -30.0_dp, -3.0_dp, &
                                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_bsp02 = reshape([ &
      -0.3_dp, 0.0_dp,  0.0_dp, &
      0.0_dp, -30.0_dp, 0.0_dp, &
      0.0_dp,  0.0_dp, -3.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_bsp03 = reshape([ &                                       ! tridiagonal with eigenvalues 1, 5, 7
      5.0_dp, 2.0_dp, 0.0_dp, &
      2.0_dp, 3.0_dp, 2.0_dp, &
      0.0_dp, 2.0_dp, 5.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_bsp03 = [5.0_dp, 3.0_dp, 5.0_dp, 2.0_dp, 2.0_dp, 0.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_bsp03 = [5.0_dp, 3.0_dp, 5.0_dp, &
                                                        2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_bsp03 = reshape([ &
      5.0_dp, 2.0_dp, 0.0_dp, &
      2.0_dp, 3.0_dp, 2.0_dp, &
      0.0_dp, 2.0_dp, 5.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_bsp04 = reshape([ &                                       ! tridiagonal with eigenvalues 1, 2, 11
      2.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 3.0_dp, 4.0_dp, &
      0.0_dp, 4.0_dp, 9.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_bsp04 = [2.0_dp, 3.0_dp, 9.0_dp, 0.0_dp, 4.0_dp, 0.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_bsp04 = [2.0_dp, 3.0_dp, 9.0_dp, &
                                                        0.0_dp, 0.0_dp, 4.0_dp, 4.0_dp, 0.0_dp, 0.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_bsp04 = reshape([ &
      2.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 3.0_dp, 4.0_dp, &
      0.0_dp, 4.0_dp, 9.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_bsp05 = reshape([ &
      0.5_dp,  1.0_dp, 2.0_dp, &
      1.0_dp, -2.0_dp, 3.0_dp, &
      2.0_dp,  3.0_dp, 8.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_bsp05 = [0.5_dp, -2.0_dp, 8.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_bsp05 = [0.5_dp, -2.0_dp, 8.0_dp, &
                                                        1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 2.0_dp, 2.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_bsp05 = reshape([ &
      0.5_dp,  1.0_dp, 2.0_dp, &
      1.0_dp, -2.0_dp, 3.0_dp, &
      2.0_dp,  3.0_dp, 8.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_bsp06 = reshape([ &
      -3.0_dp, -2.0_dp, 1.0_dp, &
      1.0_dp,  3.0_dp,  1.0_dp, &
      1.0_dp,  -1.0_dp, 5.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_bsp06 = [-3.0_dp, 3.0_dp, 5.0_dp, -2.0_dp, 1.0_dp, 1.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_bsp06 = [-3.0_dp, 3.0_dp, 5.0_dp, &
                                                        -2.0_dp, -2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_bsp06 = reshape([ &
      -3.0_dp, -2.0_dp, 1.0_dp, &
      1.0_dp,  3.0_dp,  1.0_dp, &
      1.0_dp,  -1.0_dp, 5.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_bsp07 = reshape([ &                                       ! eigenvalues -1, 1, 7
      2.0_dp, -3.0_dp, -2.0_dp, &
     -3.0_dp,  2.0_dp,  2.0_dp, &
     -2.0_dp,  2.0_dp,  3.0_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_bsp07 = [2.0_dp, 2.0_dp, 3.0_dp, -3.0_dp, 2.0_dp, -2.0_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_bsp07 = [2.0_dp, 2.0_dp, 3.0_dp, &
                                                        -3.0_dp, -3.0_dp, 2.0_dp, 2.0_dp, -2.0_dp, -2.0_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_bsp07 = reshape([ &
      2.0_dp, -3.0_dp, -2.0_dp, &
     -3.0_dp,  2.0_dp,  2.0_dp, &
     -2.0_dp,  2.0_dp,  3.0_dp], [3, 3])
#endif
#endif
#endif

   real(dp), parameter, dimension(3, 3) :: t_mat33_bsp08 = reshape([ &                                       ! eigenvalues -0.7, -0.5, -0.2
     -0.5_dp,  0.2_dp,  0.1_dp, &
      0.2_dp, -0.5_dp,  0.1_dp, &
      0.1_dp,  0.1_dp, -0.4_dp], [3, 3])
#ifdef REPR_MAT66
   real(dp), parameter, dimension(6) :: t_repr_bsp08 = [-0.5_dp, -0.5_dp, -0.4_dp, 0.2_dp, 0.1_dp, 0.1_dp]
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9) :: t_repr_bsp08 = [-0.5_dp, -0.5_dp, -0.4_dp, &
                                                         0.2_dp, 0.2_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp]
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3) :: t_repr_bsp08 = reshape([ &
     -0.5_dp,  0.2_dp,  0.1_dp, &
      0.2_dp, -0.5_dp,  0.1_dp, &
      0.1_dp,  0.1_dp, -0.4_dp], [3, 3])
#endif
#endif
#endif

#ifdef REPR_MAT66
   real(dp), parameter, dimension(6, 6) :: t_repr_tens_id = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [6, 6])
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9, 9) :: t_repr_tens_id = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [9, 9])
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3, 3, 3) :: t_repr_tens_id = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3, 3, 3])
#endif
#endif
#endif

#ifdef REPR_MAT66
   real(dp), parameter, dimension(6, 6) :: t_repr_tens_tr = reshape([ &
      1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [6, 6])
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9, 9) :: t_repr_tens_tr = reshape([ &
      1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [9, 9])
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3, 3, 3) :: t_repr_tens_tr = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3, 3, 3])
#endif
#endif
#endif

#ifdef REPR_MAT66
   real(dp), parameter, dimension(6, 6) :: t_repr_bsp401 = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp], [6, 6])
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9, 9) :: t_repr_bsp401 = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp], [9, 9])
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3, 3, 3) :: t_repr_bsp401 = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3, 3, 3])
#endif
#endif
#endif

#ifdef REPR_MAT66
   real(dp), parameter, dimension(6, 6) :: t_repr_bsp404 = reshape([ &
      17690.0_dp, 17690.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      17690.0_dp, 17690.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      17690.0_dp, 17690.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp], [6, 6])
#else
#ifdef REPR_MAT99
   real(dp), parameter, dimension(9, 9) :: t_repr_bsp404 = reshape([ &
      17690.0_dp, 17690.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      17690.0_dp, 17690.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      17690.0_dp, 17690.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, &
          0.0_dp,     0.0_dp,     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp], [9, 9])
#else
#ifdef REPR_TENS3333
   real(dp), parameter, dimension(3, 3, 3, 3) :: t_repr_bsp404 = reshape([ &
      17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 17690.0_dp, &
          0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp,     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,     0.0_dp, &
          0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp,     0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp,     0.0_dp, &
          0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp,     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,     0.0_dp, &
      17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 17690.0_dp, &
          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,     0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp,     0.0_dp, &
          0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp,     0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp,     0.0_dp, &
          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,     0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp,     0.0_dp, &
      17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 17690.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 17690.0_dp], [3, 3, 3, 3])
#endif
#endif
#endif

!    Fruit assert_equals fuer Werte, Vektoren und Matrizen:
!     - Werte:          assert_equals(var1, var2, delta)
!     - Vektoren (n):   assert_equals(var1, var2, n, delta)
!     - Matrizen (n,m): assert_equalsvar1, var2, n, m, delta)

! Some hints about FRUIT on https://www.software.ac.uk/blog/2016-09-28-look-fortran-unit-test-frameworks


   private
   public :: Math_Operations_Testing


   contains


   ! --------------------------------------------------------------- !
   pure function Sort_Vector(vec) result(sorted)
   ! --------------------------------------------------------------- !
      real(dp), dimension(:), intent(in) :: vec
      real(dp), dimension(size(vec)) :: sorted
      ! ------------------------------------------------------------ ! Might not be the most efficient implementation
      integer :: idx, pos, nel

      nel = size(vec)
      sorted = vec

      do idx = 1, nel-1
         pos = idx - 1 + minloc(sorted(idx:nel), dim=1)
         sorted([idx, pos]) = sorted([pos, idx])
      end do
   end function Sort_Vector


   ! --------------------------------------------------------------- !
   pure function Zero_Third_Dimension(vec) result(vec_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(__matrix__), intent(in) :: vec
      real(dp), dimension(__matrix__) :: vec_out
      ! ------------------------------------------------------------ !
      vec_out = vec
#ifndef REPR_TENS3333
      vec_out(3) = 0.0_dp
#else
      vec_out(3, 3) = 0.0_dp
#endif
   end function Zero_Third_Dimension


   ! --------------------------------------------------------------- !
   pure function V9_To_M(vec9) result(mat)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(__matrix__) :: mat
      ! ------------------------------------------------------------ !
#ifdef REPR_MAT66
      mat(1:3) = vec9(1:3)
      mat(4) = 0.5_dp*(vec9(4) + vec9(5))
      mat(5) = 0.5_dp*(vec9(6) + vec9(7))
      mat(6) = 0.5_dp*(vec9(8) + vec9(9))
#else
#ifdef REPR_MAT99
      mat = vec9
#else
#ifdef REPR_TENS3333
      mat(1, 1) = vec9(1)
      mat(2, 2) = vec9(2)
      mat(3, 3) = vec9(3)
      mat(1, 2) = vec9(4)
      mat(2, 1) = vec9(5)
      mat(2, 3) = vec9(6)
      mat(3, 2) = vec9(7)
      mat(1, 3) = vec9(8)
      mat(3, 1) = vec9(9)
#endif
#endif
#endif
   end function V9_To_M


   ! --------------------------------------------------------------- !
   pure function M_To_V9(mat) result(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(__matrix__), intent(in) :: mat
      real(dp), dimension(9) :: vec9
      ! ------------------------------------------------------------ !
#ifdef REPR_MAT66
      vec9(1:3) = mat(1:3)
      vec9(4:5) = mat(4)
      vec9(6:7) = mat(5)
      vec9(8:9) = mat(6)
#else
#ifdef REPR_MAT99
      vec9 = mat
#else
#ifdef REPR_TENS3333
      vec9(1) = mat(1, 1)
      vec9(2) = mat(2, 2)
      vec9(3) = mat(3, 3)
      vec9(4) = mat(1, 2)
      vec9(5) = mat(2, 1)
      vec9(6) = mat(2, 3)
      vec9(7) = mat(3, 2)
      vec9(8) = mat(1, 3)
      vec9(9) = mat(3, 1)
#endif
#endif
#endif
   end function M_To_V9


#ifdef REPR_TENS3333
   ! --------------------------------------------------------------- !
   pure function Vec_To_Mat33(vec) result(mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(__matrix__), intent(in) :: vec
      real(dp), dimension(__matrix__) :: mat33
      ! ------------------------------------------------------------ !
      mat33 = vec
   end function Vec_To_Mat33


   ! --------------------------------------------------------------- !
   pure function Mat33_To_Vec(mat33) result(vec)
   ! --------------------------------------------------------------- !
      real(dp), dimension(__matrix__), intent(in) :: mat33
      real(dp), dimension(__matrix__) :: vec
      ! ------------------------------------------------------------ !
      vec = mat33
   end function Mat33_To_Vec
#endif


   ! --------------------------------------------------------------- !
   pure function T_To_M99(tens) result(mat99)
   ! --------------------------------------------------------------- !
      real(dp), dimension(__tensor__), intent(in) :: tens
      real(dp), dimension(9, 9) :: mat99
      ! ------------------------------------------------------------ !
#ifdef REPR_MAT66
      mat99(1:3, 1:3) = tens(1:3, 1:3)

      mat99(1:3, 4) = tens(1:3, 4)
      mat99(1:3, 5) = tens(1:3, 4)
      mat99(1:3, 6) = tens(1:3, 5)
      mat99(1:3, 7) = tens(1:3, 5)
      mat99(1:3, 8) = tens(1:3, 6)
      mat99(1:3, 9) = tens(1:3, 6)

      mat99(4, 1:3) = tens(4, 1:3)
      mat99(5, 1:3) = tens(4, 1:3)
      mat99(6, 1:3) = tens(5, 1:3)
      mat99(7, 1:3) = tens(5, 1:3)
      mat99(8, 1:3) = tens(6, 1:3)
      mat99(9, 1:3) = tens(6, 1:3)
      
      mat99(4:5, 4:5) = tens(4, 4)
      mat99(4:5, 6:7) = tens(4, 5)
      mat99(4:5, 8:9) = tens(4, 6)
      
      mat99(6:7, 4:5) = tens(5, 4)
      mat99(6:7, 6:7) = tens(5, 5)
      mat99(6:7, 8:9) = tens(5, 6)

      mat99(8:9, 4:5) = tens(6, 4)
      mat99(8:9, 6:7) = tens(6, 5)
      mat99(8:9, 8:9) = tens(6, 6)
#else
#ifdef REPR_MAT99
      mat99 = tens
#else
#ifdef REPR_TENS3333
      mat99(1, 1) = tens(1, 1, 1, 1)
      mat99(2, 2) = tens(2, 2, 2, 2)
      mat99(3, 3) = tens(3, 3, 3, 3)
      mat99(1, 2) = tens(1, 1, 2, 2)
      mat99(2, 3) = tens(2, 2, 3, 3)
      mat99(1, 3) = tens(1, 1, 3, 3)
      mat99(2, 1) = tens(2, 2, 1, 1)
      mat99(3, 2) = tens(3, 3, 2, 2)
      mat99(3, 1) = tens(3, 3, 1, 1)
      
      mat99(1, 4) = tens(1, 1, 1, 2)
      mat99(1, 5) = tens(1, 1, 2, 1)
      mat99(1, 6) = tens(1, 1, 2, 3)
      mat99(1, 7) = tens(1, 1, 3, 2)
      mat99(1, 8) = tens(1, 1, 1, 3)
      mat99(1, 9) = tens(1, 1, 3, 1)
      mat99(2, 4) = tens(2, 2, 1, 2)
      mat99(2, 5) = tens(2, 2, 2, 1)
      mat99(2, 6) = tens(2, 2, 2, 3)
      mat99(2, 7) = tens(2, 2, 3, 2)
      mat99(2, 8) = tens(2, 2, 1, 3)
      mat99(2, 9) = tens(2, 2, 3, 1)
      mat99(3, 4) = tens(3, 3, 1, 2)
      mat99(3, 5) = tens(3, 3, 2, 1)
      mat99(3, 6) = tens(3, 3, 2, 3)
      mat99(3, 7) = tens(3, 3, 3, 2)
      mat99(3, 8) = tens(3, 3, 1, 3)
      mat99(3, 9) = tens(3, 3, 3, 1)

      mat99(4, 1) = tens(1, 2, 1, 1)
      mat99(5, 1) = tens(2, 1, 1, 1)
      mat99(6, 1) = tens(2, 3, 1, 1)
      mat99(7, 1) = tens(3, 2, 1, 1)
      mat99(8, 1) = tens(1, 3, 1, 1)
      mat99(9, 1) = tens(3, 1, 1, 1)
      mat99(4, 2) = tens(1, 2, 2, 2)
      mat99(5, 2) = tens(2, 1, 2, 2)
      mat99(6, 2) = tens(2, 3, 2, 2)
      mat99(7, 2) = tens(3, 2, 2, 2)
      mat99(8, 2) = tens(1, 3, 2, 2)
      mat99(9, 2) = tens(3, 1, 2, 2)
      mat99(4, 3) = tens(1, 2, 3, 3)
      mat99(5, 3) = tens(2, 1, 3, 3)
      mat99(6, 3) = tens(2, 3, 3, 3)
      mat99(7, 3) = tens(3, 2, 3, 3)
      mat99(8, 3) = tens(1, 3, 3, 3)
      mat99(9, 3) = tens(3, 1, 3, 3)

      mat99(4, 4) = tens(1, 2, 1, 2)
      mat99(4, 5) = tens(1, 2, 2, 1)
      mat99(5, 4) = tens(2, 1, 1, 2)
      mat99(5, 5) = tens(2, 1, 2, 1)
      mat99(6, 6) = tens(2, 3, 2, 3)
      mat99(6, 7) = tens(2, 3, 3, 2)
      mat99(7, 6) = tens(3, 2, 2, 3)
      mat99(7, 7) = tens(3, 2, 3, 2)
      mat99(8, 8) = tens(1, 3, 1, 3)
      mat99(8, 9) = tens(1, 3, 3, 1)
      mat99(9, 8) = tens(3, 1, 1, 3)
      mat99(9, 9) = tens(3, 1, 3, 1)

      mat99(4, 6) = tens(1, 2, 2, 3)
      mat99(4, 7) = tens(1, 2, 3, 2)
      mat99(5, 6) = tens(2, 1, 2, 3)
      mat99(5, 7) = tens(2, 1, 3, 2)
      mat99(4, 8) = tens(1, 2, 1, 3)
      mat99(4, 9) = tens(1, 2, 3, 1)
      mat99(5, 8) = tens(2, 1, 1, 3)
      mat99(5, 9) = tens(2, 1, 3, 1)
      mat99(6, 8) = tens(2, 3, 1, 3)
      mat99(6, 9) = tens(2, 3, 3, 1)
      mat99(7, 8) = tens(3, 2, 1, 3)
      mat99(7, 9) = tens(3, 2, 3, 1)

      mat99(6, 4) = tens(2, 3, 1, 2)
      mat99(7, 4) = tens(3, 2, 1, 2)
      mat99(6, 5) = tens(2, 3, 2, 1)
      mat99(7, 5) = tens(3, 2, 2, 1)
      mat99(8, 4) = tens(1, 3, 1, 2)
      mat99(9, 4) = tens(3, 1, 1, 2)
      mat99(8, 5) = tens(1, 3, 2, 1)
      mat99(9, 5) = tens(3, 1, 2, 1)
      mat99(8, 6) = tens(1, 3, 2, 3)
      mat99(9, 6) = tens(3, 1, 2, 3)
      mat99(8, 7) = tens(1, 3, 3, 2)
      mat99(9, 7) = tens(3, 1, 3, 2)
#endif
#endif
#endif
   end function T_To_M99
   


   ! --------------------------------------------------------------- !
   subroutine Math_Operations_Testing()
   ! --------------------------------------------------------------- !
      call Reset_Dimension()

      ! Call all subroutines of this module except the utilits functions
      ! Set_Dimension() and Reset_Dimension()

      call Parameters_Test()
      call Nonzero_Division_Test()
      call Trace_Test()
      call Dimensionless_Test()
      call Deviatoric_Part_Test()
      call Norm_Test()
      call Determinant_Test()
      call Squared_Test()
      call LU_Test()
      call Inverse_Tensor_Test()
      call Inverse_Internal_Test()
      call Double_Contraction42_Test()
      call Double_Contraction44_Test()
      call Double_Contraction22_Test()
      call Dyadic_Product22_Test()
      call Tensor_Partialtrace_Test()
      call Matrix_Rotation_Test()
      call Sqrt_Of_Sum_Of_Squares_Test()
      call Tridiagonal_Matrix_Test()
      call Givens_CS_Test()
      call QL_Decomposition_Test()
      call Eigendecomposition6_Test()
      call Matrix_Exponential_Test()
      call Is_Diagonal_Test()
      call Value_In_Interval_Test()
      !call Abort_If_Not_In_Interval_Test()
      call Element_In_Tensor_Test()
      call Matrix_And_Tensorpart_Test()
      call Elements_And_Matrixlist_Test()
      call Ref_Index_Test()
      call Perturbate_Test()
      call Vec_And_Mat33_Test()
      call Vec9_And_Mat_Test()
      call Tens_And_Mat99_Test()
      call Packing_Test()
      !call Import_Export_Test()
   end subroutine Math_Operations_Testing


   ! --------------------------------------------------------------- !
   function Set_Dimension(num_direct, num_shear) result(check_dim)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Check_Input_Dimensions

      integer, intent(in) :: num_direct, num_shear
      logical :: check_dim

      check_dim = Check_Input_Dimensions(num_dimensions=num_direct, num_shear=num_shear)
   end function Set_Dimension


   ! --------------------------------------------------------------- !
   subroutine Reset_Dimension()
   ! --------------------------------------------------------------- !
      logical :: check_dim

      check_dim = Set_Dimension(num_direct=3, num_shear=3)
   end subroutine Reset_Dimension


   ! --------------------------------------------------------------- !
   subroutine Parameters_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=const_root2, var2=1.4142135623730951_dp, delta=setting_epsilon, message='sqrt(2)')
      call assert_equals(var1=const_root3, var2=1.7320508075688772_dp, delta=setting_epsilon, message='sqrt(3)')
      call assert_equals(var1=const_root6, var2=2.4494897427831779_dp, delta=setting_epsilon, message='sqrt(6)')

      call assert_equals(var1=Vec_To_Mat33(const_identity2d), var2=t_mat33_id, n=3, m=3, &
         delta=setting_epsilon, message='Parameters(const_identity2d)')
      call assert_equals(var1=T_To_M99(const_identity4d_sym), var2=T_To_M99(t_repr_bsp401), n=9, m=9, &
         delta=setting_epsilon, message='Parameters(const_identity4d_sym)')
      call assert_equals(var1=T_To_M99(const_identity4d_tr), var2=T_To_M99(t_repr_tens_tr), n=9, m=9, &
         delta=setting_epsilon, message='Parameters(const_identity4d_tr)')
   end subroutine Parameters_Test


   ! --------------------------------------------------------------- !
   subroutine Nonzero_Division_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(9) :: div_mat
      real(dp) :: div_val, fac

      fac = 1.0_dp
      div_mat = M_To_V9(Nonzero_Division(val=t_repr_bsp05, fac=fac))
      call assert_equals(var1=fac*M_To_V9(t_repr_bsp05), var2=div_mat, n=9, &
         delta=setting_epsilon, message='Nonzero_Division(val=t_repr_bsp05, fac=1.0_dp)')

      fac = -0.3_dp
      div_mat = M_To_V9(Nonzero_Division(val=t_repr_bsp05, fac=fac))
      call assert_equals(var1=M_To_V9(t_repr_bsp05)/fac, var2=div_mat, n=9, &
         delta=setting_epsilon, message='Nonzero_Division(val=t_repr_bsp05, fac=-0.3_dp)')

      fac = 0.0_dp
      div_mat = M_To_V9(Nonzero_Division(val=t_repr_bsp05, fac=fac))
      call assert_equals(var1=M_To_V9(t_repr_bsp05), var2=div_mat, n=9, &
         delta=setting_epsilon, message='Nonzero_Division(val=t_repr_bsp05, fac=0.0_dp)')

      div_val = -0.00345_dp
      fac = -15922.33_dp
      call assert_equals(var1=div_val/fac, var2=Nonzero_Division(val=div_val, fac=fac), &
         delta=setting_epsilon, message='Nonzero_Division(val=-0.00345_dp, fac=-15922.33_dp)')

      div_val = 0.0_dp
      fac = -15922.33_dp
      call assert_equals(var1=div_val/fac, var2=Nonzero_Division(val=div_val, fac=fac), &
         delta=setting_epsilon, message='Nonzero_Division(val=0.0_dp, fac=-15922.33_dp)')

      div_val = 0.0_dp
      fac = 0.0_dp
      call assert_equals(var1=div_val, var2=Nonzero_Division(val=div_val, fac=fac), &
         delta=setting_epsilon, message='Nonzero_Division(val=0.0_dp, fac=0.0_dp)')
   end subroutine Nonzero_Division_Test


   ! --------------------------------------------------------------- !
   subroutine Trace_Test()
   ! --------------------------------------------------------------- !
      real(dp) :: res
      real(dp), dimension(__matrix__) :: invec
      integer :: num_direct, num_shear
      character(len=14) :: id_num

      id_num = 'ndir= , nshr= '

      do num_direct = 2, 3
         write(id_num(6:6), '(i1)') num_direct
         do num_shear = 0, 3
            if (.not. Set_Dimension(num_direct=num_direct, num_shear=num_shear)) then
               cycle
            end if
            write(id_num(14:14), '(i1)') num_shear

            invec = t_repr_id
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
               res = sum([t_mat33_id(1, 1), t_mat33_id(2, 2)])
            else
               res = sum([t_mat33_id(1, 1), t_mat33_id(2, 2), t_mat33_id(3, 3)])
            end if
            call assert_equals(var1=res, var2=Trace(invec), &
               delta=setting_epsilon, message='Trace(t_repr_id; ' // id_num // ')')

            invec = t_repr_zero
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
               res = sum([t_mat33_zero(1, 1), t_mat33_zero(2, 2)])
            else
               res = sum([t_mat33_zero(1, 1), t_mat33_zero(2, 2), t_mat33_zero(3, 3)])
            end if
            call assert_equals(var1=res, var2=Trace(invec), &
               delta=setting_epsilon, message='Trace(t_repr_zero; ' // id_num // ')')

            invec = t_repr_bsp01
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
               res = sum([t_mat33_bsp01(1, 1), t_mat33_bsp01(2, 2)])
            else
               res = sum([t_mat33_bsp01(1, 1), t_mat33_bsp01(2, 2), t_mat33_bsp01(3, 3)])
            end if
            call assert_equals(var1=res, var2=Trace(invec), &
               delta=setting_epsilon, message='Trace(t_repr_bsp01; ' // id_num // ')')

            invec = t_repr_bsp02
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
               res = sum([t_mat33_bsp02(1, 1), t_mat33_bsp02(2, 2)])
            else
               res = sum([t_mat33_bsp02(1, 1), t_mat33_bsp02(2, 2), t_mat33_bsp02(3, 3)])
            end if
            call assert_equals(var1=res, var2=Trace(invec), &
               delta=setting_epsilon, message='Trace(t_repr_bsp02; ' // id_num // ')')

            invec = t_repr_bsp03
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
               res = sum([t_mat33_bsp03(1, 1), t_mat33_bsp03(2, 2)])
            else
               res = sum([t_mat33_bsp03(1, 1), t_mat33_bsp03(2, 2), t_mat33_bsp03(3, 3)])
            end if
            call assert_equals(var1=res, var2=Trace(invec), &
               delta=setting_epsilon, message='Trace(t_repr_bsp03; ' // id_num // ')')

            invec = t_repr_bsp04
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
               res = sum([t_mat33_bsp04(1, 1), t_mat33_bsp04(2, 2)])
            else
               res = sum([t_mat33_bsp04(1, 1), t_mat33_bsp04(2, 2), t_mat33_bsp04(3, 3)])
            end if
            call assert_equals(var1=res, var2=Trace(invec), &
               delta=setting_epsilon, message='Trace(t_repr_bsp04; ' // id_num // ')')

            invec = t_repr_bsp05
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
               res = sum([t_mat33_bsp05(1, 1), t_mat33_bsp05(2, 2)])
            else
               res = sum([t_mat33_bsp05(1, 1), t_mat33_bsp05(2, 2), t_mat33_bsp05(3, 3)])
            end if
            call assert_equals(var1=res, var2=Trace(invec), &
               delta=setting_epsilon, message='Trace(t_repr_bsp05; ' // id_num // ')')
         end do
      end do

      call Reset_Dimension()
   end subroutine Trace_Test


   ! --------------------------------------------------------------- !
   subroutine Dimensionless_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(9) :: dless

      dless = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call assert_equals(var1=dless, var2=M_To_V9(Dimensionless(t_repr_id)), n=9, &
         delta=setting_epsilon, message='Dimensionless(t_repr_id)')

      dless = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call assert_equals(var1=dless, var2=M_To_V9(Dimensionless(t_repr_zero)), n=9, &
         delta=setting_epsilon, message='Dimensionless(t_repr_zero)')

      dless = [5.0_dp/55.5_dp, 0.5_dp/55.5_dp, 50.0_dp/55.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call assert_equals(var1=dless, var2=M_To_V9(Dimensionless(t_repr_bsp01)), n=9, &
         delta=setting_epsilon, message='Dimensionless(t_repr_bsp01)')

      dless = [0.3_dp/33.3_dp, 30.0_dp/33.3_dp, 3.0_dp/33.3_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call assert_equals(var1=dless, var2=M_To_V9(Dimensionless(t_repr_bsp02)), n=9, &
         delta=setting_epsilon, message='Dimensionless(t_repr_bsp02)')

      dless = [5.0_dp/13.0_dp, 3.0_dp/13.0_dp, 5.0_dp/13.0_dp, &
               2.0_dp/13.0_dp, 2.0_dp/13.0_dp, 2.0_dp/13.0_dp, 2.0_dp/13.0_dp, 0.0_dp, 0.0_dp]
      call assert_equals(var1=dless, var2=M_To_V9(Dimensionless(t_repr_bsp03)), n=9, &
         delta=setting_epsilon, message='Dimensionless(t_repr_bsp03)')

      dless = [2.0_dp/14.0_dp, 3.0_dp/14.0_dp, 9.0_dp/14.0_dp, &
               0.0_dp, 0.0_dp, 4.0_dp/14.0_dp, 4.0_dp/14.0_dp, 0.0_dp, 0.0_dp]
      call assert_equals(var1=dless, var2=M_To_V9(Dimensionless(t_repr_bsp04)), n=9, &
         delta=setting_epsilon, message='Dimensionless(t_repr_bsp04)')

      dless = [0.5_dp/6.5_dp, -2.0_dp/6.5_dp, 8.0_dp/6.5_dp, &
               1.0_dp/6.5_dp,  1.0_dp/6.5_dp, 3.0_dp/6.5_dp, 3.0_dp/6.5_dp, 2.0_dp/6.5_dp, 2.0_dp/6.5_dp]
      call assert_equals(var1=dless, var2=M_To_V9(Dimensionless(t_repr_bsp05)), n=9, &
         delta=setting_epsilon, message='Dimensionless(t_repr_bsp05)')
   end subroutine Dimensionless_Test


   ! --------------------------------------------------------------- !
   subroutine Deviatoric_Part_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(9) :: dev, dev_res
      real(dp), dimension(__matrix__) :: invec
      integer :: num_direct, num_shear
      character(len=14) :: id_num

      id_num = 'ndir= , nshr= '

      do num_direct = 2, 3
         write(id_num(6:6), '(i1)') num_direct
         do num_shear = 0, 3
            if (.not. Set_Dimension(num_direct=num_direct, num_shear=num_shear)) then
               cycle
            end if
            write(id_num(14:14), '(i1)') num_shear

            invec = t_repr_id
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
            end if
            dev = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
            dev_res = M_To_V9(Deviatoric_Part(invec))
            call assert_equals(var1=dev, var2=dev_res, n=9, &
               delta=setting_epsilon, message='Deviatoric_Part(t_repr_id; ' // id_num // ')')

            invec = t_repr_zero
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
            end if
            dev = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
            dev_res = M_To_V9(Deviatoric_Part(invec))
            call assert_equals(var1=dev, var2=dev_res, n=9, &
               delta=setting_epsilon, message='Deviatoric_Part(t_repr_zero; ' // id_num // ')')

            invec = t_repr_bsp01
            if (num_direct == 2) then
               invec = Zero_Third_Dimension(invec)
               dev = [2.25_dp, -2.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
            else
               dev = [5.0_dp-55.5_dp/3.0_dp, 0.5_dp-55.5_dp/3.0_dp, 50.0_dp-55.5_dp/3.0_dp, &
                      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
            end if
            dev_res = M_To_V9(Deviatoric_Part(invec))
            call assert_equals(var1=dev, var2=dev_res, n=9, &
               delta=setting_epsilon, message='Deviatoric_Part(t_repr_bsp01; ' // id_num // ')')

            invec = t_repr_bsp02
            if (num_direct == 2) then
               dev = [-0.3_dp+15.15_dp, -30.0_dp+15.15_dp, &
                      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
               invec = Zero_Third_Dimension(invec)
            else
               dev = [-0.3_dp+33.3_dp/3.0_dp, -30.0_dp+33.3_dp/3.0_dp, -3.0_dp+33.3_dp/3.0_dp, &
                      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
            end if
            dev_res = M_To_V9(Deviatoric_Part(invec))
            call assert_equals(var1=dev, var2=dev_res, n=9, &
               delta=setting_epsilon, message='Deviatoric_Part(t_repr_bsp02; ' // id_num // ')')

            invec = t_repr_bsp03
            if (num_direct == 2) then
               dev = [1.0_dp, -1.0_dp, 0.0_dp, 2.0_dp,  2.0_dp, 2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
               invec = Zero_Third_Dimension(invec)
            else
               dev = [5.0_dp-13.0_dp/3.0_dp, 3.0_dp-13.0_dp/3.0_dp, 5.0_dp-13.0_dp/3.0_dp, &
                      2.0_dp,  2.0_dp, 2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
            end if
            dev_res = M_To_V9(Deviatoric_Part(invec))
            call assert_equals(var1=dev, var2=dev_res, n=9, &
               delta=setting_epsilon, message='Deviatoric_Part(t_repr_bsp03; ' // id_num // ')')

            invec = t_repr_bsp04
            if (num_direct == 2) then
               dev = [-0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 4.0_dp, 0.0_dp, 0.0_dp]
               invec = Zero_Third_Dimension(invec)
            else
               dev = [2.0_dp-14.0_dp/3.0_dp, 3.0_dp-14.0_dp/3.0_dp, 9.0_dp-14.0_dp/3.0_dp, &
                      0.0_dp, 0.0_dp, 4.0_dp, 4.0_dp, 0.0_dp, 0.0_dp]
            end if
            dev_res = M_To_V9(Deviatoric_Part(invec))
            call assert_equals(var1=dev, var2=dev_res, n=9, &
               delta=setting_epsilon, message='Deviatoric_Part(t_repr_bsp04; ' // id_num // ')')

            invec = t_repr_bsp05
            if (num_direct == 2) then
               dev = [1.25_dp, -1.25_dp, 0.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 2.0_dp, 2.0_dp]
               invec = Zero_Third_Dimension(invec)
            else
               dev = [0.5_dp-6.5_dp/3.0_dp, -2.0_dp-6.5_dp/3.0_dp, 8.0_dp-6.5_dp/3.0_dp, &
                      1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 2.0_dp, 2.0_dp]
            end if
            dev_res = M_To_V9(Deviatoric_Part(invec))
            call assert_equals(var1=dev, var2=dev_res, n=9, &
               delta=setting_epsilon, message='Deviatoric_Part(t_repr_bsp05; ' // id_num // ')')
         end do
      end do

      call Reset_Dimension()
   end subroutine Deviatoric_Part_Test


   ! --------------------------------------------------------------- !
   subroutine Norm_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=sqrt(3.0_dp), var2=Norm(t_repr_id), &
         delta=setting_epsilon, message='Norm(t_repr_id)')

      call assert_equals(var1=0.0_dp, var2=Norm(t_repr_zero), &
         delta=setting_epsilon, message='Norm(t_repr_zero)')

      call assert_equals(var1=sqrt(2525.25_dp), var2=Norm(t_repr_bsp01), &
         delta=setting_epsilon, message='Norm(t_repr_bsp01)')

      call assert_equals(var1=sqrt(909.09_dp), var2=Norm(t_repr_bsp02), &
         delta=setting_epsilon, message='Norm(t_repr_bsp02)')

      call assert_equals(var1=sqrt(75.0_dp), var2=Norm(t_repr_bsp03), &
         delta=setting_epsilon, message='Norm(t_repr_bsp03)')

      call assert_equals(var1=sqrt(126.0_dp), var2=Norm(t_repr_bsp04), &
         delta=setting_epsilon, message='Norm(t_repr_bsp04)')

      call assert_equals(var1=sqrt(96.25_dp), var2=Norm(t_repr_bsp05), &
         delta=setting_epsilon, message='Norm(t_repr_bsp05)')
   end subroutine Norm_Test


   ! --------------------------------------------------------------- !
   subroutine Determinant_Test()
   ! --------------------------------------------------------------- !
      real(dp) :: det

      call assert_equals(var1=1.0_dp, var2=Determinant(t_repr_id), &
         delta=setting_epsilon, message='Determinant(t_repr_id)')

      call assert_equals(var1=0.0_dp, var2=Determinant(t_repr_zero), &
         delta=setting_epsilon, message='Determinant(t_repr_zero)')

      det = t_mat33_bsp01(1, 1)*t_mat33_bsp01(2, 2)*t_mat33_bsp01(3, 3) &
          + t_mat33_bsp01(1, 2)*t_mat33_bsp01(2, 3)*t_mat33_bsp01(3, 1) &
          + t_mat33_bsp01(1, 3)*t_mat33_bsp01(2, 1)*t_mat33_bsp01(3, 2) &
          - t_mat33_bsp01(1, 1)*t_mat33_bsp01(3, 2)*t_mat33_bsp01(2, 3) &
          - t_mat33_bsp01(1, 2)*t_mat33_bsp01(3, 3)*t_mat33_bsp01(2, 1) &
          - t_mat33_bsp01(1, 3)*t_mat33_bsp01(3, 1)*t_mat33_bsp01(2, 2)
      call assert_equals(var1=det, var2=Determinant(t_repr_bsp01), &
         delta=setting_epsilon, message='Determinant(t_repr_bsp01)')

      det = t_mat33_bsp02(1, 1)*t_mat33_bsp02(2, 2)*t_mat33_bsp02(3, 3) &
          + t_mat33_bsp02(1, 2)*t_mat33_bsp02(2, 3)*t_mat33_bsp02(3, 1) &
          + t_mat33_bsp02(1, 3)*t_mat33_bsp02(2, 1)*t_mat33_bsp02(3, 2) &
          - t_mat33_bsp02(1, 1)*t_mat33_bsp02(3, 2)*t_mat33_bsp02(2, 3) &
          - t_mat33_bsp02(1, 2)*t_mat33_bsp02(3, 3)*t_mat33_bsp02(2, 1) &
          - t_mat33_bsp02(1, 3)*t_mat33_bsp02(3, 1)*t_mat33_bsp02(2, 2)
      call assert_equals(var1=det, var2=Determinant(t_repr_bsp02), &
         delta=setting_epsilon, message='Determinant(t_repr_bsp02)')

      det = t_mat33_bsp03(1, 1)*t_mat33_bsp03(2, 2)*t_mat33_bsp03(3, 3) &
          + t_mat33_bsp03(1, 2)*t_mat33_bsp03(2, 3)*t_mat33_bsp03(3, 1) &
          + t_mat33_bsp03(1, 3)*t_mat33_bsp03(2, 1)*t_mat33_bsp03(3, 2) &
          - t_mat33_bsp03(1, 1)*t_mat33_bsp03(3, 2)*t_mat33_bsp03(2, 3) &
          - t_mat33_bsp03(1, 2)*t_mat33_bsp03(3, 3)*t_mat33_bsp03(2, 1) &
          - t_mat33_bsp03(1, 3)*t_mat33_bsp03(3, 1)*t_mat33_bsp03(2, 2)
      call assert_equals(var1=det, var2=Determinant(t_repr_bsp03), &
         delta=setting_epsilon, message='Determinant(t_repr_bsp03)')

      det = t_mat33_bsp04(1, 1)*t_mat33_bsp04(2, 2)*t_mat33_bsp04(3, 3) &
          + t_mat33_bsp04(1, 2)*t_mat33_bsp04(2, 3)*t_mat33_bsp04(3, 1) &
          + t_mat33_bsp04(1, 3)*t_mat33_bsp04(2, 1)*t_mat33_bsp04(3, 2) &
          - t_mat33_bsp04(1, 1)*t_mat33_bsp04(3, 2)*t_mat33_bsp04(2, 3) &
          - t_mat33_bsp04(1, 2)*t_mat33_bsp04(3, 3)*t_mat33_bsp04(2, 1) &
          - t_mat33_bsp04(1, 3)*t_mat33_bsp04(3, 1)*t_mat33_bsp04(2, 2)
      call assert_equals(var1=det, var2=Determinant(t_repr_bsp04), &
         delta=setting_epsilon, message='Determinant(t_repr_bsp04)')

      det = t_mat33_bsp05(1, 1)*t_mat33_bsp05(2, 2)*t_mat33_bsp05(3, 3) &
          + t_mat33_bsp05(1, 2)*t_mat33_bsp05(2, 3)*t_mat33_bsp05(3, 1) &
          + t_mat33_bsp05(1, 3)*t_mat33_bsp05(2, 1)*t_mat33_bsp05(3, 2) &
          - t_mat33_bsp05(1, 1)*t_mat33_bsp05(3, 2)*t_mat33_bsp05(2, 3) &
          - t_mat33_bsp05(1, 2)*t_mat33_bsp05(3, 3)*t_mat33_bsp05(2, 1) &
          - t_mat33_bsp05(1, 3)*t_mat33_bsp05(3, 1)*t_mat33_bsp05(2, 2)
      call assert_equals(var1=det, var2=Determinant(t_repr_bsp05), &
         delta=setting_epsilon, message='Determinant(t_repr_bsp05)')
   end subroutine Determinant_Test


   ! --------------------------------------------------------------- !
   subroutine Squared_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3) :: res

      res = matmul(t_mat33_id, t_mat33_id)
      call assert_equals(var1=res, var2=Vec_To_Mat33(Squared(t_repr_id)), n=3, m=3, &
         delta=setting_epsilon, message='Squared(t_repr_id)')

      res = matmul(t_mat33_zero, t_mat33_zero)
      call assert_equals(var1=res, var2=Vec_To_Mat33(Squared(t_repr_zero)), n=3, m=3, &
         delta=setting_epsilon, message='Squared(t_repr_zero)')

      res = matmul(t_mat33_bsp01, t_mat33_bsp01)
      call assert_equals(var1=res, var2=Vec_To_Mat33(Squared(t_repr_bsp01)), n=3, m=3, &
         delta=setting_epsilon, message='Squared(t_repr_bsp01)')

      res = matmul(t_mat33_bsp02, t_mat33_bsp02)
      call assert_equals(var1=res, var2=Vec_To_Mat33(Squared(t_repr_bsp02)), n=3, m=3, &
         delta=setting_epsilon, message='Squared(t_repr_bsp02)')

      res = matmul(t_mat33_bsp03, t_mat33_bsp03)
      call assert_equals(var1=res, var2=Vec_To_Mat33(Squared(t_repr_bsp03)), n=3, m=3, &
         delta=setting_epsilon, message='Squared(t_repr_bsp03)')

      res = matmul(t_mat33_bsp04, t_mat33_bsp04)
      call assert_equals(var1=res, var2=Vec_To_Mat33(Squared(t_repr_bsp04)), n=3, m=3, &
         delta=setting_epsilon, message='Squared(t_repr_bsp04)')

      res = matmul(t_mat33_bsp05, t_mat33_bsp05)
      call assert_equals(var1=res, var2=Vec_To_Mat33(Squared(t_repr_bsp05)), n=3, m=3, &
         delta=setting_epsilon, message='Squared(t_repr_bsp05)')
   end subroutine Squared_Test


   ! --------------------------------------------------------------- !
   subroutine LU_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3) :: matrix, p_mat, temp
      real(dp), dimension(3) :: b_vec
      logical :: is_singular

      matrix = t_mat33_id
      b_vec = [0.0_dp, 1.0_dp, 0.0_dp];
      call LU_Decomposition(matrix=matrix, nel=3, trans_mat=p_mat, is_singular=is_singular)
      call assert_equals(var1=t_mat33_id, var2=matrix, n=3, m=3, &
         delta=setting_epsilon, message='LU_Decomposition(t_mat33_id)')
      call assert_equals(var1=.False., var2=is_singular, &
         message='is_singular LU_Decomposition(t_mat33_id)?')
      if (.not. is_singular) then
         call LU_Solve(lu_mat=matrix, nel=3, trans_mat=p_mat, b_vec=b_vec)
         call assert_equals(var1=[0.0_dp, 1.0_dp, 0.0_dp], var2=b_vec, n=3, &
            delta=setting_epsilon, message='LU_Solve(t_mat33_id, [0.0_dp, 1.0_dp, 0.0_dp])')
      end if

      matrix = t_mat33_bsp06
      b_vec = [-4.0_dp, -9.0_dp, 48.0_dp];
      temp = reshape([-3.0_dp, 2.0_dp/3.0_dp,  -1.0_dp/3.0_dp, &
                      1.0_dp,  7.0_dp/3.0_dp,  4.0_dp/7.0_dp, &
                      1.0_dp,  -5.0_dp/3.0_dp, 44.0_dp/7.0_dp], [3, 3])
      call LU_Decomposition(matrix=matrix, nel=3, trans_mat=p_mat, is_singular=is_singular)
      call assert_equals(var1=temp, var2=matrix, n=3, m=3, &
         delta=setting_epsilon, message='LU_Decomposition(t_mat33_bsp06)')
      call assert_equals(var1=.False., var2=is_singular, &
         message='is_singular LU_Decomposition(t_mat33_bsp06)?')
      if (.not. is_singular) then
         call LU_Solve(lu_mat=matrix, nel=3, trans_mat=p_mat, b_vec=b_vec)
         call assert_equals(var1=[5.0_dp, 3.0_dp, 8.0_dp], var2=b_vec, n=3, &
            delta=setting_epsilon, message='LU_Solve(t_mat33_bsp06, [-4.0_dp, -9.0_dp, 48.0_dp])')
      end if

      ! Add more test cases (different nel?)
   end subroutine LU_Test


   ! --------------------------------------------------------------- !
   subroutine Inverse_Tensor_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(__tensor__) :: res
      logical :: success

      call Inverse_Tensor(tensor=t_repr_bsp401, inv_tensor=res, successful_inversion=success)
      call assert_equals(var1=T_To_M99(t_repr_bsp401), var2=T_To_M99(res), n=9, m=9, &
         delta=setting_epsilon, message='Inverse_Tensor(t_repr_bsp401)')
      call assert_equals(var1=.True., var2=success, &
         message='success Inverse_Tensor(t_repr_bsp401)?')

      call Inverse_Tensor(tensor=t_repr_bsp404, inv_tensor=res, successful_inversion=success)
      call assert_equals(var1=.False., var2=success, &
         message='success Inverse_Tensor(t_repr_bsp404)?')

      ! Add more test cases
   end subroutine Inverse_Tensor_Test


   ! --------------------------------------------------------------- !
   subroutine Inverse_Internal_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3) ::  inv_mat
      logical :: is_singular

      call Inverse_Internal(matrix=t_mat33_id, inv_matrix=inv_mat, nel=3, is_singular=is_singular)
      call assert_equals(var1=t_mat33_id, var2=inv_mat, n=3, m=3, &
         delta=setting_epsilon, message='Inverse_Internal(t_mat33_id)')
      call assert_equals(var1=.False., var2=is_singular, &
         message='is_singular Inverse_Internal(t_mat33_id)?')

      ! Add more test cases (different nel?)
   end subroutine Inverse_Internal_Test


   ! --------------------------------------------------------------- !
   subroutine Double_Contraction42_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(__tensor__) :: res

      call assert_equals(var1=M_To_V9(t_repr_id), &
         var2=M_To_V9(Double_Contraction42(t_repr_bsp401, t_repr_id)), n=9, &
         delta=setting_epsilon, message='Double_Contraction42(t_repr_bsp401, t_repr_id)')

      call assert_equals(var1=M_To_V9(t_repr_bsp04), &
         var2=M_To_V9(Double_Contraction42(t_repr_bsp401, t_repr_bsp04)), n=9, &
         delta=setting_epsilon, message='Double_Contraction42(t_repr_bsp401, t_repr_bsp04)')

      ! Add more test cases
   end subroutine Double_Contraction42_Test


   ! --------------------------------------------------------------- !
   subroutine Double_Contraction44_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(__tensor__) :: res

      call assert_equals(var1=T_To_M99(t_repr_bsp401), &
         var2=T_To_M99(Double_Contraction44(t_repr_bsp401, t_repr_bsp401)), n=9, m=9, &
         delta=setting_epsilon, message='Double_Contraction44(t_repr_bsp401, t_repr_bsp401)')

      ! Add more test cases
   end subroutine Double_Contraction44_Test

   
   ! --------------------------------------------------------------- !
   subroutine Double_Contraction22_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=3.0_dp, var2=Double_Contraction22(t_repr_id, t_repr_id), &
         delta=setting_epsilon, message='Double_Contraction22(t_repr_id, t_repr_id)')

      call assert_equals(var1=0.0_dp, var2=Double_Contraction22(t_repr_zero, t_repr_zero), &
         delta=setting_epsilon, message='Double_Contraction22(t_repr_zero, t_repr_zero)')

      call assert_equals(var1=sum(t_mat33_id*t_mat33_bsp01), &
         var2=Double_Contraction22(t_repr_id, t_repr_bsp01), &
         delta=setting_epsilon, message='Double_Contraction22(t_repr_id, t_repr_bsp01)')

      call assert_equals(var1=sum(t_mat33_bsp01*t_mat33_bsp01), &
         var2=Double_Contraction22(t_repr_bsp01, t_repr_bsp01), &
         delta=setting_epsilon, message='Double_Contraction22(t_repr_bsp01, t_repr_bsp01)')

      call assert_equals(var1=sum(t_mat33_bsp02*t_mat33_bsp03), &
         var2=Double_Contraction22(t_repr_bsp02, t_repr_bsp03), &
         delta=setting_epsilon, message='Double_Contraction22(t_repr_bsp02, t_repr_bsp03)')

      call assert_equals(var1=0.0_dp, var2=Double_Contraction22(t_repr_zero, t_repr_bsp04), &
         delta=setting_epsilon, message='Double_Contraction22(t_repr_zero, t_repr_bsp04)')
   end subroutine Double_Contraction22_Test


   ! --------------------------------------------------------------- !
   subroutine Dyadic_Product22_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=T_To_M99(t_repr_tens_tr), &
         var2=T_To_M99(Dyadic_Product22(t_repr_id, t_repr_id)), n=9, m=9, &
         delta=setting_epsilon, message='Dyadic_Product22(t_repr_id, t_repr_id)')

      ! Add more test cases
   end subroutine Dyadic_Product22_Test


   ! --------------------------------------------------------------- !
   subroutine Tensor_Partialtrace_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=3.0_dp, &
         var2=Tensor_Partialtrace(t_repr_tens_id), delta=setting_epsilon, &
         message='Tensor_Partialtrace(t_repr_tens_id)')

      call assert_equals(var1=-3.0_dp*17690.0_dp, &
         var2=Tensor_Partialtrace(-t_repr_bsp404), delta=setting_epsilon, &
         message='Tensor_Partialtrace(-t_repr_bsp404)')

      ! Add more test cases
   end subroutine Tensor_Partialtrace_Test


   ! --------------------------------------------------------------- !
   subroutine Matrix_Rotation_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3) :: res
      real(dp), dimension(__matrix__) :: res_rot

      res_rot = Matrix_Rotation(t_repr_id, rot_mat33=t_mat33_id)
      call assert_equals(var1=t_mat33_id, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_id, t_mat33_id)')

      res = matmul(t_mat33_rot01, matmul(t_mat33_id, transpose(t_mat33_rot01)))
      res_rot = Matrix_Rotation(t_repr_id, rot_mat33=t_mat33_rot01)
      call assert_equals(var1=t_mat33_id, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_id, t_mat33_rot01)')

      res = matmul(t_mat33_rot02, matmul(t_mat33_id, transpose(t_mat33_rot02)))
      res_rot = Matrix_Rotation(t_repr_id, rot_mat33=t_mat33_rot02)
      call assert_equals(var1=t_mat33_id, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_id, t_mat33_rot02)')

      res = matmul(t_mat33_rot01, matmul(t_mat33_bsp01, transpose(t_mat33_rot01)))
      res_rot = Matrix_Rotation(t_repr_bsp01, rot_mat33=t_mat33_rot01)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp01, t_mat33_rot01)')

      res = matmul(t_mat33_rot02, matmul(t_mat33_bsp01, transpose(t_mat33_rot02)))
      res_rot = Matrix_Rotation(t_repr_bsp01, rot_mat33=t_mat33_rot02)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp01, t_mat33_rot02)')

      res = matmul(t_mat33_rot01, matmul(t_mat33_bsp02, transpose(t_mat33_rot01)))
      res_rot = Matrix_Rotation(t_repr_bsp02, rot_mat33=t_mat33_rot01)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp02, t_mat33_rot01)')

      res = matmul(t_mat33_rot02, matmul(t_mat33_bsp02, transpose(t_mat33_rot02)))
      res_rot = Matrix_Rotation(t_repr_bsp02, rot_mat33=t_mat33_rot02)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp02, t_mat33_rot02)')

      res = matmul(t_mat33_rot01, matmul(t_mat33_bsp03, transpose(t_mat33_rot01)))
      res_rot = Matrix_Rotation(t_repr_bsp03, rot_mat33=t_mat33_rot01)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp03, t_mat33_rot01)')

      res = matmul(t_mat33_rot02, matmul(t_mat33_bsp03, transpose(t_mat33_rot02)))
      res_rot = Matrix_Rotation(t_repr_bsp03, rot_mat33=t_mat33_rot02)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp03, t_mat33_rot02)')

      res = matmul(t_mat33_rot01, matmul(t_mat33_bsp04, transpose(t_mat33_rot01)))
      res_rot = Matrix_Rotation(t_repr_bsp04, rot_mat33=t_mat33_rot01)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp04, t_mat33_rot01)')

      res = matmul(t_mat33_rot02, matmul(t_mat33_bsp04, transpose(t_mat33_rot02)))
      res_rot = Matrix_Rotation(t_repr_bsp04, rot_mat33=t_mat33_rot02)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp04, t_mat33_rot02)')

      res = matmul(t_mat33_rot01, matmul(t_mat33_bsp05, transpose(t_mat33_rot01)))
      res_rot = Matrix_Rotation(t_repr_bsp05, rot_mat33=t_mat33_rot01)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp05, t_mat33_rot01)')

      res = matmul(t_mat33_rot02, matmul(t_mat33_bsp05, transpose(t_mat33_rot02)))
      res_rot = Matrix_Rotation(t_repr_bsp05, rot_mat33=t_mat33_rot02)
      call assert_equals(var1=res, var2=Vec_To_Mat33(res_rot), n=3, m=3, &
         delta=setting_epsilon, message='Matrix_Rotation(t_repr_bsp05, t_mat33_rot02)')
   end subroutine Matrix_Rotation_Test


   ! --------------------------------------------------------------- !
   subroutine Sqrt_Of_Sum_Of_Squares_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=0.0_dp, var2=Sqrt_Of_Sum_Of_Squares(a=0.0_dp, b=0.0_dp), &
         delta=setting_epsilon, message='Sqrt_Of_Sum_Of_Squares(a=0.0_dp, b=0.0_dp)')
      call assert_equals(var1=5.0_dp, var2=Sqrt_Of_Sum_Of_Squares(a=4.0_dp, b=3.0_dp), &
         delta=setting_epsilon, message='Sqrt_Of_Sum_Of_Squares(a=4.0_dp, b=3.0_dp)')
      call assert_equals(var1=10.00000005_dp, var2=Sqrt_Of_Sum_Of_Squares(a=10.0_dp, b=-0.001_dp), &
         delta=setting_epsilon, message='Sqrt_Of_Sum_Of_Squares(a=10.0_dp, b=-0.001_dp)')
      call assert_equals(var1=40994768077.32698822021484_dp, &
         var2=Sqrt_Of_Sum_Of_Squares(a=-39410224789.0_dp, b=11287390832.0_dp), &
         delta=setting_epsilon, message='Sqrt_Of_Sum_Of_Squares(a=-39410224789.0_dp, b=11287390832.0_dp)')
   end subroutine Sqrt_Of_Sum_Of_Squares_Test


   ! --------------------------------------------------------------- !
   subroutine Tridiagonal_Matrix_Test()
   ! --------------------------------------------------------------- !
      use Debug, only: Formatval
      real(dp), dimension(3, 3) :: matrix, trans_mat
      real(dp), dimension(3, 3) :: res_mat, res_transmat

      matrix = t_mat33_id
      call Tridiagonal_Matrix(matrix=matrix, nel=3, trans_mat=trans_mat)
      res_mat = matmul(trans_mat, matmul(matrix, transpose(trans_mat)))
      call assert_equals(var1=t_mat33_id, var2=res_mat, n=3, m=3, &
         delta=setting_epsilon, message='Tridiagonal_Matrix(t_mat33_id)')
      call assert_equals(var1=[0.0_dp, 0.0_dp], var2=[matrix(3, 1), matrix(1, 3)], n=2, &
         delta=setting_epsilon, message='Tridiagonal_Matrix(t_mat33_id) zero_check')

      matrix = t_mat33_bsp03
      call Tridiagonal_Matrix(matrix=matrix, nel=3, trans_mat=trans_mat)
      res_mat = matmul(trans_mat, matmul(matrix, transpose(trans_mat)))
      call assert_equals(var1=t_mat33_bsp03, var2=res_mat, n=3, m=3, &
         delta=setting_epsilon, message='Tridiagonal_Matrix(t_mat33_bsp03)')
      call assert_equals(var1=[0.0_dp, 0.0_dp], var2=[matrix(3, 1), matrix(1, 3)], n=2, &
         delta=setting_epsilon, message='Tridiagonal_Matrix(t_mat33_bsp03) zero_check')

      matrix = t_mat33_bsp04
      call Tridiagonal_Matrix(matrix=matrix, nel=3, trans_mat=trans_mat)
      res_mat = matmul(trans_mat, matmul(matrix, transpose(trans_mat)))
      call assert_equals(var1=t_mat33_bsp04, var2=res_mat, n=3, m=3, &
         delta=setting_epsilon, message='Tridiagonal_Matrix(t_mat33_bsp04)')
      call assert_equals(var1=[0.0_dp, 0.0_dp], var2=[matrix(3, 1), matrix(1, 3)], n=2, &
         delta=setting_epsilon, message='Tridiagonal_Matrix(t_mat33_bsp04) zero_check')

      matrix = t_mat33_bsp05
      call Tridiagonal_Matrix(matrix=matrix, nel=3, trans_mat=trans_mat)
      res_mat = matmul(trans_mat, matmul(matrix, transpose(trans_mat)))
      call assert_equals(var1=t_mat33_bsp05, var2=res_mat, n=3, m=3, &
         delta=setting_epsilon, message='Tridiagonal_Matrix(t_mat33_bsp05)')
      call assert_equals(var1=[0.0_dp, 0.0_dp], var2=[matrix(3, 1), matrix(1, 3)], n=2, &
         delta=setting_epsilon, message='Tridiagonal_Matrix(t_mat33_bsp05) zero_check')

      ! Add more test cases (different nel?)
   end subroutine Tridiagonal_Matrix_Test


   ! --------------------------------------------------------------- !
   subroutine Givens_CS_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=t_mat33_id(1:2, 1:2), var2=Givens_CS(refval=1.0_dp, zeroval=0.0_dp), &
         n=2, m=2, delta=setting_epsilon, message='Givens_CS(refval=1.0_dp, zeroval=0.0_dp)')
      call assert_equals(var1=t_mat33_rot01(2:3, 2:3), &
         var2=Givens_CS(refval=sqrt(0.75_dp), zeroval=-0.5_dp), n=2, m=2, &
         delta=setting_epsilon, message='Givens_CS(refval=1.0_dp, zeroval=0.0_dp)')

      ! Add more test cases
   end subroutine Givens_CS_Test


   ! --------------------------------------------------------------- !
   subroutine QL_Decomposition_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3) :: matrix, eigenvec_mat, trans_mat
      real(dp), dimension(3) :: eigenval, res_eigenval
      logical :: ql_success

      ! NOTE: Sort_Vector() has to be called since the returned eigenvalues are not sorted

      matrix = t_mat33_id
      call QL_Decomposition(matrix=matrix, nel=3, eigenvec_mat=eigenvec_mat, ql_success=ql_success)
      call assert_equals(var1=.True., var2=ql_success, message='QL_Decomposition(t_mat33_id) successful')
      call assert_equals(var1=t_mat33_id, var2=matrix, n=3, m=3, &
         delta=setting_epsilon, message='QL_Decomposition(t_mat33_id)')

      matrix = t_mat33_bsp03
      eigenval = [1.0_dp, 5.0_dp, 7.0_dp]
      call QL_Decomposition(matrix=matrix, nel=3, eigenvec_mat=eigenvec_mat, ql_success=ql_success)
      call assert_equals(var1=.True., var2=ql_success, message='QL_Decomposition(t_mat33_bsp03) successful')
      res_eigenval = Sort_Vector(vec=[matrix(1, 1), matrix(2, 2), matrix(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='QL_Decomposition(t_mat33_bsp03)')

      matrix = t_mat33_bsp04
      eigenval = [1.0_dp, 2.0_dp, 11.0_dp]
      call QL_Decomposition(matrix=matrix, nel=3, eigenvec_mat=eigenvec_mat, ql_success=ql_success)
      call assert_equals(var1=.True., var2=ql_success, message='QL_Decomposition(t_mat33_bsp04) successful')
      res_eigenval = Sort_Vector(vec=[matrix(1, 1), matrix(2, 2), matrix(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='QL_Decomposition(t_mat33_bsp04)')

      matrix = t_mat33_bsp07
      eigenval = [-1.0_dp, 1.0_dp, 7.0_dp]
      call Tridiagonal_Matrix(matrix=matrix, nel=3, trans_mat=trans_mat)
      call QL_Decomposition(matrix=matrix, nel=3, eigenvec_mat=eigenvec_mat, ql_success=ql_success)
      call assert_equals(var1=.True., var2=ql_success, message='QL_Decomposition(t_mat33_bsp07) successful')
      res_eigenval = Sort_Vector(vec=[matrix(1, 1), matrix(2, 2), matrix(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='QL_Decomposition(t_mat33_bsp07)')

      matrix = t_mat33_bsp08
      eigenval = [-0.7_dp, -0.5_dp, -0.2_dp]
      call Tridiagonal_Matrix(matrix=matrix, nel=3, trans_mat=trans_mat)
      call QL_Decomposition(matrix=matrix, nel=3, eigenvec_mat=eigenvec_mat, ql_success=ql_success)
      call assert_equals(var1=.True., var2=ql_success, message='QL_Decomposition(t_mat33_bsp08) successful')
      res_eigenval = Sort_Vector(vec=[matrix(1, 1), matrix(2, 2), matrix(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='QL_Decomposition(t_mat33_bsp08)')
   end subroutine QL_Decomposition_Test


   ! --------------------------------------------------------------- !
   subroutine Eigendecomposition6_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3) :: matrix, eigenval_mat, eigenvec_mat
      real(dp), dimension(3) :: eigenval, res_eigenval
      logical :: ql_success

      ! NOTE: Sort_Vector() has to be called since the returned eigenvalues are not sorted.
      !       Consider also checking eigenvectors

      matrix = t_mat33_id
      eigenval = [1.0_dp, 1.0_dp, 1.0_dp]
      call Eigendecomposition(mat=matrix, nel=3, eigenvalues=eigenval_mat, eigenvectors=eigenvec_mat)
      call assert_equals(var1=.True., var2=ql_success, message='Eigendecomposition(t_mat33_id) successful')
      res_eigenval = Sort_Vector(vec=[eigenval_mat(1, 1), eigenval_mat(2, 2), eigenval_mat(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='eigenvalues of Eigendecomposition(t_mat33_id)')

      matrix = t_mat33_bsp03
      eigenval = [1.0_dp, 5.0_dp, 7.0_dp]
      call Eigendecomposition(mat=matrix, nel=3, eigenvalues=eigenval_mat, eigenvectors=eigenvec_mat)
      call assert_equals(var1=.True., var2=ql_success, message='Eigendecomposition(t_mat33_bsp03) successful')
      res_eigenval = Sort_Vector(vec=[eigenval_mat(1, 1), eigenval_mat(2, 2), eigenval_mat(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='eigenvalues of Eigendecomposition(t_mat33_bsp03)')

      matrix = t_mat33_bsp04
      eigenval = [1.0_dp, 2.0_dp, 11.0_dp]
      call Eigendecomposition(mat=matrix, nel=3, eigenvalues=eigenval_mat, eigenvectors=eigenvec_mat)
      call assert_equals(var1=.True., var2=ql_success, message='Eigendecomposition(t_mat33_bsp04) successful')
      res_eigenval = Sort_Vector(vec=[eigenval_mat(1, 1), eigenval_mat(2, 2), eigenval_mat(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='eigenvalues of Eigendecomposition(t_mat33_bsp04)')

      matrix = t_mat33_bsp07
      eigenval = [-1.0_dp, 1.0_dp, 7.0_dp]
      call Eigendecomposition(mat=matrix, nel=3, eigenvalues=eigenval_mat, eigenvectors=eigenvec_mat)
      call assert_equals(var1=.True., var2=ql_success, message='Eigendecomposition(t_mat33_bsp07) successful')
      res_eigenval = Sort_Vector(vec=[eigenval_mat(1, 1), eigenval_mat(2, 2), eigenval_mat(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='eigenvalues of Eigendecomposition(t_mat33_bsp07)')

      matrix = t_mat33_bsp08
      eigenval = [-0.7_dp, -0.5_dp, -0.2_dp]
      call Eigendecomposition(mat=matrix, nel=3, eigenvalues=eigenval_mat, eigenvectors=eigenvec_mat)
      call assert_equals(var1=.True., var2=ql_success, message='Eigendecomposition(t_mat33_bsp08) successful')
      res_eigenval = Sort_Vector(vec=[eigenval_mat(1, 1), eigenval_mat(2, 2), eigenval_mat(3, 3)])
      call assert_equals(var1=eigenval, var2=res_eigenval, n=3, &
         delta=setting_epsilon, message='eigenvalues of Eigendecomposition(t_mat33_bsp08)')
   end subroutine Eigendecomposition6_Test


   ! --------------------------------------------------------------- !
   subroutine Matrix_Exponential_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(9) :: res

      ! NOTE: Matrices have to be symmetric. Also positive eigenvalues decrease precision

      res = [exp(-0.3_dp), exp(-30.0_dp), exp(-3.0_dp), &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call assert_equals(var1=res, var2=M_To_V9(Matrix_Exponential(t_repr_bsp02)), n=9, &
         delta=setting_epsilon, message='Matrix_Exponential(t_repr_bsp02)')

      res = [682.033820683619_dp, 841.705848802297_dp, 10261.03461076645_dp, &
         757.193790747539_dp, 757.193790747539_dp, 2938.694736363561_dp, &
         2938.694736363561_dp, 2643.348659238216_dp, 2643.348659238216_dp]
      call assert_equals(var1=res, var2=M_To_V9(Matrix_Exponential(t_repr_bsp05)), n=9, &
         delta=10000.0_dp*setting_epsilon, message='Matrix_Exponentialt_repr_bsp05)')

      res = [366.1813728348143968_dp, 366.1813728348141694_dp, 367.3565740284584535_dp, &
         -365.8134933936428297_dp, -365.8134933936428297_dp, 364.6382921999992_dp, &
         364.6382921999992_dp, -364.6382921999993414_dp, -364.6382921999993414_dp]
      call assert_equals(var1=res, var2=M_To_V9(Matrix_Exponential(t_repr_bsp07)), n=9, &
         delta=10.0_dp*setting_epsilon, message='Matrix_Exponentialt_repr_bsp07)')

      res = [0.6222913462071376_dp, 0.6222913462071376_dp, 0.6772640241677494_dp, &
         0.1257060424157281_dp, 0.1257060424157281_dp, 0.0707333644551161_dp, &
         0.0707333644551161_dp, 0.0707333644551161_dp, 0.0707333644551161_dp]
      call assert_equals(var1=res, var2=M_To_V9(Matrix_Exponential(t_repr_bsp08)), n=9, &
         delta=setting_epsilon, message='Matrix_Exponentialt_repr_bsp08)')
   end subroutine Matrix_Exponential_Test


   ! --------------------------------------------------------------- !
   subroutine Is_Diagonal_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=.True., var2=Is_Diagonal(t_repr_id), &
         message='Is_Diagonal(t_repr_id)')
      call assert_equals(var1=.True., var2=Is_Diagonal(t_repr_zero), &
         message='Is_Diagonal(t_repr_zero)')
      call assert_equals(var1=.True., var2=Is_Diagonal(t_repr_bsp01), &
         message='Is_Diagonal(t_repr_bsp01)')
      call assert_equals(var1=.True., var2=Is_Diagonal(t_repr_bsp02), &
         message='Is_Diagonal(t_repr_bsp02)')

      call assert_equals(var1=.False., var2=Is_Diagonal(t_repr_bsp03), &
         message='Is_Diagonal(t_repr_bsp03)')
      call assert_equals(var1=.False., var2=Is_Diagonal(t_repr_bsp04), &
         message='Is_Diagonal(t_repr_bsp04)')
      call assert_equals(var1=.False., var2=Is_Diagonal(t_repr_bsp05), &
         message='Is_Diagonal(t_repr_bsp05)')
      call assert_equals(var1=.False., var2=Is_Diagonal(t_repr_bsp06), &
         message='Is_Diagonal(t_repr_bsp06)')
   end subroutine Is_Diagonal_Test


   ! --------------------------------------------------------------- !
   subroutine Value_In_Interval_Test()
   ! --------------------------------------------------------------- !
      call assert_equals(var1=.True., var2=Value_In_Interval(value=0.5_dp, limits=[0.0_dp, 1.0_dp]), &
         message='Value_In_Interval(value=0.5_dp, limits=[0.0_dp, 1.0_dp])')
      call assert_equals(var1=.True., var2=Value_In_Interval(value=-6.002_dp, limits=[-7.0_dp, -6.0_dp]), &
         message='Value_In_Interval(value=-6.002_dp, limits=[-7.0_dp, -6.0_dp])')
      call assert_equals(var1=.True., var2=Value_In_Interval(value=50.0_dp, limits=[50.0_dp, 100.0_dp]), &
         message='Value_In_Interval(value=100.0_dp, limits=[50.0_dp, 100.0_dp])')
      call assert_equals(var1=.True., var2=Value_In_Interval(value=100.0_dp, limits=[50.0_dp, 100.0_dp]), &
         message='Value_In_Interval(value=100.0_dp, limits=[50.0_dp, 100.0_dp])')
      call assert_equals(var1=.True., var2=Value_In_Interval(value=-15.0_dp, limits=[-10.0_dp, -20.0_dp]), &
         message='Value_In_Interval(value=-15.0_dp, limits=[-10.0_dp, -20.0_dp])')
      call assert_equals(var1=.False., var2=Value_In_Interval(value=15.0_dp, limits=[10.0_dp, -5.0_dp]), &
         message='Value_In_Interval(value=15.0_dp, limits=[10.0_dp, -5.0_dp])')
      call assert_equals(var1=.False., var2=Value_In_Interval(value=-90.0_dp, limits=[80.0_dp, 100.0_dp]), &
         message='Value_In_Interval(value=-90.0_dp, limits=[80.0_dp, 100.0_dp])')
      call assert_equals(var1=.False., var2=Value_In_Interval(value=90.0_dp, limits=[-80.0_dp, -100.0_dp]), &
         message='Value_In_Interval(value=90.0_dp, limits=[-80.0_dp, -100.0_dp])')
   end subroutine Value_In_Interval_Test


!   ! --------------------------------------------------------------- !
!   subroutine Abort_If_Not_In_Interval_Test()
!   ! --------------------------------------------------------------- !
!
!   end subroutine Abort_If_Not_In_Interval_Test


   ! --------------------------------------------------------------- !
   subroutine Element_In_Tensor_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(__tensor__) :: tens, ref_tens
      real(dp) :: value, elem
      integer :: ref_idx, ref_jdx, idx, jdx, kdx, ldx

      ! Checks that every position will be assigned once.
      ! Should also check that individual positions match expected positions

      tens = 0.0_dp
      ref_tens = 1.0_dp

      do ref_idx = 1, __nelmat__
         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
         do ref_jdx = 1, __nelmat__
            call Ref_Index(ref_idx=ref_jdx, idx=kdx, jdx=ldx)
            call Set_Element_In_Tensor(tens=tens, idx=idx, jdx=jdx, kdx=kdx, ldx=ldx, val=1.0_dp)
         end do
      end do

      call assert_equals(var1=T_To_M99(ref_tens), var2=T_To_M99(tens), n=9, m=9, &
         message='Set_Element_In_Tensor() (all filled)')
   end subroutine Element_In_Tensor_Test


   ! --------------------------------------------------------------- !
   subroutine Matrix_And_Tensorpart_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(__matrix__) :: mat
      real(dp), dimension(__tensor__) :: tens
      integer :: ref_idx, idx, jdx, kdx
      character(len=9) :: id_num

      id_num = 'ref_idx= '

      ! Depends on Ref_Index()
      do ref_idx = 1, __nelmat__
         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
         mat = reshape([(10.0*ref_idx+kdx, kdx = 1, __nelmat__)], [__matrix__])
         call Set_Matrix_In_Tensorpart(tens=tens, idx=idx, jdx=jdx, mat=mat)
      end do

      do ref_idx = 1, __nelmat__
         write(id_num(9:9), '(i1)') ref_idx

         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
         mat = reshape([(10.0*ref_idx+kdx, kdx = 1, __nelmat__)], [__matrix__])
         call assert_equals(var1=M_To_V9(mat), &
            var2=M_To_V9(Get_Matrix_From_Tensorpart(tens=tens, idx=idx, jdx=jdx)), &
            n=__nelmat__, message='Set_Matrix_In_Tensorpart/Get_Matrix_From_Tensorpart() (ref_idx='// id_num // ')')
      end do
   end subroutine Matrix_And_Tensorpart_Test


   ! --------------------------------------------------------------- !
   subroutine Elements_And_Matrixlist_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(1, __matrix__) :: onemat_list
      real(dp), dimension(2, __matrix__) :: twomat_list
      real(dp), dimension(4, __matrix__) :: fourmat_list
      integer :: ref_idx, idx, jdx
      character(len=9) :: id_num

      id_num = 'ref_idx= '

      ! Depends on Ref_Index()
      do ref_idx = 1, __nelmat__
         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
         call Set_Elements_In_Matrixlist(matlist=onemat_list, nel=1, &
            idx=idx, jdx=jdx, elemlist=[-1.0_dp*ref_idx])
         call Set_Elements_In_Matrixlist(matlist=twomat_list, nel=2, &
            idx=idx, jdx=jdx, elemlist=[2.0_dp*ref_idx, 2.0_dp*ref_idx+1])
         call Set_Elements_In_Matrixlist(matlist=fourmat_list, nel=4, &
            idx=idx, jdx=jdx, elemlist=[-4.0_dp*ref_idx, -4.0_dp*ref_idx-1, &
            -4.0_dp*ref_idx-2, -4.0_dp*ref_idx-3])
      end do

      do ref_idx = 1, __nelmat__
         write(id_num(9:9), '(i1)') ref_idx

         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)

         call assert_equals(var1=[-1.0_dp*ref_idx], var2=Get_Elements_From_Matrixlist(matlist=onemat_list, &
            nel=1, idx=idx, jdx=jdx), n=1, &
            message='Set_Elements_In_Matrixlist/Get_Elements_From_Matrixlist ' // &
               '(onemat_list, ref_idx=' // id_num //')')
         call assert_equals(var1=[2.0_dp*ref_idx, 2.0_dp*ref_idx+1], &
            var2=Get_Elements_From_Matrixlist(matlist=twomat_list, nel=2, idx=idx, jdx=jdx), &
            n=2, message='Set_Elements_In_Matrixlist/Get_Elements_From_Matrixlist ' // &
               '(twomat_list, ref_idx=' // id_num //')')
         call assert_equals(var1=[-4.0_dp*ref_idx, -4.0_dp*ref_idx-1, -4.0_dp*ref_idx-2, -4.0_dp*ref_idx-3], &
            var2=Get_Elements_From_Matrixlist(matlist=fourmat_list, nel=4, idx=idx, jdx=jdx), &
            n=4, message='Set_Elements_In_Matrixlist/Get_Elements_From_Matrixlist ' // &
               '(fourmat_list, ref_idx=' // id_num //')')
      end do
   end subroutine Elements_And_Matrixlist_Test


   ! --------------------------------------------------------------- !
   subroutine Ref_Index_Test()
   ! --------------------------------------------------------------- !
      integer :: idx, jdx

      call Ref_Index(ref_idx=1, idx=idx, jdx=jdx)
      call assert_equals(var1=[1, 1], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=1, ...)')

      call Ref_Index(ref_idx=2, idx=idx, jdx=jdx)
      call assert_equals(var1=[2, 2], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=2, ...)')

      call Ref_Index(ref_idx=3, idx=idx, jdx=jdx)
      call assert_equals(var1=[3, 3], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=3, ...)')

      call Ref_Index(ref_idx=4, idx=idx, jdx=jdx)
      call assert_equals(var1=[1, 2], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=4, ...)')

#ifdef REPR_MAT66
      call Ref_Index(ref_idx=5, idx=idx, jdx=jdx)
      call assert_equals(var1=[2, 3], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=5, ...)')

      call Ref_Index(ref_idx=6, idx=idx, jdx=jdx)
      call assert_equals(var1=[1, 3], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=6, ...)')
#else
      call Ref_Index(ref_idx=5, idx=idx, jdx=jdx)
      call assert_equals(var1=[2, 1], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=5, ...)')

      call Ref_Index(ref_idx=6, idx=idx, jdx=jdx)
      call assert_equals(var1=[2, 3], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=6, ...)')

      call Ref_Index(ref_idx=7, idx=idx, jdx=jdx)
      call assert_equals(var1=[3, 2], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=7, ...)')

      call Ref_Index(ref_idx=8, idx=idx, jdx=jdx)
      call assert_equals(var1=[1, 3], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=8, ...)')

      call Ref_Index(ref_idx=9, idx=idx, jdx=jdx)
      call assert_equals(var1=[3, 1], var2=[idx, jdx], n=2, &
         message='Ref_Index(ref_idx=9, ...)')
#endif
   end subroutine Ref_Index_Test


   ! --------------------------------------------------------------- !
   subroutine Perturbate_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(__matrix__) :: mat
      real(dp) :: theta, elem_before_pert, elem_after
      integer :: ref_idx, idx, jdx
      character(len=9) :: id_num

      id_num = 'ref_idx= '

      mat = t_repr_bsp06
      theta = 500.0_dp

      do ref_idx = 1, __nelmat__
         write(id_num(9:9), '(i1)') ref_idx

         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
#ifdef REPR_MAT66
         if (ref_idx > 3) then
            elem_before_pert = mat(ref_idx) + 0.5_dp*theta
         else
            elem_before_pert = mat(ref_idx) + theta
         end if
#else
#ifdef REPR_MAT99
         elem_before_pert = mat(ref_idx) + theta
#else
         elem_before_pert = mat(idx, jdx) + theta
#endif
#endif
         call Perturbate(mat=mat, idx=idx, jdx=jdx, theta=theta)

#ifndef REPR_TENS3333
         elem_after = mat(ref_idx)
#else
         elem_after = mat(idx, jdx)
#endif
         call assert_equals(var1=elem_after, var2=elem_before_pert, &
            message='Perturbate(mat=t_repr_bsp06) (ref_idx=' // id_num // ')')
      end do
   end subroutine Perturbate_Test


   ! --------------------------------------------------------------- !
   subroutine Vec_And_Mat33_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3) :: mat33
      real(dp), dimension(__matrix__) :: repr_mat

      ! Nothing to do for REPR_TENS3333
#ifndef REPR_TENS3333
      ! Transform defined examples back and forth
      repr_mat = Mat33_To_Vec(mat33=t_mat33_id)
      call assert_equals(var1=repr_mat, var2=t_repr_id, n=__nelmat__, &
         delta=setting_epsilon, message='Mat33_To_Vec(t_mat33_id)')
      mat33 = Vec_To_Mat33(repr_mat)
      call assert_equals(var1=mat33, var2=t_mat33_id, n=3, m=3, &
         delta=setting_epsilon, message='Vec_To_Mat33(t_repr_id)')
      !
      repr_mat = Mat33_To_Vec(mat33=t_mat33_zero)
      call assert_equals(var1=repr_mat, var2=t_repr_zero, n=__nelmat__, &
         delta=setting_epsilon, message='Mat33_To_Vec(t_mat33_zero)')
      mat33 = Vec_To_Mat33(repr_mat)
      call assert_equals(var1=mat33, var2=t_mat33_zero, n=3, m=3, &
         delta=setting_epsilon, message='Vec_To_Mat33(t_repr_zero)')
      !
      repr_mat = Mat33_To_Vec(mat33=t_mat33_bsp01)
      call assert_equals(var1=repr_mat, var2=t_repr_bsp01, n=__nelmat__, &
         delta=setting_epsilon, message='Mat33_To_Vec(t_mat33_bsp01)')
      mat33 = Vec_To_Mat33(t_repr_bsp01)
      call assert_equals(var1=mat33, var2=t_mat33_bsp01, n=3, m=3, &
         delta=setting_epsilon, message='Vec_To_Mat33(t_repr_bsp01)')
      !
      repr_mat = Mat33_To_Vec(mat33=t_mat33_bsp02)
      call assert_equals(var1=repr_mat, var2=t_repr_bsp02, n=__nelmat__, &
         delta=setting_epsilon, message='Mat33_To_Vec(t_mat33_bsp02)')
      mat33 = Vec_To_Mat33(t_repr_bsp02)
      call assert_equals(var1=mat33, var2=t_mat33_bsp02, n=3, m=3, &
         delta=setting_epsilon, message='Vec_To_Mat33(t_repr_bsp02)')
      !
      repr_mat = Mat33_To_Vec(mat33=t_mat33_bsp03)
      call assert_equals(var1=repr_mat, var2=t_repr_bsp03, n=__nelmat__, &
         delta=setting_epsilon, message='Mat33_To_Vec(t_mat33_bsp03)')
      mat33 = Vec_To_Mat33(t_repr_bsp03)
      call assert_equals(var1=mat33, var2=t_mat33_bsp03, n=3, m=3, &
         delta=setting_epsilon, message='Vec_To_Mat33(t_repr_bsp03)')
      !
      repr_mat = Mat33_To_Vec(mat33=t_mat33_bsp04)
      call assert_equals(var1=repr_mat, var2=t_repr_bsp04, n=__nelmat__, &
         delta=setting_epsilon, message='Mat33_To_Vec(t_mat33_bsp04)')
      mat33 = Vec_To_Mat33(t_repr_bsp04)
      call assert_equals(var1=mat33, var2=t_mat33_bsp04, n=3, m=3, &
         delta=setting_epsilon, message='Vec_To_Mat33(t_repr_bsp04)')
      !
      repr_mat = Mat33_To_Vec(mat33=t_mat33_bsp05)
      call assert_equals(var1=repr_mat, var2=t_repr_bsp05, n=__nelmat__, &
         delta=setting_epsilon, message='Mat33_To_Vec(t_mat33_bsp05)')
      mat33 = Vec_To_Mat33(t_repr_bsp05)
      call assert_equals(var1=mat33, var2=t_mat33_bsp05, n=3, m=3, &
         delta=setting_epsilon, message='Vec_To_Mat33(t_repr_bsp05)')
#endif
   end subroutine Vec_And_Mat33_Test


   ! --------------------------------------------------------------- !
   subroutine Vec9_And_Mat_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(9) :: vec_out
      real(dp), dimension(__matrix__) :: temp_mat

      temp_mat = Vec9_To_Mat(vec9=t_vec9_bsp01)
      vec_out = Mat_To_Vec9(mat=temp_mat)
      call assert_equals(var1=t_vec9_bsp01, var2=vec_out, n=9, &
         delta=setting_epsilon, message='Mat_To_Vec9(Vec9_To_Mat(t_vec9_bsp01))')

      temp_mat = Vec9_To_Mat(vec9=t_vec9_bsp02)
      vec_out = Mat_To_Vec9(mat=temp_mat)
      call assert_equals(var1=t_vec9_bsp02, var2=vec_out, n=9, &
         delta=setting_epsilon, message='Mat_To_Vec9(Vec9_To_Mat(t_vec9_bsp02))')

#ifndef REPR_MAT66
      temp_mat = Vec9_To_Mat(vec9=t_vec9_bsp03)
      vec_out = Mat_To_Vec9(mat=temp_mat)
      call assert_equals(var1=t_vec9_bsp03, var2=vec_out, n=9, &
         delta=setting_epsilon, message='Mat_To_Vec9(Vec9_To_Mat(t_vec9_bsp03))')
#endif

      ! Add more test cases
   end subroutine Vec9_And_Mat_Test


   ! --------------------------------------------------------------- !
   subroutine Tens_And_Mat99_Test()
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9) :: mat_out
      real(dp), dimension(__tensor__) :: temp_tens

      ! Nothing to do for other representations than REPR_TENS3333
#ifdef REPR_TENS3333
      temp_tens = Mat99_To_Tens(mat99=t_mat99_bsp01)
      mat_out = Tens_To_Mat99(tens3333=temp_tens)
      call assert_equals(var1=t_mat99_bsp01, var2=mat_out, n=9, m=9, &
         delta=setting_epsilon, message='Tens_To_Mat99(Mat99_To_Tens(t_mat99_bsp01))')

      temp_tens = Mat99_To_Tens(mat99=t_mat99_bsp02)
      mat_out = Tens_To_Mat99(tens3333=temp_tens)
      call assert_equals(var1=t_mat99_bsp02, var2=mat_out, n=9, m=9, &
         delta=setting_epsilon, message='Tens_To_Mat99(Mat99_To_Tens(t_mat99_bsp02))')

      ! Add more test cases
#endif
   end subroutine Tens_And_Mat99_Test


   ! --------------------------------------------------------------- !
   subroutine Packing_Test()
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_num_statevariables, setting_max_internal_states
      !
      real(dp), dimension(__matrix__) :: stress_in, stress_out
      real(dp), dimension(__tensor__) :: jac_stress_in, jac_stress_out
      real(dp), dimension(setting_num_statevariables) :: statevariables_in, statevariables_out, tmp_statevariables
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statev_in, jac_statev_out
      real(dp), dimension(setting_max_internal_states) :: packed_states
      integer :: ref_idx, idx, jdx, kdx
      character(len=9) :: id_num

      id_num = 'ref_idx= '

      ! Checking the packed vector is less important than that the packed was extracted properly

      stress_in = t_repr_bsp05
      jac_stress_in = t_repr_bsp404
      statevariables_in = [(idx, idx = 1, setting_num_statevariables)]

      ! A thorough initialization of jac_statev_in currently depends on
      ! Ref_Index() and Set_Elements_In_Matrixlist() and Get_Elements_From_Matrixlist() for checking
      jac_statev_in = 0.0_dp
      do ref_idx = 1, __nelmat__
         tmp_statevariables = [(-ref_idx+0.02*kdx, kdx = 1, setting_num_statevariables)]
         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
         call Set_Elements_In_Matrixlist(matlist=jac_statev_in, nel=setting_num_statevariables, &
            idx=idx, jdx=jdx, elemlist=tmp_statevariables)
      end do
      packed_states = Pack_States(stress=stress_in, jac_stress=jac_stress_in, &
         statevariables=statevariables_in, jac_statevariables=jac_statev_in)

      call Unpack_States(input_states=packed_states, stress=stress_out, jac_stress=jac_stress_out, &
         statevariables=statevariables_out, jac_statevariables=jac_statev_out)
      call assert_equals(var1=M_To_V9(stress_in), var2=M_To_V9(stress_out), n=9, &
         message='Pack/Unpack_States() stress')
      call assert_equals(var1=T_To_M99(jac_stress_in), var2=T_To_M99(jac_stress_out), n=9, m=9, &
         message='Pack/Unpack_States() jac_stress')
      call assert_equals(var1=statevariables_in, var2=statevariables_out, n=setting_num_statevariables, &
         message='Pack/Unpack_States() statevariables')

      do ref_idx = 1, __nelmat__
         write(id_num(9:9), '(i1)') ref_idx

         tmp_statevariables = [(-ref_idx+0.02*kdx, kdx = 1, setting_num_statevariables)]
         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
         call assert_equals(var1=tmp_statevariables, var2=Get_Elements_From_Matrixlist(matlist=jac_statev_out, &
            nel=setting_num_statevariables, idx=idx, jdx=jdx), &
            n=setting_num_statevariables, message='Pack/Unpack_States() jac_statevariables ' // &
            '(ref_idx=' // id_num // ')')
      end do
   end subroutine Packing_Test


!   ! --------------------------------------------------------------- !
!   subroutine Import_Export_Test()
!   ! --------------------------------------------------------------- !
!
!   end subroutine Import_Export_Test
end module Math_Operations_Test
