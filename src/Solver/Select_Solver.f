   ! --------------------------------------------------------------- !
   subroutine Select_Solver(new_solver, name)
   ! --------------------------------------------------------------- !
      use Euler_Explicit_Class
      use Euler_Richardson_Class
      use RK23_Simpson_Class
      use RK23_Bogacki_Shampine_Class
      use RK45_Cash_Karp_Class
      use RK45_Dormand_Prince_Class
      use RK45_Fehlberg_Class
      !
      class(Solver), allocatable, intent(out) :: new_solver
      character(len=*), intent(in) :: name
      ! ------------------------------------------------------------ !
      associate(identifier => name(1:7))
         if (identifier == 'Eul-Exp') then
            allocate(Euler_Explicit::new_solver)
            call new_solver%Initialize(name='Euler-Explicit', stepsize_fixed=.True.)
         !
         else if (identifier == 'Richard') then
            allocate(Euler_Richardson::new_solver)
            ! Similar to Fellin et al. (2009): Adaptive integration of constitutive rate equations, p. 3
            call new_solver%Initialize(name='Euler-Richardson', exp_stepgrow=-0.5_dp, exp_stepshrink=-0.5_dp, &
               max_stepgrow=2.0_dp, min_stepshrink=0.2_dp)
         !
         else if (identifier == 'RK-S 23') then
            allocate(RK23_Simpson::new_solver)
            ! Inspired by from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK2/3 (Simpson)', exp_stepgrow=-1.0_dp/3.0_dp, &
               exp_stepshrink=-0.5_dp, max_stepgrow=4.0_dp, min_stepshrink=0.2_dp)
         !
         else if (identifier == 'RK-BS23') then
            allocate(RK23_Bogacki_Shampine::new_solver)
            ! Inspired by from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK2/3 (Bogacki-Shampine)', exp_stepgrow=-1.0_dp/3.0_dp, &
               exp_stepshrink=-0.5_dp, max_stepgrow=3.0_dp, min_stepshrink=0.2_dp)
         !
         else if (identifier == 'RK-CK45') then
            allocate(RK45_Cash_Karp::new_solver)
            ! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Cash-Karp)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)
         !
         else if (identifier == 'RK-DP45') then
            allocate(RK45_Dormand_Prince::new_solver)
            ! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Dormand-Prince)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)
         !
         else if (identifier == 'RK-F 45') then
            allocate(RK45_Fehlberg::new_solver)
            ! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Fehlberg)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)
         !
         else
            call Write_Error_And_Exit('Select_Solver: Identifier >' // identifier // '< unknown')
         end if
      end associate
   end subroutine Select_Solver
