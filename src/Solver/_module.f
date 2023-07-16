#adddirectory module Baseclass


#adddirectory module Euler_Explicit


#adddirectory module Euler_Richardson


#adddirectory module RK23_Simpson


#adddirectory module RK23_Bogacki_Shampine


#adddirectory module RK45_Cash_Karp


#adddirectory module RK45_Dormand_Prince


#adddirectory module RK45_Fehlberg


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Solver_Class
   use General_Settings, only: dp, Write_Error_And_Exit
   use Solver_Baseclass
   implicit none

   private
   public :: Select_Solver, Solver


   contains


#addfile subroutine Select_Solver
end module Solver_Class
