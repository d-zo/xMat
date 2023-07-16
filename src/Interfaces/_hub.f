#ifdef PLAXIS_DLL
#adddirectory hub Plaxis
#endif


#ifdef MATLAB_CALLING
#addfile subroutine mexFunction
#endif


#addfile subroutine Abaqus/UMAT


#addfile subroutine Abaqus/VUMAT


#addfile subroutine dot_values


#addfile subroutine xmat_console


#ifdef CINTER
#addfile subroutine dot_values_cinter


#addfile subroutine xmat_console_cinter
#endif


