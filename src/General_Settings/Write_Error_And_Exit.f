   ! --------------------------------------------------------------- !
   subroutine Write_Error_And_Exit(message)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: message
      ! ------------------------------------------------------------ !
      character(len=len(identifier_log_start) + len(message)) :: ErrorMessage

      ErrorMessage = identifier_log_start // message

      ! For Abaqus calls either stdb_abqerr (Abaqus/Standard) or xplb_abqerr (Abaqus/Explicit) is
      ! used. Besides a printf-ed message string, an int (array), a real (array) and a char (array)
      ! are expected - See Abaqus user routine guide
#ifdef ABQ_STD_CALLING
      call stdb_abqerr(-2, ErrorMessage, 0, 0.0, ' ')                ! Output on Abaqus error log channel (-2)
      call xit()
#else
#ifdef ABQ_EXP_CALLING
      call xplb_abqerr(-2, ErrorMessage, 0, 0.0, ' ')                ! Output on Abaqus error log channel (-2)
      call xplb_exit()
#else
#ifdef MATLAB_CALLING
      call mexErrMsgTxt(ErrorMessage)
#else
#ifdef PLAXIS_DLL
      write(1, *) ErrorMessage                                       ! Output on Plaxis log channel (1)
      stop 'Execution aborted'
#else
#ifdef OCTAVE_CALLING
      write(0, *) ErrorMessage
      call xStopx('Execution aborted')
#else
      write(0, *) ErrorMessage                                       ! Output on default log channel (0)
      stop 'Execution aborted'
#endif
#endif
#endif
#endif
#endif
   end subroutine Write_Error_And_Exit
