   ! --------------------------------------------------------------- !
   subroutine Write_Message(message)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: message
      ! ------------------------------------------------------------ !
      character(len=len(identifier_log_start) + len(message)) :: InfoMessage

      InfoMessage = identifier_log_start // message

      ! For Abaqus calls either stdb_abqerr (Abaqus/Standard) or xplb_abqerr (Abaqus/Explicit) is
      ! used. Besides a printf-ed message string, an int (array), a real (array) and a char (array)
      ! are expected - See Abaqus user routine guide
#ifdef ABQ_STD_CALLING
      call stdb_abqerr(1, InfoMessage, 0, 0.0, ' ')                  ! Output on Abaqus message log channel (1)
#else
#ifdef ABQ_EXP_CALLING
      call xplb_abqerr(1, InfoMessage, 0, 0.0, ' ')                  ! Output on Abaqus message log channel (1)
#else
#ifdef MATLAB_CALLING
      call mexPrintf(InfoMessage)
#else
#ifdef PLAXIS_DLL
      write(1, *) InfoMessage                                        ! Output on Plaxis log channel (1)
#else
      write(0, *) InfoMessage                                        ! Output on default log channel (0)
#endif
#endif
#endif
#endif
   end subroutine Write_Message
