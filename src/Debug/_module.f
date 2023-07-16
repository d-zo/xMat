! ==================================================================================================================== !
module Debug
   use General_Settings, only: dp
   implicit none

   integer, parameter :: onenumberlength = 12
   integer, parameter :: mat_identifierlength = 6
   integer, parameter :: tens_identifierlength = 12
   character(len=6), parameter :: submatrixname = '     M'
   character(len=6), parameter :: debug_format = 'f12.5 '            ! character length of debug_format and
   character(len=6), parameter :: alternative_format = 'es12.3'      ! alternative_format must be the same

   interface Formatval
      module procedure Format_logical, Format_int_dim0, Format_dp_dim0, Format_dp_dim1, Format_dp_dim2, Format_dp_dim4
   end interface

   private
   public Formatval


   contains


#addfile function Number_To_String


#addfile function Format_logical


#addfile function Format_int_dim0


#addfile function Format_dp_dim0


#addfile function Format_dp_dim1


#addfile function Format_dp_dim2


#addfile function Format_dp_dim4
end module Debug
