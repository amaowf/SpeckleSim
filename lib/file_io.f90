module file_io

    interface operator( .f. )
      module procedure file_exists
    end interface

    contains

    function file_exists(filename) result(res)
      implicit none
      character(len=*),intent(in) :: filename
      logical                     :: res

      ! Check if the file exists
      inquire( file=trim(filename), exist=res )
    end function

end module
