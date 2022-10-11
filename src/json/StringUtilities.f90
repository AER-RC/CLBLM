module stringmethods
use json_data_types, only : MAX_BUFFER_LENGTH => MAX_STRING_LENGTH

public :: toStringInt, toStringFloat, toStringDouble, splitFilePath

interface toString
    module procedure toStringInt
    module procedure toStringFloat
    !module procedure toStringDouble
end interface

contains

    character(len=MAX_BUFFER_LENGTH) function toStringInt(iValue)
        integer, intent(in) :: iValue
        toStringInt = repeat(' ', MAX_BUFFER_LENGTH)
        write(toStringInt, *) iValue
        toStringInt = adjustl(toStringInt)
    end function toStringInt

    character(len=MAX_BUFFER_LENGTH) function toStringFloat(fValue)
        real, intent(in) :: fValue
        toStringFloat = repeat(' ', MAX_BUFFER_LENGTH)
        write(toStringFloat, *) fValue
        toStringFloat = adjustl(toStringFloat)
    end function toStringFloat

    subroutine splitFilePath(fullpath, basename, directory)
        implicit none

        character(len=*), intent(in) :: fullpath
        character(len=*), intent(out) :: basename, directory
        integer :: lastSlash

        basename = repeat(' ', len(basename))
        directory = repeat(' ', len(directory))

        lastSlash = index(fullpath, '/', BACK=.true.)
        if (lastSlash .gt. 0) directory = fullpath(1:lastSlash)
        basename = fullpath(lastSlash + 1:len(trim(fullpath)))
    end subroutine

    !character(len=MAX_BUFFER_LENGTH) function toStringDouble(dValue)
    !    real*8, intent(in) :: dValue
    !    toStringDouble = repeat(' ', MAX_BUFFER_LENGTH)
    !    write(toStringDouble, *) dValue
    !    toStringDouble = adjustl(toStringDouble)
    !end function toStringDouble

end module
