
! utils module for fortran routine. Right now just ini file.

module utils

#ifndef __GFORTRAN__
        integer iargc
        external iargc
#endif

contains

#ifdef __GFORTRAN__

  ! ===========================================================
  function iargc()
    ! ===========================================================
    integer iargc
    ! ===========================================================

    iargc=command_argument_count()

  end function iargc

  ! ===========================================================
  subroutine getarg(num, res)
    ! ===========================================================
    integer, intent(in) :: num
    character(len=*), intent(out) :: res
    integer l, err
    ! ===========================================================
    call get_command_argument(num,res,l,err)
  end subroutine getarg

#endif


  function GetParamCount()
   integer GetParamCount

    GetParamCount = iargc()

  end function GetParamCount




    function GetParam(i)

   character(LEN=512) GetParam
   integer, intent(in) :: i

   if (iargc() < i) then
     GetParam = ''
   else
    call getarg(i,GetParam)
   end if
  end function GetParam



  end module utils