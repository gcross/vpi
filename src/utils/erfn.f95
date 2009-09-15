!@+leo-ver=4-thin
!@+node:gcross.20090914154930.2019:@thin erfn.f95
!@@language fortran90

module erfn

  interface
    pure subroutine calerf(x,result,jint)
      double precision, intent(in) :: x
      integer, intent(in) :: jint
      double precision, intent(out) :: result
    end subroutine
  end interface

contains

  elemental function erf(x) result (result)
    double precision, intent(in) :: x
    double precision :: result
    call calerf(x,result,0)
  end function

  elemental function erfc(x) result (result)
    double precision, intent(in) :: x
    double precision :: result
    call calerf(x,result,1)
  end function

  elemental function erfcx(x) result (result)
    double precision, intent(in) :: x
    double precision :: result
    call calerf(x,result,2)
  end function

end module
!@-node:gcross.20090914154930.2019:@thin erfn.f95
!@-leo
