module vpi_coordinate_utils
  use vpi_defines
  implicit none

contains 

function eval_polar_angle( x, y ) result ( angle )
  real(kind=b8) ::  x, y

  real(kind=b8) :: angle

  if( x .gt. 0 ) then
    angle = atan(y/x) + M_PI_2
  else
    if( x .eq. 0 ) then
      if( y .gt. 0 ) then
        angle = M_PI_2
      else
        if( y .lt. 0 ) then
          angle = M_3PI_2
        else
          angle = -1
        end if
      end if
    else
      angle = atan(y/x) + M_3PI_2
    end if
  end if

end function eval_polar_angle

end module vpi_coordinate_utils
