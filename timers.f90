MODULE timers
 implicit none

 contains

!returns time in milliseconds from the beginning of the month
function day_timer() result(t)
 double precision t
 integer tv(8)

 call date_and_time(values=tv)
! tv[3] =  day
! tv[5] =  hour
! tv[6] =  minutes
! tv[7] =  seconds
! tv[8] =  millliseconds

 t = 60*(60*(24*tv(3) + tv(5)) + tv(6)) + tv(7) + real(tv(8))/1000

end function day_timer

!returns time in milliseconds from last minute
function sec_timer() result(t)
 double precision t
 integer tv(8)

 call date_and_time(values=tv)
! tv[7] =  seconds
! tv[8] =  millliseconds

 t = 1000*tv(7) + tv(8)

end function sec_timer

end module timers
