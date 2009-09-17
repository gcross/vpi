!@+leo-ver=4-thin
!@+node:gcross.20090805153643.1852:@thin rand_utils.f95
!@@language fortran90
module rand_utils
  use timers
  implicit none

contains 

!@+others
!@+node:gcross.20090819093822.1393:init_seed
subroutine init_seed(pid)
  integer, intent(in) :: pid
  ! Random number seed
  integer :: seed(150)
  integer :: seedsize
  integer :: i
  double precision time0

   call random_seed(size=seedsize)  ! find out the size of seed
   call random_seed(get=seed)  ! get the current seeds
   time0 = sec_timer()
   seed(1) = simple_lcg(time0 + 100.0D0*pid)
   do i = 2,seedsize
     seed(i) = simple_lcg(1.0D0*seed(i-1))
   end do
   call random_seed(put=seed)  ! initialize the random number generator
   call random_seed(size=seedsize)
   call random_seed(get=seed(1:seedsize))

end subroutine init_seed
!@-node:gcross.20090819093822.1393:init_seed
!@+node:gcross.20090819093822.1392:ru_gasdev
subroutine ru_gasdev( nu )
  double precision, dimension( : ), intent(inout):: nu

  double precision, dimension( 2 ) :: v
  double precision :: fac,rsq
  integer:: np
  integer:: i

  double precision, dimension( size(nu)+1 ):: tnu

  np = size(nu)

  do i = 1, np , 2
    rsq = 0.0
    do while ( (rsq .ge. 1) .or. (rsq .eq. 0.0) )
      call random_number( v )
      v(:) = 2.0*v(:) - 1.0
      rsq = dot_product(v,v)
    end do
    fac = sqrt(-2.0 * log(rsq) / rsq)
    tnu(i) = v(1) * fac
    tnu(i+1) = v(2) * fac
  end do

  nu(:) = tnu(1:size(nu))

end subroutine ru_gasdev
!@-node:gcross.20090819093822.1392:ru_gasdev
!@+node:gcross.20090819093822.1391:simple_lcg
pure function simple_lcg(yin) result(yout)
  double precision, intent(in) :: yin
  double precision :: yout
  double precision, parameter :: a = 273673163155d0
  double precision, parameter :: c = 13d0
  double precision, parameter :: m = 65000d0

  yout = mod((a*yin +c),m)

end function simple_lcg
!@-node:gcross.20090819093822.1391:simple_lcg
!@-others

end module rand_utils
!@-node:gcross.20090805153643.1852:@thin rand_utils.f95
!@-leo
