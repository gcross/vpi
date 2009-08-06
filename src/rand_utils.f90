!@+leo-ver=4-thin
!@+node:gcross.20090805153643.1852:@thin rand_utils.f90
!@@language fortran90
module rand_utils
  use timers
  implicit none

contains 

subroutine ru_init_seed_file(pid)
  integer :: pid
  ! Random number seed
  integer :: seed(150)
  integer :: seedsize
  integer :: i
  double precision time0
  logical :: iexist
  character (len=30) :: seed_in, seed_out

  write(seed_in,"(a10,i4.4)") ".seed.dat.",pid
  write(seed_out,"(a10,i4.4)") ".seed.dat.",pid

  inquire(file=seed_in, exist=iexist)
  if (iexist) then
     ! Seed file exists, read seed from seedfile
     open(10, file=seed_in, status="old")
     read(10, *) seedsize
     read(10, *) seed(1:seedsize)
     close(10)
     call random_seed(size=seedsize)  ! find out the size of seed
     call random_seed(put=seed(1:seedsize))
  else
     ! Seed file does not exist, generate new seed and write to file
     call random_seed(size=seedsize)  ! find out the size of seed
     print *, seedsize
     call random_seed(get=seed)  ! get the current seeds
     print *, seed
     time0 = sec_timer()
     seed(1) = simple_lcg(time0 + 100.0D0*pid)
     do i = 2,seedsize
       seed(i) = simple_lcg(1.0D0*seed(i-1))
     end do
     print *, seed
     call random_seed(put=seed)  ! initialize the random number generator
     call random_seed(size=seedsize)
     call random_seed(get=seed(1:seedsize))
     open(10, file=seed_in, status="new")
     write(10, *) seedsize
     write(10, *) seed(1:seedsize)
     close(10)
  end if
end subroutine ru_init_seed_file


! Write random number seed
subroutine ru_write_seed_file( pid )
  integer :: pid
  ! Random number seed
  integer, dimension(50) :: seed
  integer :: seedsize
  character (len=30) :: seed_in, seed_out

  write(seed_in,"(a10,i4.4)") ".seed.dat.",pid
  write(seed_out,"(a10,i4.4)") ".seed.dat.",pid

  call random_seed(size=seedsize)
  call random_seed(get=seed(1:seedsize))
  open(10, file=seed_out, status="replace")
  write(10, *) seedsize
  write(10, *) seed(1:seedsize)
  close(10)
end subroutine ru_write_seed_file

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

function simple_lcg(yin) result(yout)
  double precision :: yin,yout
  double precision, parameter :: a = 273673163155d0
  double precision, parameter :: c = 13d0
  double precision, parameter :: m = 65000d0

  yout = mod((a*yin +c),m)

end function simple_lcg


end module rand_utils
!@-node:gcross.20090805153643.1852:@thin rand_utils.f90
!@-leo
