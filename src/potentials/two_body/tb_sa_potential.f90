!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1686:@thin tb_sa_potential.f90
!@@language fortran90
module tb_sa_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1687:<< Imported modules >>
  use kinds
  use vpi_defines
  !@-node:gcross.20090624144408.1687:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1688:<< Variables >>
  !@-node:gcross.20090624144408.1688:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1689:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1690:init_tb_potential
  subroutine init_tb_potential ()
    write(*,*) "Using two-body S.A. potential."
  end subroutine init_tb_potential
  !@-node:gcross.20090624144408.1690:init_tb_potential
  !@+node:gcross.20090624144408.1691:Uij
  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim
    logical :: acc_flag

    real(kind=b8) :: Uij

    real(kind=b8) :: tmp, t6,t8

    integer :: j,k

    acc_flag = .true.

    Uij = 0

    do j = 1, np
      if(j == ip) cycle
      if( (p_lat_r0(ip,4) == 2) .and. (p_lat_r0(j,4) == 2) ) then
        k = 3
      else 
        if( (p_lat_r0(ip,4) == 2) .or.  (p_lat_r0(j,4) == 2) ) then
          k = 2
        else 
          k = 1
        endif
      endif
      t6 = xij2(slice,ip,j)**(-3)
      t8 = xij2(slice,ip,j)**(-4)
      tmp = exp(-p_sa_a(k)*sqrt(xij2(slice,ip,j)))
      Uij = Uij + p_sa_b(k)*tmp - p_sa_c(k)*t6 - p_sa_d(k)*t8
    end do

  end function Uij_func
  !@-node:gcross.20090624144408.1691:Uij
  !@+node:gcross.20090624144408.1692:gUij
  ! f := b*exp(-a*r) - c/r**6 - d/r**8;
  ! grad(f,r) = 
  !                                            c        d
  !                       -b a exp(-a r) + 6 ---- + 8 ----
  !                                            7        9
  !                                           r        r
  function gUij_func( x, xij2, slice, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    real(kind = b8), dimension ( ndim ) :: t_gUij
    real(kind = b8) :: rij, tij
    integer :: i,j,k

    gUij = 0

    do i = 1, np
      do j = i+1, np
        if( (p_lat_r0(i,4) == 2) .and. (p_lat_r0(j,4) == 2) ) then
            k = 3
        else 
          if( (p_lat_r0(i,4) == 2) .or.  (p_lat_r0(j,4) == 2) ) then
            k = 2
          else 
            k = 1
          endif
        endif
        rij =  sqrt(xij2(slice,i,j))
        tij = -p_sa_b(k)*p_sa_a(k)*exp(-p_sa_a(k)*rij) + 6.0_b8*p_sa_c(k)*rij**(-7) + 8.0_b8*p_sa_d(k)*rij**(-9)
        t_gUij(:) = ( x(slice,i,:) - x(slice,j,:) ) * tij
        gUij(i,:) = gUij(i,:) + t_gUij(:)
        gUij(j,:) = gUij(j,:) - t_gUij(:)
      end do
    end do

  end function gUij_func
  !@-node:gcross.20090624144408.1692:gUij
  !@-others
  !@-node:gcross.20090624144408.1689:<< Subroutines >>
  !@nl

end module tb_sa_potential
!@-node:gcross.20090624144408.1686:@thin tb_sa_potential.f90
!@-leo
