module vpi_aziz
  use vpi_defines
  implicit none

  ! He-He potential parameters, Aziz 1992 (Kelvin and Angstroms)
  real (kind=b8), parameter :: &
       Astar=192215.29_b8, alphastar=10.73520708_b8, betastar=-1.89296514_b8, &
       Daziz=1.4135_b8, c6aziz=1.34920045_b8, c8aziz=0.41365922_b8, &
       c10aziz=0.17078164_b8, eps=10.94_b8, rm_inv=0.3367003367_b8
contains 

!*************************************************************************
! Aziz(1992) potential for He
! HW 11/1/06
!*************************************************************************

! He-He potential, in Kelvin
!function vpi_Uij_aziz( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
function vpi_Uij_aziz(x,r2,slice,ip,nslice,np,ndim,acc_flag) result ( Uij )
!

  integer :: nslice, np, ndim, ip, j, slice
  real (kind=b8), dimension(nslice,np,ndim) :: x
  real (kind=b8), dimension ( nslice, np , np ) :: r2
  real (kind=b8) :: Uij
  real (kind=b8) :: r_inv, damp, r

  logical :: acc_flag

  acc_flag = .true.

  Uij=0.0_b8

  do j=1,ip-1
        r=sqrt(r2(slice,ip,j))
!        r= sqrt((x(slice,ip,1)-x(slice,j,1))**2+(x(slice,ip,2)-x(slice,j,2))**2+ &
!             (x(slice,ip,3)-x(slice,j,3))**2)
        r = r*rm_inv
        r_inv = 1.0_b8/r

        if (r < Daziz) then
           damp = exp(-((Daziz*r_inv-1.0_b8)**2))
        else
           damp = 1.0_b8
        end if

        Uij = Uij + eps*(Astar*exp(-alphastar*r+betastar*r*r) - &
             damp*(c6aziz*(r_inv**6) + c8aziz*(r_inv**8) + &
             c10aziz*(r_inv**10)))

!  write(6,*) 'aziz', sqrt(r2(slice,ip,j)), eps*(Astar*exp(-alphastar*r+betastar*r*r) - &
!             damp*(c6aziz*(r_inv**6) + c8aziz*(r_inv**8) + &
!             c10aziz*(r_inv**10))), Uij

  enddo

     do j=ip+1,np 
        r=sqrt(r2(slice,ip,j))
!        r= sqrt((x(slice,ip,1)-x(slice,j,1))**2+(x(slice,ip,2)-x(slice,j,2))**2+ &
!             (x(slice,ip,3)-x(slice,j,3))**2)
        r = r*rm_inv
        r_inv = 1.0_b8/r
        
        if (r < Daziz) then
           damp = exp(-((Daziz*r_inv-1.0_b8)**2))
        else
           damp = 1.0_b8
        end if

        Uij = Uij + eps*(Astar*exp(-alphastar*r+betastar*r*r) - &
             damp*(c6aziz*(r_inv**6) + c8aziz*(r_inv**8) + &
             c10aziz*(r_inv**10)))


!  write(6,*) 'aziz', sqrt(r2(slice,ip,j)), eps*(Astar*exp(-alphastar*r+betastar*r*r) - &
!             damp*(c6aziz*(r_inv**6) + c8aziz*(r_inv**8) + &
!             c10aziz*(r_inv**10))), Uij
     end do

     Uij=Uij/2.0_b8


end function vpi_Uij_aziz


function vpi_gUij_aziz(x,r2,slice,nslice,np,ndim) result ( gUij )
!

  integer :: nslice, np, ndim, ip, j, slice, m
  real (kind=b8), dimension(nslice,np,ndim) :: x
  real (kind=b8), dimension(nslice,np,np) :: r2
!  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  real (kind=b8), dimension(np,ndim) :: gUij
  real (kind=b8) :: r_inv, damp, r, damp2

  logical :: acc_flag

  acc_flag = .true.


  gUij(1:np,1:ndim)=0.0_b8

  do ip=1,np
     do j=ip+1,np 
        r=sqrt(r2(slice,ip,j))
!        r= sqrt((x(slice,ip,1)-x(slice,j,1))**2+(x(slice,ip,2)-x(slice,j,2))**2+ &
!             (x(slice,ip,3)-x(slice,j,3))**2)
        r = r*rm_inv
        r_inv = 1.0_b8/r
        
        if (r < Daziz) then
           damp = exp(-((Daziz*r_inv-1.0_b8)**2))
           damp2 = damp*(2.0_b8*Daziz*r_inv**2*(Daziz*r_inv-1.0_b8))
        else
           damp = 1.0_b8
           damp2= 0.0_b8
        end if

        gUij(ip,1:3) = gUij(ip,1:3) + (eps*(Astar*(-alphastar*2.0_b8*betastar*r)* &
             exp(-alphastar*r+betastar*r*r) - &
             damp2*(c6aziz*(r_inv**6) + c8aziz*(r_inv**8) + &
             c10aziz*(r_inv**10)) + &
             damp*(6.0_b8*c6aziz*(r_inv**7)+8.0_b8*c8aziz*(r_inv**9)+10.0_b8* &
             c10aziz*(r_inv**11))))*(x(slice,ip,1:3)-x(slice,j,1:3))/sqrt(r2(slice,ip,j))
        gUij(j,1:3) = gUij(j,1:3) + (eps*(Astar*(-alphastar*2.0_b8*betastar*r)* &
             exp(-alphastar*r+betastar*r*r) - &
             damp2*(c6aziz*(r_inv**6) + c8aziz*(r_inv**8) + &
             c10aziz*(r_inv**10)) + &
             damp*(6.0_b8*c6aziz*(r_inv**7)+8.0_b8*c8aziz*(r_inv**9)+10.0_b8* &
             c10aziz*(r_inv**11))))*(x(slice,j,1:3)-x(slice,ip,1:3))/sqrt(r2(slice,ip,j))

     end do
  end do

end function vpi_gUij_aziz


end module vpi_aziz
