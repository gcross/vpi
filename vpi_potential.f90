module vpi_potential
  use vpi_defines
  use vpi_xij
  use vpi_aziz
  implicit none

contains 

function vpi_Usp_NULL( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  Usp = 0.0_b8

end function vpi_Usp_NULL

function vpi_gUsp_NULL( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp = 0.0_b8

end function vpi_gUsp_NULL


function vpi_Usp_box( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp,r

  Usp = 0
!  r = sqrt(dot_product(x(slice,ip,:),x(slice,ip,:))) 
  if ( (abs(x(slice,ip,1)) > abox) .or. (abs(x(slice,ip,2)) > abox) .or. (abs(x(slice,ip,3)) > abox)  ) then
    Usp = realbignumber
  end if

end function vpi_Usp_box

function vpi_gUsp_box( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp = 0

end function vpi_gUsp_box

function vpi_Usp_cylinder( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp,r

  Usp = 0
  r = sqrt( x(slice,ip,1)**2 + x(slice,ip,2)**2 ) 
  if(r > r_cylinder) then
    Usp = realbignumber*r**2
  end if

end function vpi_Usp_cylinder

function vpi_gUsp_cylinder( x, slice, nslice, np, ndim ) result ( gUsp )
  integer :: slice, nslice, np, ndim
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  real(kind=b8), dimension ( np, ndim ) :: gUsp
  real(kind=b8), dimension ( np ) :: r

  gUsp = 0

end function vpi_gUsp_cylinder


function vpi_Usp_wlink( x, slice, ip, nslice, np, ndim ) result( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8) :: Usp

  real(kind=b8) :: term1, term2, r2

  r2 = sum(x(slice,ip,:)**2)
  term1 = exp(-(x(slice,ip,3)/a1wlink)**12)
  term2 = exp(-((x(slice,ip,1)**2+x(slice,ip,2)**2)/a2wlink**2)**6)
  Usp = vbwlink*term1*(1.0_b8-term2) + w0wlink*r2
end function vpi_Usp_wlink

function vpi_gUsp_wlink( x, ip, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: ip, slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp
  real(kind=b8) :: term1, term2, r2

  r2 = sum(x(slice,ip,:)**2)
  term1 = exp(-(x(slice,ip,3)/a1wlink)**12)
  term2 = exp(-((x(slice,ip,1)**2+x(slice,ip,2)**2)/a2wlink**2)**6)

  gUsp = 0

  gUsp(:,1) = 2.0_b8*x(slice,:,1)
  gUsp(:,2) = 2.0_b8*x(slice,:,2)
  gUsp(:,3) = 2.0_b8*lam_ho*x(slice,:,3)

end function vpi_gUsp_wlink

function vpi_Usp_HO( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  Usp = ( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2) / 2.0_b8

end function vpi_Usp_HO

function vpi_gUsp_HO( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp(:,1) = x(slice,:,1)
  gUsp(:,2) = x(slice,:,2)
  gUsp(:,3) = lam_ho*x(slice,:,3)

end function vpi_gUsp_HO

function vpi_Usp_HOz( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  Usp = lam_ho*x(slice,ip,3)**2 

end function vpi_Usp_HOz

function vpi_gUsp_HOz( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp(:,1) = 0.0_b8
  gUsp(:,2) = 0.0_b8
  gUsp(:,3) = lam_ho*x(slice,:,3)

end function vpi_gUsp_HOz

function vpi_Usp_HO_overlap( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  if(slice .gt. CSLICE) then
    Usp = domega*( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2)
  else
    Usp = ( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2)
  end if

end function vpi_Usp_HO_overlap

function vpi_gUsp_HO_overlap( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  if(slice .gt. CSLICE) then
    gUsp(:,1) = 2.0_b8*domega*x(slice,:,1)
    gUsp(:,2) = 2.0_b8*domega*x(slice,:,2)
    gUsp(:,3) = 2.0_b8*domega*lam_ho*x(slice,:,3)
  else
    gUsp(:,1) = 2.0_b8*x(slice,:,1)
    gUsp(:,2) = 2.0_b8*x(slice,:,2)
    gUsp(:,3) = 2.0_b8*lam_ho*x(slice,:,3)
  end if

end function vpi_gUsp_HO_overlap

!dimensionless 3D morse oscilator H = -0.5T +De[exp(-2a(r-r0))-2*exp(-a(r-r0))]
function vpi_Usp_Morse( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp
  real(kind=b8) :: r,t1,t2

  r = sqrt(x(slice,ip,1)**2 + x(slice,ip,2)**2 + x(slice,ip,3)**2)

  if(slice .le. CSLICE) then
    t1 = exp(-2.0_b8*p_MO_a*(r-p_MO_r0))
    t2 = exp(-p_MO_a*(r-p_MO_r0))
    Usp = p_MO_De*(t1-2.0_b8*t2)
  else
    t1 = exp(-2.0_b8*p_MO_ap*(r-p_MO_r0))
    t2 = exp(-p_MO_ap*(r-p_MO_r0))
    Usp = p_MO_Dep*(t1-2.0_b8*t2)
  end if

end function vpi_Usp_Morse

!                                   2    2    2 1/2
!-2 p_MO_De p_MO_a x exp(p_MO_a (-(x  + y  + z )    + p_MO_r0))
!
!                    2    2    2 1/2                    /   2    2    2 1/2
!    (exp(p_MO_a (-(x  + y  + z )    + p_MO_r0)) - 1)  /  (x  + y  + z )

function vpi_gUsp_Morse( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  real(kind=b8), dimension ( np ) :: r,t1,t2

  r(:) = sqrt(x(slice,:,1)**2 + x(slice,:,2)**2 + x(slice,:,3)**2)
  if(slice .le. CSLICE) then
    t1(:) = exp(-2.0_b8*p_MO_a*(r(:)-p_MO_r0))
    t2(:) = exp(-p_MO_a*(r(:)-p_MO_r0))
    gUsp(:,1) = 2.0_b8*p_MO_De*p_MO_a*(t2(:) - t1(:))*x(slice,:,1)/r(:)
    gUsp(:,2) = 2.0_b8*p_MO_De*p_MO_a*(t2(:) - t1(:))*x(slice,:,2)/r(:)
    gUsp(:,3) = 2.0_b8*p_MO_De*p_MO_a*(t2(:) - t1(:))*x(slice,:,3)/r(:)
  else
    t1(:) = exp(-2.0_b8*p_MO_ap*(r(:)-p_MO_r0))
    t2(:) = exp(-p_MO_ap*(r(:)-p_MO_r0))
    gUsp(:,1) = 2.0_b8*p_MO_Dep*p_MO_ap*(t2(:) - t1(:))*x(slice,:,1)/r(:)
    gUsp(:,2) = 2.0_b8*p_MO_Dep*p_MO_ap*(t2(:) - t1(:))*x(slice,:,2)/r(:)
    gUsp(:,3) = 2.0_b8*p_MO_Dep*p_MO_ap*(t2(:) - t1(:))*x(slice,:,3)/r(:)
  end if

end function vpi_gUsp_Morse

function vpi_Usp_Gwell( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  Usp = ( x(slice,ip,1)**2 + x(slice,ip,2)**2  + lam_ho*x(slice,ip,3)**2 ) + &
        gw_e_norm*exp(-(x(slice,ip,3)-gw_z0)**2/gw_2asq)

end function vpi_Usp_Gwell

function vpi_gUsp_Gwell( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp(:,1) = 2.0_b8*x(slice,:,1)
  gUsp(:,2) = 2.0_b8*x(slice,:,2)
  gUsp(:,3) = 2.0_b8*lam_ho*x(slice,:,3) - &
              (x(slice,:,3)-gw_z0)*gw_e_norm*exp(-(x(slice,:,3)-gw_z0)**2/gw_2asq)/gw_asq

end function vpi_gUsp_Gwell

function vpi_Usp_Dwell( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  Usp = ( x(slice,ip,1)**2 + x(slice,ip,2)**2 )/2.0_b8 + ep_dw*((x(slice,ip,3)/a_dw)**2 - 1.0_b8)**2

end function vpi_Usp_Dwell

function vpi_gUsp_Dwell( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp(:,1) = x(slice,:,1)
  gUsp(:,2) = x(slice,:,2)
  gUsp(:,3) = 4.0_b8*ep_dw*((x(slice,:,3)/a_dw)**2 - 1.0_b8)*(x(slice,:,3)/a_dw2) 

end function vpi_gUsp_Dwell

function vpi_Usp_Annulus( x, slice, ip, nslice, np, ndim ) result ( Usp )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  integer :: slice, ip, nslice, np, ndim
  real(kind=b8) :: Usp

  Usp = ( x(slice,ip,2)**2 + ap_lam*( x(slice,ip,1)**2  + x(slice,ip,3)**2) ) + &
    ap_e_norm*exp(-( x(slice,ip,1)**2 + x(slice,ip,3)**2 )/ap_2asq)

  if ( x(slice,ip,1) .gt. ab_x0 ) then 
    Usp = Usp + ab_e_norm*exp(-( x(slice,ip,3)**2 )/ab_2asq)
  end if

end function vpi_Usp_Annulus

function vpi_gUsp_Annulus( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  integer :: i

  gUsp(:,1) = 2.0_b8*ap_lam*x(slice,:,1) - &
    x(slice,:,1)*ap_e_norm*exp(-(x(slice,:,1)**2+x(slice,:,3)**2)/ap_2asq)/ap_asq 
  gUsp(:,2) = 2.0_b8*x(slice,:,2)
  gUsp(:,3) = 2.0_b8*ap_lam*x(slice,:,3) - &
    x(slice,:,3)*ap_e_norm*exp(-(x(slice,:,2)**2+x(slice,:,3)**2)/ap_2asq)/ap_asq 

  do i = 1, np
    if ( x(slice,i,1) .gt. ab_x0 ) then 
      gUsp(i,3) = gUsp(i,3) + x(slice,i,3)*ab_e_norm*exp(-(x(slice,i,3)**2)/ab_2asq)/ab_asq 
    end if
  end do

end function vpi_gUsp_Annulus

function vpi_Usp_Annulus2( x, slice, ip, nslice, np, ndim ) result ( Usp )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  integer :: slice, ip, nslice, np, np2, ndim
  real(kind=b8) :: Usp

  Usp = ( x(slice,ip,2)**2 + ap_lam*( x(slice,ip,1)**2  + x(slice,ip,3)**2) ) + &
        ap_e_norm*exp(-( x(slice,ip,1)**2 + x(slice,ip,3)**2 )/ap_2asq)

  if ( ip .gt. np-N_PARTICLE2 ) then 
    if ( x(slice,ip,1) .gt. ab_x0 ) then 
      Usp = Usp + ab_e_norm*exp(-( x(slice,ip,3)**2 )/ab_2asq)
    end if
  end if

end function vpi_Usp_Annulus2

function vpi_gUsp_Annulus2( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  integer :: i

  gUsp(:,1) = 2.0_b8*ap_lam*x(slice,:,1) - &
    x(slice,:,1)*ap_e_norm*exp(-(x(slice,:,1)**2+x(slice,:,3)**2)/ap_2asq)/ap_asq 
  gUsp(:,2) = 2.0_b8*x(slice,:,2)
  gUsp(:,3) = 2.0_b8*ap_lam*x(slice,:,3) - &
    x(slice,:,3)*ap_e_norm*exp(-(x(slice,:,2)**2+x(slice,:,3)**2)/ap_2asq)/ap_asq 

  do i = (np-N_PARTICLE2)+1, np
    if ( x(slice,i,1) .gt. ab_x0 ) then 
      gUsp(i,3) = gUsp(i,3) + x(slice,i,3)*ab_e_norm*exp(-(x(slice,i,3)**2)/ab_2asq)/ab_asq 
    end if
  end do

end function vpi_gUsp_Annulus2

function vpi_Usp_Annulus3( x, slice, ip, nslice, np, ndim ) result ( Usp )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  integer :: slice, ip, nslice, np, np2, ndim
  real(kind=b8) :: Usp

  Usp = ( x(slice,ip,2)**2 + ap_lam*( x(slice,ip,1)**2  + x(slice,ip,3)**2) ) + &
        ap_e_norm*exp(-( x(slice,ip,1)**2 + x(slice,ip,3)**2 )/ap_2asq)

  if ( x(slice,ip,1) .gt. ab_x0 ) then 
    Usp = Usp + ab_e_norm*exp(-( x(slice,ip,3)**2 )/ab_2asq)
  end if

  if ( ip .gt. np-N_PARTICLE2 ) then 
    Usp = Usp + ae3*x(slice,ip,3)**2 
  end if

end function vpi_Usp_Annulus3

function vpi_gUsp_Annulus3( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  integer :: i

  gUsp(:,1) = 2.0_b8*ap_lam*x(slice,:,1) - &
    x(slice,:,1)*ap_e_norm*exp(-(x(slice,:,1)**2 + &
    x(slice,:,3)**2)/ap_2asq)/ap_asq 
  gUsp(np-N_PARTICLE2+1:np,1) = gUsp(np-N_PARTICLE2+1:np,3) + &
    2.0_b8*ae3*x(slice,np-N_PARTICLE2+1:np,3)
  gUsp(:,2) = 2.0_b8*x(slice,:,2)
  gUsp(:,3) = 2.0_b8*ap_lam*x(slice,:,3) - &
    x(slice,:,3)*ap_e_norm*exp(-(x(slice,:,2)**2+x(slice,:,3)**2)/ap_2asq)/ap_asq 

  do i = 1, np
    if ( x(slice,i,1) .gt. ab_x0 ) then 
      gUsp(i,3) = gUsp(i,3) + x(slice,i,3)*ab_e_norm*exp(-(x(slice,i,3)**2)/ab_2asq)/ab_asq 
    end if
  end do

end function vpi_gUsp_Annulus3

function vpi_Usp_Atom( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8) :: Usp

  real(kind=b8) :: ri

  ri = dot_product(x(slice,ip,:),x(slice,ip,:))**(-0.5_b8)

  Usp = atom_qn * ri

end function vpi_Usp_Atom

function vpi_gUsp_Atom( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  real(kind=b8) :: r3i
  integer :: i

  do i = 1, np
    r3i = dot_product(x(slice,i,:),x(slice,i,:))**(-3.0_b8/2)
    gUsp(i,1) = -atom_qn*x(slice,i,1)*r3i
    gUsp(i,2) = -atom_qn*x(slice,i,2)*r3i
    gUsp(i,3) = -atom_qn*x(slice,i,3)*r3i
  end do

end function vpi_gUsp_Atom

function vpi_Usp_Atom_overlap( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8) :: Usp

  real(kind=b8) :: ri

  ri = sqrt( dot_product(x(slice,ip,:),x(slice,ip,:)) )

  if(slice .gt. CSLICE) then
    Usp = atom_qn_p * ri
  else
    Usp = atom_qn * ri
  end if

end function vpi_Usp_Atom_overlap

function vpi_gUsp_Atom_overlap( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  real(kind=b8) :: r3i
  integer :: i

  if(slice .gt. CSLICE) then
    do i = 1, np
      r3i = dot_product(x(slice,i,:),x(slice,i,:))**(-3.0/2)
      gUsp(i,1) = -atom_qn_p*x(slice,i,1)*r3i
      gUsp(i,2) = -atom_qn_p*x(slice,i,2)*r3i
      gUsp(i,3) = -atom_qn_p*x(slice,i,3)*r3i
    end do
  else
    do i = 1, np
      r3i = dot_product(x(slice,i,:),x(slice,i,:))**(-3.0/2)
      gUsp(i,1) = -atom_qn*x(slice,i,1)*r3i
      gUsp(i,2) = -atom_qn*x(slice,i,2)*r3i
      gUsp(i,3) = -atom_qn*x(slice,i,3)*r3i
    end do
  end if

end function vpi_gUsp_Atom_overlap

function vpi_Uij_Charge( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  logical :: acc_flag

  real(kind=b8) :: Uij

  acc_flag = .true.

  Uij = atom_qe * ( sum(xij2(slice,ip,1:ip-1)**(-0.5)) + sum(xij2(slice,ip,ip+1:np)**(-0.5)) )

end function vpi_Uij_Charge

function vpi_gUij_Charge( x, xij2, slice, nslice, np, ndim ) result ( gUij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, nslice, np, ndim

  real(kind = b8), dimension ( np , ndim ) :: gUij

  real(kind = b8), dimension ( ndim ) :: tg_ee
  real(kind = b8), dimension ( np , ndim ) :: g_ee
  real(kind=b8) :: rij3i
  integer :: i,j

  g_ee = 0

  do i = 1, np
    do j = i+1, np
      rij3i =  xij2(slice,i,j)**(-3.0/2)
      tg_ee(:) = ( x(slice,i,:) - x(slice,j,:) ) * rij3i
      g_ee(i,:) = g_ee(i,:) + tg_ee(:)
      g_ee(j,:) = g_ee(j,:) - tg_ee(:)
    end do
  end do

  gUij(:,1) = -atom_qe*g_ee(:,1)
  gUij(:,2) = -atom_qe*g_ee(:,2)
  gUij(:,3) = -atom_qe*g_ee(:,3)
  
end function vpi_gUij_Charge

function vpi_Uij_LJ( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  logical :: acc_flag

  real(kind=b8) :: Uij

  real(kind=b8) :: t3,t6

  acc_flag = .true.

  t6 = a_lj2**6 * (sum( xij2(slice,ip,1:ip-1)**(-6) ) + sum( xij2(slice,ip,ip+1:np)**(-6) ))
  t3 = a_lj2**3 * (sum( xij2(slice,ip,1:ip-1)**(-3) ) + sum( xij2(slice,ip,ip+1:np)**(-3) ))

  Uij = e_lj * ( t6 - t3 )/2.0_b8

end function vpi_Uij_LJ

!> f := e*((r/s)^(-12) - (r/s)^(-6))/2;
!> v := [r, theta, phi]:
!> simplify(grad(f,v,coords=spherical)[1]/r);
!                                              6      6    6
!                                         3 e s  (-2 s  + r )
!                                         -------------------
!                                                  14
!                                                 r


function vpi_gUij_LJ( x, xij2, slice, nslice, np, ndim ) result ( gUij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, nslice, np, ndim

  real(kind = b8), dimension ( np , ndim ) :: gUij

  real(kind = b8), dimension ( ndim ) :: t_gUij
  real(kind = b8) :: rij, tij
  real(kind = b8) :: gc1,gc2
  integer :: i,j

  gc1 =  -6.0_b8*a_lj2**6
  gc2  =  3.0_b8*a_lj2**3

  gUij = 0

  do i = 1, np
    do j = i+1, np
      rij =  xij2(slice,i,j)
      tij = gc1*rij**(-7) + gc2*rij**(-4)
      t_gUij(:) = ( x(slice,i,:) - x(slice,j,:) ) * tij
      gUij(i,:) = gUij(i,:) + e_lj*t_gUij(:)
      gUij(j,:) = gUij(j,:) - e_lj*t_gUij(:)
    end do
  end do
  
end function vpi_gUij_LJ

function vpi_Uij_SA( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
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

end function vpi_Uij_SA

! f := b*exp(-a*r) - c/r**6 - d/r**8;
! grad(f,r) = 
!                                            c        d
!                       -b a exp(-a r) + 6 ---- + 8 ----
!                                            7        9
!                                           r        r
function vpi_gUij_SA( x, xij2, slice, nslice, np, ndim ) result ( gUij )
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
  
end function vpi_gUij_SA

function vpi_Uij_Sc( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  logical :: acc_flag

  integer k
  real(kind=b8) :: Uij

  acc_flag = .true.

  Uij = 0

  do k = 1,ip-1
    Uij = Uij + p_sc_a * exp(-p_sc_w*xij2(slice,ip,k))
  end do
  do k = ip+1,np
    Uij = Uij + p_sc_a * exp(-p_sc_w*xij2(slice,ip,k))
  end do

end function vpi_Uij_Sc

function vpi_gUij_Sc( x, xij2, slice, nslice, np, ndim ) result ( gUij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, nslice, np, ndim

  real(kind = b8), dimension ( np , ndim ) :: gUij

  real(kind = b8), dimension ( ndim ) :: t_gUij
  real(kind = b8) :: Ui
  real(kind = b8) :: gc1,gc2
  integer :: i,j,k

  gUij = 0

  do i = 1, np

    Ui = 0
    do k = 1,i-1
      Ui = Ui + p_sc_a * exp(-p_sc_w*xij2(slice,i,k))
    end do
    do k = i+1,np
      Ui = Ui + p_sc_a * exp(-p_sc_w*xij2(slice,i,k))
    end do

    do j = i+1, np
      t_gUij(:) = ( x(slice,i,:) - x(slice,j,:) ) * Ui 
      gUij(i,:) = gUij(i,:) - 2.0*p_sc_w*t_gUij(:)
      gUij(j,:) = gUij(j,:) + 2.0*p_sc_w*t_gUij(:)
    end do
  end do
  
end function vpi_gUij_Sc

function vpi_Uij_Hs( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  logical, intent(inout):: acc_flag

  integer k
  real(kind=b8) :: Uij

  Uij = 0.0_b8

  do k=1,ip-1
    if(xij2(slice, ip, k) .le. a_hs2) then
      acc_flag = .false.
      Uij = realbignumber
      goto 111
    end if
  end do

  do k=ip+1,np
    if(xij2(slice, ip, k) .le. a_hs2) then
      acc_flag = .false.
      Uij = realbignumber
      goto 111
    end if
  end do

111 return

end function vpi_Uij_Hs

function vpi_gUij_Hs( x, xij2, slice, nslice, np, ndim ) result ( gUij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, nslice, np, ndim

  real(kind = b8), dimension ( np , ndim ) :: gUij

  gUij = 0

end function vpi_gUij_Hs

function vpi_Usp_Nwell( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  Usp = ( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2) + p_nw_vb*cos(x(slice,ip,3)*p_nw_l)

end function vpi_Usp_Nwell

function vpi_gUsp_Nwell( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp(:,1) = 2.0_b8*x(slice,:,1)
  gUsp(:,2) = 2.0_b8*x(slice,:,2)
  gUsp(:,3) = 2.0_b8*lam_ho*x(slice,:,3) - p_nw_vb*p_nw_l*sin(x(slice,:,3)*p_nw_l)

end function vpi_gUsp_Nwell

function vpi_Usp_lattice( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  Usp = 0.0_b8

  Usp = ( x(slice,ip,1)**2 + lam_ho*x(slice,ip,2)**2 + x(slice,ip,3)**2) &
       + p_lattice_vb*cos(x(slice,ip,1)*p_lattice_ax)**2 &
       + p_lattice_vb*cos(x(slice,ip,3)*p_lattice_az)**2

end function vpi_Usp_lattice

function vpi_gUsp_lattice( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp(:,1) = 2.0_b8*x(slice,:,1) - 2.0_b8*p_lattice_vb*p_lattice_ax*cos(x(slice,:,1)*p_lattice_ax)*sin(x(slice,:,1)*p_lattice_ax)
  gUsp(:,2) = 2.0_b8*lam_ho*x(slice,:,2)
  gUsp(:,3) = 2.0_b8*x(slice,:,3) - 2.0_b8*p_lattice_vb*p_lattice_ax*cos(x(slice,:,3)*p_lattice_az)*sin(x(slice,:,3)*p_lattice_az)

end function vpi_gUsp_lattice

function vpi_Usp_3dlattice( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x

  real(kind=b8) :: Usp

  Usp = 0.0_b8

  Usp = 0.5_b8*( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2) +  &
     p_lattice_vb*cos(x(slice,ip,1)*p_lattice_ax+p_lattice_phase_x)**2 + &
     p_lattice_vb*cos(x(slice,ip,2)*p_lattice_ay+p_lattice_phase_y)**2 + & 
     p_lattice_vb*cos(x(slice,ip,3)*p_lattice_az+p_lattice_phase_z)**2

end function vpi_Usp_3dlattice

function vpi_gUsp_3dlattice( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim

  real(kind=b8), dimension ( np, ndim ) :: gUsp

  gUsp(:,1) = x(slice,:,1) - &
     2.0_b8*p_lattice_vb*p_lattice_ax*cos(x(slice,:,1)*p_lattice_ax + p_lattice_phase_x)* &
     sin(x(slice,:,1)*p_lattice_ax+p_lattice_phase_x)
  gUsp(:,2) = lam_ho*x(slice,:,2) - &
     2.0_b8*p_lattice_vb*p_lattice_ay*cos(x(slice,:,2)*p_lattice_ay+p_lattice_phase_y)* &
     sin(x(slice,:,2)*p_lattice_ay+p_lattice_phase_y)
  gUsp(:,3) = x(slice,:,3) - &
     2.0_b8*p_lattice_vb*p_lattice_az*cos(x(slice,:,3)*p_lattice_az+p_lattice_phase_z)* &
     sin(x(slice,:,3)*p_lattice_az+p_lattice_phase_z)

end function vpi_gUsp_3dlattice

function vpi_Uij_dimer( x, x_rot, xij2, slice, ip, nslice, np, ndim ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim

  real(kind=b8) :: Uij

  real(kind=b8) :: h1,h2
  real(kind=b8) :: t3,t6,costh
  real(kind=b8), dimension (3) :: dr
  integer :: i

  t3 = 0d0
  do i=1,ip-1
    costh = dot_product(x_rot(slice,ip,:),x_rot(slice,i,:))
    dr = x(slice,ip,:)-x(slice,i,:)
    h1 = dot_product(dr,x_rot(slice,i ,:))
    h2 = dot_product(dr,x_rot(slice,ip,:))
    t3 = t3 + (costh-3d0*h1*h2/xij2(slice,ip,i))*xij2(slice,ip,i)**(-1.5d0)
  end do
  do i=ip+1,np
    costh = dot_product(x_rot(slice,ip,:),x_rot(slice,i,:))
    dr = x(slice,ip,:)-x(slice,i,:)
    h1 = dot_product(dr,x_rot(slice,i ,:))
    h2 = dot_product(dr,x_rot(slice,ip,:))
    t3 = t3 + (costh-3d0*h1*h2/xij2(slice,ip,i))*xij2(slice,ip,i)**(-1.5d0)
  end do

  Uij = e_dimer*t3

end function vpi_Uij_dimer

!function vpi_Uij_polarized_dimer( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
!  real(kind=b8), dimension ( nslice, np , ndim ) :: x
!  real(kind=b8), dimension ( nslice, np , np ) :: xij2
!  integer :: slice, ip, nslice, np, ndim
!  logical, intent(inout):: acc_flag
!
!  real(kind=b8) :: Uij
!
!  real(kind=b8) :: t3,t6
!
!  integer :: i
!
!  if(acc_flag == .true.) then
!    t3 = 0
!    do i=1,ip-1
!      t3 = t3 + xij2(slice,ip,i)**(-1.5)*(1.0_b8-3.0_b8*dot_product(x_rot(slice,ip,:),x_rot(slice,i,:))**2)
!    end do
!    do i=ip+1,np
!      t3 = t3 + xij2(slice,ip,i)**(-1.5)*(1.0_b8-3.0_b8*dot_product(x_rot(slice,ip,:),x_rot(slice,i,:))**2)
!    end do
!
!    Uij = e_dimer*t3 
!  end if
!
!end function vpi_Uij_polarized_dimer

function vpi_Uij_z_polarized_dimer_hs( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  logical :: acc_flag

  real(kind=b8) :: Uij
  real(kind=b8) :: t3,t6

  integer :: i,k

  do k=1,ip-1
    if(xij2(slice, ip, k) .le. a_hs2) then
      acc_flag = .false.
      Uij = realbignumber
      return
    end if
  end do

  do k=ip+1,np
    if(xij2(slice, ip, k) .le. a_hs2) then
      acc_flag = .false.
      Uij = realbignumber
      return
    end if
  end do

  if(acc_flag .eqv. .true.) then
    t3 = 0
    do i=1,ip-1
      t3 = t3 + xij2(slice,ip,i)**(-1.5)*(1.0_b8-3.0_b8*xij2(slice,ip,i)**(-1)*(x(slice,ip,3)-x(slice,i,3))**2)
    end do
    do i=ip+1,np
      t3 = t3 + xij2(slice,ip,i)**(-1.5)*(1.0_b8-3.0_b8*xij2(slice,ip,i)**(-1)*(x(slice,ip,3)-x(slice,i,3))**2)
    end do

    Uij = e_dimer*t3 
  end if
end function vpi_Uij_z_polarized_dimer_hs

function vpi_gUij_z_polarized_dimer_hs( x, xij2, sl, nslice, np, ndim ) result ( gUij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: sl, nslice, np, ndim

  real(kind = b8), dimension ( np , ndim ) :: gUij

  real(kind=b8), dimension( nslice, np, ndim ) :: tx
  real(kind=b8), dimension( nslice, np, np ) :: txij2
  real(kind=b8), dimension( ndim ) :: fhi,flo
  integer :: i, k
  logical :: acc_flag
  
  gUij = 0.0
  do i = 1, np
    tx(sl,:,:) = x(sl,:,:)
    txij2(sl,:,:) = xij2(sl,:,:)
    do k = 1, ndim
      tx(sl,i,k) = x(sl,i,k) + ntol_eps
      if(use_pbc) then
        call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
      else
        call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
      end if
      fhi(k) = vpi_Uij_z_polarized_dimer( tx, txij2, sl, i, nslice, np, ndim, acc_flag )

      tx(sl,i,k) = x(sl,i,k) - ntol_eps
      if(use_pbc) then
        call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
      else
        call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
      end if
      flo(k) =  vpi_Uij_z_polarized_dimer( tx, txij2, sl, i, nslice, np, ndim, acc_flag )
      tx(sl,i,k) = x(sl,i,k)
    end do
    gUij(i,:) = fhi(:) - flo(:)
  end do
  gUij(:,:) = gUij(:,:)/(2.0_b8*ntol_eps)

end function vpi_gUij_z_polarized_dimer_hs

function vpi_Uij_z_polarized_dimer( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  logical :: acc_flag

  real(kind=b8) :: Uij
  real(kind=b8) :: t3,tmp

  integer :: i,k

  acc_flag = .true.
  t3 = 0.0_b8
  do i=1,ip-1
    tmp =  a_lj*xij2(slice,ip,i)**(-6) &
         + xij2(slice,ip,i)**(-1.5)*(1.0_b8-3.0_b8*xij2(slice,ip,i)**(-1)*(x(slice,ip,3)-x(slice,i,3))**2)
    t3 = t3 + exp(-(xij2(slice,ip,i)/(p_pbc_L*.99))**50)*tmp
  end do

  do i=ip+1,np
    tmp =  a_lj*xij2(slice,ip,i)**(-6) &
         + xij2(slice,ip,i)**(-1.5)*(1.0_b8-3.0_b8*xij2(slice,ip,i)**(-1)*(x(slice,ip,3)-x(slice,i,3))**2)
    t3 = t3 + exp(-(xij2(slice,ip,i)/(p_pbc_L*.99))**50)*tmp
  end do

  Uij = t3 
end function vpi_Uij_z_polarized_dimer

function vpi_gUij_z_polarized_dimer( x, xij2, sl, nslice, np, ndim ) result ( gUij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: sl, nslice, np, ndim

  real(kind = b8), dimension ( np , ndim ) :: gUij

  real(kind=b8), dimension( nslice, np, ndim ) :: tx
  real(kind=b8), dimension( nslice, np, np ) :: txij2
  real(kind=b8), dimension( ndim ) :: fhi,flo
  integer :: i, k
  logical :: acc_flag
  
  gUij = 0.0
  do i = 1, np
    tx(sl,:,:) = x(sl,:,:)
    do k = 1, ndim
      txij2(sl,:,:) = xij2(sl,:,:)
      tx(sl,i,k) = x(sl,i,k) + ntol_eps
      if(use_pbc) then
        call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
      else
        call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
      end if
      fhi(k) = vpi_Uij_z_polarized_dimer( tx, txij2, sl, i, nslice, np, ndim, acc_flag )

      tx(sl,i,k) = x(sl,i,k) - ntol_eps
      if(use_pbc) then
        call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
      else
        call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
      end if
      flo(k) =  vpi_Uij_z_polarized_dimer( tx, txij2, sl, i, nslice, np, ndim, acc_flag )
      tx(sl,i,k) = x(sl,i,k)
    end do
    gUij(i,:) = fhi(:) - flo(:)
  end do
  gUij(:,:) = gUij(:,:)/(2.0_b8*ntol_eps)

end function vpi_gUij_z_polarized_dimer

function vpi_gUij_dimer( x, x_rot, xij2, slice, nslice, np, ndim ) result ( gUij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, nslice, np, ndim

  real(kind = b8), dimension ( np , ndim ) :: gUij

  gUij = 0
end function vpi_gUij_dimer

end module vpi_potential
