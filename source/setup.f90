subroutine setup
use commondata
use recipdata
use boxdata
implicit none
!
! declare local variables...
!
integer :: i,j,l
double precision :: vol3
!
! allocate arrays...
!
allocate (x(num),y(num),z(num))
allocate (xmu(num),ymu(num),zmu(num))
allocate (dxmu(num,3),dymu(num,3),dzmu(num,3))
allocate (elecx(num),elecy(num),elecz(num))
allocate (delecx(num,3),delecy(num,3),delecz(num,3))
allocate (erfc(num,num))
allocate (dxsav(num,num),dysav(num,num),dzsav(num,num))
!allocate (inertie(6,num),inertiec(6,num))
!
l = 0
!
!
rsqmax = rcut**2.0d0
!
eta = etainpt/(2.0d0*rcut)
etasq = eta*eta
etapi = 2.0d0*eta/sqrpi
etaconst = -1.0d0/(4.0d0*eta*eta)
!
! calculate various quantities from cell matrix h...
!
call dcell(h,b)
!
vol3 = boxlenx*boxleny*boxlenz*b(10)
hlab3 = h
hlab2 = h
h3(:,1) = h(:,1)*boxlenx/(vol3**(1.0d0/3.0d0))
h3(:,2) = h(:,2)*boxleny/(vol3**(1.0d0/3.0d0))
h3(:,3) = h(:,3)*boxlenz/(vol3**(1.0d0/3.0d0))
!
call invert(h3,hi3)
call boxreset
call dcell(fullh,b)
!
halfboxminsq = (0.5d0*dmin1(b(7),b(8),b(9)))**2
rsqmax = rcut**2
!
call kset
!
allocate (elcall(num,0:kmaxx+1),elsall(num,0:kmaxx+1))
allocate (emcall(num,0:kmaxy+1),emsall(num,0:kmaxy+1))
allocate (encall(num,0:kmaxz+1),ensall(num,0:kmaxz+1))
!

!
!open (unit=48,file='dipole_fields.out',status='unknown')
open (unit=49,file='alpha_au.out',status='unknown')
write(49,*)'# alpha(1,1) (2,2) (3,3) (1,2) (1,3) (2,1) (2,3) (3,1) &
 (3,2)'! inertie tot xx xy xz yy yz zz inertie sans spread xx xy ...'
open (unit=50,file='alpha_ang.out',status='unknown')
write(50,*)'# alpha(1,1) (2,2) (3,3) (1,2) (1,3) (2,1) (2,3) (3,1) &
 (3,2)' ! inertie tot xx xy xz yy yz zz inertie sans spread xx xy ...'
!open (unit=51,file='dipole_change.out',status='unknown')
!open (unit=52,file='total_field.out',status='unknown')
open (unit=53,file='polarizability.out',status='unknown')
write(53,*)'# trace(pol)' !,trace(inertie),trace(inertie sans spread)'
open (unit=54,file='staticdipole_au.out',status='unknown')
write(54,*)'# abs(mu) (sans champ)'
open (unit=55,file='staticdipolecomponents_au.out',status='unknown')
write(55,*)'# mux, muy, muz (sans champ)'
!open (unit=55,file='radius.out',status='unknown')
!open (unit=56,file='radius2.out',status='unknown')
!write(55,*)'# avec spread'
!write(56,*)'# sans spread'
return
end subroutine setup
