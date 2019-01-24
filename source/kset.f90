subroutine kset
use commondata
use boxdata
use recipdata
implicit none
!
! declare local variables...
!
integer :: kmax
!
double precision :: rkmin
!
fullhit = transpose(fullhi)
!
call dcell(fullhit,bh)
!
rkmin = (twopi*dmin1(dabs(bh(7)),dabs(bh(8)),dabs(bh(9))))
rksqmax = -log(conv1)*4.0d0*etasq/(rkmin*rkmin)
kmax = int(dsqrt(rksqmax))
!
kmaxx=int((rkmin/(twopi*bh(7)))*float(kmax))+1
kmaxy=int((rkmin/(twopi*bh(8)))*float(kmax))+1
kmaxz=int((rkmin/(twopi*bh(9)))*float(kmax))+1
!
rksqmax = rksqmax*convfac
!
return
end subroutine kset
