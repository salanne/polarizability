subroutine boxreset
use commondata
use boxdata
implicit none
!
! sets/resets all box length dependent variables.
!
! declare local variables...
!
integer :: i
!
double precision, dimension(3)  :: lengths
double precision, dimension(10) :: dcellinfo
!
lengths(1) = boxlenx
lengths(2) = boxleny
lengths(3) = boxlenz
!
call dcell(h,dcellinfo)
!
cellvol = dcellinfo(10)*boxlenx*boxleny*boxlenz
halfboxx = boxlenx/2.0d0
halfboxxrec = 1.0d0/halfboxx
halfboxy = boxleny/2.0d0
halfboxyrec = 1.0d0/halfboxy
halfboxz = boxlenz/2.0d0
halfboxzrec = 1.0d0/halfboxz
twopiboxx = twopi/boxlenx
twopiboxy = twopi/boxleny
twopiboxz = twopi/boxlenz
fourpicell = fourpi/cellvol
cellvol3rec = onethird/cellvol
!
! now invert the cell matrix for use with the force calculations
! and the ewald summation...
!
call invert(h,hi)
!
do i = 1,3
   fullh(i,:) = h(i,:)*lengths(:)
end do
!
call invert(fullh,fullhi)
!
! calculate factor for reciprocal space forces...
!
fac = eightpi/cellvol
!
return
end subroutine boxreset
