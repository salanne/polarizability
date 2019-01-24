subroutine dcell(hloc,bloc)
implicit none
!
! calculate the properies of a cell specified by the
! matrix hloc(3,3)
!
! the results are returned in the array bloc, with :
!
! bloc(1->3) - lengths of cell vectors
! bloc(4->6) - cosines of cell angles
! bloc(7->9) - perpendicular cell widths
! bloc(10)   - cell volume
!
! declare local variables...
!
double precision :: axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3
!
double precision, dimension(3,3), intent(in) :: hloc
!
double precision, dimension(10), intent(out) :: bloc
!
! calculate lengths of cell vectors...
!
bloc(1:3) = dsqrt(hloc(1,1:3)*hloc(1,1:3)+hloc(2,1:3)*hloc(2,1:3)+hloc(3,1:3)*hloc(3,1:3))
!
! calculate the cosines of cell angles...
!
bloc(4) = (hloc(1,1)*hloc(1,2)+hloc(2,1)*hloc(2,2)+hloc(3,1)*hloc(3,2))/(bloc(1)*bloc(2))
bloc(5) = (hloc(1,1)*hloc(1,3)+hloc(2,1)*hloc(2,3)+hloc(3,1)*hloc(3,3))/(bloc(1)*bloc(3))
bloc(6) = (hloc(1,2)*hloc(1,3)+hloc(2,2)*hloc(2,3)+hloc(3,2)*hloc(3,3))/(bloc(2)*bloc(3))
!
! calculate the cofactors of matrix hloc...
!
axb1 = hloc(2,1)*hloc(3,2)-hloc(3,1)*hloc(2,2)
axb2 = hloc(3,1)*hloc(1,2)-hloc(1,1)*hloc(3,2)
axb3 = hloc(1,1)*hloc(2,2)-hloc(2,1)*hloc(1,2)
bxc1 = hloc(2,2)*hloc(3,3)-hloc(3,2)*hloc(2,3)
bxc2 = hloc(3,2)*hloc(1,3)-hloc(1,2)*hloc(3,3)
bxc3 = hloc(1,2)*hloc(2,3)-hloc(2,2)*hloc(1,3)
cxa1 = hloc(2,3)*hloc(3,1)-hloc(2,1)*hloc(3,3)
cxa2 = hloc(1,1)*hloc(3,3)-hloc(3,1)*hloc(1,3)
cxa3 = hloc(2,1)*hloc(1,3)-hloc(1,1)*hloc(2,3)
!
! calculate the volume of the cell...
!
bloc(10) = abs(hloc(1,1)*bxc1+hloc(2,1)*bxc2+hloc(3,1)*bxc3)
!
! calculate the perpendicular widths...
!
bloc(7) = bloc(10)/dsqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
bloc(8) = bloc(10)/dsqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
bloc(9) = bloc(10)/dsqrt(axb1*axb1+axb2*axb2+axb3*axb3)
!
return
end subroutine dcell
