subroutine invert(hloc,hiloc)
implicit none
!
! define local variables...
!
double precision :: d,r
!
double precision, dimension(3,3), intent(in)  :: hloc
double precision, dimension(3,3), intent(out) :: hiloc
!
! calculate adjoint matrix.
!
hiloc(1,1) = hloc(2,2)*hloc(3,3)-hloc(3,2)*hloc(2,3)
hiloc(2,1) = hloc(3,1)*hloc(2,3)-hloc(2,1)*hloc(3,3)
hiloc(3,1) = hloc(2,1)*hloc(3,2)-hloc(3,1)*hloc(2,2)
hiloc(1,2) = hloc(3,2)*hloc(1,3)-hloc(1,2)*hloc(3,3)
hiloc(2,2) = hloc(1,1)*hloc(3,3)-hloc(3,1)*hloc(1,3)
hiloc(3,2) = hloc(3,1)*hloc(1,2)-hloc(1,1)*hloc(3,2)
hiloc(1,3) = hloc(1,2)*hloc(2,3)-hloc(2,2)*hloc(1,3)
hiloc(2,3) = hloc(2,1)*hloc(1,3)-hloc(1,1)*hloc(2,3)
hiloc(3,3) = hloc(1,1)*hloc(2,2)-hloc(2,1)*hloc(1,2)
!
! calculate determinant.
!
d = hloc(1,1)*hiloc(1,1)+hloc(1,2)*hiloc(2,1)+hloc(1,3)*hiloc(3,1)
r = 0.0d0
if (abs(d)>0.0d0) r = 1.0d0/d
!
! complete inverse matrix.
!
hiloc = r*hiloc
!
return
end subroutine invert
