subroutine readin
use commondata
use boxdata
use recipdata
implicit none
!
! declare local variables...
!
integer :: i
character*80 :: wannierf1,wannierf2,wannierf3,wannierf4,wannierf5

!
! read the runtime.inpt file...
!
open (unit=10,file='runtime.inpt',status='old')
!
read (10,*) nstep 
read (10,*) num
read (10,*) nwannier
!
read (10,*) etainpt
read (10,*) rcut
read (10,*) conv1
read (10,*) convfac
read (10,*) fldstr
read (10,*) h(1,1),h(1,2),h(1,3)
read (10,*) h(2,1),h(2,2),h(2,3)
read (10,*) h(3,1),h(3,2),h(3,3)
read (10,*) boxlenx
read (10,*) boxleny
read (10,*) boxlenz
read (10,*) rcutwannier
read (10,'(a)') wannierf1
read (10,'(a)') wannierf2
read (10,'(a)') wannierf3
read (10,'(a)') wannierf4
!read (10,'(a)') wannierf5
open(11,file=wannierf1)
open(12,file=wannierf2)
open(13,file=wannierf3)
open(14,file=wannierf4)
!open(15,file=wannierf5)
!
close (unit=10)
!
return
end subroutine readin
