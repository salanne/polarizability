program main
use commondata
use boxdata
use recipdata
implicit none
!
! declare local variables...
!
integer :: i,j,k
!
! read input files and set up accordingly...
!
call readin
call setup

do k=1,nstep

   call wannier2dipoles
!
! calculate the fields...
!
   do i = 1,3
      xmu = dxmu(:,i)
      ymu = dymu(:,i)
      zmu = dzmu(:,i)
      call conjgrad
      call recipe_dippim
      call reale_dippim
      delecx(:,i) = elecx
      delecy(:,i) = elecy
      delecz(:,i) = elecz
   end do
!
! output the fields on each ion...
!
!   do i = 1,3
!      do j = 1,num
!         write (48,*) delecx(j,i),delecy(j,i),delecz(j,i)
!      end do
!   end do
!
   call polarizability
enddo

close(11)
close(12)
close(13)
close(14)
!close(15)
!close(48)
close(49)
close(50)
!close(51)
!close(52)
close(53)
close(54)
close(55)
!
stop
end program main
