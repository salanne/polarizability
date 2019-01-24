subroutine polarizability
use commondata
implicit none
!
! declare local variables
!
integer :: i,j,k
!
double precision :: tempconv,tempconvin,polartot
double precision :: trinertie,trinertiec
!
double precision, dimension(3,3) :: ddip,dfld,dfldinv,alpha
!
tempconv = ao**3.0
tempconvin = ao**2.0
!
! extract the polarizability tensor for each ion...
!
!
!
do i = 1,num
   ddip(1,:) = dxmu(i,:)  
   ddip(2,:) = dymu(i,:)  
   ddip(3,:) = dzmu(i,:)  
   dfld(1,:) = delecx(i,:)  
   dfld(2,:) = delecy(i,:)  
   dfld(3,:) = delecz(i,:)  
   do j = 1,3
      dfld(j,j) = dfld(j,j)+fldstr
   end do
   call invert(dfld,dfldinv)
   alpha = matmul(ddip,dfldinv)
   write(49,'(9(2x,f10.4))')alpha(1,1),alpha(2,2),alpha(3,3), &
   alpha(1,2),alpha(1,3),alpha(2,1),alpha(2,3),alpha(3,1),alpha(3,2)
!   (inertie(j,i),j=1,6),(inertiec(j,i),j=1,6)
!   do j = 1,3
!      write (49,*) (alpha(j,k),k=1,3)
!   end do
!   do j = 1,3
!      write (50,*) (tempconv*alpha(j,k),k=1,3)
!   end do
   alpha=alpha*tempconv
!  inertie(:,i)=inertie(:,i)*tempconvin
!  inertiec(:,i)=inertiec(:,i)*tempconvin

   write(50,'(9(2x,f10.4))')alpha(1,1),alpha(2,2),alpha(3,3), &
   alpha(1,2),alpha(1,3),alpha(2,1),alpha(2,3),alpha(3,1),alpha(3,2)
!  (inertie(j,i),j=1,6),(inertiec(j,i),j=1,6)

   polartot=(alpha(1,1)+alpha(2,2)+alpha(3,3))/3.0
!  trinertie=(inertie(1,i)+inertie(4,i)+inertie(6,i))/3.0
!  trinertiec=(inertiec(1,i)+inertiec(4,i)+inertiec(6,i))/3.0
   write(53,*)polartot!,trinertie,trinertiec
!  do j = 1,3
!     write (51,*) (ddip(j,k),k=1,3)
!  end do
!  do j = 1,3
!     write (52,*) (dfld(j,k),k=1,3)
!  end do
end do
!
return
end subroutine polarizability
