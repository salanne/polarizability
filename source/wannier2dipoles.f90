subroutine wannier2dipoles

use commondata
use boxdata
use recipdata

implicit none

integer :: i,j,k,nw

double precision :: halflboxrecx,halflboxrecy,halflboxrecz
double precision,dimension(nwannier) :: xw,yw,zw,spreadw
double precision :: dx,dy,dz,dr,vol3,du
double precision :: dx2,dy2,dz2

double precision, dimension(4,num) :: dipx,dipy,dipz
double precision, dimension(6,nwannier) :: quad
!double precision, dimension(num) :: radius,radius2

character*2 :: atom

halflboxrecx=2.0d0/boxlenx
halflboxrecy=2.0d0/boxleny
halflboxrecz=2.0d0/boxlenz

!do i=1,nwannier
!   read(15,*)du,(quad(j,i),j=1,6)
!enddo

!inertie=0.d0
!inertiec=0.d0

read(11,*)
read(11,*)
do i=1,num
   read(11,*)atom,x(i),y(i),z(i)
   x(i)=x(i)/ao
   y(i)=y(i)/ao
   z(i)=z(i)/ao
   do while (x(i).lt.0)
      x(i)=x(i)+boxlenx
   enddo
   do while (x(i).ge.boxlenx)
      x(i)=x(i)-boxlenx
   enddo
   do while (y(i).lt.0)
      y(i)=y(i)+boxleny
   enddo
   do while (y(i).ge.boxleny)
      y(i)=y(i)-boxleny
   enddo
   do while (z(i).lt.0)
      z(i)=z(i)+boxlenz
   enddo
   do while (z(i).ge.boxlenz)
      z(i)=z(i)-boxlenz
   enddo
enddo
do i=1,nwannier
   read(11,*)atom,xw(i),yw(i),zw(i)
   xw(i)=xw(i)/ao
   yw(i)=yw(i)/ao
   zw(i)=zw(i)/ao
   do while (xw(i).lt.0)
      xw(i)=xw(i)+boxlenx
   enddo
   do while (xw(i).ge.boxlenx)
      xw(i)=xw(i)-boxlenx
   enddo
   do while (yw(i).lt.0)
      yw(i)=yw(i)+boxleny
   enddo
   do while (yw(i).ge.boxleny)
      yw(i)=yw(i)-boxleny
   enddo
   do while (zw(i).lt.0)
      zw(i)=zw(i)+boxlenz
   enddo
   do while (zw(i).ge.boxlenz)
      zw(i)=zw(i)-boxlenz
   enddo
!   read(15,*)du,du,du,du,spreadw(i)
enddo  

dipx=0.d0
dipy=0.d0
dipz=0.d0
!radius=0.d0
!radius2=0.d0

write(6,*)' Verification:'
write(6,*)' -sans champ'
do i=1,num
   nw=0
   do k=1,nwannier
      dx=xw(k)-x(i)
      dy=yw(k)-y(i)
      dz=zw(k)-z(i)
      dx=dx-boxlenx*int(dx*halflboxrecx)
      dy=dy-boxleny*int(dy*halflboxrecy)
      dz=dz-boxlenz*int(dz*halflboxrecz)
      dr=dsqrt(dx**2+dy**2+dz**2)
      if(dr.lt.rcutwannier)then
         nw=nw+1
         dipx(1,i)=dipx(1,i)-2.0*dx
         dipy(1,i)=dipy(1,i)-2.0*dy
         dipz(1,i)=dipz(1,i)-2.0*dz
!        inertiec(1,i)=inertiec(1,i)+dx*dx
!        inertiec(2,i)=inertiec(2,i)+dx*dy
!        inertiec(3,i)=inertiec(3,i)+dx*dz
!        inertiec(4,i)=inertiec(4,i)+dy*dy
!        inertiec(5,i)=inertiec(5,i)+dy*dz
!        inertiec(6,i)=inertiec(6,i)+dz*dz
!        inertie(1,i)=inertie(1,i)+dx*dx+quad(1,k)
!        inertie(2,i)=inertie(2,i)+dx*dy+quad(2,k)
!        inertie(3,i)=inertie(3,i)+dx*dz+quad(3,k)
!        inertie(4,i)=inertie(4,i)+dy*dy+quad(4,k)
!        inertie(5,i)=inertie(5,i)+dy*dz+quad(5,k)
!        inertie(6,i)=inertie(6,i)+dz*dz+quad(6,k)
      endif
   enddo
   write(6,*)i,nw
   write(54,*)dsqrt(dipx(1,i)**2+dipy(1,i)**2+dipz(1,i)**2)
   write(55,*)dipx(1,i),dipy(1,i),dipz(1,i)

   
enddo

read(12,*)
read(12,*)
do i=1,num
   read(12,*)
enddo
do i=1,nwannier
   read(12,*)atom,xw(i),yw(i),zw(i)
   xw(i)=xw(i)/ao
   yw(i)=yw(i)/ao
   zw(i)=zw(i)/ao
   do while (xw(i).lt.0)
      xw(i)=xw(i)+boxlenx
   enddo
   do while (xw(i).ge.boxlenx)
      xw(i)=xw(i)-boxlenx
   enddo
   do while (yw(i).lt.0)
      yw(i)=yw(i)+boxleny
   enddo
   do while (yw(i).ge.boxleny)
      yw(i)=yw(i)-boxleny
   enddo
   do while (zw(i).lt.0)
      zw(i)=zw(i)+boxlenz
   enddo
   do while (zw(i).ge.boxlenz)
      zw(i)=zw(i)-boxlenz
   enddo
enddo  

write(6,*)' -champ x'
do i=1,num
   nw=0
   do k=1,nwannier
      dx=xw(k)-x(i)
      dy=yw(k)-y(i)
      dz=zw(k)-z(i)
      dx=dx-boxlenx*int(dx*halflboxrecx)
      dy=dy-boxleny*int(dy*halflboxrecy)
      dz=dz-boxlenz*int(dz*halflboxrecz)
      dr=dsqrt(dx**2+dy**2+dz**2)
      if(dr.lt.rcutwannier)then
         nw=nw+1
         dipx(2,i)=dipx(2,i)-2.0*dx
         dipy(2,i)=dipy(2,i)-2.0*dy
         dipz(2,i)=dipz(2,i)-2.0*dz
      endif
   enddo
   write(6,*)i,nw
enddo

read(13,*)
read(13,*)
do i=1,num
   read(13,*)
enddo
do i=1,nwannier
   read(13,*)atom,xw(i),yw(i),zw(i)
   xw(i)=xw(i)/ao
   yw(i)=yw(i)/ao
   zw(i)=zw(i)/ao
   do while (xw(i).lt.0)
      xw(i)=xw(i)+boxlenx
   enddo
   do while (xw(i).ge.boxlenx)
      xw(i)=xw(i)-boxlenx
   enddo
   do while (yw(i).lt.0)
      yw(i)=yw(i)+boxleny
   enddo
   do while (yw(i).ge.boxleny)
      yw(i)=yw(i)-boxleny
   enddo
   do while (zw(i).lt.0)
      zw(i)=zw(i)+boxlenz
   enddo
   do while (zw(i).ge.boxlenz)
      zw(i)=zw(i)-boxlenz
   enddo
enddo  

write(6,*)' -champ y'
do i=1,num
   nw=0
   do k=1,nwannier
      dx=xw(k)-x(i)
      dy=yw(k)-y(i)
      dz=zw(k)-z(i)
      dx=dx-boxlenx*int(dx*halflboxrecx)
      dy=dy-boxleny*int(dy*halflboxrecy)
      dz=dz-boxlenz*int(dz*halflboxrecz)
      dr=dsqrt(dx**2+dy**2+dz**2)
      if(dr.lt.rcutwannier)then
         nw=nw+1
         dipx(3,i)=dipx(3,i)-2.0*dx
         dipy(3,i)=dipy(3,i)-2.0*dy
         dipz(3,i)=dipz(3,i)-2.0*dz
      endif
   enddo
   write(6,*)i,nw
enddo

read(14,*)
read(14,*)
do i=1,num
   read(14,*)
enddo
do i=1,nwannier
   read(14,*)atom,xw(i),yw(i),zw(i)
   xw(i)=xw(i)/ao
   yw(i)=yw(i)/ao
   zw(i)=zw(i)/ao
   do while (xw(i).lt.0)
      xw(i)=xw(i)+boxlenx
   enddo
   do while (xw(i).ge.boxlenx)
      xw(i)=xw(i)-boxlenx
   enddo
   do while (yw(i).lt.0)
      yw(i)=yw(i)+boxleny
   enddo
   do while (yw(i).ge.boxleny)
      yw(i)=yw(i)-boxleny
   enddo
   do while (zw(i).lt.0)
      zw(i)=zw(i)+boxlenz
   enddo
   do while (zw(i).ge.boxlenz)
      zw(i)=zw(i)-boxlenz
   enddo
enddo  

write(6,*)' -champ z'
do i=1,num
   nw=0
   do k=1,nwannier
      dx=xw(k)-x(i)
      dy=yw(k)-y(i)
      dz=zw(k)-z(i)
      dx=dx-boxlenx*int(dx*halflboxrecx)
      dy=dy-boxleny*int(dy*halflboxrecy)
      dz=dz-boxlenz*int(dz*halflboxrecz)
      dr=dsqrt(dx**2+dy**2+dz**2)
      if(dr.lt.rcutwannier)then
         nw=nw+1
         dipx(4,i)=dipx(4,i)-2.0*dx
         dipy(4,i)=dipy(4,i)-2.0*dy
         dipz(4,i)=dipz(4,i)-2.0*dz
      endif
   enddo
   write(6,*)i,nw
enddo

do i=2,4
   do j=1,num
      dxmu(j,i-1)=dipx(i,j)-dipx(1,j)
      dymu(j,i-1)=dipy(i,j)-dipy(1,j)
      dzmu(j,i-1)=dipz(i,j)-dipz(1,j)
   enddo
enddo

!  
!
return

end subroutine
