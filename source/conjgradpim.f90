SUBROUTINE conjgrad

USE commondata, ONLY: num,x,y,z,dxsav,dysav,dzsav,erfc,eta,xmu,ymu,zmu
USE recipdata, ONLY: elcall,elsall,emcall,emsall,encall,ensall,kmaxx,kmaxy, &
                     kmaxz
USE boxdata, ONLY: twopiboxx,twopiboxy,twopiboxz,boxlenx,boxleny,boxlenz, &
                   halfboxxrec,halfboxyrec,halfboxzrec,hlab2

IMPLICIT NONE

DOUBLE PRECISION :: erfunc
DOUBLE PRECISION :: dxcf,dycf,dzcf,drsq,dr
INTEGER :: i,l,j,ipoint
!
! Setup cos and sin arrays for addition rule.
!
 elcall(:,0)=1.0d0
 emcall(:,0)=1.0d0
 encall(:,0)=1.0d0
 elsall(:,0)=0.0d0
 emsall(:,0)=0.0d0
 ensall(:,0)=0.0d0

 elcall(:,1)=cos(twopiboxx*x)
 emcall(:,1)=cos(twopiboxy*y)
 encall(:,1)=cos(twopiboxz*z)
 elsall(:,1)=sin(twopiboxx*x)
 emsall(:,1)=sin(twopiboxy*y)
 ensall(:,1)=sin(twopiboxz*z)
!
! Calculate all cosines and sines
!
do l=2,kmaxx
   elcall(:,l)=elcall(:,l-1)*elcall(:,1)-elsall(:,l-1)*elsall(:,1)
   elsall(:,l)=elsall(:,l-1)*elcall(:,1)+elcall(:,l-1)*elsall(:,1)
enddo   
do l=2,kmaxy
   emcall(:,l)=emcall(:,l-1)*emcall(:,1)-emsall(:,l-1)*emsall(:,1)
   emsall(:,l)=emsall(:,l-1)*emcall(:,1)+emcall(:,l-1)*emsall(:,1)
enddo   
do l=2,kmaxz
   encall(:,l)=encall(:,l-1)*encall(:,1)-ensall(:,l-1)*ensall(:,1)
   ensall(:,l)=ensall(:,l-1)*encall(:,1)+encall(:,l-1)*ensall(:,1)
enddo   

do j=2,num
   do i=1,j-1
      dxcf=x(i)-x(j)
      dycf=y(i)-y(j)
      dzcf=z(i)-z(j)
      dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
      dycf=dycf-boxleny*int(dycf*halfboxyrec)
      dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)
      dxsav(i,j)=hlab2(1,1)*dxcf+hlab2(1,2)*dycf+hlab2(1,3)*dzcf
      dysav(i,j)=hlab2(2,1)*dxcf+hlab2(2,2)*dycf+hlab2(2,3)*dzcf
      dzsav(i,j)=hlab2(3,1)*dxcf+hlab2(3,2)*dycf+hlab2(3,3)*dzcf

      drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)
      dr=dsqrt(drsq)
      erfc(i,j)=erfunc(eta*dr)
   enddo   
enddo   


return
END SUBROUTINE
