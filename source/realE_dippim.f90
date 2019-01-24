SUBROUTINE realE_dippim

USE commondata, ONLY: num,q,x,y,z,xmu,ymu,zmu,elecx,elecy,elecz, &
                      twopi,rsqmax,dxsav,dysav,dzsav,etapi,etasq,erfc,  &
                      onethird
USE boxdata, ONLY: hlab2,boxlenx,boxleny,boxlenz

IMPLICIT NONE

INTEGER :: i,j,kk,m,n
DOUBLE PRECISION ::frrxdipqij,frrxdipqji, &
  frrydipqij,frrydipqji,frrzdipqij,frrzdipqji,xmuidotrij, &
  xmujdotrij,expewld
DOUBLE PRECISION :: efact,efact2, &
  chgdipengij,chgdipengji
DOUBLE PRECISION :: drsq,drsqrec,dr,dr3rec, &
  qirec,qjrec,drrec, qdotmuxij,qdotmuxji,qdotmuyij, &
  qdotmuyji,qdotmuzij,qdotmuzji

do j=2,num
!   qjrec=1.0d0/q(j)
   do i=1,j-1
!      qirec=1.0d0/q(i)

      drsq=dxsav(i,j)*dxsav(i,j)+dysav(i,j)*dysav(i,j)+dzsav(i,j)*dzsav(i,j)

      if (drsq.ge.rsqmax) CYCLE

      drsqrec=1.0d0/drsq
      dr=dsqrt(drsq)
      drrec=1.0d0/dr
      dr3rec=drrec*drsqrec

      expewld=exp(-etasq*drsq)


      xmuidotrij=xmu(i)*dxsav(i,j)+ymu(i)*dysav(i,j)+zmu(i)*dzsav(i,j)
      xmujdotrij=xmu(j)*dxsav(i,j)+ymu(j)*dysav(i,j)+zmu(j)*dzsav(i,j)

      efact=etapi*expewld*dr+erfc(i,j)
      chgdipengij=(xmuidotrij)*dr3rec
      chgdipengji=-(xmujdotrij)*dr3rec

      qdotmuxij=xmu(i)
      qdotmuxji=-xmu(j)
      qdotmuyij=ymu(i)
      qdotmuyji=-ymu(j)
      qdotmuzij=zmu(i)
      qdotmuzji=-zmu(j)
      efact2=2.0d0*etasq*drsqrec*etapi*expewld

      frrxdipqij=(qdotmuxij*dr3rec &
                 -3.0d0*dxsav(i,j)*drsqrec*chgdipengij)*efact &
                 -dxsav(i,j)*xmuidotrij*efact2
      frrxdipqji=(qdotmuxji*dr3rec &
                 -3.0d0*dxsav(i,j)*drsqrec*chgdipengji)*efact &
                 +dxsav(i,j)*xmujdotrij*efact2
      frrydipqij=(qdotmuyij*dr3rec &
                 -3.0d0*dysav(i,j)*drsqrec*chgdipengij)*efact &
                 -dysav(i,j)*xmuidotrij*efact2
      frrydipqji=(qdotmuyji*dr3rec &
                 -3.0d0*dysav(i,j)*drsqrec*chgdipengji)*efact &
                 +dysav(i,j)*xmujdotrij*efact2
      frrzdipqij=(qdotmuzij*dr3rec &
                 -3.0d0*dzsav(i,j)*drsqrec*chgdipengij)*efact &
                 -dzsav(i,j)*xmuidotrij*efact2
      frrzdipqji=(qdotmuzji*dr3rec &
                 -3.0d0*dzsav(i,j)*drsqrec*chgdipengji)*efact &
                 +dzsav(i,j)*xmujdotrij*efact2

      elecx(j)=elecx(j)-frrxdipqij
      elecy(j)=elecy(j)-frrydipqij
      elecz(j)=elecz(j)-frrzdipqij
      elecx(i)=elecx(i)+frrxdipqji
      elecy(i)=elecy(i)+frrydipqji
      elecz(i)=elecz(i)+frrzdipqji
   enddo   
enddo   

return
END SUBROUTINE
