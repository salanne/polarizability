SUBROUTINE recipE_dippim

USE commondata, ONLY: num,elecx,elecy,elecz,xmu,ymu,zmu,q,x,y,z, &
                      twopi,sqrpi,pithreehalf,etasq, &
                      eta,etaconst,onethird
USE recipdata, ONLY: kmaxx,kmaxy,kmaxz,elcall,elsall,emcall,emsall,encall, &
                     ensall,rksqmax
USE boxdata, ONLY: fullhi,cellvol,fourpicell,fac

IMPLICIT NONE

INTEGER :: i,j,l,m,n,ll,mm,nn,nmin,mmin,kbin,lk,ilower,iupper
DOUBLE PRECISION :: rl,rkx1,rky1,rkz1,rm,rkx2,rky2, &
                 rkz2,rn,rkx3,rky3,rkz3,xkk,ckcs,ckss,cx,cy,cz,sx,sy,sz, &
                 akk,arl,arm,arn,clx,cly, &
                 clz,slx,sly,slz,rdotmus,rdotmuc, &
                 fac2
DOUBLE PRECISION :: sqrxkk,erfunc
DOUBLE PRECISION, DIMENSION(num) :: clmall,slmall,ckcnoqall,cksnoqall, &
                                    temp,temp2


elecx=0.d0
elecy=0.d0
elecz=0.d0



mmin=0
nmin=1

do ll=0,kmaxx
   l=iabs(ll)
   rl=dble(ll)*twopi

   rkx1=rl*fullhi(1,1)
   rky1=rl*fullhi(1,2)
   rkz1=rl*fullhi(1,3)

   do mm=mmin,kmaxy
      m=iabs(mm)
      rm=dble(mm)*twopi

      rkx2=rkx1+rm*fullhi(2,1)
      rky2=rky1+rm*fullhi(2,2)
      rkz2=rkz1+rm*fullhi(2,3)

      if(mm.ge.0) then
         clmall(:)=elcall(:,l)*emcall(:,m)-elsall(:,l)*emsall(:,m)
         slmall(:)=elsall(:,l)*emcall(:,m)+emsall(:,m)*elcall(:,l)
      else
         clmall(:)=elcall(:,l)*emcall(:,m)+elsall(:,l)*emsall(:,m)
         slmall(:)=elsall(:,l)*emcall(:,m)-emsall(:,m)*elcall(:,l)
      endif

      do nn=nmin,kmaxz
         n=iabs(nn)
         rn=dble(nn)*twopi

         rkx3=rkx2+rn*fullhi(3,1)
         rky3=rky2+rn*fullhi(3,2)
         rkz3=rkz2+rn*fullhi(3,3)

         xkk=rkx3*rkx3+rky3*rky3+rkz3*rkz3
         sqrxkk=dsqrt(xkk)

         if(xkk.le.rksqmax) then
            if(nn.ge.0) then
               ckcnoqall(:)=clmall(:)*encall(:,n)-slmall(:)*ensall(:,n)
               cksnoqall(:)=slmall(:)*encall(:,n)+clmall(:)*ensall(:,n)
            else
               ckcnoqall(:)=clmall(:)*encall(:,n)+slmall(:)*ensall(:,n)
               cksnoqall(:)=slmall(:)*encall(:,n)-clmall(:)*ensall(:,n)
            endif

!            ckcs=SUM(ckcnoqall*q)
!            ckss=SUM(cksnoqall*q)

            cx=SUM(ckcnoqall*xmu)
            cy=SUM(ckcnoqall*ymu)
            cz=SUM(ckcnoqall*zmu)

            sx=SUM(cksnoqall*xmu)
            sy=SUM(cksnoqall*ymu)
            sz=SUM(cksnoqall*zmu)


            akk=exp(etaconst*xkk)/xkk

            arl=akk*rkx3
            arm=akk*rky3
            arn=akk*rkz3

            clx=cx*rkx3
            cly=cy*rky3
            clz=cz*rkz3
            slx=sx*rkx3
            sly=sy*rky3
            slz=sz*rkz3

            rdotmuc=clx+cly+clz
            rdotmus=slx+sly+slz
!  next line is just the charge contribution to field

!            temp=cksnoqall*ckcs-ckcnoqall*ckss


!            elecx=elecx+arl*temp
!            elecy=elecy+arm*temp
!            elecz=elecz+arn*temp

            temp=-ckcnoqall*rdotmuc-cksnoqall*rdotmus

            elecx=elecx+arl*temp
            elecy=elecy+arm*temp
            elecz=elecz+arn*temp

         endif
      enddo   
      nmin=-kmaxz
   enddo   
   mmin=-kmaxy
enddo   


elecx=elecx*fac
elecy=elecy*fac
elecz=elecz*fac

fac2=4.0d0*(eta**3.0d0)/(3.0d0*sqrpi)

elecx=elecx+(fac2*xmu)
elecy=elecy+(fac2*ymu)
elecz=elecz+(fac2*zmu)

return
END SUBROUTINE
