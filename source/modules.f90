module commondata
implicit none
save
!
! constants...
!
double precision, parameter :: pi=3.14159265358979323d0
double precision, parameter :: pisq=pi*pi
double precision, parameter :: sqrpi=1.7724538509055160249d0
double precision, parameter :: pithreehalf=sqrpi*pi
double precision, parameter :: twopi=2.d0*pi 
double precision, parameter :: fourpi=4.d0*pi
double precision, parameter :: eightpi=8.d0*pi
double precision, parameter :: fourpisq=4.d0*pisq
double precision, parameter :: onethird=1.d0/3.d0
double precision, parameter :: ao=0.529177d0
!
! integer variables...
!
integer :: num,nanion,ncation,nskip,nwannier,nstep
!
! double precision variables...
!
double precision :: etainpt,eta,etapi,etaconst,etasq ! ewald smearing parameters
double precision :: rsqmax,rsqmaxsr,rcut,rcutsr      ! several cutoffs
double precision :: fldstr                           ! applied field strength
double precision :: rcutwannier
!
double precision, allocatable, dimension(:) :: x,y,z                ! coordinates
double precision, allocatable, dimension(:) :: q                    ! charges
double precision, allocatable, dimension(:) :: xmu,ymu,zmu          ! dipoles
double precision, allocatable, dimension(:) :: elecx,elecy,elecz    ! electric fields
!
double precision, allocatable, dimension(:,:) :: dxmu,dymu,dzmu       ! dipole change
double precision, allocatable, dimension(:,:) :: delecx,delecy,delecz ! electric field change
double precision, allocatable, dimension(:,:) :: dxsav,dysav,dzsav    ! interatomic distances
double precision, allocatable, dimension(:,:) :: erfc                 ! error function
!
double precision, allocatable, dimension(:,:) :: inertie,inertiec
!character*80 :: wannierfile
end module commondata

module boxdata
implicit none
save
!
! double precision variables...
!
double precision :: boxlenx,boxleny,boxlenz,halfboxx,halfboxy,halfboxz
double precision :: halfboxxrec,halfboxyrec,halfboxzrec,halfboxminsq
double precision :: twopiboxx,twopiboxy,twopiboxz
double precision :: a0,b0,c0,cellvol,fourpicell,cellvol3rec,fac
!
double precision, dimension(10) :: bh,b,bee
!
double precision, dimension(3,3) :: h,hi,fullhi,fullhit,fullh ! cell matrices
double precision, dimension(3,3) :: h2,hi2,h3,hi3,vg2,vg3,hlab2,hlab3       
!
end module boxdata

module recipdata
implicit none
save
!
! integer variables...
!
integer :: kmaxx,kmaxy,kmaxz,ksqmax ! reciprocal space maximum vectors...
!
! double precision variables...
!
double precision :: conv1,convfac ! convergence factors for reciprocal space summations...
double precision :: rksqmax       ! reciprocal space cutoff
!
double precision, allocatable, dimension(:,:) :: elcall,emcall,encall,elsall,emsall,ensall
!
end module recipdata
