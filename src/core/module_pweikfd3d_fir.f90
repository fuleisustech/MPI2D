module module_pweikfd3d_fir
use module_array
use module_utility

implicit none
integer,allocatable,private :: mx(:),my(:),mz(:),mxwk(:),mywk(:),mzwk(:)
real,pointer,private :: s(:,:,:),t(:,:,:)
real,private :: tmin,h,hh
real,allocatable,target,private :: s_work(:,:,:),t_work(:,:,:)
integer,private :: ntwk,nx,ny,nz,ncp
real,private :: t1(4,4,4),s1(4,4,4)
integer,private :: nx1,nx2,ny1,ny2,nz1,nz2

private :: xsou, qcksrt, arr2, eband, treat1, treat2, treat3

contains

! remember s0 is (nx,ny,nz)
!subroutine eikfd3d_fir(s0,t0,zs0,ys0,xs0,nz0,ny0,nx0,ncp0,h0,vmax0)
subroutine pweikfd3d_fir(s0,t0,t_init0,pz0,nz0,ny0,nx0,h0,vmax0,vmin0)
! inclue "eik3d.pa1"
real,intent(in) :: h0, vmax0, vmin0, pz0, s0(:,:,:), t_init0(:,:)
integer,intent(in) :: nz0, ny0, nx0
real,intent(out):: t0(:,:,:)
integer,allocatable :: mxwk1(:),mywk1(:),mzwk1(:)
integer :: ntt, nt, i, j, jj1, i1, j1, k1
real    :: pz, vmax, vmin, ddx, dt, t00
real,allocatable :: temp3(:,:,:),temp2(:,:)
logical :: isPad
isPad=.false.
call copy_variable(pz0,pz,h0,h,vmax0,vmax)
call copy_variable(nz0,nz,ny0,ny,nx0,nx)
call allocate_and_initial(temp2,nx,ny)
temp2=transpose(t_init0)
t00=0.01
temp2=temp2+t00
if (pz0<=5.0*h) then
   isPad=.true.
   call allocate_and_initial(temp3,nx,ny,nz)
   call permute(s0,temp3,3,2,1)
   nx=nx+10
   ny=ny+10
   nz=nz+10
   pz=pz+5.0*h
endif   
call allocate_and_initial(s_work,t_work,nx,ny,nz)
if (isPad) then
   call Padding(temp3,nx-10,5,5,ny-10,5,5,nz-10,5,5,1.0/vmin,s_work)
   call deallocate_and_free(temp3)
else
   call permute(s0,s_work,3,2,1)
endif
s=>s_work
t=>t_work
!t=0.0
ncp=2*(nz*ny+nz*nx+ny*nx)
call allocate_and_initial(mx,my,mz,ncp)
call allocate_and_initial(mxwk,mywk,mzwk,ncp)
call allocate_and_initial(mxwk1,mywk1,mzwk1,ncp)
!s=transpose(s0)
hh=h*h
ddx=1.0
!c%%%% calculate tolerance parameter for expanding wave front:
dt=ddx*h/vmax
tmin=t00
ntwk=0
!call xsou_old(xs,ys,zs,t00)
call pw_init(isPad,pz,temp2)
call deallocate_and_free(temp2)

nt=ntwk
do i=1,nt
   mx(i)=mxwk(i)
   my(i)=mywk(i)
   mz(i)=mzwk(i)
enddo
    
!c%%% Sorting the traveltimes along the edge of timed region
call qcksrt(nt,t,mx,my,mz)
101 ntwk=0
tmin=t(mx(1),my(1),mz(1))
!write(*,*)'nt=',nt
do i=1,nt
   i1=mx(i)
   j1=my(i)
   k1=mz(i)
!c*To check whether the traveltime of this point is too big.
   if (t(i1,j1,k1).gt.tmin+dt) then
      if (nt-i.gt.100) then
         do j=i,nt
            jj1=j-i+1
            mxwk1(jj1)=mx(j)
            mywk1(jj1)=my(j)
            mzwk1(jj1)=mz(j)
 !           write(*,*)mxwk1(jj1),mywk1(jj1),mzwk1(jj1)
         enddo
      else
         do j=i,nt
            ntwk=ntwk+1
            mxwk(ntwk)=mx(j)
            mywk(ntwk)=my(j)
            mzwk(ntwk)=mz(j)
         enddo
      endif
      goto 11
   endif
!c%%%%
!  Expand the wave front 
   call eband(i)
enddo
 
11 continue
ntt=jj1
jj1=0

!c %%%  The wave front has been expaned to the outmost edge of the model:
if (ntwk+ntt.le.0) goto 102
!c
!c  Sorting traveltimes of the new edge points of the solution region
!c
call qcksrt(ntwk,t,mxwk,mywk,mzwk)
call arr2(t,nt,ntwk,ntt,mxwk,mywk,mzwk,mxwk1,&
           mywk1,mzwk1,mx,my,mz)
goto 101

!c****
!c* Turn to the real arrival times.
!c****
102 continue
t=t-t00
if (isPad) then
   call permute(t_work(6:nx-5,6:ny-5,6:nz-5),t0,3,2,1)
else
   call permute(t_work,t0,3,2,1)
endif
call deallocate_and_free(t_work,s_work)

end subroutine pweikfd3d_fir

subroutine pw_init(isPad,pz,t_init)
real,intent(in) :: pz,t_init(:,:)
logical, intent(in) :: isPad
real, allocatable :: vt_init(:)
integer, allocatable :: indx(:)
integer              :: n1,n2, i, i1, i2, ipz
! put the initial traveltime
if (isPad) then
   ipz=int(pz/h)+6
   call copy_array(t_init,t_work(6:nx-5,6:ny-5,ipz))
else
   ipz=int(pz/h)+1
   call copy_array(t_init,t_work(:,:,ipz))
endif
! sort them
n1=size(t_init,1)
n2=size(t_init,2)
call allocate_and_initial(vt_init,n1*n2)
call allocate_and_initial(indx,n1*n2)
vt_init=arr2vec(t_init,n1,n2)
call indexx(n1*n2,vt_init,indx)
ntwk=0
do i=1,n1*n2 
   ntwk=ntwk+1
   call ind2sub(n1,n2,indx(i),i1,i2)
   if (isPad) then
      mxwk(ntwk)=i1+5
      mywk(ntwk)=i2+5
   else
      mxwk(ntwk)=i1
      mywk(ntwk)=i2
   endif
   mzwk(ntwk)=ipz
enddo
call deallocate_and_free(vt_init)
call deallocate_and_free(indx)
end subroutine pw_init

subroutine xsou(xs,ys,zs,t00)
real, intent(in) :: xs, ys, zs, t00
real, allocatable :: tsou(:,:,:),vtsou(:)
integer, allocatable :: indx(:)
integer :: isx, isy, isz, ix1, ix2, iy1, iy2, iz1, iz2, ns,nxs, nys, nzs
integer :: i,i1,i2,i3
real    :: s_ave
isx=nint(xs/h)+1
isy=nint(ys/h)+1
isz=nint(zs/h)+1
ix1=max(isx-2,1)
ix2=min(isx+2,nx)
iy1=max(isy-2,1)
iy2=min(isy+2,ny)
iz1=max(isz-2,1)
iz2=min(isz+2,nz)
nxs=ix2-ix1+1
nys=iy2-iy1+1
nzs=iz2-iz1+1
ns=nxs*nys*nzs
do i3=ix1,ix2
   do i2=iy1,iy2
      do i1=iz1,iz2
         s_ave=sum(s(min(i3,isx):max(i3,isx),min(i2,isy):max(i2,isy),min(i1,isz):max(i1,isz))) &
              /real((max(i1,isz)-min(i1,isz)+1)*(max(i2,isy)-min(i2,isy)+1)*(max(i3,isx)-min(i3,isx)+1))
         t(i3,i2,i1)=sqrt(((i3-1)*h-xs)**2.0+((i2-1)*h-ys)**2.0+((i1-1)*h-zs)**2.0) * s_ave + t00
         write(*,*)i1,i2,i3,t(i3,i2,i1)
      enddo
   enddo
enddo
call allocate_and_initial(tsou,ix2-ix1+1,iy2-iy1+1,iz2-iz1+1)
call allocate_and_initial(vtsou,nxs*nys*nzs)
call allocate_and_initial(indx,nxs*nys*nzs)
call copy_array(t(ix1:ix2,iy1:iy2,iz1:iz2),tsou)
vtsou=arr2vec(tsou,nxs,nys,nzs)
call indexx(nxs*nys*nzs,vtsou,indx)
ntwk=0
do i=1,ns
   ntwk=ntwk+1
   call ind2sub(nxs,nys,nzs,indx(i),i1,i2,i3)
   mxwk(ntwk)=i1+ix1-1  
   mywk(ntwk)=i2+iy1-1  
   mzwk(ntwk)=i3+iz1-1
enddo
call deallocate_and_free(tsou)  
call deallocate_and_free(vtsou)  
call deallocate_and_free(indx)  
end subroutine xsou

!=====================================================
!c      ******************************************
!c      *                                        *
!c      * Expanding Wave Front By Timing         *
!c      * The Grid Points Nearest To The Source  *
!c      * In A Volume Of 5x5x5 (if the source is *
!c      * located on a grid) or 4x4x4 (if the    *
!c      * source is located in a cell)           *
!c      * That Surround The Source               *
!c      *                                        *
!c      ******************************************
subroutine xsou_old(xs,ys,zs,t00)
!       include 'eik3d.pa1'
real,intent(in)   :: xs, ys, zs, t00
real :: c1, c2, c3, h0, xl, s0
integer :: n1, j1, k1
integer :: i, j, k
c1=abs(xs/h-int(xs/h))
c2=abs(ys/h-int(ys/h))
c3=abs(zs/h-int(zs/h))
h0=0.001*h
n1=int(xs/h)+1 
j1=int(ys/h)+1 
k1=int(zs/h)+1
if (c1.lt.h0.and.c2.lt.h0.and.c3.lt.h0) then
   t(n1,j1,k1)=t00
   do 10 i=n1-2,n1+2
      do 10 j=j1-2,j1+2
         do 10 k=k1-2,k1+2
            xl=sqrt( real((i-n1)**2.+(j-j1)**2.+(k-k1)**2.) )
            t(i,j,k)=t(n1,j1,k1)+xl*h*(s(i,j,k)+s(n1,j1,k1))*0.5
            if (i==1.or.i==nx.or.j==1.or.j==ny.or.k==1.or.k==nz) goto 10 
            if (xl.lt.1.5) goto 10
            ntwk=ntwk+1
            mxwk(ntwk)=i
            mywk(ntwk)=j
            mzwk(ntwk)=k
10 continue
else
   s0=average(s(n1,j1,k1),s(n1+1,j1,k1),s(n1,j1+1,k1),s(n1+1,j1+1,k1),&
              s(n1,j1,k1+1),s(n1+1,j1,k1+1),s(n1,j1+1,k1+1),s(n1+1,j1+1,k1+1))
   do 20 i=n1,n1+1
      do 20 j=j1,j1+1
         do 20 k=k1,k1+1
            xl=sqrt(((i-1)*h-xs)**2.+((j-1)*h-ys)**2.+((k-1)*h-zs)**2.)
            t(i,j,k)=t00+xl*s0
            if (i==1.or.i==nx.or.j==1.or.j==ny.or.k==1.or.k==nz) goto 20
            ntwk=ntwk+1  
            mxwk(ntwk)=i
            mywk(ntwk)=j
            mzwk(ntwk)=k
20 continue
   do 30 i=n1-1,n1+2,3
      do 30 j=j1-1,j1+2,3
         do 30 k=k1-1,k1+2,3
            xl=sqrt(((i-1)*h-xs)**2.+((j-1)*h-ys)**2.+((k-1)*h-zs)**2.)
            t(i,j,k)=t00+xl*(s0+s(i,j,k))*.5
            if (i==1.or.i==nx.or.j==1.or.j==ny.or.k==1.or.k==nz) goto 30
            ntwk=ntwk+1
            mxwk(ntwk)=i
            mywk(ntwk)=j
            mzwk(ntwk)=k
30 continue
endif
end subroutine xsou_old      
 
!c      $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c      $                                       $
!c      $  Quick Sorting of The Traveltimes on  $
!c      $  The Edge of Having Timed Region      $
!c      $                                       $
!c      $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


subroutine qcksrt(n,arr,mx,my,mz)
!c******************************************
!         include 'eik3d.pa3'
integer, intent(in)    :: n
real,    intent(in)    :: arr(:,:,:)
integer, intent(inout) :: mx(:),my(:),mz(:)
integer :: jstack,l,ir,ia,ib,ic,i,j,iq,l1,l2,l3,m1,m2,m3,n1,n2,n3 
integer, parameter :: m=7, nstack=50
real, parameter :: fm=7875.0, fa=211.0, fc=1663.0, fmi=1.0/fm
real    :: fx,a
integer :: istack(nstack)
jstack=0
l=1
ir=n
fx=0.
10 continue
if (ir-l.lt.m) then
   do j=l+1,ir
      a=arr(mx(j),my(j),mz(j))
      ia=mx(j)
      ib=my(j)
      ic=mz(j)
      do 11 i=j-1,1,-1
         if (arr(mx(i),my(i),mz(i)).le.a) goto 12
         n1=mx(i+1)
         n2=my(i+1)
         n3=mz(i+1)
         mx(i+1)=mx(i)
         my(i+1)=my(i)
         mz(i+1)=mz(i) 
         mx(i)=n1
         my(i)=n2
         mz(i)=n3
11    continue
      i=0
12    continue   
      mx(i+1)=ia
      my(i+1)=ib
      mz(i+1)=ic
   enddo
   if (jstack==0) return
   ir=istack(jstack)
   l=istack(jstack-1)
   jstack=jstack-2
else
   i=l
   j=ir
   fx=mod(fx*fa+fc,fm)
   iq=l+(ir-l+1)*(fx*fmi)
   a=arr(mx(iq),my(iq),mz(iq))
   ia=mx(iq)
   ib=my(iq)
   ic=mz(iq)
   mx(iq)=mx(l)
   my(iq)=my(l)
   mz(iq)=mz(l)
   mx(l)=ia
   my(l)=ib
   mz(l)=ic
20 continue
21 if (j.gt.0) then
      if (a.lt.arr(mx(j),my(j),mz(j))) then
         j=j-1
         goto 21
      endif
   endif
   if (j.le.i) then
      mx(i)=ia
      my(i)=ib
      mz(i)=ic
      goto 30
   endif
   m1=mx(i)
   m2=my(i)
   m3=mz(i)
   mx(i)=mx(j)
   my(i)=my(j)
   mz(i)=mz(j)
   mx(j)=m1
   my(j)=m2
   mz(j)=m3
   i=i+1
22 if (i.le.n) then
      if (a.gt.arr(mx(i),my(i),mz(i))) then
         i=i+1
         goto 22
      endif
   endif
   if (j.le.i) then
      mx(j)=ia
      my(j)=ib
      mz(j)=ic
      i=j
      goto 30
   endif
   l1=mx(j)
   l2=my(j)
   l3=mz(j)
   mx(j)=mx(i)
   my(j)=my(i)
   mz(j)=mz(i)
   mx(i)=l1
   my(i)=l2
   mz(i)=l3
   j=j-1
   goto 20
30 jstack=jstack+2
   if (jstack.gt.nstack) pause 'nstack must be made larger.'
   if (ir-i.ge.i-l) then
      istack(jstack)=ir
      istack(jstack-1)=i+1
      ir=i-1
   else
      istack(jstack)=i-1
      istack(jstack-1)=l
      l=i+1
   endif
endif
goto 10
return
end subroutine qcksrt

!c     *********************************            
!c     * Sorting the unexpanded points * 
!c     * and the new points expaned    *
!c     * in the last step              *
!c     *********************************
subroutine arr2(t,nt,ntwk,ntwk1,mxwk,mywk,mzwk,mxwk1,&
         mywk1,mzwk1,mx,my,mz)
!        include 'eik3d.pa2'
real, intent(in) :: t(:,:,:)
integer,intent(in) :: ntwk, ntwk1, mxwk(:), mywk(:), &
   mzwk(:), mxwk1(:),mywk1(:), mzwk1(:)
integer, intent(inout) :: nt, mx(:), my(:), mz(:)
integer :: i1, i2, j
i1=1
i2=1
nt=0
if (ntwk1.eq.0) then
   mx(1:ntwk)=mxwk(1:ntwk)
   my(1:ntwk)=mywk(1:ntwk)
   mz(1:ntwk)=mzwk(1:ntwk)
   nt=ntwk
   goto 100
endif
30 continue
nt=nt+1
if (t(mxwk(i1),mywk(i1),mzwk(i1)).lt. &
    t(mxwk1(i2),mywk1(i2),mzwk1(i2))) then
   mx(nt)=mxwk(i1)
   my(nt)=mywk(i1)
   mz(nt)=mzwk(i1)
   i1=i1+1
   if (i1.gt.ntwk) then
      do j=i2,ntwk1
         nt=nt+1
         mx(nt)=mxwk1(j)
         my(nt)=mywk1(j)
         mz(nt)=mzwk1(j)
      enddo
      goto 100
   endif
else
   mx(nt)=mxwk1(i2)
   my(nt)=mywk1(i2)
   mz(nt)=mzwk1(i2)
   i2=i2+1
   if (i2.gt.ntwk1) then
      do j=i1,ntwk
        nt=nt+1
        mx(nt)=mxwk(j)
        my(nt)=mywk(j)
        mz(nt)=mzwk(j)
      enddo
      goto 100
   endif
endif
goto 30
100 return
end subroutine arr2
        
!c      **************************** 
!c      *  Expanding a wave front  *
!c      ****************************
subroutine eband(mi) 
!include 'eik3d.pa1'
integer, intent(in) :: mi
integer :: i,j,k,l,n1,j1,k1,ntwk0
real    :: a, b
!common/sub/ t1(4,4,4),s1(4,4,4),nx1,nx2,ny1,ny2,nz1,nz2

!c**********************************************************************
!c* This part wants to expand the wavefront boundary from the minimum  *
!c* arriving time point.                                               *
!c**********************************************************************
n1=mx(mi)
k1=mz(mi) 
j1=my(mi) 
!c**** upper half points are known      
if (t(n1-1,j1,k1)*t(n1-1,j1-1,k1)*t(n1-1,j1+1,k1)*t(n1+1,j1,k1) &
   *t(n1+1,j1-1,k1)*t(n1+1,j1+1,k1)*t(n1,j1-1,k1)*t(n1,j1+1,k1) &
   *t(n1-1,j1,k1-1)*t(n1-1,j1-1,k1-1)*t(n1-1,j1+1,k1-1)&
   *t(n1+1,j1,k1-1)*t(n1+1,j1-1,k1-1)*t(n1+1,j1+1,k1-1)&
   *t(n1,j1,k1-1)*t(n1,j1-1,k1-1)*t(n1,j1+1,k1-1).ne.0.) then   
   forall (k=k1-1:k1+1,j=j1-1:j1+1,i=n1-1:n1+1)
      s1(i-n1+2,j-j1+2,k-k1+2)=s(i,j,k)
      t1(i-n1+2,j-j1+2,k-k1+2)=t(i,j,k)
   end forall
   nx1=3-n1
   ny1=3-j1
   nz1=3-k1
   nx2=nx-n1+2
   ny2=ny-j1+2
   nz2=nz-k1+2
   ntwk0=ntwk
   call treat1
   do l=ntwk0+1,ntwk
      mxwk(l)=mxwk(l)+n1-2
      mywk(l)=mywk(l)+j1-2
      mzwk(l)=mzwk(l)+k1-2
   enddo
   forall (k=k1-1:k1+1,j=j1-1:j1+1,i=n1-1:n1+1)
      t(i,j,k)=t1(i-n1+2,j-j1+2,k-k1+2)
   end forall
   return
endif

if (t(n1-1,j1+1,k1)*t(n1-1,j1,k1)*t(n1-1,j1-1,k1)*t(n1,j1-1,k1) &
   *t(n1,j1+1,k1)*t(n1+1,j1-1,k1)*t(n1+1,j1,k1)*t(n1+1,j1+1,k1) &
   *t(n1-1,j1+1,k1+1)*t(n1-1,j1,k1+1)*t(n1-1,j1-1,k1+1)&
   *t(n1,j1-1,k1+1)*t(n1,j1+1,k1+1)*t(n1+1,j1-1,k1+1)&
   *t(n1+1,j1,k1+1)*t(n1+1,j1+1,k1+1)*t(n1,j1,k1+1).ne.0.) then 
   do i=n1-1,n1+1
      do j=j1+1,j1-1,-1
         do k=k1+1,k1-1,-1
            s1(i-n1+2,j1+2-j,k1+2-k)=s(i,j,k) 
            t1(i-n1+2,j1+2-j,k1+2-k)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-n1 
   ny1=j1+2-ny
   nz1=k1+2-nz
   nx2=nx-n1+2 
   ny2=j1+2-1
   nz2=k1+2-1
   ntwk0=ntwk
   call treat1 
   do l=ntwk0+1,ntwk
      mxwk(l)=mxwk(l)+n1-2 
      mywk(l)=j1+2-mywk(l) 
      mzwk(l)=k1+2-mzwk(l)
   enddo
   do i=n1-1,n1+1
      do j=j1+1,j1-1,-1
         do k=k1+1,k1-1,-1
            t(i,j,k)=t1(i-n1+2,j1+2-j,k1+2-k)
         enddo
      enddo
   enddo
   return
endif

if (t(n1,j1-1,k1-1)*t(n1,j1,k1-1)*t(n1,j1+1,k1-1)*t(n1,j1-1,k1)&
   *t(n1,j1+1,k1)*t(n1,j1-1,k1+1)*t(n1,j1,k1+1)*t(n1,j1+1,k1+1)&
   *t(n1-1,j1-1,k1-1)*t(n1-1,j1,k1-1)*t(n1-1,j1+1,k1-1)&
   *t(n1-1,j1-1,k1)*t(n1-1,j1+1,k1)*t(n1-1,j1-1,k1+1)&
   *t(n1-1,j1,k1+1)*t(n1-1,j1+1,k1+1)*t(n1-1,j1,k1).ne.0.) then
   do i=n1-1,n1+1
       do j=j1-1,j1+1
          do k=k1+1,k1-1,-1
             s1(k1+2-k,j-j1+2,i-n1+2)=s(i,j,k)
             t1(k1+2-k,j-j1+2,i-n1+2)=t(i,j,k)  
          enddo
       enddo
   enddo
   nx1=k1+2-nz
   ny1=3-j1
   nz1=3-n1
   nx2=k1+2-1
   ny2=ny-j1+2
   nz2=nx-n1+2
   ntwk0=ntwk
   call treat1
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=mzwk(l)+n1-2
      mywk(l)=mywk(l)+j1-2
      mzwk(l)=k1+2-a
   enddo
   do i=n1-1,n1+1
      do j=j1-1,j1+1
         do k=k1+1,k1-1,-1
            t(i,j,k)=t1(k1+2-k,j-j1+2,i-n1+2)
         enddo
      enddo
   enddo
   return
endif

if (t(n1,j1-1,k1-1)*t(n1,j1,k1-1)*t(n1,j1+1,k1-1)*t(n1,j1-1,k1)&
   *t(n1,j1+1,k1)*t(n1,j1-1,k1+1)*t(n1,j1,k1+1)*t(n1,j1+1,k1+1)&
   *t(n1+1,j1-1,k1-1)*t(n1+1,j1,k1-1)*t(n1+1,j1+1,k1-1)&
   *t(n1+1,j1-1,k1)*t(n1+1,j1+1,k1)*t(n1+1,j1-1,k1+1)&
   *t(n1+1,j1,k1+1)*t(n1+1,j1+1,k1+1)*t(n1+1,j1,k1).ne.0.) then
   do i=n1+1,n1-1,-1
      do j=j1-1,j1+1
         do k=k1-1,k1+1
            s1(k-k1+2,j-j1+2,n1+2-i)=s(i,j,k)
            t1(k-k1+2,j-j1+2,n1+2-i)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-k1
   ny1=3-j1
   nz1=n1+2-nx
   nx2=nz-k1+2
   ny2=ny-j1+2
   nz2=n1+2-1
   ntwk0=ntwk
   call treat1
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=n1+2-mzwk(l)
      mywk(l)=mywk(l)+j1-2
      mzwk(l)=a+k1-2
   enddo
   do i=n1+1,n1-1,-1
      do j=j1-1,j1+1
         do k=k1-1,k1+1
            t(i,j,k)=t1(k-k1+2,j-j1+2,n1+2-i)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1,k1-1)*t(n1,j1,k1-1)*t(n1+1,j1,k1-1)*t(n1-1,j1,k1)&
   *t(n1+1,j1,k1)*t(n1-1,j1,k1+1)*t(n1,j1,k1+1)*t(n1+1,j1,k1+1)&
   *t(n1-1,j1+1,k1-1)*t(n1,j1+1,k1-1)*t(n1+1,j1+1,k1-1)&
   *t(n1-1,j1+1,k1)*t(n1,j1+1,k1)*t(n1+1,j1+1,k1)&
   *t(n1-1,j1+1,k1+1)*t(n1,j1+1,k1+1)*t(n1+1,j1+1,k1+1).ne.0.) then
   do i=n1-1,n1+1
      do j=j1+1,j1-1,-1
         do k=k1-1,k1+1
            s1(i-n1+2,k-k1+2,j1+2-j)=s(i,j,k)
            t1(i-n1+2,k-k1+2,j1+2-j)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-n1
   ny1=3-k1
   nz1=j1+2-ny
   nx2=nx-n1+2
   ny2=nz-k1+2
   nz2=j1+2-1
   ntwk0=ntwk
   call treat1
   do l=ntwk0+1,ntwk
      b=mywk(l)
      mxwk(l)=mxwk(l)+n1-2
      mywk(l)=j1+2-mzwk(l)
      mzwk(l)=b+k1-2
   enddo
   do i=n1-1,n1+1
      do j=j1+1,j1-1,-1
         do k=k1-1,k1+1
            t(i,j,k)=t1(i-n1+2,k-k1+2,j1+2-j)
         enddo
      enddo
   enddo
   return
endif
        
if (t(n1-1,j1,k1-1)*t(n1,j1,k1-1)*t(n1+1,j1,k1-1)*t(n1-1,j1,k1)&
   *t(n1+1,j1,k1)*t(n1-1,j1,k1+1)*t(n1,j1,k1+1)*t(n1+1,j1,k1+1)&
   *t(n1-1,j1-1,k1-1)*t(n1,j1-1,k1-1)*t(n1+1,j1-1,k1-1)&
   *t(n1-1,j1-1,k1)*t(n1,j1-1,k1)*t(n1+1,j1-1,k1)&
   *t(n1-1,j1-1,k1+1)*t(n1,j1-1,k1+1)*t(n1+1,j1-1,k1+1).ne.0.) then
   do i=n1-1,n1+1
      do k=k1+1,k1-1,-1
         do j=j1-1,j1+1  
            s1(i-n1+2,k1+2-k,j-j1+2)=s(i,j,k)
            t1(i-n1+2,k1+2-k,j-j1+2)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-n1
   ny1=k1+2-nz
   nz1=3-j1
   nx2=nx-n1+2
   ny2=k1+2-1
   nz2=ny-j1+2
   ntwk0=ntwk
   call treat1
   do l=ntwk0+1,ntwk
      b=mywk(l)
      mxwk(l)=mxwk(l)+n1-2
      mywk(l)=mzwk(l)+j1-2
      mzwk(l)=k1+2-b
   enddo
   do i=n1-1,n1+1
      do k=k1+1,k1-1,-1
         do j=j1-1,j1+1
            t(i,j,k)=t1(i-n1+2,k1+2-k,j-j1+2)
         enddo
      enddo
   enddo
   return
endif

!c%%%%%% 1/4 is known
if (t(n1-1,j1+1,k1-1)*t(n1-1,j1,k1-1)*t(n1,j1+1,k1-1)&
   *t(n1,j1,k1-1)*t(n1+1,j1+1,k1-1)*t(n1+1,j1,k1-1)&
   *t(n1-1,j1+1,k1)*t(n1-1,j1,k1)*t(n1,j1+1,k1)&
   *t(n1+1,j1+1,k1)*t(n1+1,j1,k1).ne.0.)then
   forall (i=n1-1:n1+1,j=j1-1:j1+2,k=k1-2:k1+1)
      s1(i-n1+2,j-j1+2,k-k1+3)=s(i,j,k)
      t1(i-n1+2,j-j1+2,k-k1+3)=t(i,j,k)
   end forall
   nx1=3-n1
   nx2=nx-n1+2
   ny1=3-j1
   ny2=ny-j1+2
   nz1=4-k1
   nz2=nz-k1+3
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      mxwk(l)=mxwk(l)+n1-2
      mywk(l)=mywk(l)+j1-2
      mzwk(l)=mzwk(l)+k1-3
   enddo
   forall (i=n1-1:n1+1,j=j1-1:j1+2,k=k1-2:k1+1)
      t(i,j,k)=t1(i-n1+2,j-j1+2,k-k1+3)
   end forall
   return
endif

if (t(n1+1,j1-1,k1-1)*t(n1,j1-1,k1-1)*t(n1+1,j1,k1-1)&
   *t(n1,j1,k1-1)*t(n1+1,j1+1,k1-1)*t(n1,j1+1,k1-1)&
   *t(n1+1,j1-1,k1)*t(n1,j1-1,k1)*t(n1+1,j1,k1)&
   *t(n1+1,j1+1,k1)*t(n1,j1+1,k1).ne.0.) then
   do i=n1-1,n1+2
      do j=j1+1,j1-1,-1
         do k=k1-2,k1+1
            s1(j1+2-j,i-n1+2,k-k1+3)=s(i,j,k)
            t1(j1+2-j,i-n1+2,k-k1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=j1+2-ny
   ny1=3-n1
   nz1=4-k1
   nx2=j1+2-1
   ny2=nx-n1+2
   nz2=nz-k1+3
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=mywk(l)+n1-2
      mywk(l)=j1+2-a
      mzwk(l)=mzwk(l)+k1-3
   enddo
   do i=n1-1,n1+2
      do j=j1+1,j1-1,-1
         do k=k1-2,k1+1
            t(i,j,k)=t1(j1+2-j,i-n1+2,k-k1+3)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1-1,k1-1)*t(n1-1,j1,k1-1)*t(n1,j1-1,k1-1)*t(n1,j1,k1-1)&
   *t(n1+1,j1-1,k1-1)*t(n1+1,j1,k1-1)*t(n1-1,j1-1,k1)*t(n1-1,j1,k1)&
   *t(n1,j1-1,k1)*t(n1+1,j1-1,k1)*t(n1+1,j1,k1).ne.0.) then
   do i=n1+1,n1-1,-1
      do j=j1+1,j1-2,-1
         do k=k1-2,k1+1
            s1(n1+2-i,j1+2-j,k-k1+3)=s(i,j,k)
            t1(n1+2-i,j1+2-j,k-k1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=n1+2-nx
   ny1=j1+2-ny
   nz1=4-k1
   nx2=n1+2-1
   ny2=j1+2-1
   nz2=nz-k1+3
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      mxwk(l)=n1+2-mxwk(l)
      mywk(l)=j1+2-mywk(l)
      mzwk(l)=mzwk(l)+k1-3
   enddo
   do i=n1+1,n1-1,-1
      do j=j1+1,j1-2,-1
         do k=k1-2,k1+1
            t(i,j,k)=t1(n1+2-i,j1+2-j,k-k1+3)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1+1,k1-1)*t(n1-1,j1,k1-1)*t(n1-1,j1-1,k1-1)&
   *t(n1,j1+1,k1-1)*t(n1,j1,k1-1)*t(n1,j1-1,k1-1)*t(n1-1,j1+1,k1)&
   *t(n1-1,j1,k1)*t(n1-1,j1-1,k1)*t(n1,j1+1,k1)*t(n1,j1-1,k1).ne.0.) then
   do i=n1+1,n1-2,-1
      do j=j1-1,j1+1
         do k=k1-2,k1+1
            s1(j-j1+2,n1+2-i,k-k1+3)=s(i,j,k)
            t1(j-j1+2,n1+2-i,k-k1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-j1
   ny1=n1+2-nx
   nz1=4-k1
   nx2=ny-j1+2
   ny2=n1+2-1
   nz2=nz-k1+3
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=n1+2-mywk(l)
      mywk(l)=a+j1-2
      mzwk(l)=mzwk(l)+k1-3
   enddo
   do i=n1+1,n1-2,-1
      do j=j1-1,j1+1
         do k=k1-2,k1+1
            t(i,j,k)=t1(j-j1+2,n1+2-i,k-k1+3)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1+1,k1+1)*t(n1-1,j1,k1+1)*t(n1,j1+1,k1+1)*t(n1,j1,k1+1)&
   *t(n1+1,j1+1,k1+1)*t(n1+1,j1,k1+1)*t(n1-1,j1+1,k1)*t(n1-1,j1,k1)&
   *t(n1,j1+1,k1)*t(n1+1,j1+1,k1)*t(n1+1,j1,k1).ne.0.) then
   do i=n1-1,n1+1
      do j=j1+2,j1-1,-1
         do k=k1-1,k1+2
            s1(i-n1+2,k-k1+2,j1+3-j)=s(i,j,k)
            t1(i-n1+2,k-k1+2,j1+3-j)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-n1
   ny1=3-k1
   nz1=j1+3-ny
   nx2=nx-n1+2
   ny2=nz-k1+2
   nz2=j1+3-1
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      b=mywk(l)
      mxwk(l)=mxwk(l)+n1-2
      mywk(l)=j1+3-mzwk(l)
      mzwk(l)=b+k1-2
   enddo
   do i=n1-1,n1+1
      do j=j1+2,j1-1,-1
         do k=k1-1,k1+2 
            t(i,j,k)=t1(i-n1+2,k-k1+2,j1+3-j)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1-1,k1+1)*t(n1-1,j1,k1+1)*t(n1,j1-1,k1+1)&
   *t(n1,j1,k1+1)*t(n1+1,j1-1,k1+1)*t(n1+1,j1,k1+1)*t(n1-1,j1-1,k1)&
   *t(n1-1,j1,k1)*t(n1,j1-1,k1)*t(n1+1,j1-1,k1)*t(n1+1,j1,k1).ne.0.) then
   do i=n1-1,n1+1
      do j=j1+1,j1-2,-1
         do k=k1+2,k1-1,-1
            s1(i-n1+2,j1+2-j,k1+3-k)=s(i,j,k)
            t1(i-n1+2,j1+2-j,k1+3-k)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-n1
   ny1=j1+2-ny
   nz1=k1+3-nz
   nx2=nx-n1+2
   ny2=j1+2-1
   nz2=k1+3-1
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      mxwk(l)=mxwk(l)+n1-2
      mywk(l)=j1+2-mywk(l)
      mzwk(l)=k1+3-mzwk(l)
   enddo
   do i=n1-1,n1+1
      do j=j1+1,j1-2,-1
         do k=k1+2,k1-1,-1
            t(i,j,k)=t1(i-n1+2,j1+2-j,k1+3-k)
         enddo
      enddo
   enddo
   return
endif       

if (t(n1+1,j1-1,k1+1)*t(n1,j1-1,k1+1)*t(n1+1,j1,k1+1)&
   *t(n1,j1,k1+1)*t(n1+1,j1+1,k1+1)*t(n1,j1+1,k1+1)&
   *t(n1+1,j1-1,k1)*t(n1,j1-1,k1)*t(n1+1,j1,k1)&
   *t(n1+1,j1+1,k1)*t(n1,j1+1,k1).ne.0.) then 
   do j=j1+1,j1-1,-1
      do k=k1-1,k1+2
         do i=n1+2,n1-1,-1
            s1(j1+2-j,k-k1+2,n1+3-i)=s(i,j,k)
            t1(j1+2-j,k-k1+2,n1+3-i)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=j1+2-ny
   ny1=3-k1
   nz1=n1+3-nx
   nx2=j1+2-1
   ny2=nz-k1+2
   nz2=n1+3-1
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      b=mywk(l)
      mxwk(l)=n1+3-mzwk(l)
      mywk(l)=j1+2-a
      mzwk(l)=b+k1-2
   enddo
   do j=j1+1,j1-1,-1
      do k=k1-1,k1+2
         do i=n1+2,n1-1,-1
            t(i,j,k)=t1(j1+2-j,k-k1+2,n1+3-i)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1+1,k1+1)*t(n1-1,j1,k1+1)*t(n1-1,j1-1,k1+1)&
   *t(n1,j1+1,k1+1)*t(n1,j1,k1+1)*t(n1,j1-1,k1+1)&
   *t(n1-1,j1+1,k1)*t(n1-1,j1,k1)*t(n1-1,j1-1,k1)&
   *t(n1,j1+1,k1)*t(n1,j1-1,k1).ne.0.) then
   do j=j1-1,j1+1
      do k=k1-1,k1+2
         do i=n1-2,n1+1
            s1(j-j1+2,k-k1+2,i-n1+3)=s(i,j,k)
            t1(j-j1+2,k-k1+2,i-n1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-j1
   ny1=3-k1
   nz1=4-n1
   nx2=ny-j1+2
   ny2=nz-k1+2
   nz2=nx-n1+3
   ntwk0=ntwk
   call treat2 
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      b=mywk(l) 
      mxwk(l)=mzwk(l)+n1-3
      mywk(l)=a+j1-2
      mzwk(l)=b+k1-2
   enddo
   do j=j1-1,j1+1
      do k=k1-1,k1+2
         do i=n1-2,n1+1
            t(i,j,k)=t1(j-j1+2,k-k1+2,i-n1+3)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1,k1-1)*t(n1-1,j1+1,k1-1)*t(n1,j1,k1-1)*t(n1,j1+1,k1-1)&
   *t(n1-1,j1,k1)*t(n1-1,j1+1,k1)*t(n1,j1+1,k1)*t(n1-1,j1,k1+1)&
   *t(n1-1,j1+1,k1+1)*t(n1,j1,k1+1)*t(n1,j1+1,k1+1).ne.0.) then
   do k=k1-1,k1+1
      do i=n1+1,n1-2,-1
         do j=j1+2,j1-1,-1
            s1(k-k1+2,n1+2-i,j1+3-j)=s(i,j,k)
            t1(k-k1+2,n1+2-i,j1+3-j)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=3-k1
   ny1=n1+2-nx
   nz1=j1+3-ny
   nx2=nz-k1+2
   ny2=n1+2-1
   nz2=j1+3-1
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=n1+2-mywk(l)
      mywk(l)=j1+3-mzwk(l)
      mzwk(l)=a+k1-2
   enddo
   do k=k1-1,k1+1 
      do i=n1+1,n1-2,-1
         do j=j1+2,j1-1,-1
            t(i,j,k)=t1(k-k1+2,n1+2-i,j1+3-j)
         enddo
      enddo
   enddo
   return
endif

if (t(n1+1,j1,k1-1)*t(n1+1,j1+1,k1-1)*t(n1,j1,k1-1)&
   *t(n1,j1+1,k1-1)*t(n1+1,j1,k1)*t(n1+1,j1+1,k1)&
   *t(n1,j1+1,k1)*t(n1+1,j1,k1+1)*t(n1+1,j1+1,k1+1)&
   *t(n1,j1,k1+1)*t(n1,j1+1,k1+1).ne.0.) then
   do k=k1+1,k1-1,-1
      do i=n1-1,n1+2
         do j=j1+2,j1-1,-1
            s1(k1+2-k,i-n1+2,j1+3-j)=s(i,j,k)
            t1(k1+2-k,i-n1+2,j1+3-j)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=k1+2-nz
   ny1=3-n1
   nz1=j1+3-ny
   nx2=k1+2-1
   ny2=nx-n1+2
   nz2=j1+3-1
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=mywk(l)+n1-2
      mywk(l)=j1+3-mzwk(l)
      mzwk(l)=k1+2-a
   enddo
   do k=k1+1,k1-1,-1
      do i=n1-1,n1+2
         do j=j1+2,j1-1,-1
            t(i,j,k)=t1(k1+2-k,i-n1+2,j1+3-j)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1,k1-1)*t(n1-1,j1-1,k1-1)*t(n1,j1,k1-1)*t(n1,j1-1,k1-1)&
   *t(n1-1,j1,k1)*t(n1-1,j1-1,k1)*t(n1,j1-1,k1)*t(n1-1,j1,k1+1)&
   *t(n1-1,j1-1,k1+1)*t(n1,j1,k1+1)*t(n1,j1-1,k1+1).ne.0.) then
   do k=k1+1,k1-1,-1
      do i=n1+1,n1-2,-1  
         do j=j1-2,j1+1
            s1(k1+2-k,n1+2-i,j-j1+3)=s(i,j,k)
            t1(k1+2-k,n1+2-i,j-j1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=k1+2-nz
   ny1=n1+2-nx
   nz1=4-j1
   nx2=k1+2-1
   ny2=n1+2-1
   nz2=ny-j1+3
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=n1+2-mywk(l)
      mywk(l)=mzwk(l)+j1-3
      mzwk(l)=k1+2-a
   enddo
   do k=k1+1,k1-1,-1
      do i=n1+1,n1-2,-1
         do j=j1-2,j1+1
            t(i,j,k)=t1(k1+2-k,n1+2-i,j-j1+3)
         enddo
      enddo
   enddo
   return
endif

if (t(n1+1,j1,k1-1)*t(n1+1,j1-1,k1-1)*t(n1,j1,k1-1)&
   *t(n1,j1-1,k1-1)*t(n1+1,j1,k1)*t(n1+1,j1-1,k1)&
   *t(n1,j1-1,k1)*t(n1+1,j1,k1+1)*t(n1+1,j1-1,k1+1)&
   *t(n1,j1,k1+1)*t(n1,j1-1,k1+1).ne.0.) then 
   forall (k=k1-1:k1+1,i=n1-1:n1+2,j=j1-2:j1+1)
      s1(k-k1+2,i-n1+2,j-j1+3)=s(i,j,k)
      t1(k-k1+2,i-n1+2,j-j1+3)=t(i,j,k)
   end forall
   nx1=3-k1
   ny1=3-n1
   nz1=4-j1
   nx2=nz-k1+2
   ny2=nx-n1+2
   nz2=ny-j1+3
   ntwk0=ntwk
   call treat2
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=mywk(l)+n1-2
      mywk(l)=mzwk(l)+j1-3
      mzwk(l)=a+k1-2
   enddo
   forall (k=k1-1:k1+1,i=n1-1:n1+2,j=j1-2:j1+1)
      t(i,j,k)=t1(k-k1+2,i-n1+2,j-j1+3)
   end forall       
   return
endif

! 1/8 is known %%         
if (t(n1-1,j1-1,k1-1)*t(n1,j1-1,k1-1)*t(n1-1,j1,k1-1)*t(n1,j1,k1-1)&
   *t(n1-1,j1-1,k1)*t(n1,j1-1,k1)*t(n1-1,j1,k1).ne.0.) then
   forall (i=n1-2:n1+1,j=j1-2:j1+1,k=k1-2:k1+1)
      s1(i-n1+3,j-j1+3,k-k1+3)=s(i,j,k)
      t1(i-n1+3,j-j1+3,k-k1+3)=t(i,j,k)
   end forall
   nx1=4-n1
   ny1=4-j1
   nz1=4-k1
   nx2=nx-n1+3
   ny2=ny-j1+3
   nz2=nz-k1+3
   ntwk0=ntwk
   call treat3
   do l=ntwk0+1,ntwk
      mxwk(l)=mxwk(l)+n1-3
      mywk(l)=mywk(l)+j1-3
      mzwk(l)=mzwk(l)+k1-3
   enddo
   forall (i=n1-2:n1+1,j=j1-2:j1+1,k=k1-2:k1+1)
      t(i,j,k)=t1(i-n1+3,j-j1+3,k-k1+3)
   end forall
   return
endif

if (t(n1-1,j1+1,k1-1)*t(n1,j1+1,k1-1)*t(n1-1,j1,k1-1)*t(n1,j1,k1-1)&
   *t(n1-1,j1+1,k1)*t(n1,j1+1,k1)*t(n1-1,j1,k1).ne.0.) then
   do j=j1+2,j1-1,-1
      do i=n1-2,n1+1
         do k=k1-2,k1+1
            s1(j1+3-j,i-n1+3,k-k1+3)=s(i,j,k)
            t1(j1+3-j,i-n1+3,k-k1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=j1+3-ny
   ny1=4-n1
   nz1=4-k1
   nx2=j1+3-1
   ny2=nx-n1+3
   nz2=nz-k1+3
   ntwk0=ntwk
   call treat3
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=mywk(l)+n1-3
      mywk(l)=j1+3-a
      mzwk(l)=mzwk(l)+k1-3
   enddo
   do j=j1+2,j1-1,-1
      do i=n1-2,n1+1
         do k=k1-2,k1+1
            t(i,j,k)=t1(j1+3-j,i-n1+3,k-k1+3)
         enddo
      enddo
   enddo   
   return
endif
        
if (t(n1+1,j1-1,k1-1)*t(n1,j1-1,k1-1)*t(n1+1,j1,k1-1)&
   *t(n1,j1,k1-1)*t(n1+1,j1-1,k1)*t(n1,j1-1,k1)&
   *t(n1+1,j1,k1).ne.0.) then
   do j=j1-2,j1+1
      do i=n1+2,n1-1,-1
         do k=k1-2,k1+1
            s1(j-j1+3,n1+3-i,k-k1+3)=s(i,j,k)
            t1(j-j1+3,n1+3-i,k-k1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=4-j1
   ny1=n1+3-nx
   nz1=4-k1
   nx2=ny-j1+3
   ny2=n1+3-1
   nz2=nz-k1+3
   ntwk0=ntwk
   call treat3
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=n1+3-mywk(l)
      mywk(l)=a+j1-3
      mzwk(l)=mzwk(l)+k1-3
   enddo
   do j=j1-2,j1+1
      do i=n1+2,n1-1,-1
         do k=k1-2,k1+1
            t(i,j,k)=t1(j-j1+3,n1+3-i,k-k1+3)
         enddo
      enddo
   enddo
   return
endif

if (t(n1+1,j1+1,k1-1)*t(n1,j1+1,k1-1)*t(n1+1,j1,k1-1)*t(n1,j1,k1-1)&
   *t(n1+1,j1+1,k1)*t(n1,j1+1,k1)*t(n1+1,j1,k1).ne.0.) then
   do i=n1+2,n1-1,-1
      do j=j1+2,j1-1,-1
         do k=k1-2,k1+1
            s1(n1+3-i,j1+3-j,k-k1+3)=s(i,j,k)
            t1(n1+3-i,j1+3-j,k-k1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=n1+3-nx
   ny1=j1+3-ny
   nz1=4-k1
   nx2=n1+3-1
   ny2=j1+3-1
   nz2=nz-k1+3
   ntwk0=ntwk
   call treat3
   do l=ntwk0+1,ntwk
      mxwk(l)=n1+3-mxwk(l)
      mywk(l)=j1+3-mywk(l)
      mzwk(l)=mzwk(l)+k1-3
   enddo
   do i=n1+2,n1-1,-1
      do j=j1+2,j1-1,-1
         do k=k1-2,k1+1
            t(i,j,k)=t1(n1+3-i,j1+3-j,k-k1+3)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1-1,k1+1)*t(n1,j1-1,k1+1)*t(n1-1,j1,k1+1)*t(n1,j1,k1+1)&
   *t(n1-1,j1-1,k1)*t(n1,j1-1,k1)*t(n1-1,j1,k1).ne.0.) then
   do k=k1+2,k1-1,-1
      do j=j1-2,j1+1
         do i=n1-2,n1+1
            s1(k1+3-k,j-j1+3,i-n1+3)=s(i,j,k)
            t1(k1+3-k,j-j1+3,i-n1+3)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=k1+3-nz
   ny1=4-j1
   nz1=4-n1
   nx2=k1+3-1
   ny2=ny-j1+3
   nz2=nx-n1+3
   ntwk0=ntwk  
   call treat3
   do l=ntwk0+1,ntwk
      a=mxwk(l)
      mxwk(l)=mzwk(l)+n1-3
      mywk(l)=mywk(l)+j1-3
      mzwk(l)=k1+3-a
   enddo
   do k=k1+2,k1-1,-1
      do j=j1-2,j1+1
         do i=n1-2,n1+1
            t(i,j,k)=t1(k1+3-k,j-j1+3,i-n1+3)
         enddo
      enddo
   enddo
   return
endif

if (t(n1-1,j1+1,k1+1)*t(n1,j1+1,k1+1)*t(n1-1,j1,k1+1)*t(n1,j1,k1+1)&
   *t(n1-1,j1+1,k1)*t(n1,j1+1,k1)*t(n1-1,j1,k1).ne.0.) then
   do i=n1-2,n1+1
      do j=j1+2,j1-1,-1
         do k=k1+2,k1-1,-1
            s1(i-n1+3,j1+3-j,k1+3-k)=s(i,j,k)
            t1(i-n1+3,j1+3-j,k1+3-k)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=4-n1
   ny1=j1+3-ny
   nz1=k1+3-nz
   nx2=nx-n1+3
   ny2=j1+3-1
   nz2=k1+3-1
   ntwk0=ntwk  
   call treat3
   do l=ntwk0+1,ntwk
      mxwk(l)=mxwk(l)+n1-3
      mywk(l)=j1+3-mywk(l)
      mzwk(l)=k1+3-mzwk(l)
   enddo
   do i=n1-2,n1+1
      do j=j1+2,j1-1,-1
         do k=k1+2,k1-1,-1
            t(i,j,k)=t1(i-n1+3,j1+3-j,k1+3-k)
         enddo
      enddo
   enddo
   return
endif

if (t(n1+1,j1+1,k1+1)*t(n1,j1+1,k1+1)*t(n1+1,j1,k1+1)*t(n1,j1,k1+1)&
   *t(n1+1,j1+1,k1)*t(n1,j1+1,k1)*t(n1+1,j1,k1).ne.0.) then   
   do i=n1+2,n1-1,-1
      do k=k1+2,k1-1,-1
         do j=j1+2,j1-1,-1
            s1(n1+3-i,k1+3-k,j1+3-j)=s(i,j,k)
            t1(n1+3-i,k1+3-k,j1+3-j)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=n1+3-nx
   ny1=k1+3-nz
   nz1=j1+3-ny
   nx2=n1+3-1
   ny2=k1+3-1
   nz2=j1+3-1
   ntwk0=ntwk
   call treat3
   do l=ntwk0+1,ntwk
      b=mywk(l)
      mxwk(l)=n1+3-mxwk(l)
      mywk(l)=j1+3-mzwk(l)
      mzwk(l)=k1+3-b
   enddo
   do i=n1+2,n1-1,-1
      do k=k1+2,k1-1,-1
         do j=j1+2,j1-1,-1
            t(i,j,k)=t1(n1+3-i,k1+3-k,j1+3-j)
         enddo
      enddo
   enddo
   return
endif

if (t(n1+1,j1-1,k1+1)*t(n1,j1-1,k1+1)*t(n1+1,j1,k1+1)*t(n1,j1,k1+1)&
   *t(n1+1,j1-1,k1)*t(n1,j1-1,k1)*t(n1+1,j1,k1).ne.0.) then
   do i=n1+2,n1-1,-1
      do j=j1-2,j1+1
         do k=k1+2,k1-1,-1
            s1(n1+3-i,j-j1+3,k1+3-k)=s(i,j,k)
            t1(n1+3-i,j-j1+3,k1+3-k)=t(i,j,k)
         enddo
      enddo
   enddo
   nx1=n1+3-nx
   ny1=4-j1
   nz1=k1+3-nz
   nx2=n1+3-1
   ny2=ny-j1+3
   nz2=k1+3-1
   ntwk0=ntwk 
   call treat3 
   do l=ntwk0+1,ntwk
      mxwk(l)=n1+3-mxwk(l)
      mywk(l)=mywk(l)+j1-3
      mzwk(l)=k1+3-mzwk(l)
   enddo
   do i=n1+2,n1-1,-1 
      do j=j1-2,j1+1
         do k=k1+2,k1-1,-1
            t(i,j,k)=t1(n1+3-i,j-j1+3,k1+3-k)
         enddo
      enddo
   enddo
   return
endif

return
end subroutine eband         


! 3 subroutines for the cases of 1/2, 1/4 and 1/8 of a cube are known

!  1/2 of a cube (3x3x3 grid points) are known
subroutine treat1
!include 'eik3d.pa1'
!******************************************************
! upper half of the cube is known
! calculate t(1,1,3),t(1,2,3),t(1,3,3),
!           t(2,1,3),t(2,2,3),t(2,3,3),
!           t(3,1,3),t(3,2,3),t(3,3,3)
!*************************************
!common/sub/ t1(4,4,4),s1(4,4,4),nx1,nx2,ny1,ny2,nz1,nz2
real :: sx,sq
if (t1(2,2,3).eq.0.) then
   sx=average(average(s1(1,2,2),s1(2,2,2),s1(3,2,2),&
                      s1(2,1,2),s1(2,3,2)),s1(2,2,3))
   sq=hh*sx*sx-.25*((t1(1,2,2)-t1(3,2,2))**2.+(t1(2,1,2)-t1(2,3,2))**2.)
   if (sq.lt.0.) sq=hh*sx*sx
   t1(2,2,3)=t1(2,2,2)+sqrt(sq)
   if (nz2.le.3) goto 1
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=2
   mzwk(ntwk)=3
1  continue
endif

if (t1(1,2,3).eq.0.) then
   sx=average(average(s1(2,1,2),s1(2,2,2),s1(2,3,2)),&
      s1(1,2,3),s1(2,2,3),s1(1,2,2))
   sq=2.*hh*sx*sx-.5*(t1(2,1,2)-t1(2,3,2))**2.-(t1(2,2,3)-t1(1,2,2))**2.
   if(sq.lt.0.0) sq=2.*hh*sx*sx
   t1(1,2,3)=t1(2,2,2)+sqrt(sq)
   if (nx1.ge.1.or.nz2.le.3) goto 2
   ntwk=ntwk+1
   mxwk(ntwk)=1
   mywk(ntwk)=2
   mzwk(ntwk)=3
2  continue
endif

if (t1(3,2,3).eq.0.) then
   sx=average(average(s1(2,1,2),s1(2,2,2),s1(2,3,2)),&
              s1(2,2,3),s1(3,2,3),s1(3,2,2))
   sq=2.*hh*sx*sx-.5*(t1(2,1,2)-t1(2,3,2))**2.-(t1(2,2,3)-t1(3,2,2))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(3,2,3)=t1(2,2,2)+sqrt(sq)
   if (nx2.le.3.or.nz2.le.3) goto 3
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=2
   mzwk(ntwk)=3
3  continue
endif

if (t1(2,3,3).eq.0.) then
   sx=average(average(s1(1,2,2),s1(2,2,2),s1(3,2,2)),&
      s1(2,3,3),s1(2,2,3),s1(2,3,2))
   sq=2.*hh*sx*sx-.5*(t1(1,2,2)-t1(3,2,2))**2.-(t1(2,2,3)-t1(2,3,2))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(2,3,3)=t1(2,2,2)+sqrt(sq)
   if (ny2.le.3.or.nz2.le.3) goto 4
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=3
   mzwk(ntwk)=3
4  continue
endif

if (t1(2,1,3).eq.0.) then
   sx=average(average(s1(1,2,2),s1(2,2,2),s1(3,2,2)),&
              s1(2,1,3),s1(2,2,3),s1(2,1,2))
   sq=2.*hh*sx*sx-.5*(t1(1,2,2)-t1(3,2,2))**2.-(t1(2,2,3)-t1(2,1,2))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(2,1,3)=t1(2,2,2)+sqrt(sq)
   if (ny1.ge.1.or.nz2.le.3) goto 5
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=1
   mzwk(ntwk)=3
5  continue
endif

if (t1(1,3,3).eq.0.) then
   sx=average(s1(1,3,3),s1(1,3,2),s1(2,2,3),s1(2,2,2),&
              s1(2,3,3),s1(2,3,2),s1(1,2,3),s1(1,2,2))
   sq=3.*hh*sx*sx-.5*((t1(2,3,2)-t1(1,2,2))**2.+(t1(1,2,2)-t1(2,2,3))**2.&
                     +(t1(2,2,3)-t1(2,3,2))**2.+(t1(1,3,2)-t1(2,3,3))**2.&
                     +(t1(2,3,3)-t1(1,2,3))**2.+(t1(1,2,3)-t1(1,3,2))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(1,3,3)=t1(2,2,2)+sqrt(sq)
   if (nx1.ge.1.or.ny2.le.3.or.nz2.le.3) goto 6
   ntwk=ntwk+1
   mxwk(ntwk)=1
   mywk(ntwk)=3
   mzwk(ntwk)=3
6  continue
endif

if (t1(3,3,3).eq.0.) then
   sx=average(s1(3,3,3),s1(3,3,2),s1(2,2,3),s1(2,2,2),&
              s1(2,3,3),s1(2,3,2),s1(3,2,3),s1(3,2,2))
   sq=3.*hh*sx*sx-.5*( (t1(2,3,2)-t1(3,2,2))**2.+(t1(3,2,2)-t1(2,2,3))**2.&
                      +(t1(2,2,3)-t1(2,3,2))**2.+(t1(3,3,2)-t1(2,3,3))**2.&
                      +(t1(2,3,3)-t1(3,2,3))**2.+(t1(3,2,3)-t1(3,3,2))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(3,3,3)=t1(2,2,2)+sqrt(sq)
   if (nx2.le.3.or.ny2.le.3.or.nz2.le.3) goto 7
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=3
   mzwk(ntwk)=3
7  continue
endif

if (t1(3,1,3).eq.0.) then
   sx=average(s1(3,1,3),s1(3,2,3),s1(2,2,3),s1(2,1,3),&
              s1(3,1,2),s1(3,2,2),s1(2,2,2),s1(2,1,2))
   sq=3.*hh*sx*sx-.5*((t1(3,2,2)-t1(2,1,2))**2.+(t1(2,1,2)-t1(2,2,3))**2.&
                     +(t1(2,2,3)-t1(3,2,2))**2.+(t1(3,1,2)-t1(3,2,3))**2.&
                     +(t1(3,2,3)-t1(2,1,3))**2.+(t1(2,1,3)-t1(3,1,2))**2)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(3,1,3)=t1(2,2,2)+sqrt(sq)
   if (nx2.le.3.or.ny1.ge.1.or.nz2.le.3) goto 8
      ntwk=ntwk+1
      mxwk(ntwk)=3
      mywk(ntwk)=1
      mzwk(ntwk)=3
8  continue
endif

if (t1(1,1,3).eq.0.) then
   sx=average(s1(1,1,3),s1(1,2,3),s1(2,2,3),s1(2,1,3),&
              s1(1,1,2),s1(1,2,2),s1(2,2,2),s1(2,1,2))
   sq=3.*hh*sx*sx-.5*((t1(1,2,2)-t1(2,1,2))**2.+(t1(2,1,2)-t1(2,2,3))**2.&
                     +(t1(2,2,3)-t1(1,2,2))**2.+(t1(1,1,2)-t1(1,2,3))**2.&
                     +(t1(1,2,3)-t1(2,1,3))**2.+(t1(2,1,3)-t1(1,1,2))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(1,1,3)=t1(2,2,2)+sqrt(sq)
   if (nx1.ge.1.or.ny1.ge.1.or.nz2.le.3) goto 9
   ntwk=ntwk+1
   mxwk(ntwk)=1
   mywk(ntwk)=1
   mzwk(ntwk)=3
9  continue
endif
return
end subroutine treat1

subroutine treat2
!include 'eik3d.pa1'
!*********************************
! t1(1,2,1),t1(2,2,1),t1(3,2,1),t1(1,3,1),t1(2,3,1),t1(3,3,1), 
!c t1(1,2,2),t1(2,2,2),t1(3,2,2),t1(1,3,2),t1(2,3,2),t1(3,3,2) are known 
! common/sub/ t1(4,4,4),s1(4,4,4),nx1,nx2,ny1,ny2,nz1,nz2
real :: sx, sq
if (t1(2,1,2).eq.0.) then
   sx=average(average(s1(2,2,2),s1(1,2,2),s1(3,2,2),s1(2,2,3),s1(2,2,1)),s1(2,1,2))
   sq=hh*sx*sx-.25*((t1(2,2,1)-t1(2,2,3))**2.+(t1(1,2,2)-t1(3,2,2))**2.)
   if (sq.lt.0.0) sq=hh*sx*sx
   t1(2,1,2)=t1(2,2,2)+sqrt(sq)
   if (ny1.ge.1.or.nz1.ge.2) goto 1
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=1
   mzwk(ntwk)=2
1  continue
endif

if (t1(2,1,3).eq.0.) then
   sx=average(average(s1(1,2,2),s1(2,2,2),s1(3,2,2)),&
              s1(2,2,3),s1(2,1,3),s1(2,1,2))
   sq=2.*hh*sx*sx-.5*(t1(1,2,2)-t1(3,2,2))**2.-(t1(2,2,3)-t1(2,1,2))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(2,1,3)=t1(2,2,2)+sqrt(sq)
   if(ny1.ge.1)goto 2
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=1
   mzwk(ntwk)=3
2  continue
endif

if (t1(2,2,4).eq.0.) then
   sx=average(average(s1(2,2,3),s1(1,2,3),s1(3,2,3),s1(2,1,3),s1(2,3,3)),s1(2,2,4))
   sq=hh*sx*sx-.25*((t1(1,2,3)-t1(3,2,3))**2.+(t1(2,1,3)-t1(2,3,3))**2.)
   if (sq.lt.0.0) sq=hh*sx*sx
   t1(2,2,4)=t1(2,2,3)+sqrt(sq)
   if (nz2.le.4) goto 3
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=2
   mzwk(ntwk)=4
3  continue
endif

if (t1(2,3,4).eq.0.) then
   sx=average(average(s1(1,2,3),s1(2,2,3),s1(3,2,3)),s1(2,3,3),s1(2,3,4),s1(2,2,4))
   sq=2.*hh*sx*sx-.5*(t1(1,2,3)-t1(3,2,3))**2.-(t1(2,3,3)-t1(2,2,4))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(2,3,4)=t1(2,2,3)+sqrt(sq)
   if (ny2.le.3.or.nz2.le.4) goto 4
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=3
   mzwk(ntwk)=4
4  continue
endif

if (t1(2,1,4).eq.0.) then
   sx=average(average(s1(1,2,3),s1(2,2,3),s1(3,2,3)),s1(2,2,4),s1(2,1,3),s1(2,1,4))
   sq=2.*hh*sx*sx-.5*(t1(1,2,3)-t1(3,2,3))**2.-(t1(2,1,3)-t1(2,2,4))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(2,1,4)=t1(2,2,3)+sqrt(sq)
   if (ny1.ge.1.or.nz2.le.4)goto 5
      ntwk=ntwk+1
      mxwk(ntwk)=2
      mywk(ntwk)=1
      mzwk(ntwk)=4
5  continue
endif

if (t1(1,2,4).eq.0.) then
   sx=average(average(s1(2,1,3),s1(2,2,3),s1(2,3,3)),s1(2,2,4),s1(1,2,3),s1(1,2,4))
   sq=2.*hh*sx*sx-.5*(t1(2,1,3)-t1(2,3,3))**2.-(t1(1,2,3)-t1(2,2,4))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(1,2,4)=t1(2,2,3)+sqrt(sq)
   if (nx1.ge.1.or.nz2.le.4) goto 6
      ntwk=ntwk+1
      mxwk(ntwk)=1
      mywk(ntwk)=2
      mzwk(ntwk)=4
6  continue
endif

if (t1(3,2,4).eq.0.) then
   sx=average(average(s1(2,1,3),s1(2,2,3),s1(2,3,3)),s1(2,2,4),s1(3,2,3),s1(3,2,4))
   sq=2.*hh*sx*sx-.5*(t1(2,1,3)-t1(2,3,3))**2.-(t1(3,2,3)-t1(2,2,4))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(3,2,4)=t1(2,2,3)+sqrt(sq)
   if (nx2.le.3.or.nz2.le.4) goto 7
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=2
   mzwk(ntwk)=4
7  continue
endif

if (t1(1,1,3).eq.0.) then
   sx=average(average(s1(2,2,2),s1(2,2,3),s1(2,2,4)),s1(1,2,3),s1(2,1,3),s1(1,1,3))
   sq=2.*hh*sx*sx-.5*(t1(2,2,2)-t1(2,2,4))**2.-(t1(1,2,3)-t1(2,1,3))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(1,1,3)=t1(2,2,3)+sqrt(sq)
   if (nx1.ge.1.or.ny1.ge.1) goto 8
   ntwk=ntwk+1
   mxwk(ntwk)=1
   mywk(ntwk)=1
   mzwk(ntwk)=3
8  continue
endif

if (t1(3,1,3).eq.0.) then
   sx=average(average(s1(2,2,2),s1(2,2,3),s1(2,2,4)),s1(3,2,3),s1(2,1,3),s1(3,1,3))
   sq=2.*hh*sx*sx-.5*(t1(2,2,2)-t1(2,2,4))**2.-(t1(3,2,3)-t1(2,1,3))**2.
   if (sq.lt.0.0) sq=2.*hh*sx*sx
   t1(3,1,3)=t1(2,2,3)+sqrt(sq)
   if (nx2.le.3.or.ny1.ge.1) goto 9
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=1
   mzwk(ntwk)=3
9  continue
endif

if (t1(1,1,2).eq.0.) then
   sx=average(s1(1,1,2),s1(1,2,2),s1(2,1,2),s1(2,2,2),&
              s1(1,1,3),s1(1,2,3),s1(2,1,3),s1(2,2,3))
   sq=3.*hh*sx*sx-.5*((t1(1,2,3)-t1(2,1,3))**2.+(t1(2,1,3)-t1(2,2,2))**2.&
                     +(t1(2,2,2)-t1(1,2,3))**2.+(t1(1,1,3)-t1(1,2,2))**2.&
                     +(t1(1,2,2)-t1(2,1,2))**2.+(t1(2,1,2)-t1(1,1,3))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(1,1,2)=t1(2,2,3)+sqrt(sq)
   if (nx1.ge.1.or.ny1.ge.1.or.nz1.ge.2) goto 10
   ntwk=ntwk+1
   mxwk(ntwk)=1
   mywk(ntwk)=1
   mzwk(ntwk)=2
10 continue
endif

if (t1(3,1,2).eq.0.) then
   sx=average(s1(3,1,2),s1(3,2,2),s1(2,1,2),s1(2,2,2),&
              s1(3,1,3),s1(3,2,3),s1(2,1,3),s1(2,2,3))
   sq=3.*hh*sx*sx-.5*((t1(3,2,3)-t1(2,1,3))**2.+(t1(2,1,3)-t1(2,2,2))**2.&
                     +(t1(2,2,2)-t1(3,2,3))**2.+(t1(3,1,3)-t1(3,2,2))**2.&
                     +(t1(3,2,2)-t1(2,1,2))**2.+(t1(2,1,2)-t1(3,1,3))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(3,1,2)=t1(2,2,3)+sqrt(sq)
   if (nx2.le.3.or.ny1.ge.1.or.nz1.ge.2) goto 11
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=1
   mzwk(ntwk)=2
11 continue
endif

if (t1(1,1,4).eq.0.) then
   sx=average(s1(1,1,4),s1(1,2,4),s1(2,1,4),s1(2,2,4),&
              s1(1,1,3),s1(1,2,3),s1(2,1,3),s1(2,2,3))
   sq=3.*hh*sx*sx-.5*((t1(1,2,3)-t1(2,1,3))**2.+(t1(2,1,3)-t1(2,2,4))**2.&
                     +(t1(2,2,4)-t1(1,2,3))**2.+(t1(1,1,3)-t1(1,2,4))**2.&
                     +(t1(1,2,4)-t1(2,1,4))**2.+(t1(2,1,4)-t1(1,1,3))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(1,1,4)=t1(2,2,3)+sqrt(sq)
   if (nx1.ge.1.or.ny1.ge.1.or.nz2.le.4) goto 12
   ntwk=ntwk+1
   mxwk(ntwk)=1
   mywk(ntwk)=1
   mzwk(ntwk)=4
12 continue
endif

if (t1(3,1,4).eq.0.) then
   sx=average(s1(3,1,4),s1(3,2,4),s1(2,1,4),s1(2,2,4),&
              s1(3,1,3),s1(3,2,3),s1(2,1,3),s1(2,2,3))
   sq=3.*hh*sx*sx-.5*((t1(3,2,3)-t1(2,1,3))**2.+(t1(2,1,3)-t1(2,2,4))**2.&
                     +(t1(2,2,4)-t1(3,2,3))**2.+(t1(3,1,3)-t1(3,2,4))**2.&
                     +(t1(3,2,4)-t1(2,1,4))**2.+(t1(2,1,4)-t1(3,1,3))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(3,1,4)=t1(2,2,3)+sqrt(sq)
   if (nx2.le.3.or.ny1.ge.1.or.nz2.le.4) goto 13
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=1
   mzwk(ntwk)=4
13 continue
endif

if (t1(1,3,4).eq.0.) then
   sx=average(s1(1,3,4),s1(1,2,4),s1(2,3,4),s1(2,2,4),&
              s1(1,3,3),s1(1,2,3),s1(2,3,3),s1(2,2,3))
   sq=3.*hh*sx*sx-.5*((t1(2,3,3)-t1(1,2,3))**2.+(t1(1,2,3)-t1(2,2,4))**2.&
                     +(t1(2,2,4)-t1(2,3,3))**2.+(t1(1,3,3)-t1(2,3,4))**2.&
                     +(t1(2,3,4)-t1(1,2,4))**2.+(t1(1,2,4)-t1(1,3,3))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(1,3,4)=t1(2,2,3)+sqrt(sq)
   if (nx1.ge.1.or.ny2.le.3.or.nz2.le.4) goto 14
   ntwk=ntwk+1
   mxwk(ntwk)=1
   mywk(ntwk)=3
   mzwk(ntwk)=4
14 continue
endif

if (t1(3,3,4).eq.0.) then
   sx=average(s1(3,3,4),s1(3,2,4),s1(2,3,4),s1(2,2,4),&
              s1(3,3,3)+s1(3,2,3)+s1(2,3,3)+s1(2,2,3))
   sq=3.*hh*sx*sx-.5*((t1(2,3,3)-t1(3,2,3))**2.+(t1(3,2,3)-t1(2,2,4))**2.&
                     +(t1(2,2,4)-t1(2,3,3))**2.+(t1(3,3,3)-t1(2,3,4))**2.&
                     +(t1(2,3,4)-t1(3,2,4))**2.+(t1(3,2,4)-t1(3,3,3))**2.)
   if (sq.lt.0.0) sq=3.*hh*sx*sx
   t1(3,3,4)=t1(2,2,3)+sqrt(sq)
   if (nx2.le.3.or.ny2.le.3.or.nz2.le.4) goto 15
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=3
   mzwk(ntwk)=4
15 continue
endif
return
end subroutine treat2

subroutine treat3
!include 'eik3d.pa1'
!    common/sub/ t1(4,4,4),s1(4,4,4),nx1,nx2,ny1,ny2,nz1,nz2
!*************************************
! t(n1-1,j1-1,k1-1),t(n1-1,j1,k1-1),t(n1,j1-1,k1-1),   
! t(n1,j1,k1-1),t(n1-1,j1-1,k1),t(n1-1,j1,k1),t(n1,j1-1,k1) are known 
real :: sx, sq
if (t1(4,2,2).eq.0.) then
   sx=average(average(s1(3,2,2),s1(3,2,1),s1(3,1,2),s1(3,3,2),s1(3,2,3)),s1(4,2,2))
   sq=hh*sx*sx-.25*((t1(3,1,2)-t1(3,3,2))**2.+(t1(3,2,1)-t1(3,2,3))**2.)
   if (sq.lt.0.) sq=hh*sx*sx
   t1(4,2,2)=t1(3,2,2)+sqrt(sq)
   if (nx2.le.4.or.ny1.ge.2.or.nz1.ge.2) goto 1
   ntwk=ntwk+1
   mxwk(ntwk)=4
   mywk(ntwk)=2
   mzwk(ntwk)=2
1  continue
endif

if (t1(2,4,2).eq.0.) then
   sx=average(average(s1(2,3,2),s1(1,3,2),s1(3,3,2),s1(2,3,1),s1(2,3,3)),s1(2,4,2))
   sq=hh*sx*sx-.25*((t1(2,3,1)-t1(2,3,3))**2.+(t1(1,3,2)-t1(3,3,2))**2.)
   if (sq.lt.0.) sq=hh*sx*sx
   t1(2,4,2)=t1(2,3,2)+sqrt(sq)
   if (nx1.ge.2.or.ny2.le.4.or.nz1.ge.2) goto 2
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=4
   mzwk(ntwk)=2
2  continue
endif

if (t1(4,3,2).eq.0.) then
   sx=average(average(s1(3,2,3),s1(3,2,2),s1(3,2,1)),s1(3,3,2),s1(4,2,2),s1(4,3,2))
   sq=2.*hh*sx*sx-.5*(t1(3,2,1)-t1(3,2,3))**2.-(t1(4,2,2)-t1(3,3,2))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(4,3,2)=t1(3,2,2)+sqrt(sq)
   if (nx2.le.4.or.nz1.ge.2) goto 3
      ntwk=ntwk+1
      mxwk(ntwk)=4
      mywk(ntwk)=3
      mzwk(ntwk)=2
3  continue
endif

if (t1(3,4,2).eq.0.) then
   sx=average(average(s1(2,3,3),s1(2,3,2),s1(2,3,1)),s1(3,3,2),s1(2,4,2),s1(3,4,2))
   sq=2.*hh*sx*sx-.5*(t1(2,3,1)-t1(2,3,3))**2.-(t1(2,4,2)-t1(3,3,2))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(3,4,2)=t1(2,3,2)+sqrt(sq)
   if (ny2.le.4.or.nz1.ge.2) goto 4
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=4
   mzwk(ntwk)=2
4  continue
endif

if (t1(4,4,2).eq.0.) then
   sx=average(average(s1(3,3,1),s1(3,3,2),s1(3,3,3)),s1(3,4,2),s1(4,3,2),s1(4,4,2))
   sq=2.*hh*sx*sx-.5*(t1(3,3,1)-t1(3,3,3))**2.-(t1(4,3,2)-t1(3,4,2))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(4,4,2)=t1(3,3,2)+sqrt(sq)
   if (nx2.le.4.or.ny2.le.4.or.nz1.ge.2) goto 5
   ntwk=ntwk+1
   mxwk(ntwk)=4
   mywk(ntwk)=4
   mzwk(ntwk)=2
5  continue
endif

if (t1(4,3,3).eq.0.) then
   sx=average(average(s1(3,4,2),s1(3,3,2),s1(3,2,2)),s1(3,3,3),s1(4,3,3),s1(4,3,2))
   sq=2.*hh*sx*sx-.5*(t1(3,4,2)-t1(3,2,2))**2.-(t1(4,3,2)-t1(3,3,3))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(4,3,3)=t1(3,3,2)+sqrt(sq)
   if (nx2.le.4) goto 6
   ntwk=ntwk+1
   mxwk(ntwk)=4
   mywk(ntwk)=3
   mzwk(ntwk)=3
6  continue
endif

if (t1(3,4,3).eq.0.) then
   sx=average(average(s1(4,3,2),s1(3,3,2),s1(2,3,2)),s1(3,3,3),s1(3,4,3),s1(3,4,2))
   sq=2.*hh*sx*sx-.5*(t1(4,3,2)-t1(2,3,2))**2.-(t1(3,4,2)-t1(3,3,3))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(3,4,3)=t1(3,3,2)+sqrt(sq)
   if (ny2.le.4) goto 7
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=4
   mzwk(ntwk)=3
7  continue
endif

if (t1(4,2,3).eq.0.) then
   sx=average(s1(4,2,3),s1(4,3,3),s1(3,3,3),s1(3,2,3),&
              s1(3,3,2),s1(4,3,2),s1(4,2,2),s1(3,2,2))
   sq=3.*hh*sx*sx-.5*((t1(4,3,3)-t1(4,2,2))**2.+(t1(4,2,2)-t1(3,2,3))**2.&
                     +(t1(3,2,3)-t1(4,3,3))**2.+(t1(4,3,2)-t1(3,2,2))**2.&
                     +(t1(3,2,2)-t1(3,3,3))**2.+(t1(3,3,3)-t1(4,3,2))**2.)
   if (sq.lt.0.) sq=3.*hh*sx*sx
   t1(4,2,3)=t1(3,3,2)+sqrt(sq)
   if (nx2.le.4.or.ny1.ge.2) goto 8
   ntwk=ntwk+1
   mxwk(ntwk)=4
   mywk(ntwk)=2
   mzwk(ntwk)=3
8  continue
endif

if (t1(4,4,3).eq.0.) then
   sx=average(s1(4,4,3),s1(4,3,3),s1(3,3,3),s1(3,4,3),&
              s1(3,3,2),s1(4,3,2),s1(4,4,2),s1(3,4,2))
   sq=3.*hh*sx*sx-.5*((t1(4,3,3)-t1(4,4,2))**2.+(t1(4,4,2)-t1(3,4,3))**2.&
                     +(t1(3,4,3)-t1(4,3,3))**2.+(t1(4,3,2)-t1(3,4,2))**2.&
                     +(t1(3,4,2)-t1(3,3,3))**2.+(t1(3,3,3)-t1(4,3,2))**2.)
   if (sq.lt.0.) sq=3.*hh*sx*sx
   t1(4,4,3)=t1(3,3,2)+sqrt(sq)
   if (nx2.le.4.or.ny2.le.4) goto 9
   ntwk=ntwk+1
   mxwk(ntwk)=4
   mywk(ntwk)=4
   mzwk(ntwk)=3
9  continue
endif

if (t1(2,4,3).eq.0.) then
   sx=average(s1(2,4,3),s1(2,3,3),s1(3,3,3),s1(3,4,3),&
              s1(3,3,2),s1(2,3,2),s1(2,4,2),s1(3,4,2))
   sq=3.*hh*sx*sx-.5*((t1(2,3,3)-t1(2,4,2))**2.+(t1(2,4,2)-t1(3,4,3))**2.&
                     +(t1(3,4,3)-t1(2,3,3))**2.+(t1(2,3,2)-t1(3,4,2))**2.&
                     +(t1(3,4,2)-t1(3,3,3))**2.+(t1(3,3,3)-t1(2,3,2))**2.)
   if (sq.lt.0.) sq=3.*hh*sx*sx
   t1(2,4,3)=t1(3,3,2)+sqrt(sq)
   if (nx1.ge.2.or.ny2.le.4) goto 10
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=4
   mzwk(ntwk)=3
10 continue
endif

if (t1(3,3,4).eq.0.) then
   sx=average(average(s1(4,3,3),s1(2,3,3),s1(3,3,3),s1(3,4,3),s1(3,2,3)),s1(3,3,4))
   sq=hh*sx*sx-.25*((t1(2,3,3)-t1(4,3,3))**2.+(t1(3,4,3)-t1(3,2,3))**2.)
   if (sq.lt.0.) sq=hh*sx*sx
   t1(3,3,4)=t1(3,3,3)+sqrt(sq)
   if (nz2.le.4) goto 11
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=3
   mzwk(ntwk)=4
11 continue
endif

if (t1(3,4,4).eq.0.) then
   sx=average(average(s1(3,3,3),s1(4,3,3),s1(2,3,3)),s1(3,4,3),s1(3,3,4),s1(3,4,4))
   sq=2.*hh*sx*sx-.5*(t1(4,3,3)-t1(2,3,3))**2.-(t1(3,4,3)-t1(3,3,4))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(3,4,4)=t1(3,3,3)+sqrt(sq)
   if (ny2.le.4.or.nz2.le.4) goto 12
      ntwk=ntwk+1
      mxwk(ntwk)=3
      mywk(ntwk)=4
      mzwk(ntwk)=4
12 continue
endif

if (t1(3,2,4).eq.0.) then
   sx=average(average(s1(3,3,3),s1(4,3,3),s1(2,3,3)),s1(3,2,3),s1(3,3,4),s1(3,2,4))
   sq=2.*hh*sx*sx-.5*(t1(4,3,3)-t1(2,3,3))**2.-(t1(3,2,3)-t1(3,3,4))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(3,2,4)=t1(3,3,3)+sqrt(sq)
   if (ny1.ge.2.or.nz2.le.4) goto 13
   ntwk=ntwk+1
   mxwk(ntwk)=3
   mywk(ntwk)=2
   mzwk(ntwk)=4
13 continue
endif

if (t1(4,3,4).eq.0.) then
   sx=average(average(s1(3,3,3),s1(3,4,3),s1(3,2,3)),s1(4,3,3),s1(3,3,4),s1(4,3,4))
   sq=2.*hh*sx*sx-.5*(t1(3,4,3)-t1(3,2,3))**2.-(t1(4,3,3)-t1(3,3,4))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(4,3,4)=t1(3,3,3)+sqrt(sq)
   if (nx2.le.4.or.nz2.le.4) goto 14
   ntwk=ntwk+1
   mxwk(ntwk)=4
   mywk(ntwk)=3
   mzwk(ntwk)=4
14 continue
endif

if (t1(2,3,4).eq.0.) then
   sx=average(average(s1(3,3,3),s1(3,4,3),s1(3,2,3)),s1(2,3,3),s1(3,3,4),s1(2,3,4))
   sq=2.*hh*sx*sx-.5*(t1(3,4,3)-t1(3,2,3))**2.-(t1(2,3,3)-t1(3,3,4))**2.
   if (sq.lt.0.) sq=2.*hh*sx*sx
   t1(2,3,4)=t1(3,3,3)+sqrt(sq)
   if (nx1.ge.2.or.nz2.le.4) goto 15
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=3
   mzwk(ntwk)=4
15 continue
endif

if (t1(4,2,4).eq.0.) then
   sx=average(s1(4,2,4),s1(4,3,4),s1(3,3,4),s1(3,2,4),&
              s1(3,3,3),s1(4,3,3),s1(4,2,3),s1(3,2,3))
   sq=3.*hh*sx*sx-.5*((t1(4,3,4)-t1(4,2,3))**2.+(t1(4,2,3)-t1(3,2,4))**2.&
                     +(t1(3,2,4)-t1(4,3,4))**2.+(t1(4,3,3)-t1(3,2,3))**2.&
                     +(t1(3,2,3)-t1(3,3,4))**2.+(t1(3,3,4)-t1(4,3,3))**2.)
   if (sq.lt.0.) sq=3.*hh*sx*sx
   t1(4,2,4)=t1(3,3,3)+sqrt(sq)
   if (nx2.le.4.or.ny1.ge.2.or.nz2.le.4) goto 16
   ntwk=ntwk+1
   mxwk(ntwk)=4
   mywk(ntwk)=2
   mzwk(ntwk)=4
16 continue
endif
 
if (t1(4,4,4).eq.0.) then
   sx=average(s1(4,4,4),s1(4,3,4),s1(3,3,4),s1(3,4,4),&
              s1(3,3,3),s1(4,3,3),s1(4,4,3),s1(3,4,3))
   sq=3.*hh*sx*sx-.5*((t1(4,3,4)-t1(4,4,3))**2.+(t1(4,4,3)-t1(3,4,4))**2.&
                     +(t1(3,4,4)-t1(4,3,4))**2.+(t1(4,3,3)-t1(3,4,3))**2.&
                     +(t1(3,4,3)-t1(3,3,4))**2.+(t1(3,3,4)-t1(4,3,3))**2.)
   if (sq.lt.0.) sq=3.*hh*sx*sx
   t1(4,4,4)=t1(3,3,3)+sqrt(sq)
   if (nx2.le.4.or.ny2.le.4.or.nz2.le.4) goto 17
      ntwk=ntwk+1
      mxwk(ntwk)=4
      mywk(ntwk)=4
      mzwk(ntwk)=4
17 continue
endif

if (t1(2,4,4).eq.0.) then
   sx=average(s1(2,4,4),s1(2,3,4),s1(3,3,4),s1(3,4,4),&
              s1(3,3,3),s1(2,3,3),s1(2,4,3),s1(3,4,3))
   sq=3.*hh*sx*sx-.5*((t1(2,3,4)-t1(2,4,3))**2.+(t1(2,4,3)-t1(3,4,4))**2.&
                     +(t1(3,4,4)-t1(2,3,4))**2.+(t1(2,3,3)-t1(3,4,3))**2.&
                     +(t1(3,4,3)-t1(3,3,4))**2.+(t1(3,3,4)-t1(2,3,3))**2.)
   if (sq.lt.0.) sq=3.*hh*sx*sx
   t1(2,4,4)=t1(3,3,3)+sqrt(sq)
   if (nx1.ge.2.or.ny2.le.4.or.nz2.le.4) goto 18
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=4
   mzwk(ntwk)=4
18 continue
endif

if (t1(2,2,4).eq.0.) then
   sx=average(s1(2,2,4),s1(2,3,4),s1(3,3,4),s1(3,2,4),&
              s1(3,3,3),s1(2,3,3),s1(2,2,3),s1(3,2,3))
   sq=3.*hh*sx*sx-.5*((t1(2,3,4)-t1(2,2,3))**2.+(t1(2,2,3)-t1(3,2,4))**2.&
                     +(t1(3,2,4)-t1(2,3,4))**2.+(t1(2,3,3)-t1(3,2,3))**2.&
                     +(t1(3,2,3)-t1(3,3,4))**2.+(t1(3,3,4)-t1(2,3,3))**2.)
   if (sq.lt.0.) sq=3.*hh*sx*sx
   t1(2,2,4)=t1(3,3,3)+sqrt(sq)
   if (nx1.ge.2.or.ny1.ge.2.or.nz2.le.4) goto 19
   ntwk=ntwk+1
   mxwk(ntwk)=2
   mywk(ntwk)=2
   mzwk(ntwk)=4
19 continue
endif
return
end subroutine treat3

end module module_pweikfd3d_fir
