module module_eikfd3d_dir
use module_utility
use module_io
use module_string

implicit none
integer, private    :: nz,ny,nx, isx,isy,isz, iTop, iBottom, iLeft, &
   iRight, iFront, iBack
integer, private    :: n1, n2, n
logical, private    :: isTop, isBottom, isLeft, isRight, isFront, isBack
real,    private    :: sx,sy,sz,h,hh,tm
real,pointer,private:: s(:,:,:),t(:,:,:)
real,private,parameter :: c1=1.0/(2.0**0.5)
real, private, allocatable :: s1(:,:),s2(:,:),t1(:,:),t2(:,:),t1sort(:)
integer, private, allocatable :: indx(:)
!real,private,target,allocatable :: ss(:,:,:),tt(:,:,:)
real, private :: s11(2,2),s12(2,2),t11(2,2),t12(2,2)
real, private :: s21(3,2),s22(3,2),t21(3,2),t22(3,2)
real, private :: s31(3,3),s32(3,3),t31(3,3),t32(3,3)
logical, private :: is_t2near(4), is_t2far(4)
real, private, allocatable :: es1(:),es2(:),es3(:),es4(:)
real, private, allocatable :: et1(:),et2(:),et3(:),et4(:)

private :: aroundSource, cubic_expand, surface_expand, eadge_expand, &
   treat1, treat2, treat3, treat1_ur, treat1_dr, treat1_dl, treat1_ul,&
   treat2_u, treat2_r, treat2_d, treat2_l, treat3_c, is_near, is_far, count_me

contains

!======================================================================
! the main subroutine for 3D FD solver of direct arrival of eikonal eqn
!======================================================================
subroutine eikfd3d_dir(s0,t0,sz0,sy0,sx0,nz0,ny0,nx0,h0)
implicit none
real,target,intent(in)   :: s0(:,:,:)
real,intent(in)          :: sx0,sy0,sz0,h0
integer, intent(in)      :: nz0,ny0,nx0
real,target,intent(out)  :: t0(:,:,:)
integer                  :: nmax
call copy_variable(sz0,sz,sy0,sy,sx0,sx,h0,h)
call copy_variable(nz0,nz,ny0,ny,nx0,nx)
tm=0.0001
s=>s0
t=>t0
hh=h*h
nmax=max(nz,ny,nx)
call allocate_and_initial(s1,s2,t1,t2,nmax,nmax)
call allocate_and_initial(es1,es2,es3,es4,nmax)
call allocate_and_initial(et1,et2,et3,et4,nmax)
call allocate_and_initial(t1sort,nz*ny*nx/min(nz,ny,nx))
call allocate_and_initial(indx,nz*ny*nx/min(nz,ny,nx))
isx=nint(sx/h)+1
isy=nint(sy/h)+1
isz=nint(sz/h)+1
call aroundSource()
!t=t+tm+0.00001
!goto 999
call cubic_expand()
999 continue
call deallocate_and_free(s1,s2,t1,t2)
call deallocate_and_free(es1,es2,es3,es4)
call deallocate_and_free(et1,et2,et3,et4)
end subroutine eikfd3d_dir

!----------------------------------------------------------------------
! This subroutine is used to calculate the traveltime around the source
!----------------------------------------------------------------------
subroutine aroundSource
implicit none
integer               :: ix1,ix2,iy1,iy2,iz1,iz2
real                  :: s_ave
integer               :: i1,i2,i3
integer :: ix, iy,iz
ix1=max(isx-2,1)
ix2=min(isx+2,nx)
iy1=max(isy-2,1)
iy2=min(isy+2,ny)
iz1=max(isz-2,1)
iz2=min(isz+2,nz)
iTop=iz1
iBottom=iz2
iLeft=ix1
iRight=ix2
iBack=iy1
iFront=iy2
do i3=ix1,ix2
   do i2=iy1,iy2
      do i1=iz1,iz2
         s_ave=average(s(min(i1,isz):max(i1,isz),min(i2,isy):max(i2,isy),&
                         min(i3,isx):max(i3,isx))) 
         t(i1,i2,i3)=sqrt(((i3-1)*h-sx)**2.0+((i2-1)*h-sy)**2.0 &
                      +((i1-1)*h-sz)**2.0) * s_ave
      enddo
   enddo
enddo

end subroutine aroundSource

!-----------------------------------------------------------------
! This subroutine is used to expand the traveltime table by square 
!-----------------------------------------------------------------
subroutine cubic_expand()
implicit none
integer :: i_outer,i_inner
integer :: j1,j2,k1,k2
integer :: i

do i=3,max(isx-1,isy-1,isz-1,nx-isx,ny-isy,nz-isz)
!  deal with the top face
   if (isz-i>=1) then
      isTop=.true.
      iTop=iTop-1
      i_outer=isz-i
      i_inner=i_outer+1
      j1=max(1,isy-i+1)
      j2=min(isy+i-1,ny)
      k1=max(1,isx-i+1)
      k2=min(isx+i-1,nx)
      n1=j2-j1+1
      n2=k2-k1+1
      s1(1:n1,1:n2)=s(i_inner,j1:j2,k1:k2)
      t1(1:n1,1:n2)=t(i_inner,j1:j2,k1:k2)
      s2(1:n1,1:n2)=s(i_outer,j1:j2,k1:k2)
      t2(1:n1,1:n2)=0.0
      call surface_expand
      t(i_outer,j1:j2,k1:k2)=t2(1:n1,1:n2)
   else
      isTop=.false.
   endif
! the end of the dealing with the top face

! deal with the right face
   if (isx+i<=nx) then
      isRight=.true.
      iRight=iRight+1
      i_outer=isx+i
      i_inner=i_outer-1
      j1=max(1,isz-i+1)
      j2=min(isz+i-1,nz)
      k1=max(1,isy-i+1)
      k2=min(isy+i-1,ny)
      n1=j2-j1+1
      n2=k2-k1+1
      s1(1:n1,1:n2)=s(j1:j2,k1:k2,i_inner)
      t1(1:n1,1:n2)=t(j1:j2,k1:k2,i_inner)
      s2(1:n1,1:n2)=s(j1:j2,k1:k2,i_outer)
      t2(1:n1,1:n2)=0.0
      call surface_expand()
      t(j1:j2,k1:k2,i_outer)=t2(1:n1,1:n2)
   else
      isRight=.false.
   endif    
! the end of the dealing with the right face

! deal with the back face
   if (isy-i>=1) then
      isBack=.true.
      iBack=iBack-1
      i_outer=isy-i
      i_inner=i_outer+1
      j1=max(1,isz-i+1)
      j2=min(isz+i-1,nz)
      k1=max(1,isx-i+1)
      k2=min(isx+i-1,nx)
      n1=j2-j1+1
      n2=k2-k1+1
      s1(1:n1,1:n2)=s(j1:j2,i_inner,k1:k2)
      t1(1:n1,1:n2)=t(j1:j2,i_inner,k1:k2)
      s2(1:n1,1:n2)=s(j1:j2,i_outer,k1:k2)
      t2(1:n1,1:n2)=0.0
      call surface_expand()
      t(j1:j2,i_outer,k1:k2)=t2(1:n1,1:n2)
   else
      isBack=.false.
   endif
! the end of the dealing with the back face

! deal with the left face
   if (isx-i>=1) then
      isLeft=.true.
      iLeft=iLeft-1
      i_outer=isx-i
      i_inner=i_outer+1
      j1=max(1,isz-i+1)
      j2=min(isz+i-1,nz)
      k1=max(1,isy-i+1)
      k2=min(isy+i-1,ny)
      n1=j2-j1+1
      n2=k2-k1+1
      s1(1:n1,1:n2)=s(j1:j2,k1:k2,i_inner)
      t1(1:n1,1:n2)=t(j1:j2,k1:k2,i_inner)
      s2(1:n1,1:n2)=s(j1:j2,k1:k2,i_outer)
      t2(1:n1,1:n2)=0.0
      call surface_expand()
      t(j1:j2,k1:k2,i_outer)=t2(1:n1,1:n2)
   else
      isLeft=.false.
   endif
! the end of the dealing with the left face

! deal with the front face
   if (isy+i<=ny) then
      isFront=.true.
      iFront=iFront+1
      i_outer=isy+i
      i_inner=i_outer-1
      j1=max(1,isz-i+1)
      j2=min(isz+i-1,nz)
      k1=max(1,isx-i+1)
      k2=min(isx+i-1,nx)
      n1=j2-j1+1
      n2=k2-k1+1
      s1(1:n1,1:n2)=s(j1:j2,i_inner,k1:k2)
      t1(1:n1,1:n2)=t(j1:j2,i_inner,k1:k2)
      s2(1:n1,1:n2)=s(j1:j2,i_outer,k1:k2)
      t2(1:n1,1:n2)=0.0
      call surface_expand()
      t(j1:j2,i_outer,k1:k2)=t2(1:n1,1:n2)
   else
      isFront=.false.
   endif
! end deling with the front face

! deal with the bottom face
   if (isz+i<=nz) then
      isBottom=.true.
      iBottom=iBottom+1
      i_outer=isz+i
      i_inner=i_outer-1
      j1=max(1,isy-i+1)
      j2=min(isy+i-1,ny)
      k1=max(1,isx-i+1)
      k2=min(isx+i-1,nx)
      n1=j2-j1+1
      n2=k2-k1+1
      s1(1:n1,1:n2)=s(i_inner,j1:j2,k1:k2)
      t1(1:n1,1:n2)=t(i_inner,j1:j2,k1:k2)
      s2(1:n1,1:n2)=s(i_outer,j1:j2,k1:k2)
      t2(1:n1,1:n2)=0.0
      call surface_expand()
      t(i_outer,j1:j2,k1:k2)=t2(1:n1,1:n2)
   else
      isBottom=.false.
   endif
! end dealing with the bottom face

! dealing with eadges
   ! eadge between top and left
   !if (isTop.or.isLeft) then
   if (isTop.and.isLeft) then
      n=iFront-iBack+1
      es1(1:n)=s(iTop+1,iBack:iFront,iLeft+1)      
      es2(1:n)=s(iTop+1,iBack:iFront,iLeft)      
      es3(1:n)=s(iTop,iBack:iFront,iLeft+1)      
      es4(1:n)=s(iTop,iBack:iFront,iLeft)      
      et1(1:n)=t(iTop+1,iBack:iFront,iLeft+1)      
      et2(1:n)=t(iTop+1,iBack:iFront,iLeft)      
      et3(1:n)=t(iTop,iBack:iFront,iLeft+1)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iTop,iBack+1:iFront-1,iLeft)=et4(2:n-1)
   endif

   ! eadge between top and right
   !if (isTop.or.isRight) then
   if (isTop.and.isRight) then
      n=iFront-iBack+1
      es1(1:n)=s(iTop+1,iBack:iFront,iRight-1)      
      es2(1:n)=s(iTop+1,iBack:iFront,iRight)      
      es3(1:n)=s(iTop,iBack:iFront,iRight-1)      
      es4(1:n)=s(iTop,iBack:iFront,iRight)      
      et1(1:n)=t(iTop+1,iBack:iFront,iRight-1)      
      et2(1:n)=t(iTop+1,iBack:iFront,iRight)      
      et3(1:n)=t(iTop,iBack:iFront,iRight-1)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iTop,iBack+1:iFront-1,iRight)=et4(2:n-1)
   endif

   ! eadge between top and front
   !if (isTop.or.isFront) then
   if (isTop.and.isFront) then
      n=iRight-iLeft+1
      es1(1:n)=s(iTop+1,iFront-1,iLeft:iRight)      
      es2(1:n)=s(iTop+1,iFront,iLeft:iRight)      
      es3(1:n)=s(iTop,iFront-1,iLeft:iRight)      
      es4(1:n)=s(iTop,iFront,iLeft:iRight)      
      et1(1:n)=t(iTop+1,iFront-1,iLeft:iRight)      
      et2(1:n)=t(iTop+1,iFront,iLeft:iRight)      
      et3(1:n)=t(iTop,iFront-1,iLeft:iRight)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iTop,iFront,iLeft+1:iRight-1)=et4(2:n-1)
   endif

   ! eadge between top and back
   !if (isTop.or.isBack) then
   if (isTop.and.isBack) then
      n=iRight-iLeft+1
      es1(1:n)=s(iTop+1,iBack+1,iLeft:iRight)      
      es2(1:n)=s(iTop+1,iBack,iLeft:iRight)      
      es3(1:n)=s(iTop,iBack+1,iLeft:iRight)      
      es4(1:n)=s(iTop,iBack,iLeft:iRight)      
      et1(1:n)=t(iTop+1,iBack+1,iLeft:iRight)      
      et2(1:n)=t(iTop+1,iBack,iLeft:iRight)      
      et3(1:n)=t(iTop,iBack+1,iLeft:iRight)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iTop,iBack,iLeft+1:iRight-1)=et4(2:n-1)
   endif

   ! eadge between bottom and left
   !if (isBottom.or.isLeft) then
   if (isBottom.and.isLeft) then
      n=iFront-iBack+1
      es1(1:n)=s(iBottom-1,iBack:iFront,iLeft+1)      
      es2(1:n)=s(iBottom-1,iBack:iFront,iLeft)      
      es3(1:n)=s(iBottom,iBack:iFront,iLeft+1)      
      es4(1:n)=s(iBottom,iBack:iFront,iLeft)      
      et1(1:n)=t(iBottom-1,iBack:iFront,iLeft+1)      
      et2(1:n)=t(iBottom-1,iBack:iFront,iLeft)      
      et3(1:n)=t(iBottom,iBack:iFront,iLeft+1)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iBottom,iBack+1:iFront-1,iLeft)=et4(2:n-1)
   endif

   ! eadge between bottom and right
   !if (isBottom.or.isRight) then
   if (isBottom.and.isRight) then
      n=iFront-iBack+1
      es1(1:n)=s(iBottom-1,iBack:iFront,iRight-1)      
      es2(1:n)=s(iBottom-1,iBack:iFront,iRight)      
      es3(1:n)=s(iBottom,iBack:iFront,iRight-1)      
      es4(1:n)=s(iBottom,iBack:iFront,iRight)      
      et1(1:n)=t(iBottom-1,iBack:iFront,iRight-1)      
      et2(1:n)=t(iBottom-1,iBack:iFront,iRight)      
      et3(1:n)=t(iBottom,iBack:iFront,iRight-1)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iBottom,iBack+1:iFront-1,iRight)=et4(2:n-1)
   endif

   ! eadge between bottom and front
   !if (isBottom.or.isFront) then
   if (isBottom.and.isFront) then
      n=iRight-iLeft+1
      es1(1:n)=s(iBottom-1,iFront-1,iLeft:iRight)      
      es2(1:n)=s(iBottom-1,iFront,iLeft:iRight)      
      es3(1:n)=s(iBottom,iFront-1,iLeft:iRight)      
      es4(1:n)=s(iBottom,iFront,iLeft:iRight)      
      et1(1:n)=t(iBottom-1,iFront-1,iLeft:iRight)      
      et2(1:n)=t(iBottom-1,iFront,iLeft:iRight)      
      et3(1:n)=t(iBottom,iFront-1,iLeft:iRight)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iBottom,iFront,iLeft+1:iRight-1)=et4(2:n-1)
   endif

   ! eadge between bottom and back
   !if (isBottom.or.isBack) then
   if (isBottom.and.isBack) then
      n=iRight-iLeft+1
      es1(1:n)=s(iBottom-1,iBack+1,iLeft:iRight)      
      es2(1:n)=s(iBottom-1,iBack,iLeft:iRight)      
      es3(1:n)=s(iBottom,iBack+1,iLeft:iRight)      
      es4(1:n)=s(iBottom,iBack,iLeft:iRight)      
      et1(1:n)=t(iBottom-1,iBack+1,iLeft:iRight)      
      et2(1:n)=t(iBottom-1,iBack,iLeft:iRight)      
      et3(1:n)=t(iBottom,iBack+1,iLeft:iRight)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iBottom,iBack,iLeft+1:iRight-1)=et4(2:n-1)
   endif

   ! eadge between front and left
   !if (isFront.or.isLeft) then
   if (isFront.and.isLeft) then
      n=iBottom-iTop+1
      es1(1:n)=s(iTop:iBottom,iFront-1,iLeft+1)      
      es2(1:n)=s(iTop:iBottom,iFront-1,iLeft)      
      es3(1:n)=s(iTop:iBottom,iFront,iLeft+1)      
      es4(1:n)=s(iTop:iBottom,iFront,iLeft)      
      et1(1:n)=t(iTop:iBottom,iFront-1,iLeft+1)      
      et2(1:n)=t(iTop:iBottom,iFront-1,iLeft)      
      et3(1:n)=t(iTop:iBottom,iFront,iLeft+1)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iTop+1:iBottom-1,iFront,iLeft)=et4(2:n-1)
   endif

   ! eadge between front and right
   !if (isFront.or.isRight) then
   if (isFront.and.isRight) then
      n=iBottom-iTop+1
      es1(1:n)=s(iTop:iBottom,iFront-1,iRight-1)      
      es2(1:n)=s(iTop:iBottom,iFront-1,iRight)      
      es3(1:n)=s(iTop:iBottom,iFront,iRight-1)      
      es4(1:n)=s(iTop:iBottom,iFront,iRight)      
      et1(1:n)=t(iTop:iBottom,iFront-1,iRight-1)      
      et2(1:n)=t(iTop:iBottom,iFront-1,iRight)      
      et3(1:n)=t(iTop:iBottom,iFront,iRight-1)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iTop+1:iBottom-1,iFront,iRight)=et4(2:n-1)
   endif

   ! eadge between back and left
   !if (isBack.or.isLeft) then
   if (isBack.and.isLeft) then
      n=iBottom-iTop+1
      es1(1:n)=s(iTop:iBottom,iBack+1,iLeft+1)      
      es2(1:n)=s(iTop:iBottom,iBack+1,iLeft)      
      es3(1:n)=s(iTop:iBottom,iBack,iLeft+1)      
      es4(1:n)=s(iTop:iBottom,iBack,iLeft)      
      et1(1:n)=t(iTop:iBottom,iBack+1,iLeft+1)      
      et2(1:n)=t(iTop:iBottom,iBack+1,iLeft)      
      et3(1:n)=t(iTop:iBottom,iBack,iLeft+1)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iTop+1:iBottom-1,iBack,iLeft)=et4(2:n-1)
   endif

   ! eadge between back and right
   !if (isBack.or.isRight) then
   if (isBack.and.isRight) then
      n=iBottom-iTop+1
      es1(1:n)=s(iTop:iBottom,iBack+1,iRight-1)      
      es2(1:n)=s(iTop:iBottom,iBack+1,iRight)      
      es3(1:n)=s(iTop:iBottom,iBack,iRight-1)      
      es4(1:n)=s(iTop:iBottom,iBack,iRight)      
      et1(1:n)=t(iTop:iBottom,iBack+1,iRight-1)      
      et2(1:n)=t(iTop:iBottom,iBack+1,iRight)      
      et3(1:n)=t(iTop:iBottom,iBack,iRight-1)      
      et4(1:n)=0.0      
      call eadge_expand
      t(iTop+1:iBottom-1,iBack,iRight)=et4(2:n-1)
   endif

   ! coners
   ! coner of top, front and left
   if (isTop.or.isFront.or.isLeft) then
   !if (isTop.and.isFront.and.isLeft) then
      s11=s(iTop+1,iFront-1:iFront+1,iLeft+1:iLeft:-1)
      s12=s(iTop,iFront-1:iFront+1,iLeft+1:iLeft:-1)
      t11=t(iTop+1,iFront-1:iFront+1,iLeft+1:iLeft:-1)
      t12=t(iTop,iFront-1:iFront+1,iLeft+1:iLeft:-1)
      t(iTop,iFront,iLeft)=treat1()
   endif

   ! coner of top, front and right
   if (isTop.or.isFront.or.isRight) then
   !if (isTop.and.isFront.and.isRight) then
      s11=s(iTop+1,iFront-1:iFront+1,iRight-1:iRight)
      s12=s(iTop,iFront-1:iFront+1,iRight-1:iRight)
      t11=t(iTop+1,iFront-1:iFront+1,iRight-1:iRight)
      t12=t(iTop,iFront-1:iFront+1,iRight-1:iRight)
      t(iTop,iFront,iRight)=treat1()
   endif

   ! coner of top, back and left
   if (isTop.or.isBack.or.isLeft) then
   !if (isTop.and.isBack.and.isLeft) then
      s11=s(iTop+1,iBack+1:iBack:-1,iLeft+1:iLeft:-1)
      s12=s(iTop,iBack+1:iBack:-1,iLeft+1:iLeft:-1)
      t11=t(iTop+1,iBack+1:iBack:-1,iLeft+1:iLeft:-1)
      t12=t(iTop,iBack+1:iBack:-1,iLeft+1:iLeft:-1)
      t(iTop,iBack,iLeft)=treat1()
   endif

   ! coner of top, back and right
   if (isTop.or.isBack.or.isLeft) then
   !if (isTop.and.isBack.and.isLeft) then
      s11=s(iTop+1,iBack+1:iBack:-1,iRight-1:iRight)
      s12=s(iTop,iBack+1:iBack:-1,iRight-1:iRight)
      t11=t(iTop+1,iBack+1:iBack:-1,iRight-1:iRight)
      t12=t(iTop,iBack+1:iBack:-1,iRight-1:iRight)
      t(iTop,iBack,iRight)=treat1()
   endif

   ! coner of bottom, front and left
   if (isBottom.or.isFront.or.isLeft) then
   !if (isBottom.and.isFront.and.isLeft) then
      s11=s(iBottom-1,iFront-1:iFront+1,iLeft+1:iLeft:-1)
      s12=s(iBottom,iFront-1:iFront+1,iLeft+1:iLeft:-1)
      t11=t(iBottom-1,iFront-1:iFront+1,iLeft+1:iLeft:-1)
      t12=t(iBottom,iFront-1:iFront+1,iLeft+1:iLeft:-1)
      t(iBottom,iFront,iLeft)=treat1()
   endif

   ! coner of bottom, front and right
   if (isBottom.or.isFront.or.isRight) then
   !if (isBottom.and.isFront.and.isRight) then
      s11=s(iBottom-1,iFront-1:iFront+1,iRight-1:iRight)
      s12=s(iBottom,iFront-1:iFront+1,iRight-1:iRight)
      t11=t(iBottom-1,iFront-1:iFront+1,iRight-1:iRight)
      t12=t(iBottom,iFront-1:iFront+1,iRight-1:iRight)
      t(iBottom,iFront,iRight)=treat1()
   endif

   ! coner of bottom, back and left
   if (isBottom.or.isBack.or.isLeft) then
   !if (isBottom.and.isBack.and.isLeft) then
      s11=s(iBottom-1,iBack+1:iBack:-1,iLeft+1:iLeft:-1)
      s12=s(iBottom,iBack+1:iBack:-1,iLeft+1:iLeft:-1)
      t11=t(iBottom-1,iBack+1:iBack:-1,iLeft+1:iLeft:-1)
      t12=t(iBottom,iBack+1:iBack:-1,iLeft+1:iLeft:-1)
      t(iBottom,iBack,iLeft)=treat1()
   endif

   ! coner of bottom, back and right
   if (isBottom.or.isBack.or.isLeft) then
   !if (isBottom.and.isBack.and.isLeft) then
      s11=s(iBottom-1,iBack+1:iBack:-1,iRight-1:iRight)
      s12=s(iBottom,iBack+1:iBack:-1,iRight-1:iRight)
      t11=t(iBottom-1,iBack+1:iBack:-1,iRight-1:iRight)
      t12=t(iBottom,iBack+1:iBack:-1,iRight-1:iRight)
      t(iBottom,iBack,iRight)=treat1()
   endif
enddo
!write(*,*)isTop,isBottom,isLeft,isRight,isFront,isBack
999 return
end subroutine cubic_expand

! total area has been divided into 9 parts
!                n2
!   o----o-----------------o----o   expand order: 1
!   | 9  |        2        | 6  |          -> 2,3,4 and 5
!   o----o-----------------o----o          -> 6,7,8 and 9
!   |    |                 |    |
!   |    |                 |    |
!n1 | 5  |        1        | 3  |
!   |    |                 |    |
!   |    |                 |    |
!   |    |                 |    |
!   o----o-----------------o----o
!   | 8  |        4        | 7  |
!   o----o-----------------o----o
!
subroutine surface_expand()
implicit none
!integer, external :: count
integer :: nsort, i, i1, i2, count_is_t2near, count_is_t2far
! ----------------- dealing with area 1
nsort=(n1-2)*(n2-2)
t1sort(1:nsort)=arr2vec(t1(2:n1-1,2:n2-1),n1-2,n2-2)
call indexx(nsort,t1sort(1:nsort),indx(1:nsort))
do i=1,nsort
   call ind2sub(n1-2,n2-2,indx(i),i1,i2)
   i1=i1+1
   i2=i2+1
   is_t2near=is_near(t2,i1,i2)
   count_is_t2near=count_me(is_t2near)
   if (count_is_t2near.eq.0) then ! none is know on t2
      t2(i1,i2)=treat3_c(i1,i2)
   elseif (count_is_t2near.eq.1) then
      if (is_t2near(1)) then 
         t2(i1,i2)=treat2_u(i1,i2)
      elseif (is_t2near(2)) then
         t2(i1,i2)=treat2_r(i1,i2)
      elseif (is_t2near(3)) then
         t2(i1,i2)=treat2_d(i1,i2)
      else   
         t2(i1,i2)=treat2_l(i1,i2)
      endif
   else
      is_t2far=is_far(t2,i1,i2)
      count_is_t2far=count_me(is_t2far)
      if (count_is_t2near.eq.2) then
         if (is_t2near(1)) then
            if (is_t2near(2)) then
               if (is_t2far(1)) then
                  t2(i1,i2)=treat1_ur(i1,i2)
               else
                  t2(i1,i2)=min(treat2_u(i1,i2),treat2_l(i1,i2))
               endif
            elseif (is_t2near(3)) then
                  t2(i1,i2)=min(treat2_u(i1,i2),treat2_d(i1,i2))
            elseif (is_t2near(4)) then
               if (is_t2far(4)) then
                  t2(i1,i2)=treat1_ul(i1,i2)
               else
                  t2(i1,i2)=min(treat2_u(i1,i2),treat2_r(i1,i2))
               endif
            endif
         else
            if (is_t2near(2)) then
               if (is_t2near(3)) then
                  if (is_t2far(2)) then
                     t2(i1,i2)=treat1_dr(i1,i2)
                  else
                     t2(i1,i2)=min(treat2_r(i1,i2),treat2_d(i1,i2))
                  endif
               elseif (is_t2near(4)) then
                  t2(i1,i2)=min(treat2_l(i1,i2),treat2_r(i1,i2))
               endif
            else
               if (is_t2far(3)) then
                  t2(i1,i2)=treat1_dl(i1,i2)
               else
                  t2(i1,i2)=min(treat2_d(i1,i2),treat2_r(i1,i2))
               endif
            endif
         endif
      elseif (count_is_t2near.eq.3) then
         if (is_t2near(1).eq..false.) then ! is_t2near(2,3,4) known
            if (is_t2far(2)) then 
               if (is_t2far(3)) then ! is_t2far(2,3) known
                  t2(i1,i2)=min(treat1_dr(i1,i2),treat1_dl(i1,i2))
               else ! is_t2far(2) known
                  t2(i1,i2)=treat1_dr(i1,i2)
               endif
            else
               if (is_t2far(3)) then ! is_t2far(3) known
                  t2(i1,i2)=treat1_dl(i1,i2)
               else ! is_t2far all unknown
                  t2(i1,i2)=min(treat2_l(i1,i2),treat2_d(i1,i2),treat2_r(i1,i2))
               endif
            endif
         elseif (is_t2near(2).eq..false.) then ! is_t2near(1,3,4) known
            if (is_t2far(3)) then
               if (is_t2far(4)) then ! is_t2far(3,4) known
                  t2(i1,i2)=min(treat1_dl(i1,i2),treat1_ul(i1,i2))
               else ! is_t2far(3) known
                  t2(i1,i2)=treat1_dl(i1,i2)
               endif
            else
               if (is_t2far(4)) then ! is_t2far(4) known
                  t2(i1,i2)=treat1_ul(i1,i2)
               else ! is_t2far all unknown
                  t2(i1,i2)=min(treat2_d(i1,i2),treat2_r(i1,i2),treat2_u(i1,i2))
               endif
            endif
         elseif (is_t2near(3).eq..false.) then ! is_t2near(1,2,4) known
            if (is_t2far(4)) then
               if (is_t2far(1)) then ! is_t2far(1,4) known
                  t2(i1,i2)=min(treat1_ul(i1,i2),treat1_ur(i1,i2))
               else ! is_t2far(4) known
                  t2(i1,i2)=treat1_ul(i1,i2)
               endif
            else 
               if (is_t2far(1)) then ! is_t2far(1) known
                  t2(i1,i2)=treat1_ur(i1,i2)
               else ! is_t2far all unknown
                  t2(i1,i2)=min(treat2_l(i1,i2),treat2_u(i1,i2),treat2_r(i1,i2))
               endif
            endif
         else ! is_t2near(1,2,3) known 
            if (is_t2far(1)) then
               if (is_t2far(2)) then ! is_t2far(1,2) known
                  t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2))
               else  ! is_t2far(1) known
                  t2(i1,i2)=treat1_ur(i1,i2)
               endif
            else
               if (is_t2far(2)) then ! is_t2far(2) known
                  t2(i1,i2)=treat1_dr(i1,i2)
               else ! is_t2far all unknown
                  t2(i1,i2)=min(treat2_u(i1,i2),treat2_r(i1,i2),treat2_d(i1,i2))
               endif
            endif
         endif
      elseif (count_is_t2near.eq.4) then
         if (is_t2far(1)) then
            if (is_t2far(2)) then
               if (is_t2far(3)) then
                  if (is_t2far(4)) then ! is_t2far(1,2,3,4) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2),&
                                   treat1_dl(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(1,2,3) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2),treat1_dl(i1,i2))
                  endif
               else
                  if (is_t2far(4)) then ! is_t2far(1,2,4) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(1,2) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2))
                  endif
               endif
            else
               if (is_t2far(3)) then 
                  if (is_t2far(4)) then ! is_t2far(1,3,4) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dl(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(1,3) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dl(i1,i2))
                  endif
               else
                  if (is_t2far(4)) then ! is_t2far(1,4) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_ul(i1,i2))
                  else  ! is_t2far(1) known
                     t2(i1,i2)=treat1_ur(i1,i2)
                  endif
               endif
            endif
         else
            if (is_t2far(2)) then
               if (is_t2far(3)) then
                  if (is_t2far(4)) then ! is_t2far(2,3,4) known
                     t2(i1,i2)=min(treat1_dr(i1,i2),treat1_dl(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(2,3) known
                     t2(i1,i2)=min(treat1_dr(i1,i2),treat1_dl(i1,i2))
                  endif
               else
                  if (is_t2far(4)) then ! is_t2far(2,4) known
                     t2(i1,i2)=min(treat1_dr(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(2) known
                     t2(i1,i2)=treat1_dr(i1,i2)
                  endif
               endif
            else
               if (is_t2far(3)) then
                  if (is_t2far(4)) then ! is_t2far(3,4) known
                     t2(i1,i2)=min(treat1_dl(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(3) known
                     t2(i1,i2)=treat1_dl(i1,i2)
                  endif
               else 
                  if (is_t2far(4)) then ! is_t2far(4) known
                     t2(i1,i2)=treat1_ul(i1,i2)
                  else 
                     t2(i1,i2)=min(treat2_u(i1,i2),treat2_r(i1,i2),&
                                   treat2_d(i1,i2),treat2_l(i1,i2))
                  endif
               endif
            endif
         endif
      endif
   endif
enddo
! --------------end dealing with area 1

! ----------------- dealing with area 2
n=n2
es1(1:n)=s1(2,1:n)
es2(1:n)=s1(1,1:n)
es3(1:n)=s2(2,1:n)
es4(1:n)=s2(1,1:n)
et1(1:n)=t1(2,1:n)
et2(1:n)=t1(1,1:n)
et3(1:n)=t2(2,1:n)
et4(1:n)=0.0
call eadge_expand()
t2(1,2:n-1)=et4(2:n-1)
! --------------end dealing with area 2

! ----------------- dealing with area 3
n=n1
es1(1:n)=s1(1:n,n2-1)
es2(1:n)=s1(1:n,n2)
es3(1:n)=s2(1:n,n2-1)
es4(1:n)=s2(1:n,n2)
et1(1:n)=t1(1:n,n2-1)
et2(1:n)=t1(1:n,n2)
et3(1:n)=t2(1:n,n2-1)
et4(1:n)=0.0
call eadge_expand()
t2(2:n-1,n2)=et4(2:n-1)
! --------------end dealing with area 3

! ----------------- dealing with area 4
n=n2
es1(1:n)=s1(n1-1,1:n)
es2(1:n)=s1(n1,1:n)
es3(1:n)=s2(n1-1,1:n)
es4(1:n)=s2(n1,1:n)
et1(1:n)=t1(n1-1,1:n)
et2(1:n)=t1(n1,1:n)
et3(1:n)=t2(n1-1,1:n)
et4(1:n)=0.0
call eadge_expand()
t2(n1,2:n-1)=et4(2:n-1)
! --------------end dealing with area 4

! ----------------- dealing with area 5
n=n1
es1(1:n)=s1(1:n,2)
es2(1:n)=s1(1:n,1)
es3(1:n)=s2(1:n,2)
es4(1:n)=s2(1:n,1)
et1(1:n)=t1(1:n,2)
et2(1:n)=t1(1:n,1)
et3(1:n)=t2(1:n,2)
et4(1:n)=0.0
call eadge_expand()
t2(2:n-1,1)=et4(2:n-1)
! --------------end dealing with area 5

! ----------------- dealing with area 6
s11=s1(2:1:-1,n2-1:n2)
s12=s2(2:1:-1,n2-1:n2)
t11=t1(2:1:-1,n2-1:n2)
t12=t2(2:1:-1,n2-1:n2)
t2(1,n2)=treat1()
! --------------end dealing with area 6

! ----------------- dealing with area 7
s11=s1(n1-1:n1,n2-1:n2)
s12=s2(n1-1:n1,n2-1:n2)
t11=t1(n1-1:n1,n2-1:n2)
t12=t2(n1-1:n1,n2-1:n2)
t2(n1,n2)=treat1()
! --------------end dealing with area 7

! ----------------- dealing with area 8
s11=s1(n1-1:n1,2:1:-1)
s12=s2(n1-1:n1,2:1:-1)
t11=t1(n1-1:n1,2:1:-1)
t12=t2(n1-1:n1,2:1:-1)
t2(n1,1)=treat1()
! --------------end dealing with area 8

! ----------------- dealing with area 9
s11=s1(2:1:-1,2:1:-1)
s12=s2(2:1:-1,2:1:-1)
t11=t1(2:1:-1,2:1:-1)
t12=t2(2:1:-1,2:1:-1)
t2(1,1)=treat1()
! --------------end dealing with area 9

end subroutine surface_expand

! treat 1, most accurate type stencel
!     call ind2sub(n1-2,n2-2,indx(i),i1,i2)
!                      A (1,1) --
!       E--------F     B (1,2)   | t11 in plane
!      /|       /|     C (2,1)   |
!     / |      / |     D (2,2) --
!    /  |     /  |     
!   G--------H   |     E (1,1) --
!   |   A----|---B     F (1,2)   |
!   |  /     |  /      G (2,1)   | t12 out plane
!   | /      | /       H (2,2) --
!   |/       |/
!   C--------D         unknow is H, all other are known
! H  = A + 1/(2^0.5) * sq^0.5
! s_ave = the average of all 8 points
! sq = 6 * h^2 * s_ave^2 -(B-C)^2-(E-C)^2-(B-E)^2-(G-D)^2-(D-F)^2-(G-F)^2
real function treat1
implicit none
real :: s_ave, sq1, sq
s_ave=average(average(s11),average(s12))
sq1=6.0*hh*s_ave**2.0
sq=sq1-(t11(1,2)-t11(2,1))**2.0-(t12(1,1)-t11(2,1))**2.0&
      -(t11(1,2)-t12(1,1))**2.0-(t12(2,1)-t11(2,2))**2.0&
      -(t11(2,2)-t12(1,2))**2.0-(t12(2,1)-t12(1,2))**2.0
if (sq<0.0) then
   sq=sq1
endif
treat1=t11(1,1)+c1*sqrt(sq)
!if (sq<0.0) then
!   treat1=min(t11(2,2),t12(1,2),t12(2,1))+h*s_ave
!else
!   treat1=t11(1,1)+c1*sqrt(sq)
!endif
end function treat1

! treat2, second accurate type stencel
!
!             G---------H     A (1,1) ---
!            /|        /|     B (1,2)   |
!           / |       / |     C (2,1)   |_ t21 in plane
!          /  |      /  |     D (2,2)   |
!         I---------J   |     E (3,1)   |
!        /|  A/----/|---/B    F (3,2) ---  
!       / |  /    / |  /
!      /  | /    /  | /       G (1,1) ---
!     K---------L   |/        H (1,2)   |
!     |   C-----|---/D        I (2,1)   |_ t22 out plane
!     |  /      |  /          J (2,2)   |
!     | /       | /           K (3,1)   |
!     |/        |/            L (3,2) ---
!     E---------F         J is unknown    
! J = C + sq^0.5
! sq=2*h^2*s_ave^2-0.5*(A-E)^2-(I-D)^2
! s_ave= the average of all 12 points
real function treat2
implicit none
real :: s_ave, sq1, sq
!s_ave=average(average(s21),average(s22))
s_ave=average(average(s21(1:3,1)),s21(2,2),s22(2,1),s22(2,2))
sq1=2.0*hh*s_ave**2.0
sq=sq1-0.5*(t21(1,1)-t21(3,1))**2.0-(t22(2,1)-t21(2,2))**2.0
if (sq<0.0) then
   sq=sq1
endif
treat2=t21(2,1)+sqrt(sq)
!if (sq<0.0) then
!   treat2=min(t21(2,2),t22(2,1))+h*s_ave
!else
!   treat2=t21(2,1)+sqrt(sq)
!endif
end function treat2

! treat3, third type stencel            t1 in plane :
! 
!             J---------K---------L     A (1,1), B (1,2), C (1,3)   
!            /|        /|        /|     D (2,1), E (2,2), F (2,3)
!           / |       / |       / |     G (3,1), H (3,2), I (3,3)
!          /  |      /  |      /  |
!         M---A-----N---B-----O---C     t2 out plane : 
!        /|  /     /|  /     /|  /      J (1,1), K (1,2), L (1,3)
!       / | /     / | /     / | /       M (2,1), N (2,2), O (2,3)
!      /  |/     /  |/     /  |/        P (3,1), Q (3,2), R (3,3)
!     P---D-----Q---E-----R---F 
!     |  /      |  /      |  /          N is unknown
!     | /       | /       | /
!     |/        |/        |/
!     G---------H---------I   
! N = E + sq^0.5  
! sq=h^2*s_ave^2-0.25*[(D-F)^2+(B-H)^2]     
! s_ave= the average of 18 points    
real function treat3
implicit none
real :: s_ave, sq1, sq
!s_ave=average(average(s31),average(s32))
s_ave=average(average(s31(1,1),s31(2,1),s31(3,1),s31(1,2),s31(1,3)),s32(2,2))
sq1=hh*s_ave**2.0
sq=sq1-0.25*((t31(2,1)-t31(2,3))**2.0+(t31(1,2)-t31(3,2))**2.0)
if (sq<0.0) then
   sq=sq1
endif
treat3=t31(2,2)+sqrt(sq)
!if (sq<0.0) then
!   treat3=t31(2,2)+h*s_ave
!else
!   treat3=t31(2,2)+sqrt(sq)
!endif
end function treat3

! f4-----n1------f1    near(n1,n2,n3,n4)
!  |      |      |   
!  |      |      |     far(f1,f2,f3,f4)
! n4------o------n2
!  |      |      |
!  |      |      |
! f3-----n3------n4 
function is_near(t2,i1,i2)
implicit none
logical :: is_near(4)
real,   intent(in) :: t2(:,:)
integer,intent(in) :: i1, i2
is_near=.false.
if (t2(i1-1,i2)>tm) is_near(1)=.true.
if (t2(i1,i2+1)>tm) is_near(2)=.true.
if (t2(i1+1,i2)>tm) is_near(3)=.true.
if (t2(i1,i2-1)>tm) is_near(4)=.true.
end function is_near

function is_far(t2,i1,i2)
implicit none
logical :: is_far(4)
real,   intent(in) :: t2(:,:)
integer,intent(in) :: i1, i2
is_far=.false.
if (t2(i1-1,i2+1)>tm) is_far(1)=.true.
if (t2(i1+1,i2+1)>tm) is_far(2)=.true.
if (t2(i1+1,i2-1)>tm) is_far(3)=.true.
if (t2(i1-1,i2-1)>tm) is_far(4)=.true.
end function is_far

! f4-----n1------f1    near(n1,n2,n3,n4)
!  |      |      |   
!  |      |      |     far(f1,f2,f3,f4)
! n4------o------n2
!  |      |      |
!  |      |      |
! f3-----n3------f2 

real function treat1_ur(i1,i2)
implicit none
integer, intent(in) :: i1, i2
s11=s1(i1-1:i1,i2+1:i2:-1) 
s12=s2(i1-1:i1,i2+1:i2:-1) 
t11=t1(i1-1:i1,i2+1:i2:-1) 
t12=t2(i1-1:i1,i2+1:i2:-1) 
treat1_ur=treat1()
end function treat1_ur

real function treat1_dr(i1,i2)
implicit none
integer, intent(in) :: i1,i2
s11=s1(i1+1:i1:-1,i2+1:i2:-1)
s12=s2(i1+1:i1:-1,i2+1:i2:-1)
t11=t1(i1+1:i1:-1,i2+1:i2:-1)
t12=t2(i1+1:i1:-1,i2+1:i2:-1)
treat1_dr=treat1()
end function treat1_dr

real function treat1_dl(i1,i2)
implicit none
integer, intent(in) :: i1,i2
s11=s1(i1+1:i1:-1,i2-1:i2)
s12=s2(i1+1:i1:-1,i2-1:i2)
t11=t1(i1+1:i1:-1,i2-1:i2)
t12=t2(i1+1:i1:-1,i2-1:i2)
treat1_dl=treat1()
end function treat1_dl

real function treat1_ul(i1,i2)
implicit none
integer, intent(in) :: i1, i2
s11=s1(i1-1:i1,i2-1:i2)
s12=s2(i1-1:i1,i2-1:i2)
t11=t1(i1-1:i1,i2-1:i2)
t12=t2(i1-1:i1,i2-1:i2)
treat1_ul=treat1()
end function treat1_ul

real function treat2_u(i1,i2)
implicit none
integer, intent(in) :: i1, i2
s21=transpose(s1(i1-1:i1,i2-1:i2+1))
s22=transpose(s2(i1-1:i1,i2-1:i2+1))
t21=transpose(t1(i1-1:i1,i2-1:i2+1))
t22=transpose(t2(i1-1:i1,i2-1:i2+1))
treat2_u=treat2()
end function treat2_u

real function treat2_r(i1,i2)
implicit none
integer, intent(in) :: i1, i2
s21=s1(i1-1:i1+1,i2+1:i2:-1)
s22=s2(i1-1:i1+1,i2+1:i2:-1)
t21=t1(i1-1:i1+1,i2+1:i2:-1)
t22=t2(i1-1:i1+1,i2+1:i2:-1)
treat2_r=treat2()
end function treat2_r

real function treat2_d(i1,i2)
implicit none
integer, intent(in) :: i1, i2
s21=transpose(s1(i1+1:i1:-1,i2-1:i2+1))
s22=transpose(s2(i1+1:i1:-1,i2-1:i2+1))
t21=transpose(t1(i1+1:i1:-1,i2-1:i2+1))
t22=transpose(t2(i1+1:i1:-1,i2-1:i2+1))
treat2_d=treat2()
end function treat2_d

real function treat2_l(i1,i2)
implicit none
integer, intent(in) :: i1, i2
s21=s1(i1-1:i1+1,i2-1:i2)
s22=s2(i1-1:i1+1,i2-1:i2)
t21=t1(i1-1:i1+1,i2-1:i2)
t22=t2(i1-1:i1+1,i2-1:i2)
treat2_l=treat2()
end function treat2_l

real function treat3_c(i1,i2)
implicit none
integer, intent(in) :: i1,i2
s31=s1(i1-1:i1+1,i2-1:i2+1)
s32=s2(i1-1:i1+1,i2-1:i2+1)
t31=t1(i1-1:i1+1,i2-1:i2+1)
treat3_c=treat3()
end function treat3_c

! eadge
!
!           K-1---------L-1      (I-1,I,I+1) -- et1,es1 
!            /|        /|        (J-1,J,J+1) -- et2,es2
!           / |       / |        (K-1,K,K+1) -- et3,es3
!          /  |      /  |        (L-1,L,L+1) -- et4,es4
!         K---------L   |     
!        /| I-1/---/|---/J-1     I(et1(i) is the sorting point  
!       / |  /    / |  /         L(et4(i) is unknown
!      /  | /    /  | /       
!   K+1---------L+1 |/        
!     |   I-----|---/J        
!     |  /      |  /          
!     | /       | /           
!     |/        |/            
!   I+1---------J+1             
subroutine eadge_expand()
implicit none
integer :: nsort, ii, i
real :: et4i_tmp
nsort=n-2
t1sort(1:nsort)=et1(2:n-1)
call indexx(nsort,t1sort(1:nsort),indx(1:nsort))
do ii=1,nsort
   i=indx(ii)+1
   if (et4(i-1)>tm) then
      s11=compose(es1(i-1),es2(i-1),es1(i),es2(i))
      s12=compose(es3(i-1),es4(i-1),es3(i),es4(i))
      t11=compose(et1(i-1),et2(i-1),et1(i),et2(i))
      t12=compose(et3(i-1),et4(i-1),et3(i),et4(i))
      et4i_tmp=treat1()
      if (et4(i+1)>tm) then
         s11=compose(es1(i+1),es2(i+1),es1(i),es2(i))
         s12=compose(es3(i+1),es4(i+1),es3(i),es4(i))
         t11=compose(et1(i+1),et2(i+1),et1(i),et2(i))
         t12=compose(et3(i+1),et4(i+1),et3(i),et4(i))
         et4(i)=min(et4i_tmp,treat1())     
      else
         et4(i)=et4i_tmp
      endif
   else
      if (et4(i+1)>tm) then
         s11=compose(es1(i+1),es2(i+1),es1(i),es2(i))
         s12=compose(es3(i+1),es4(i+1),es3(i),es4(i))
         t11=compose(et1(i+1),et2(i+1),et1(i),et2(i))
         t12=compose(et3(i+1),et4(i+1),et3(i),et4(i))
         et4(i)=treat1()     
      else
         s21=compose(es1(i-1:i+1),es2(i-1:i+1),3)
         s22=compose(es3(i-1:i+1),es4(i-1:i+1),3)
         t21=compose(et1(i-1:i+1),et2(i-1:i+1),3)
         t22=compose(et3(i-1:i+1),et4(i-1:i+1),3)
         et4(i)=treat2()
      endif
   endif
enddo
end subroutine eadge_expand

integer function count_me(d_in)
implicit none
logical, intent(in) :: d_in(4)
integer :: i=4
count_me=0
do i=1,4
   if (d_in(i)) then
      count_me=count_me+1
   endif
enddo
end function count_me

subroutine expand_out_of_local_min
implicit none
integer :: i
do i=1,min_number
   s31=s1(i1min_sort(i)-1:i1min_sort(i)+1,i2min_sort(i)-1:i2min_sort(i)+1)
   s32=s2(i1min_sort(i)-1:i1min_sort(i)+1,i2min_sort(i)-1:i2min_sort(i)+1)
   t31=t1(i1min_sort(i)-1:i1min_sort(i)+1,i2min_sort(i)-1:i2min_sort(i)+1)
   t2(i1min_sort(i),i2min_sort(i))=treat3()
enddo
wf_up(1:4,min_number)=.true.
wf_dn(1:4,min_number)=.true.
wf_lt(1:4,min_number)=.true.
wf_rt(1:4,min_number)=.true.
do i=1,min_number
   wf_num(i)=4
   t2(i1min_sort(i)-1,i2min_sort(i))=treat2_u(i1,i2)
   i1_wf(1,i)=i1min_sort(i)-1
   i2_wf(1,i)=i2min_sort(i)
   wf_dn(1,i)=.false.
   t2(i1min_sort(i),i2min_sort(i)+1)=treat2_r(i1,i2)
   i1_wf(2,i)=i1min_sort(i)
   i2_wf(2,i)=i2min_sort(i)+1
   wf_lt(2,i)=.false.
   t2(i1min_sort(i)+1,i2min_sort(i))=treat2_d(i1,i2)
   i1_wf(3,i)=i1min_sort(i)+1
   i2_wf(3,i)=i2min_sort(i)
   wf_up(3,i)=.false.
   t2(i1min_sort(i),i2min_sort(i)-1)=treat2_l(i1,i2)
   i1_wf(4,i)=i1min_sort(i)
   i2_wf(4,i)=i2min_sort(i)-1
   wf_rt(4,i)=.false.
   isStop(i)=.false.
enddo
end subroutine expand_out_of_local_min

subroutine expand_from_local_min
implicit none
integer :: i
do i=1,min_number
   if (isStop(i).eq..false.) then
      do j=1,wf_num(i)
         if (wf_up(j,i)) then
                        
         

subroutine localMin()
implicit none
integer         :: i1,i2
min_number=0
min_loc=0
dist=0.0
isOutMin=.false.
do i2=2,n2-1
   do i1=2,n1-1
      if (t1(i1,i2)<=t1(i1-1,i2).and.t1(i1,i2)<=t1(i1+1,i2)&
          t1(i1,i2)<=t1(i1,i2-1).and.t1(i1,i2)<=t1(i1,i2+1))
          call local_min_count(i1,i2)
      endif
   enddo
enddo
call indexx(min_number,t1min(1:min_number),indx(1:min_number)
do i=1:min_number
   i1min_sort(i)=i1min(indx(i))
   i2min_sort(i)=i2min(indx(i))
enddo
end subroutine localMin

subroutine localMax()
implicit none
integer         :: i
max_number=0
max_loc=0
do i2=2,n2-1
   do i1=2,n1-1
      if (t1(i1,i2)>=t1(i1-1,i2).and.t1(i1,i2)>=t1(i1+1,i2)&
          t1(i1,i2)>=t1(i1,i2-1).and.t1(i1,i2)>=t1(i1,i2+1))
          call local_max_count(i1,i2)
      endif
   enddo
enddo
end subroutine localMax

subroutine local_min_count(i1,i2)
implicit none
integer, intent(in) :: i1, i2
min_number=min_number+1
t1min(min_number)=t1(i1,i2)
i1min(min_number)=i1
i2min(min_number)=i2
min_loc(1,min_number)=i1
min_loc(2,min_number)=i2
dist(min_number)=(i1-s1_index)**2+(i2-s2_index)**2
isOutMin(i1,i2)=.true.
end subroutine local_min_count

subroutine local_max_count(i1,i2)
implicit none
integer, intent(in) :: i1, i2
max_number=max_number+1
max_loc(1,max_number)=i1
max_loc(2,max_number)=i2
end subroutine local_max_count



end module module_eikfd3d_dir
