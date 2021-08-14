module module_ref3d
implicit none

contains

!---------------------------------------------------------------------------------------
! This subroutine is used to get the accurate traveltime for 
! a special position (zr,yr,zr)
subroutine accurate_time(tr,t,xr,yr,zr,h)
real, intent(in)  :: t(:,:,:),zr,yr,xr,h
real, intent(out) :: tr
real              :: x0,y0,z0,za,ya,xa
integer           :: iz,iy,ix
!  (iz,ix)  *-------------------* (iz,ix+1)
!  (z0,x0)  |       |           |
!           |       | (za)      |
!           | (xa)  |           |
!           |-------*  (zr,xr)  | 
!           |                   |
!           |                   |
!           |                   |
!           |                   |
! (iz+1,ix) *-------------------* (iz+1,ix+1)  
iz=int(zr/h)+1
iy=int(yr/h)+1
ix=int(xr/h)+1
z0=(iz-1)*h
y0=(iy-1)*h
x0=(ix-1)*h
za=(zr-z0)/h
ya=(yr-y0)/h
xa=(xr-x0)/h
tr=(1.-za)*(1.-ya)*(1.-xa)*t(iz,iy,ix)&
  +(1.-za)*(1.-ya)*    xa *t(iz,iy,ix+1) &
  +(1.-za)*    ya *(1.-xa)*t(iz,iy+1,ix) &
  +(1.-za)*    ya *    xa *t(iz,iy+1,ix+1) &
  +    za *(1.-ya)*(1.-xa)*t(iz+1,iy,ix)&
  +    za *(1.-ya)*    xa *t(iz+1,iy,ix+1) &
  +    za *    ya *(1.-xa)*t(iz+1,iy+1,ix) &
  +    za *    ya *    xa *t(iz+1,iy+1,ix+1)
end subroutine accurate_time

!                         ************
!                         * Sub ray3d*
!                         ************
!
subroutine ray3d(nzp,nyp,nxp,xr,yr,zr,xs,ys,zs,h,t1,ra,pre,path,a1,tt)
integer, intent(in) :: nzp, nyp, nxp
real,    intent(in) :: xr, yr, zr, xs, ys, zs, h, t(:,:,:)
real,    intent(inout) :: ra(:,:,:), pre(:,:,:)
integer, intent(inout) :: path(:,:,:)
      dimension t(nxp,nyp,nzp),g(nxp,nyp,nzp),pre(nxp,nyp,nzp)
      common/ray/pre
tol=1.5*h
kn=1
nmax=nxp+nyp+nzp
x0=xr
y0=yr
z0=zr
a=a1/sqrt((x0-xs)*(x0-xs)+(y0-ys)*(y0-ys))
nxr=nint(xr/h)+1
nyr=nint(yr/h)+1
nzr=nint(zr/h)+1
nxs=nint(xs/h)+1
nys=nint(ys/h)+1
nzs=nint(zs/h)+1
t0=tt
i2i=(sqrt(sqrt((nxr-nxs)*(nxr-nxs)*1.0+(nyr-nys)*(nyr-nys)*1.0))*rayfat)
ii1=nzr
ii2=nyr
ii3=nxr
path(ii1,ii2,ii3)=path(ii1,ii2,ii3)+1
g(ii1,ii2,ii3)=a+g(ii1,ii2,ii3)
pre(ii1,ii2,ii3)=pre(ii1,ii2,ii3)+1.

! Find the grid that contains the point

500 continue
kn=kn+1
if (kn.gt.nmax) goto 1000
nx0=nint(x0/h)+1
ny0=nint(y0/h)+1
nz0=nint(z0/h)+1

! Find the gradient at this point

vx=(t(nz0,ny0,nx0+1)-t(nz0,ny0,nx0-1))
vy=(t(nz0,ny0+1,nx0)-t(nz0,ny0-1,nx0))
vz=(t(nz0+1,ny0,nx0)-t(nz0-1,ny0,nx0))
va=sqrt(vx*vx+vy*vy+vz*vz)

! Along the gradient direction, forward one grid length
vx=-h*vx/va
vy=-h*vy/va
vz=-h*vz/va

x0=x0+vx
y0=y0+vy
z0=z0+vz

ii1=nint(z0/h+1)
ii2=nint(y0/h+1)
ii3=nint(x0/h+1)
do 989 i=max(1,ii1-1),min(nx,ii1+1)
do 989 j=max(1,ii2-1),min(ny,ii2+1)
do 989 k=max(1,ii3-i2i),min(nz,ii3+i2i)
g(i,j,k)=a+g(i,j,k)
989     pre(i,j,k)=pre(i,j,k)+1.

if (sqrt((x0-xs)*(x0-xs)+(y0-ys)*(y0-ys)+(z0-zs)*(z0-zs)).gt.tol) goto 500

1000 continue
ii1=nzs
ii2=nys
ii3=nxs
g(ii1,ii2,ii3)=a+g(ii1,ii2,ii3)
pre(ii1,ii2,ii3)=pre(ii1,ii2,ii3)+1.

return
end



end module module_ref3d
