module module_ref2d
use module_io
implicit none

contains

!---------------------------------------------------------------------------------------
! This subroutine is used to get the accurate traveltime for 
! a special position (xr,zr)
subroutine accurate_time(tr,t,xr,zr,h)
real, intent(in)  :: t(:,:),zr,xr,h
real, intent(out) :: tr
real              :: x0,z0,za,xa
integer           :: iz,ix
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
ix=int(xr/h)+1
z0=(iz-1)*h
x0=(ix-1)*h
za=(zr-z0)/h
xa=(xr-x0)/h
tr=(1.-xa)*(1.-za)*t(iz,ix)+xa*(1.-za)*t(iz,ix+1) &
  +xa*za*t(iz+1,ix+1)+(1.-xa)*za*t(iz+1,ix)
end subroutine accurate_time

subroutine ray2(nz,nx,xr,zr,xs,zs,h,ttt,ra,pre,tdiff)
implicit none
integer, intent(in) :: nz,nx
real,    intent(in) :: xr,zr,xs,zs,h,tdiff
real,    intent(in) :: ttt(:,:)
real,    intent(out):: ra(:,:),pre(:,:)
integer :: nxr,nzr,nxs,nzs,nmax,kn,ii1,ii2,ii10,ii20,nx0,nz0
real    :: x0,z0,tol,dist,vx,vz,va
ra=0.0
pre=0.0
dist=0.0

tol=(1.5*h)**2.0
kn=1
nmax=2*(nz+nx)

x0=xr
z0=zr
nxr=nint(xr/h)+1
nzr=nint(zr/h)+1
nxs=nint(xs/h)+1
nzs=nint(zs/h)+1

dist=((xr-(nxr-1)*h)**2.0+(zr-(nzr-1)*h)**2.0)**0.5

ii1=nzr
ii2=nxr
pre(ii1,ii2)=pre(ii1,ii2)+1.0
ii10=ii1
ii20=ii2

do while ((x0-xs)**2.0+(z0-zs)**2.0 > tol .and. kn<=nmax)
    kn=kn+1
    nx0=nint(x0/h)+1
    nz0=nint(z0/h)+1
    vx=ttt(nz0,nx0+1)-ttt(nz0,nx0-1)
    vz=ttt(nz0+1,nx0)-ttt(nz0-1,nx0)
    va=(vx*vx+vz*vz)**0.5
    vx=-h*vx/va
    vz=-h*vz/va
    x0=x0+vx
    z0=z0+vz
    ii1=nint(z0/h)+1
    ii2=nint(x0/h)+1
    ii1=min(ii1,nz)
    ii1=max(ii1,1)
    ii2=min(ii2,nx)
    ii2=max(ii2,1)
    if (pre(ii1,ii2)<0.5) then
        pre(ii1,ii2)=pre(ii1,ii2)+1.0
        dist=dist+((ii1-ii10)**2.0+(ii2-ii20)**2.0)**0.5*h
    endif
    ii10=ii1
    ii20=ii2
enddo
ii1=nzs
ii2=nxs
if (pre(ii1,ii2)<0.5) then
    pre(ii1,ii2)=pre(ii1,ii2)+1.0
    dist=dist+((ii1-ii10)**2.0+(ii2-ii20)**2.0)**0.5*h
endif
!write(*,*)dist,tdiff/dist
if (dist>0.001) then
    ra=tdiff/dist/sum(pre)*pre
endif
end subroutine ray2

end module module_ref2d
