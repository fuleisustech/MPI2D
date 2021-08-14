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

subroutine ray3(nz,ny,nx,xr,yr,zr,xs,ys,zs,h,ttt,ra,pre,tdiff)
implicit none
integer, intent(in) :: nz,ny,nx
real,    intent(in) :: xr,yr,zr,xs,ys,zs,h,tdiff
real,    intent(in) :: ttt(:,:,:)
real,    intent(out):: ra(:,:,:),pre(:,:,:)
integer :: nxr,nyr,nzr,nxs,nys,nzs,nmax,kn,ii1,ii2,ii3,&
   ii10,ii20,ii30,nx0,ny0,nz0
real    :: x0,y0,z0,tol,dist,vx,vy,vz,va

ra=0.0
pre=0.0
dist=0.0

tol=(1.5*h)**2.0
kn=1
nmax=nz+ny+nx

x0=xr
y0=yr
z0=zr
nxr=nint(xr/h)+1
nyr=nint(yr/h)+1
nzr=nint(zr/h)+1
nxs=nint(xs/h)+1
nys=nint(ys/h)+1
nzs=nint(zs/h)+1

dist=((xr-(nxr-1)*h)**2.0+(yr-(nyr-1)*h)**2.0+(xr-(nxr-1)*h)**2.0)**0.5

ii1=nzr
ii2=nyr
ii3=nxr
pre(ii1,ii2,ii3)=pre(ii1,ii2,ii3)+1.0
ii10=ii1
ii20=ii2
ii30=ii3

do while ((x0-xs)**2.0+(y0-ys)**2.0+(z0-zs)**2.0 > tol .and. kn<=nmax)
    kn=kn+1
    nx0=nint(x0/h)+1
    ny0=nint(y0/h)+1
    nz0=nint(z0/h)+1
    vx=ttt(nz0,ny0,nx0+1)-ttt(nz0,ny0,nx0-1)
    vy=ttt(nz0,ny0+1,nx0)-ttt(nz0,ny0-1,nx0)
    vz=ttt(nz0+1,ny0,nx0)-ttt(nz0-1,ny0,nx0)
    va=(vx*vx+vy*vy+vz*vz)**0.5
    vx=-h*vx/va
    vy=-h*vy/va
    vz=-h*vz/va
    x0=x0+vx
    y0=y0+vy
    z0=z0+vz
    ii1=nint(z0/h)+1
    ii2=nint(y0/h)+1
    ii3=nint(x0/h)+1

    ii1=min(ii1,nz)
    ii1=max(ii1,1)
    ii2=min(ii2,ny)
    ii2=max(ii2,1)
    ii3=min(ii3,nx)
    ii3=max(ii3,1)

    if (pre(ii1,ii2,ii3)<0.5) then
        pre(ii1,ii2,ii3)=pre(ii1,ii2,ii3)+1.0
        dist=dist+((ii1-ii10)**2.0+(ii2-ii20)**2.0+(ii3-ii30)**2.0)**0.5*h
    endif
    ii10=ii1
    ii20=ii2
    ii30=ii3
enddo
ii1=nzs
ii2=nys
ii3=nxs
if (pre(ii1,ii2,ii3)<0.5) then
    pre(ii1,ii2,ii3)=pre(ii1,ii2,ii3)+1.0
    dist=dist+((ii1-ii10)**2.0+(ii2-ii20)**2.0+(ii3-ii30)**2.0)**0.5*h
endif
if (dist>0.001) then
    ra=tdiff/dist/sum(pre)*pre
endif

end subroutine ray3

end module module_ref3d
