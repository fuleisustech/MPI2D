module module_aperture3d
use module_datatype
use module_array
use module_utility
use module_source
use module_fftw
implicit none
real,    private :: xmin, xmax, ymin,ymax, zmin,zmax
integer, private :: ix1,ix2,iy1,iy2,iz1,iz2,npml,nx,ny,nz,nzw,nx_pml,ny_pml,nz_pml,nzw_pml

interface aperture_init
   module procedure init_shot_fdm
   module procedure init_shot_img_fdm
   module procedure init_ms_shot_fdm
   module procedure init_ms_shot_img_fdm
   module procedure init_eik
   module procedure init_ref
   module procedure init_pw1_shot
   module procedure init_pw2_shot
!   module procedure init_eik_shot
!   module procedure init_eik_shot_img
   module procedure init_shot_kmm
   module procedure init_shot_img_kmm
   module procedure init_pw1_shot_kmm
   module procedure init_pw1_shot_img_kmm
   module procedure init_pw2_shot_kmm
   module procedure init_pw2_shot_img_kmm
end interface aperture_init

interface aperture_finalize
   module procedure finalize_shot_fdm
   module procedure finalize_shot_img_fdm
   module procedure finalize_ms_shot_fdm
   module procedure finalize_ms_shot_img_fdm
   module procedure finalize_eik
!   module procedure finalize_eik_shot
!   module procedure finalize_eik_shot_img
   module procedure finalize_shot_kmm
   module procedure finalize_shot_img_kmm
   module procedure finalize_pw1_shot_kmm
   module procedure finalize_pw1_shot_img_kmm
   module procedure finalize_pw2_shot_kmm
   module procedure finalize_pw2_shot_img_kmm
!   module procedure finalize_pwk
end interface aperture_finalize

interface aperture_image_l2g
   module procedure aperture_image_l2g_real
   module procedure aperture_image_l2g_int
end interface aperture_image_l2g

private :: nz_ny_nx

contains

subroutine init_shot_fdm(c,a,m,f,s,o,is,sh,fdm)
type(coord3d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(fd3d),           intent(in) :: f
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot3d),         intent(out) :: sh
type(fdmod3d),        intent(out) :: fdm
call init_shot(c,a,m,o,is,sh)
call init_fdmod(c,a,m,f,s,is,fdm)
end subroutine init_shot_fdm

subroutine finalize_shot_fdm(sh,fdm)
type(shot3d),         intent(inout) :: sh
type(fdmod3d),        intent(inout) :: fdm
call finalize_shot(sh)
call finalize_fdmod(fdm)
end subroutine finalize_shot_fdm

subroutine init_shot_img_fdm(c,a,m,f,s,o,is,sh,i,fdm)
type(coord3d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(fd3d),           intent(in) :: f
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot3d),         intent(out) :: sh
type(image3d),        intent(out) :: i
type(fdmod3d),        intent(out) :: fdm
call init_shot(c,a,m,o,is,sh)
call init_fdmod(c,a,m,f,s,is,fdm)
call init_image(c,a,m,o,is,i)
end subroutine init_shot_img_fdm

subroutine finalize_shot_img_fdm(sh,i,fdm)
type(shot3d),         intent(inout) :: sh
type(image3d),        intent(inout) :: i
type(fdmod3d),        intent(inout) :: fdm
call finalize_shot(sh)
call finalize_image(i)
call finalize_fdmod(fdm)
end subroutine finalize_shot_img_fdm

subroutine init_ms_shot_fdm(c,a,m,f,s,o,is,sh,fdm)
type(ms_coord3d),     intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(fd3d),           intent(in) :: f
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(ms_shot3d),      intent(out) :: sh
type(fdmod3d),        intent(out) :: fdm
call init_ms_shot(c,a,m,o,is,sh)
call init_ms_fdmod(c,a,m,f,s,is,fdm)
end subroutine init_ms_shot_fdm

subroutine finalize_ms_shot_fdm(sh,fdm)
type(ms_shot3d),      intent(inout) :: sh
type(fdmod3d),        intent(inout) :: fdm
call finalize_ms_shot(sh)
call finalize_fdmod(fdm)
end subroutine finalize_ms_shot_fdm

subroutine init_ms_shot_img_fdm(c,a,m,f,s,o,is,sh,i,fdm)
type(ms_coord3d),     intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(fd3d),           intent(in) :: f
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(ms_shot3d),      intent(out) :: sh
type(image3d),        intent(out) :: i
type(fdmod3d),        intent(out) :: fdm
call init_ms_shot(c,a,m,o,is,sh)
call init_ms_fdmod(c,a,m,f,s,is,fdm)
call init_ms_image(c,a,m,o,is,i)
end subroutine init_ms_shot_img_fdm

subroutine finalize_ms_shot_img_fdm(sh,i,fdm)
type(ms_shot3d),      intent(inout) :: sh
type(image3d),        intent(inout) :: i
type(fdmod3d),        intent(inout) :: fdm
call finalize_ms_shot(sh)
call finalize_image(i)
call finalize_fdmod(fdm)
end subroutine finalize_ms_shot_img_fdm

subroutine init_shot_kmm(c,a,m,k,s,o,is,sh,kmm)
type(coord3d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(km3d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot3d),         intent(out) :: sh
type(kmmod3d),        intent(out) :: kmm
call init_shot(c,a,m,o,is,sh)
call init_kmmod(c,a,m,k,s,is,kmm)
end subroutine init_shot_kmm

subroutine finalize_shot_kmm(sh,kmm)
type(shot3d),         intent(inout) :: sh
type(kmmod3d),        intent(inout) :: kmm
call finalize_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_shot_kmm

!subroutine init_eik_shot_img(c,a,m,o,is,e,sh,i)
!type(coord3d),        intent(in) :: c
!type(aperture),       intent(in) :: a
!type(model3d),        intent(in) :: m
!type(other),          intent(in) :: o
!integer,              intent(in) :: is
!type(eik3d),          intent(out) :: e
!type(shot3d),         intent(out) :: sh
!type(image3d),        intent(out) :: i
!call init_eik(c,a,m,is,e)
!call init_shot(c,a,m,o,is,sh)
!call init_image(c,a,m,o,is,i)
!end subroutine init_eik_shot_img

!subroutine finalize_eik_shot_img(e,sh,i)
!type(eik3d),          intent(inout) :: e
!type(shot3d),         intent(inout) :: sh
!type(image3d),        intent(inout) :: i
!call finalize_eik(e)
!call finalize_shot(sh)
!call finalize_image(i)
!end subroutine finalize_eik_shot_img

!subroutine init_eik_shot(c,a,m,o,is,e,sh)
!type(coord3d),        intent(in) :: c
!type(aperture),       intent(in) :: a
!type(model3d),        intent(in) :: m
!type(other),          intent(in) :: o
!integer,              intent(in) :: is
!type(eik3d),          intent(out) :: e
!type(shot3d),         intent(out) :: sh
!call init_eik(c,a,m,is,e)
!call init_shot(c,a,m,o,is,sh)
!end subroutine init_eik_shot

!subroutine finalize_eik_shot(e,sh)
!type(eik3d),          intent(inout) :: e
!type(shot3d),         intent(inout) :: sh
!call finalize_eik(e)
!call finalize_shot(sh)
!end subroutine finalize_eik_shot

subroutine init_shot_img_kmm(c,a,m,k,s,o,is,sh,i,kmm)
type(coord3d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(km3d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot3d),         intent(out) :: sh
type(image3d),        intent(out) :: i
type(kmmod3d),        intent(out) :: kmm
call init_shot(c,a,m,o,is,sh)
call init_image(c,a,m,o,is,i)
call init_kmmod(c,a,m,k,s,is,kmm)
end subroutine init_shot_img_kmm

subroutine finalize_shot_img_kmm(sh,i,kmm)
type(shot3d),         intent(inout) :: sh
type(image3d),        intent(inout) :: i
type(kmmod3d),        intent(inout) :: kmm
call finalize_image(i)
!deallocate(i%img)
call finalize_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_shot_img_kmm

subroutine init_pw1_shot_kmm(c,a,m,k,s,o,is,sh,kmm)
type(pw1_coord3d),    intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(km3d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(pw1_shot3d),         intent(out) :: sh
type(kmmod3d),        intent(out) :: kmm
call init_pw1_shot(c,a,m,o,is,sh)
call init_pw1_kmmod(c,a,m,k,s,is,kmm)
end subroutine init_pw1_shot_kmm

subroutine finalize_pw1_shot_kmm(sh,kmm)
type(pw1_shot3d),         intent(inout) :: sh
type(kmmod3d),        intent(inout) :: kmm
call finalize_pw1_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_pw1_shot_kmm

subroutine init_pw1_shot_img_kmm(c,a,m,k,s,o,is,sh,i,kmm)
type(pw1_coord3d),    intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(km3d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(pw1_shot3d),         intent(out) :: sh
type(image3d),        intent(out) :: i
type(kmmod3d),        intent(out) :: kmm
call init_pw1_shot(c,a,m,o,is,sh)
call init_pw1_image(c,a,m,o,is,i)
call init_pw1_kmmod(c,a,m,k,s,is,kmm)
end subroutine init_pw1_shot_img_kmm

subroutine finalize_pw1_shot_img_kmm(sh,i,kmm)
type(pw1_shot3d),         intent(inout) :: sh
type(image3d),        intent(inout) :: i
type(kmmod3d),        intent(inout) :: kmm
call finalize_image(i)
!deallocate(i%img)
call finalize_pw1_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_pw1_shot_img_kmm

subroutine init_pw2_shot_kmm(c,a,m,k,s,o,is,sh,kmm)
type(pw2_coord3d),    intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(km3d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(pw2_shot3d),     intent(out) :: sh
type(kmmod3d),        intent(out) :: kmm
call init_pw2_shot(c,a,m,o,is,sh)
call init_pw2_kmmod(a,m,k,s,kmm)
end subroutine init_pw2_shot_kmm

subroutine finalize_pw2_shot_kmm(sh,kmm)
type(pw2_shot3d),         intent(inout) :: sh
type(kmmod3d),        intent(inout) :: kmm
call finalize_pw2_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_pw2_shot_kmm

subroutine init_pw2_shot_img_kmm(c,a,m,k,s,o,is,sh,i,kmm)
type(pw2_coord3d),    intent(in) :: c
type(aperture),       intent(in) :: a
type(model3d),        intent(in) :: m
type(km3d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(pw2_shot3d),     intent(out) :: sh
type(image3d),        intent(out) :: i
type(kmmod3d),        intent(out) :: kmm
call init_pw2_shot(c,a,m,o,is,sh)
call init_pw2_image(a,m,o,i)
call init_pw2_kmmod(a,m,k,s,kmm)
end subroutine init_pw2_shot_img_kmm

subroutine finalize_pw2_shot_img_kmm(sh,i,kmm)
type(pw2_shot3d),         intent(inout) :: sh
type(image3d),        intent(inout) :: i
type(kmmod3d),        intent(inout) :: kmm
call finalize_image(i)
call finalize_pw2_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_pw2_shot_img_kmm

!=====================================================================

subroutine coord_min_max(c,a,m,is)
type(coord3d),   intent(in)    :: c
type(aperture),  intent(in)    :: a
type(model3d),   intent(in)    :: m
integer, intent(in) :: is
! get the xmin and xmax for the receiver
xmax = maxval(c%xg(1:c%ng(is),is)); xmin = minval(c%xg(1:c%ng(is),is))
! get the xmin and xmax for source and receiver
xmax = max(xmax, maxval(c%xs(1:c%ng(is),is)))
xmin = min(xmin, minval(c%xs(1:c%ng(is),is)))
! get the ymin and ymax for the receiver
ymax = maxval(c%yg(1:c%ng(is),is)); ymin = minval(c%yg(1:c%ng(is),is))
! get the ymin and ymax for source and receiver
ymax = max(ymax, maxval(c%ys(1:c%ng(is),is)))
ymin = min(ymin, minval(c%ys(1:c%ng(is),is)))
zmin=a%z0; zmax=a%z0+(m%nz-1)*m%dz
call nz_ny_nx(nz,ny,nx,a%aper_z,a%aper_y,a%aper_x,a%z0,a%y0,a%x0,m%dz,m%dy,m%dx)
end subroutine coord_min_max

subroutine ms_coord_min_max(c,a,m,isg)
type(ms_coord3d),intent(in)    :: c
type(aperture),  intent(in)    :: a
type(model3d),   intent(in)    :: m
integer, intent(in) :: isg
! get the xmin and xmax for the receiver
xmax = maxval(c%xg(1:c%ng(isg),isg)); xmin = minval(c%xg(1:c%ng(isg),isg))
! get the xmin and xmax for source and receiver
xmax = max(xmax, maxval(c%xs(1:c%ns(isg),isg)))
xmin = min(xmin, minval(c%xs(1:c%ns(isg),isg)))
! get the ymin and ymax for the receiver
ymax = maxval(c%yg(1:c%ng(isg),isg)); ymin = minval(c%yg(1:c%ng(isg),isg))
! get the ymin and ymax for source and receiver
ymax = max(ymax, maxval(c%ys(1:c%ns(isg),isg)))
ymin = min(ymin, minval(c%ys(1:c%ns(isg),isg)))
zmin=a%z0; zmax=a%z0+(m%nz-1)*m%dz
call nz_ny_nx(nz,ny,nx,a%aper_z,a%aper_y,a%aper_x,a%z0,a%y0,a%x0,m%dz,m%dy,m%dx)
end subroutine ms_coord_min_max

!=====================================================================

subroutine pw1_coord_min_max(c,a,m,is)
type(pw1_coord3d),intent(in)   :: c
type(aperture),  intent(in)    :: a
type(model3d),   intent(in)    :: m
integer, intent(in) :: is
! get the xmin and xmax for the receiver
xmax = maxval(c%xg(1:c%ng(is),is)); xmin = minval(c%xg(1:c%ng(is),is))
! get the ymin and ymax for the receiver
ymax = maxval(c%yg(1:c%ng(is),is)); ymin = minval(c%yg(1:c%ng(is),is))
zmin=a%z0; zmax=a%z0+(m%nz-1)*m%dz
call nz_ny_nx(nz,ny,nx,a%aper_z,a%aper_y,a%aper_x,a%z0,a%y0,a%x0,m%dz,m%dy,m%dx)
end subroutine pw1_coord_min_max

!=====================================================================

subroutine pw2_coord_min_max(a,m)
type(aperture),  intent(in)    :: a
type(model3d),   intent(in)    :: m
! get the xmin and xmax for the receiver
xmin = a%x0; xmax = xmin+(m%nx-1)*m%dx
ymin = a%y0; ymax = ymin+(m%ny-1)*m%dy
zmin = a%z0; zmax = zmin+(m%nz-1)*m%dz
call nz_ny_nx(nz,ny,nx,a%aper_z,a%aper_y,a%aper_x,a%z0,a%y0,a%x0,m%dz,m%dy,m%dx)
end subroutine pw2_coord_min_max

!=====================================================================

subroutine eik_coord_min_max(a,m)
type(aperture),   intent(in)    :: a
type(model3d),    intent(in)    :: m
xmin = a%x0;   xmax = xmin + (m%nx-1.0)*m%dx
ymin = a%y0;   ymax = ymin + (m%ny-1.0)*m%dy
zmin = a%z0;   zmax = zmin+(m%nz-1)*m%dz
call nz_ny_nx(nz,ny,nx,a%aper_z,a%aper_y,a%aper_x,a%z0,a%y0,a%x0,m%dz,m%dy,m%dx)
end subroutine eik_coord_min_max

!=====================================================================

subroutine init_eik(c,a,m,is,e)
type(eik_coord3d), intent(in) :: c
type(aperture),    intent(in) :: a
type(model3d),     intent(in) :: m
integer,           intent(in) :: is
type(eik3d),       intent(out):: e
call copytype(m,e%m)
if (a%isAper) then
   call eik_coord_min_max(a,m);   e%m%nz=nz;   e%m%ny=ny;   e%m%nx=nx
   call allocate_and_initial(e%m%v,e%ttt,nz,ny,nx)
   call cut(m%v,m%nz,iz1,iz2,nz,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,e%m%v)
else
   xmin=a%x0; ymin=a%y0; zmin=a%z0; nz=m%nz; ny=m%ny; nx=m%nx
   call allocate_and_initial(e%m%v,e%ttt,nz,ny,nx)
   call copy_array(m%v,e%m%v)
endif
e%xs=c%xs(is)-xmin; e%ys=c%ys(is)-ymin; e%zs=c%zs(is)-zmin
end subroutine init_eik

!=====================================================================

subroutine finalize_eik(e)
type(eik3d),       intent(inout) :: e
call deallocate_and_free(e%m%v)
if (allocated(e%ttt)) call deallocate_and_free(e%ttt)
end subroutine finalize_eik

!=====================================================================

subroutine init_ref(c,a,m,i,fs,is,r)
type(coord3d),     intent(in) :: c
type(aperture),    intent(in) :: a
type(model3d),     intent(in) :: m
type(inv),         intent(in) :: i
integer,           intent(in) :: is,fs(:,:)
type(ref3d),       intent(out):: r
real                          :: ssmax,diff
integer                       :: i2, i3
ssmax=1.0/i%vmin
call copy_variable(m%dx,r%dx,m%dy,r%dz,m%dz,r%dz)
!call copytype(m,r%m)
if (a%isAper) then
   !call init_xmin_xmax_zmin_zmax(c,a,m,is,e%m%nz,e%m%nx)
   !call allocate_and_initial(e%m%v,e%ttt,e%m%nz,e%m%nx)
   !call cut(m%v,m%nz,iz1,iz2,e%m%nz,m%nx,ix1,ix2,e%m%nx,e%m%v)
else
   xmin=a%x0-3*m%dx;   ymin=a%y0-3*m%dx;   zmin=a%z0-3*m%dz
   r%nx=m%nx+6;   r%ny=m%ny+6;   r%nz=m%nz+6
   call allocate_and_initial(r%v,r%ttt,r%nz,r%ny,r%nx)
   call padarray(m%v,m%nz,3,3,m%ny,3,3,m%nx,3,3,r%v,"coefficient",1.5)
   call allocate_and_initial(r%fs,r%ny,r%nx)
   call padarray(fs,m%ny,3,3,m%nx,3,3,r%fs,"replicate")
   r%fs=r%fs+3
   do i3=1,r%nx
      do i2=1,r%ny
         diff=(r%v(r%fs(i2,i3),i2,i3)-ssmax)/3.
         r%v(1:r%fs(i2,i3)-3,i2,i3)=ssmax
         r%v(r%fs(i2,i3)-2,i2,i3)=ssmax+diff
         r%v(r%fs(i2,i3)-1,i2,i3)=ssmax+2.*diff
      enddo
   enddo
   call allocate_and_initial(r%xs,r%ys,r%zs,r%xg,r%yg,r%zg,c%ng(is))
   r%ng=c%ng(is)
   r%xs=c%xs(:,is)-xmin;   r%ys=c%ys(:,is)-ymin;   r%zs=c%zs(:,is)-zmin
   r%xg=c%xg(:,is)-xmin;   r%yg=c%yg(:,is)-ymin;   r%zg=c%zg(:,is)-zmin
   r%ix1=4;   r%iy1=4;   r%iz1=4
   r%ix2=3+m%nx;   r%iy2=3+m%ny;   r%iz2=3+m%nz
endif
end subroutine init_ref

!=====================================================================

subroutine pt_init(v,xref,yref,zp,theta_x,theta_y,nx,ny,dx,t_init)
real,    intent(in) :: v(:,:,:),xref,yref,zp,theta_x,theta_y,dx
integer, intent(in) :: nx,ny
real,    intent(out):: t_init(:,:)
integer             :: ix, iy, izp
real                :: sin_thetax, sin_thetay, x, y
izp=nint(zp/dx)+1
sin_thetax=sin(theta_x/180.0*PI);   sin_thetay=sin(theta_y/180.0*PI)
do ix=1,nx
   x=(ix-1)*dx-xref
   do iy=1,ny
      y=(iy-1)*dx-yref
      t_init(iy,ix)=x*sin_thetax/v(izp,iy,ix)+y*sin_thetay/v(izp,iy,ix)
   enddo
enddo
end subroutine pt_init

!=====================================================================

subroutine init_shot(c,a,m,o,is,s)
type(coord3d),      intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model3d),      intent(in)    :: m
type(other),        intent(in)    :: o
integer ,           intent(in)    :: is
type(shot3d),       intent(inout) :: s
call copy_variable(a%z0,zmin,o%dt_out,s%dt)
s%nt=o%nt_out
if (a%isAper) then
   call coord_min_max(c,a,m,is)
else
   xmin=a%x0;   ymin=a%y0;   zmin=a%z0
endif
s%ng=c%ng(is)
call allocate_and_initial(s%xs,s%ys,s%zs,s%ng)
s%xs=c%xs(1:s%ng,is)-xmin; s%ys=c%ys(1:s%ng,is)-ymin; s%zs=c%zs(1:s%ng,is)-zmin
call allocate_and_initial(s%xg,s%yg,s%zg,s%ng)
s%xg=c%xg(1:s%ng,is)-xmin; s%yg=c%yg(1:s%ng,is)-ymin; s%zg=c%zg(1:s%ng,is)-zmin
call allocate_and_initial(s%sid,s%gid,s%ng)
s%sid=c%sid(1:s%ng,is); s%gid=c%gid(1:s%ng,is)
end subroutine init_shot

!=====================================================================

subroutine finalize_shot(s)
type(shot3d), intent(inout) :: s
deallocate(s%xs,s%ys,s%zs,s%xg,s%yg,s%zg,s%sid,s%gid)
if (allocated(s%seis)) deallocate(s%seis)
if (allocated(s%seis_u)) deallocate(s%seis_u)  ! modified by Bowen Guo
if (allocated(s%seis_v)) deallocate(s%seis_v)  ! modified by Bowen Guo
if (allocated(s%seis_w)) deallocate(s%seis_w)  ! modified by Bowen Guo

end subroutine finalize_shot

!====================================================================
subroutine init_ms_shot(c,a,m,o,isg,s)
type(ms_coord3d),intent(in)    :: c
type(aperture),  intent(in)    :: a
type(model3d),   intent(in)    :: m
type(other),     intent(in)    :: o
integer,         intent(in)    :: isg
type(ms_shot3d), intent(inout) :: s
call copy_variable(a%z0,zmin,o%dt_out,s%dt)
s%nt=o%nt_out
if (a%isAper) then
   call ms_coord_min_max(c,a,m,isg)
else
   xmin=a%x0;   ymin=a%y0;   ymin=a%z0
endif
s%ns=c%ns(isg);   s%ng=c%ng(isg)
call allocate_and_initial(s%xs,s%ys,s%zs,s%d,s%p,s%ns)
call allocate_and_initial(s%xg,s%yg,s%zg,s%ng)
s%d=c%d(1:s%ns,isg);   s%p=c%p(1:s%ns,isg)
s%xs=c%xs(1:s%ns,isg)-xmin; s%ys=c%ys(1:s%ns,isg)-ymin; s%zs=c%zs(1:s%ns,isg)-zmin
s%xg=c%xg(1:s%ng,isg)-xmin; s%yg=c%yg(1:s%ng,isg)-ymin; s%zg=c%zg(1:s%ng,isg)-zmin
end subroutine init_ms_shot

!=====================================================================

subroutine finalize_ms_shot(s)
type(ms_shot3d), intent(inout) :: s
deallocate(s%xs,s%ys,s%zs,s%xg,s%yg,s%zg)
if (allocated(s%seis)) deallocate(s%seis)
end subroutine finalize_ms_shot

!=====================================================================

subroutine init_pw1_shot(c,a,m,o,is,s)
type(pw1_coord3d),  intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model3d),      intent(in)    :: m
type(other),        intent(in)    :: o
integer ,           intent(in)    :: is
type(pw1_shot3d),   intent(inout) :: s
call copy_variable(a%z0,zmin,o%dt_out,s%dt)
call copy_variable(c%zp1(is),s%zp1)
call copy_variable(c%theta_x1(is),s%theta_x1,c%xref_1(is),s%xref_1)
call copy_variable(c%theta_y1(is),s%theta_y1,c%yref_1(is),s%yref_1)
call copy_variable(o%nt_out,s%nt,c%proid(is),s%proid)
if (allocated(c%ng)) call copy_variable(c%ng(is),s%ng)
if (a%isAper) then
   call pw1_coord_min_max(c,a,m,is)
else
   xmin=a%x0;   ymin=a%y0;   zmin=a%z0
endif
if (allocated(c%ng)) then
   call allocate_and_initial(s%xg,s%yg,s%zg,s%ng)
   s%xg=c%xg(1:s%ng,is)-xmin; s%yg=c%yg(1:s%ng,is)-ymin; s%zg=c%zg(1:s%ng,is)-zmin
   call allocate_and_initial(s%gid,s%ng)
   s%gid=c%gid(1:s%ng,is)
endif
end subroutine init_pw1_shot

!=====================================================================

subroutine finalize_pw1_shot(s)
type(pw1_shot3d), intent(inout) :: s
deallocate(s%xg,s%yg,s%zg,s%gid)
if (allocated(s%seis)) deallocate (s%seis)
end subroutine finalize_pw1_shot

!=====================================================================

subroutine init_pw2_shot(c,a,m,o,is,s)
type(pw2_coord3d),  intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model3d),      intent(in)    :: m
type(other),        intent(in)    :: o
integer,            intent(in)    :: is
type(pw2_shot3d),   intent(inout) :: s
call copy_variable(a%z0,zmin,o%dt_out,s%dt)
call copy_variable(c%theta_x1(is),s%theta_x1,c%xref_1(is),s%xref_1)
call copy_variable(c%theta_y1(is),s%theta_y1,c%yref_1(is),s%yref_1,c%zp1(is),s%zp1)
call copy_variable(o%nt_out,s%nt,c%proid(is),s%proid)
if (allocated(c%ng)) call copy_variable(c%ng(is),s%ng)
if (a%isAper) then
   call pw2_coord_min_max(a,m)
else
   xmin=a%x0;   ymin=a%y0;   zmin=a%z0
endif
if (allocated(c%ng)) then
   call allocate_and_initial(s%gid,s%ng)
   call allocate_and_initial(s%theta_x2,s%xref_2,s%theta_y2,s%yref_2,s%zp2,s%ng)
   s%gid=c%gid(1:s%ng,is)
   s%theta_x2=c%theta_x2(1:s%ng,is)
   s%theta_y2=c%theta_y2(1:s%ng,is)
   s%xref_2=c%xref_2(1:s%ng,is)
   s%yref_2=c%yref_2(1:s%ng,is)
   s%zp2=c%zp2(1:s%ng,is)
endif
end subroutine init_pw2_shot

!=====================================================================

subroutine finalize_pw2_shot(s)
type(pw2_shot3d), intent(inout) :: s
call deallocate_and_free(s%gid)
call deallocate_and_free(s%theta_x2,s%theta_y2,s%xref_2,s%yref_2,s%zp2)
if (allocated(s%seis)) deallocate (s%seis)
end subroutine finalize_pw2_shot

!=====================================================================

subroutine init_image(c,a,m,o,is,i)
type(coord3d),      intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model3d),      intent(in)    :: m
type(other),        intent(in)    :: o
integer,            intent(in)    :: is
type(image3d),      intent(inout) :: i
call copy_variable(a%z0,zmin,o%illum_order,i%illum_order)
call copy_variable(o%isIllum,i%isIllum,o%isIllumReceiver,i%isIllumReceiver,&
                   o%isPreStackImg,i%isPreStackImg)
if (a%isAper) then
   call coord_min_max(c,a,m,is)
   i%nz=nz; i%ny=ny; i%nx=nx; 
   i%iz1=iz1; i%iz2=iz2; i%iy1=iy1; i%iy2=iy2; i%ix1=ix1; i%ix2=ix2
else
   i%nz=m%nz; i%ny=m%ny; i%nx=m%nx; nz=m%nz; ny=m%ny; nx=m%nx
   i%iz1=1; i%iz2=nz; i%iy1=1; i%iy2=ny; i%ix1=1; i%ix2=nx
endif
call allocate_and_initial(i%img,nz,ny,nx)
if (o%isIllum) call allocate_and_initial(i%illum,nz,ny,nx)
end subroutine init_image

!=====================================================================

subroutine init_ms_image(c,a,m,o,isg,i)
type(ms_coord3d),  intent(in)    :: c
type(aperture),    intent(in)    :: a
type(model3d),     intent(in)    :: m
type(other),       intent(in)    :: o
integer,           intent(in)    :: isg
type(image3d),     intent(inout) :: i
call copy_variable(a%z0,zmin,o%illum_order,i%illum_order)
call copy_variable(o%isIllum,i%isIllum,o%isIllumReceiver,i%isIllumReceiver,&
                   o%isPreStackImg,i%isPreStackImg)
if (a%isAper) then
   call ms_coord_min_max(c,a,m,isg)
   i%nz=nz; i%ny=ny; i%nx=nx; 
   i%iz1=iz1; i%iz2=iz2; i%iy1=iy1; i%iy2=iy2; i%ix1=ix1; i%ix2=ix2
else
   i%nz=m%nz; i%ny=m%ny; i%nx=m%nx; nz=m%nz; ny=m%ny; nx=m%nx
   i%iz1=1; i%iz2=nz; i%iy1=1; i%iy2=ny; i%ix1=1; i%ix2=nx
endif
call allocate_and_initial(i%img,nz,ny,nx)
if (o%isIllum) call allocate_and_initial(i%illum,nz,ny,nx)
end subroutine init_ms_image

!=====================================================================

subroutine init_pw1_image(c,a,m,o,is,i)
type(pw1_coord3d),  intent(in)    :: c
type(aperture),    intent(in)    :: a
type(model3d),     intent(in)    :: m
type(other),       intent(in)    :: o
integer,           intent(in)    :: is
type(image3d),     intent(inout) :: i
call copy_variable(a%z0,zmin,o%illum_order,i%illum_order)
call copy_variable(o%isIllum,i%isIllum,o%isIllumReceiver,i%isIllumReceiver,&
                   o%isPreStackImg,i%isPreStackImg)
i%nz=m%nz
if (a%isAper) then
   call pw1_coord_min_max(c,a,m,is)
   i%nz=nz; i%ny=ny; i%nx=nx; 
   i%iz1=iz1; i%iz2=iz2; i%iy1=iy1; i%iy2=iy2; i%ix1=ix1; i%ix2=ix2
else
   i%nz=m%nz; i%ny=m%ny; i%nx=m%nx; nz=m%nz; ny=m%ny; nx=m%nx
   i%iz1=1; i%iz2=nz; i%iy1=1; i%iy2=ny; i%ix1=1; i%ix2=nx
endif
call allocate_and_initial(i%img,nz,ny,nx)
if (o%isIllum) call allocate_and_initial(i%illum,nz,ny,nx)
end subroutine init_pw1_image

!=====================================================================

subroutine init_pw2_image(a,m,o,i)
type(aperture),    intent(in)    :: a
type(model3d),     intent(in)    :: m
type(other),       intent(in)    :: o
type(image3d),     intent(inout) :: i
call copy_variable(a%z0,zmin,o%illum_order,i%illum_order)
call copy_variable(o%isIllum,i%isIllum,o%isIllumReceiver,i%isIllumReceiver,&
                   o%isPreStackImg,i%isPreStackImg)
i%nz=m%nz
if (a%isAper) then
   call pw2_coord_min_max(a,m)
   i%nz=nz; i%ny=ny; i%nx=nx; 
   i%iz1=iz1; i%iz2=iz2; i%iy1=iy1; i%iy2=iy2; i%ix1=ix1; i%ix2=ix2
else
   i%nz=m%nz; i%ny=m%ny; i%nx=m%nx; nz=m%nz; ny=m%ny; nx=m%nx
   i%iz1=1; i%iz2=nz; i%iy1=1; i%iy2=ny; i%ix1=1; i%ix2=nx
endif
call allocate_and_initial(i%img,nz,ny,nx)
if (o%isIllum) call allocate_and_initial(i%illum,nz,ny,nx)
end subroutine init_pw2_image

!====================================================================

subroutine finalize_image(i)
type(image3d), intent(inout) :: i
deallocate(i%img)
if (allocated(i%illum)) deallocate(i%illum)
end subroutine finalize_image

!=====================================================================

subroutine init_wf(c,a,m,f,is,wf)
type(coord3d),       intent(in)    :: c
type(aperture),      intent(in)    :: a
type(model3d),       intent(in)    :: m
type(fd3d),          intent(in)    :: f
integer,             intent(in)    :: is
type(wf3d),          intent(inout) :: wf
! setup normal parameters
wf%nt=f%nt
! setup parameters with aperture
if (a%isAper) then
   call coord_min_max(c,a,m,is)
else
   nx=m%nx;   ny=m%ny;   nz=m%nz
endif
wf%nx_2=nx+2;   wf%ny_2=ny+2;   wf%nz_2=nz+2
call allocate_and_initial(wf%pwf,wf%nz_2,wf%ny_2,wf%nx_2,wf%nt)
end subroutine init_wf

!=====================================================================

subroutine finalize_wf(wf)
type(wf3d),          intent(inout) :: wf
deallocate(wf%pwf)
end subroutine finalize_wf

!=====================================================================

subroutine init_fdmod(c,a,m,f,s,is,fdm)
type(coord3d),   intent(in) :: c
type(aperture),  intent(in) :: a
type(model3d),   intent(in) :: m
type(fd3d),      intent(in) :: f
type(source),    intent(in) :: s
integer,         intent(in) :: is
type(fdmod3d),   intent(inout) :: fdm
call copytype(f,fdm%f); call copytype(m,fdm%m); npml=f%npml
! setup source
call copytype(s,fdm%s)
call allocate_and_initial(fdm%s%sou,fdm%f%nt)
call wavelet_expand(s%sou,fdm%s%sou,fdm%f%nt)
! setup parameters with aperture
if (a%isAper) then
   call coord_min_max(c,a,m,is)
   fdm%m%nz=nz; fdm%m%ny=ny; fdm%m%nx=nx; 
   nx_pml=nx + 2*npml; ny_pml=ny + 2*npml; nz_pml=nz + 2*npml
   call allocate_and_initial(fdm%f%fs,ny_pml,nx_pml)
   call cut_and_pad(f%fs,m%ny,m%nx,iy1,iy2,ix1,ix2,npml,ny,nx,fdm%f%fs)
   if (allocated(m%v)) then
      call allocate_and_initial(fdm%m%v,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%v,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(fdm%m%refl,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%refl,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%refl)
   endif
   if (allocated(m%vp)) then
      call allocate_and_initial(fdm%m%vp,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%vp,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%vp)
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(fdm%m%vs,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%vs,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%vs)
   endif
   if (allocated(m%lambda)) then
      call allocate_and_initial(fdm%m%lambda,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%lambda,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%lambda)
   endif
   if (allocated(m%mu)) then
      call allocate_and_initial(fdm%m%mu,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%mu,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%mu)
   endif    
   if (allocated(m%den)) then
      call allocate_and_initial(fdm%m%den,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%den,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%den)
   endif
   
   
else
   nz=m%nz; ny=m%ny; nx=m%nx; nx_pml=nx+2*npml; ny_pml=ny+2*npml; nz_pml=nz+2*npml
   call allocate_and_initial(fdm%f%fs,ny_pml,nx_pml)
   call padarray(f%fs,ny,npml,npml,nx,npml,npml,fdm%f%fs,"replicate")
   if (allocated(m%v)) then
      call allocate_and_initial(fdm%m%v,nz_pml,ny_pml,nx_pml)
      call padarray(m%v,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%v,"replicate")
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(fdm%m%refl,nz_pml,ny_pml,nx_pml)
      call padarray(m%refl,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%refl,"replicate")
   endif
   if (allocated(m%vp)) then
      call allocate_and_initial(fdm%m%vp,nz_pml,ny_pml,nx_pml)
      call padarray(m%vp,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%vp,"replicate")
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(fdm%m%vs,nz_pml,ny_pml,nx_pml)
      call padarray(m%vs,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%vs,"replicate")
   endif
   if (allocated(m%lambda)) then
      call allocate_and_initial(fdm%m%lambda,nz_pml,ny_pml,nx_pml)
      call padarray(m%lambda,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%lambda,"replicate")
   endif
   if (allocated(m%mu)) then
      call allocate_and_initial(fdm%m%mu,nz_pml,ny_pml,nx_pml)
      call padarray(m%mu,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%mu,"replicate")
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(fdm%m%den,nz_pml,ny_pml,nx_pml)
      call padarray(m%den,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%den,"replicate")
   endif
endif
end subroutine init_fdmod

subroutine init_ms_fdmod(c,a,m,f,s,isg,fdm)
type(ms_coord3d),intent(in) :: c
type(aperture),  intent(in) :: a
type(model3d),   intent(in) :: m
type(fd3d),      intent(in) :: f
type(source),    intent(in) :: s
integer,         intent(in) :: isg
type(fdmod3d),   intent(inout) :: fdm
call copytype(f,fdm%f);   call copytype(m,fdm%m);  npml=f%npml
! setup source
call copytype(s,fdm%s)
call allocate_and_initial(fdm%s%sou,fdm%f%nt)
call wavelet_expand(s%sou,fdm%s%sou,fdm%f%nt)
! setup parameters with aperture
if (a%isAper) then
   call ms_coord_min_max(c,a,m,isg); fdm%m%nz=nz; fdm%m%ny=ny; fdm%m%nx=nx 
   nx_pml=nx + 2*npml; ny_pml=ny + 2*npml; nz_pml=nz + 2*npml
   call allocate_and_initial(fdm%f%fs,ny_pml,nx_pml)
   call cut_and_pad(f%fs,m%ny,m%nx,iy1,iy2,ix1,ix2,npml,ny,nx,fdm%f%fs)
   if (allocated(m%v)) then
      call allocate_and_initial(fdm%m%v,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%v,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(fdm%m%refl,nz_pml,ny_pml,nx_pml)
      call cut_and_pad(m%refl,m%nz,m%ny,m%nx,iz1,iz2,iy1,iy2,ix1,ix2,npml,&
                       nz,ny,nx,fdm%m%refl)
   endif
else  
   nz=m%nz; ny=m%ny; nx=m%nx; nx_pml=nx+2*npml; ny_pml=ny+2*npml; nz_pml=nz+2*npml
   call allocate_and_initial(fdm%f%fs,ny_pml,nx_pml)
   call padarray(f%fs,ny,npml,npml,nx,npml,npml,fdm%f%fs,"replicate")
   if (allocated(m%v)) then
      call allocate_and_initial(fdm%m%v,nz_pml,ny_pml,nx_pml)
      call padarray(m%v,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%v,"replicate")
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(fdm%m%refl,nz_pml,ny_pml,nx_pml)
      call padarray(m%refl,nz,npml,npml,ny,npml,npml,nx,npml,npml,fdm%m%refl,"replicate")
   endif
endif
end subroutine init_ms_fdmod

!====================================================================

subroutine finalize_fdmod(fdm)
type(fdmod3d), intent(inout) :: fdm
! modified by Bowen Guo
if (allocated(fdm%f%fs))deallocate(fdm%f%fs)
if (allocated(fdm%s%sou)) deallocate(fdm%s%sou)
if (allocated(fdm%m%v)) deallocate(fdm%m%v)
if (allocated(fdm%m%lambda)) deallocate(fdm%m%lambda)
if (allocated(fdm%m%mu)) deallocate(fdm%m%mu)
if (allocated(fdm%m%vp)) deallocate(fdm%m%vp)
if (allocated(fdm%m%vs)) deallocate(fdm%m%vs)
if (allocated(fdm%m%den)) deallocate(fdm%m%den)
! modified by Bowen Guo
if (allocated(fdm%m%refl)) deallocate(fdm%m%refl)
end subroutine finalize_fdmod

!=====================================================================

!subroutine init_water_fdmod3d(c,a,m,f,s,is,fdm)
!type(coord3d),   intent(in) :: c
!type(aperture),  intent(in) :: a
!type(water_model3d),   intent(in) :: m
!type(fd3d),      intent(in) :: f
!type(source),    intent(in) :: s
!integer,         intent(in) :: is
!type(water_fdmod3d),   intent(inout) :: fdm
!call copytype(f,fdm%f)
!call copytype(m,fdm%m)
!nz_pml = fdm%m%nz + 2*fdm%f%npml
!nzw_pml = fdm%m%nzw + 2*fdm%f%npml
!! setup source
!call copytype(s,fdm%s)
!call allocate_and_initial(fdm%s%sou,fdm%f%nt)
!call wavelet_expand(s%sou,fdm%s%sou,fdm%f%nt)
!! setup parameters with aperture
!if (a%isAper) then
!   call init_xmin_xmax_ymin_ymax_zmin_zmax(c,a,m%m,is,fdm%m%nx,fdm%m%ny)
!   nx_pml=fdm%m%nx + 2*fdm%f%npml
!   ny_pml=fdm%m%ny + 2*fdm%f%npml
!   call allocate_and_initial(fdm%f%fs,ny_pml,nx_pml)
!   call cut_and_pad(f%fs,m%nx,m%ny,ix1,ix2,iy1,iy2,f%npml,fdm%m%nx,fdm%m%ny,fdm%f%fs)
!   call allocate_and_initial(fdm%m%v,nz_pml,ny_pml,nx_pml)
!   call cut_and_pad(m%v,m%nx,m%ny,fdm%m%nz,ix1,ix2,iy1,iy2,fdm%f%npml,&
!                    fdm%m%nx,fdm%m%ny,fdm%m%v)
!   call allocate_and_initial(fdm%m%v_w,nzw_pml,ny_pml,nx_pml)
!   call cut_and_pad(m%v_w,m%nx,m%ny,fdm%m%nzw,ix1,ix2,iy1,iy2,fdm%f%npml,&
!                    fdm%m%nx,fdm%m%ny,fdm%m%v_w)
!   if (allocated(m%v_h)) then
!      call allocate_and_initial(fdm%m%v_h,nzw_pml,ny_pml,nx_pml)
!      call cut_and_pad(m%v_h,m%nx,m%ny,fdm%m%nzw,ix1,ix2,iy1,iy2,fdm%f%npml,&
!                       fdm%m%nx,fdm%m%ny,fdm%m%v_h)
!   endif
!   if (allocated(m%refl)) then
!      call allocate_and_initial(fdm%m%refl,nz_pml,ny_pml,nx_pml)
!      call cut_and_pad(m%refl,m%nx,m%ny,fdm%m%nz,ix1,ix2,iy1,iy2,fdm%f%npml,&
!                       fdm%m%nx,fdm%m%ny,fdm%m%refl)
!   endif
!   if (allocated(m%refl_fs)) then
!      call allocate_and_initial(fdm%m%refl_fs,nzw_pml,ny_pml,nx_pml)
!      call cut_and_pad(m%refl_fs,m%nx,m%ny,fdm%m%nzw,ix1,ix2,iy1,iy2,fdm%f%npml,&
!                       fdm%m%nx,fdm%m%ny,fdm%m%refl_fs)
!   endif
!   if (allocated(m%refl_ob)) then
!      call allocate_and_initial(fdm%m%refl_ob,nzw_pml,ny_pml,nx_pml)
!      call cut_and_pad(m%refl_ob,m%nx,m%ny,fdm%m%nzw,ix1,ix2,iy1,iy2,fdm%f%npml,&
!                       fdm%m%nx,fdm%m%ny,fdm%m%refl_ob)
!   endif
!else
!   nx_pml=fdm%m%nx+2*fdm%f%npml
!   ny_pml=fdm%m%nx+2*fdm%f%npml
!   call allocate_and_initial(fdm%f%fs,ny_pml,nx_pml)
!   call Padding(f%fs,fdm%m%ny,f%npml,f%npml,fdm%m%nx,f%npml,f%npml,fdm%f%fs)
!   call allocate_and_initial(fdm%m%v,nz_pml,ny_pml,nx_pml)
!   call Padding(m%v,m%nz,fdm%f%npml,fdm%f%npml,fdm%m%ny,fdm%f%npml,fdm%f%npml,&
!                fdm%m%nx,fdm%f%npml,fdm%f%npml,fdm%m%v)
!   call allocate_and_initial(fdm%m%v_w,nz_pml,ny_pml,nx_pml)
!   call Padding(m%v_w,m%nzw,fdm%f%npml,fdm%f%npml,fdm%m%ny,fdm%f%npml,fdm%f%npml,&
!                fdm%m%nx,fdm%f%npml,fdm%f%npml,fdm%m%v_w)
!   if (allocated(m%refl)) then
!      call allocate_and_initial(fdm%m%refl,nz_pml,ny_pml,nx_pml)
!      call Padding(m%refl,m%nz,fdm%f%npml,fdm%f%npml,fdm%m%ny,fdm%f%npml,fdm%f%npml,&
!                   fdm%m%nx,fdm%f%npml,fdm%f%npml,fdm%m%refl)
!   endif
!   if (allocated(m%v_h)) then
!      call allocate_and_initial(fdm%m%v_h,nzw_pml,ny_pml,nx_pml)
!      call Padding(m%v_h,m%nzw,fdm%f%npml,fdm%f%npml,fdm%m%ny,fdm%f%npml,fdm%f%npml,&
!                   fdm%m%nx,fdm%f%npml,fdm%f%npml,fdm%m%v_h)
!   endif
!   if (allocated(m%refl_fs)) then
!      call allocate_and_initial(fdm%m%refl_fs,nzw_pml,ny_pml,nx_pml)
!      call Padding(m%refl_fs,m%nzw,fdm%f%npml,fdm%f%npml,fdm%m%ny,fdm%f%npml,fdm%f%npml,&
!                   fdm%m%nx,fdm%f%npml,fdm%f%npml,fdm%m%refl_fs)
!   endif
!   if (allocated(m%refl_ob)) then
!      call allocate_and_initial(fdm%m%refl_ob,nzw_pml,ny_pml,nx_pml)
!      call Padding(m%refl_ob,m%nzw,fdm%f%npml,fdm%f%npml,fdm%m%ny,fdm%f%npml,fdm%f%npml,&
!                   fdm%m%nx,fdm%f%npml,fdm%f%npml,fdm%m%refl_ob)
!   endif
!endif
!end subroutine init_water_fdmod3d

!====================================================================

!subroutine finalize_water_fdmod3d(fdm)
!type(water_fdmod3d), intent(inout) :: fdm
!deallocate(fdm%f%fs,fdm%s%sou,fdm%m%v,fdm%m%v_w)
!if (allocated(fdm%m%refl)) deallocate(fdm%m%v_h)
!if (allocated(fdm%m%refl)) deallocate(fdm%m%refl)
!if (allocated(fdm%m%refl)) deallocate(fdm%m%refl_fs)
!if (allocated(fdm%m%refl)) deallocate(fdm%m%refl_ob)
!end subroutine finalize_water_fdmod3d

!=====================================================================

subroutine init_kmmod(c,a,m,k,s,is,kmm)
type(coord3d),   intent(in) :: c
type(aperture),  intent(in) :: a
type(model3d),   intent(in) :: m
type(km3d),      intent(in) :: k
type(source),    intent(in) :: s
integer,         intent(in) :: is
type(kmmod3d),   intent(inout) :: kmm
call copytype(m,kmm%m);   call copytype(k,kmm%k)
! setup source
call copytype(s,kmm%s)
kmm%k%nt_conv=kmm%k%nt+s%nw-1
call allocate_and_initial(kmm%k%wfft,kmm%k%nt_conv)
call fft_wavelet(s%sou,kmm%k%wfft,s%nw,kmm%k%nt_conv)
! setup parameters with aperture
if (a%isAper) then
   call coord_min_max(c,a,m,is); kmm%m%nz=nz; kmm%m%ny=ny; kmm%m%nx=nx
   call allocate_and_initial(kmm%k%z0_mig,ny,nx)
   call cut(k%z0_mig,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,ny,nx)
      call cut(m%v,m%nz,iz1,iz2,nz,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,ny,nx)
      call cut(m%refl,m%nz,iz1,iz2,nz,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%m%refl)
   endif
else
   nz=m%nz; ny=m%ny; nx=m%nx
   call allocate_and_initial(kmm%k%z0_mig,ny,nx)
   call copy_array(k%z0_mig,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,ny,nx)
      call copy_array(m%v,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,ny,nx)
      call copy_array(m%refl,kmm%m%refl)
   endif
endif
end subroutine init_kmmod

!====================================================================

subroutine init_pw1_kmmod(c,a,m,k,s,is,kmm)
type(pw1_coord3d),   intent(in) :: c
type(aperture),  intent(in) :: a
type(model3d),   intent(in) :: m
type(km3d),      intent(in) :: k
type(source),    intent(in) :: s
integer,         intent(in) :: is
type(kmmod3d),   intent(inout) :: kmm
call copytype(m,kmm%m);   call copytype(k,kmm%k)
! setup source
call copytype(s,kmm%s)
kmm%k%nt_conv=kmm%k%nt+s%nw-1
call allocate_and_initial(kmm%k%wfft,kmm%k%nt_conv)
call fft_wavelet(s%sou,kmm%k%wfft,s%nw,kmm%k%nt_conv)
! setup parameters with aperture
if (a%isAper) then
   call pw1_coord_min_max(c,a,m,is); kmm%m%nz=nz; kmm%m%ny=ny; kmm%m%nx=nx
   call allocate_and_initial(kmm%k%z0_mig,ny,nx)
   call cut(k%z0_mig,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,ny,nx)
      call cut(m%v,m%nz,iz1,iz2,nz,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,ny,nx)
      call cut(m%refl,m%nz,iz1,iz2,nz,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%m%refl)
   endif
else
   nz=m%nz; ny=m%ny; nx=m%nx
   call allocate_and_initial(kmm%k%z0_mig,ny,nx)
   call copy_array(k%z0_mig,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,ny,nx)
      call copy_array(m%v,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,ny,nx)
      call copy_array(m%refl,kmm%m%refl)
   endif
endif
end subroutine init_pw1_kmmod

subroutine init_pw2_kmmod(a,m,k,s,kmm)
type(aperture),  intent(in) :: a
type(model3d),   intent(in) :: m
type(km3d),      intent(in) :: k
type(source),    intent(in) :: s
type(kmmod3d),   intent(inout) :: kmm
call copytype(m,kmm%m);   call copytype(k,kmm%k)
! setup source
call copytype(s,kmm%s)
kmm%k%nt_conv=kmm%k%nt+s%nw-1
call allocate_and_initial(kmm%k%wfft,kmm%k%nt_conv)
call fft_wavelet(s%sou,kmm%k%wfft,s%nw,kmm%k%nt_conv)
! setup parameters with aperture
if (a%isAper) then
   call pw2_coord_min_max(a,m); kmm%m%nz=nz; kmm%m%ny=ny; kmm%m%nx=nx
   call allocate_and_initial(kmm%k%z0_mig,ny,nx)
   call cut(k%z0_mig,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,ny,nx)
      call cut(m%v,m%nz,iz1,iz2,nz,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,ny,nx)
      call cut(m%refl,m%nz,iz1,iz2,nz,m%ny,iy1,iy2,ny,m%nx,ix1,ix2,nx,kmm%m%refl)
   endif
else
   nz=m%nz; ny=m%ny; nx=m%nx
   call allocate_and_initial(kmm%k%z0_mig,ny,nx)
   call copy_array(k%z0_mig,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,ny,nx)
      call copy_array(m%v,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,ny,nx)
      call copy_array(m%refl,kmm%m%refl)
   endif
endif
end subroutine init_pw2_kmmod

!=====================================================================

subroutine finalize_kmmod(kmm)
type(kmmod3d), intent(inout) :: kmm
deallocate(kmm%k%wfft,kmm%k%z0_mig)
if (allocated(kmm%m%v)) deallocate(kmm%m%v)
if (allocated(kmm%m%refl)) deallocate(kmm%m%refl)
end subroutine finalize_kmmod

!=====================================================================

subroutine aperture_image_l2g_real(isAper,img,iz1,iz2,iy1,iy2,ix1,ix2,img_out,nz,ny,nx)
logical,optional,   intent(in) :: isAper
real,               intent(in) :: img(:,:,:)
integer,            intent(in) :: iz1,iz2,iy1,iy2,ix1,ix2,nz,ny,nx
real,               intent(out) :: img_out(:,:,:)
integer ::ix1_in,ix2_in,iy1_in,iy2_in,iz1_in,iz2_in,&
          ix1_out,ix2_out,iy1_out,iy2_out,iz1_out,iz2_out
if (present(isAper).and. isAper) then
   img_out=0.0
   call l2g_i1i2(ix1,ix2,nx,ix1_in,ix2_in,ix1_out,ix2_out)
   call l2g_i1i2(iy1,iy2,ny,iy1_in,iy2_in,iy1_out,iy2_out)
   call l2g_i1i2(iz1,iz2,nz,iz1_in,iz2_in,iz1_out,iz2_out)
   img_out(iz1_out:iz2_out,iy1_out:iy2_out,ix1_out:ix2_out) &
   =img(iz1_in:iz2_in,iy1_in:iy2_in,ix1_in:ix2_in)
else
   img_out=img
endif
end subroutine aperture_image_l2g_real

subroutine aperture_image_l2g_int(isAper,img,iz1,iz2,iy1,iy2,ix1,ix2,img_out,nz,ny,nx)
logical,optional,   intent(in) :: isAper
integer,            intent(in) :: img(:,:,:)
integer,            intent(in) :: iz1,iz2,iy1,iy2,ix1,ix2,nz,ny,nx
integer,            intent(out) :: img_out(:,:,:)
integer ::ix1_in,ix2_in,iy1_in,iy2_in,iz1_in,iz2_in,&
          ix1_out,ix2_out,iy1_out,iy2_out,iz1_out,iz2_out
if (present(isAper).and.isAper) then
   img_out=0.0
   call l2g_i1i2(ix1,ix2,nx,ix1_in,ix2_in,ix1_out,ix2_out)
   call l2g_i1i2(iy1,iy2,ny,iy1_in,iy2_in,iy1_out,iy2_out)
   call l2g_i1i2(iz1,iz2,nz,iz1_in,iz2_in,iz1_out,iz2_out)
   img_out(iz1_out:iz2_out,iy1_out:iy2_out,ix1_out:ix2_out) &
   =img(iz1_in:iz2_in,iy1_in:iy2_in,ix1_in:ix2_in)
else
   img_out=img
endif
end subroutine aperture_image_l2g_int

subroutine l2g_i1i2(i1,i2,n,i1_in,i2_in,i1_out,i2_out)
integer, intent(in) :: i1,i2,n
integer, intent(out):: i1_in,i2_in,i1_out,i2_out
if (i1<=0) then
   i1_in=-i1+2;   i1_out=1
   if (i2<=n) then
      i2_in=-i1+1+i2;   i2_out=i2
   else
      i2_in=-i1+1+n;    i2_out=n
   endif
else
   i1_in=1;   i1_out=i1
   if (i2<=n) then
      i2_in=-i1+1+i2;   i2_out=i2
   else
      i2_in=-i1+1+n;   i2_out=n
   endif
endif
end subroutine l2g_i1i2

subroutine nz_ny_nx(nz,ny,nx,aper_z,aper_y,aper_x,z0,y0,x0,dz,dy,dx)
integer, intent(out) :: nz,ny,nx
real,    intent(in)  :: aper_z,aper_y,aper_x,z0,y0,x0,dz,dy,dx
xmax = xmax + aper_x - x0; xmin = xmin - aper_x - x0
ix1 = int(xmin/dx) + 1; ix2 = nint(xmax/dx) + 1
ymax = ymax + aper_y - y0; ymin = ymin - aper_y - y0
iy1 = int(ymin/dy) + 1; iy2 = nint(ymax/dy) + 1
zmax = zmax + aper_z - z0;   zmin = zmin - aper_z - z0
iz1 = int(zmin/dz) + 1;   iz2 = nint(zmax/dz) + 1
nz = iz2 - iz1 + 1; ny = iy2 - iy1 + 1; nx = ix2 - ix1 + 1
end subroutine nz_ny_nx

!====================================================================

end module module_aperture3d
