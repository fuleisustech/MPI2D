module module_aperture2d
use module_datatype
use module_array
use module_utility
use module_source
use module_fftw
implicit none
real,    private :: xmin, xmax, zmin, zmax
integer, private :: ix1,ix2,iz1,iz2,npml,nx,nz,nzw,nx_pml,nz_pml,nzw_pml

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
!   module procedure init_pweik2d
   module procedure init_shot_kmm
   module procedure init_shot_img_kmm
   module procedure init_pw1_shot_kmm
   module procedure init_pw1_shot_img_kmm
   module procedure init_pw2_shot_kmm
   module procedure init_pw2_shot_img_kmm
   module procedure init_shot_gdmod_v1
   module procedure init_shot_gdmod_v2
end interface aperture_init

interface aperture_finalize
   module procedure finalize_kmmod
   module procedure finalize_shot
   module procedure finalize_shot_fdm
   module procedure finalize_shot_img_fdm
   module procedure finalize_ms_shot_fdm
   module procedure finalize_ms_shot_img_fdm
   module procedure finalize_eik
!   module procedure finalize_eik_shot
!   module procedure finalize_eik_shot_img
!   module procedure finalize_pweik2d
   module procedure finalize_shot_kmm
   module procedure finalize_shot_img_kmm
   module procedure finalize_pw1_shot_kmm
   module procedure finalize_pw1_shot_img_kmm
   module procedure finalize_pw2_shot_kmm
   module procedure finalize_pw2_shot_img_kmm
   module procedure finalize_shot_gdmod
end interface aperture_finalize

interface aperture_image_l2g
   module procedure aperture_image_l2g_real
   module procedure aperture_image_l2g_int
end interface aperture_image_l2g

private :: nz_nx

contains

subroutine init_shot_fdm(c,a,m,f,s,o,is,sh,fdm)
type(coord2d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(fd2d),           intent(in) :: f
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot2d),         intent(out) :: sh
type(fdmod2d),        intent(out) :: fdm
call init_shot(c,a,m,o,is,sh)
call init_fdmod(c,a,m,f,s,is,fdm)
end subroutine init_shot_fdm

subroutine finalize_shot_fdm(sh,fdm)
type(shot2d),         intent(inout) :: sh
type(fdmod2d),        intent(inout) :: fdm
call finalize_shot(sh)
call finalize_fdmod(fdm)
end subroutine finalize_shot_fdm

subroutine init_shot_img_fdm(c,a,m,f,s,o,is,sh,i,fdm)
type(coord2d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(fd2d),           intent(in) :: f
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot2d),         intent(out) :: sh
type(image2d),        intent(out) :: i
type(fdmod2d),        intent(out) :: fdm
call init_shot(c,a,m,o,is,sh)
call init_fdmod(c,a,m,f,s,is,fdm)
call init_image(c,a,m,o,is,i)
end subroutine init_shot_img_fdm

subroutine finalize_shot_img_fdm(sh,i,fdm)
type(shot2d),         intent(inout) :: sh
type(image2d),        intent(inout) :: i
type(fdmod2d),        intent(inout) :: fdm
call finalize_shot(sh)
call finalize_image(i)
call finalize_fdmod(fdm)
end subroutine finalize_shot_img_fdm

subroutine init_ms_shot_fdm(mc,a,m,f,s,o,isg,msh,fdm)
type(ms_coord2d),     intent(in) :: mc
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(fd2d),           intent(in) :: f
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: isg
type(ms_shot2d),      intent(out) :: msh
type(fdmod2d),        intent(out) :: fdm
call init_ms_shot(mc,a,m,o,isg,msh)
call init_ms_fdmod(mc,a,m,f,s,isg,fdm)
end subroutine init_ms_shot_fdm

subroutine finalize_ms_shot_fdm(msh,fdm)
type(ms_shot2d),      intent(inout) :: msh
type(fdmod2d),        intent(inout) :: fdm
call finalize_ms_shot(msh)
call finalize_fdmod(fdm)
end subroutine finalize_ms_shot_fdm

subroutine init_ms_shot_img_fdm(mc,a,m,f,s,o,isg,msh,i,fdm)
type(ms_coord2d),     intent(in) :: mc
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(fd2d),           intent(in) :: f
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: isg
type(ms_shot2d),      intent(out) :: msh
type(image2d),        intent(out) :: i
type(fdmod2d),        intent(out) :: fdm
call init_ms_shot(mc,a,m,o,isg,msh)
call init_ms_fdmod(mc,a,m,f,s,isg,fdm)
call init_ms_image(mc,a,m,o,isg,i)
end subroutine init_ms_shot_img_fdm

subroutine finalize_ms_shot_img_fdm(msh,i,fdm)
type(ms_shot2d),         intent(inout) :: msh
type(image2d),        intent(inout) :: i
type(fdmod2d),        intent(inout) :: fdm
call finalize_ms_shot(msh)
call finalize_image(i)
call finalize_fdmod(fdm)
end subroutine finalize_ms_shot_img_fdm

subroutine init_shot_kmm(c,a,m,k,s,o,is,sh,kmm)
type(coord2d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(km2d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot2d),         intent(out) :: sh
type(kmmod2d),        intent(out) :: kmm
call init_shot(c,a,m,o,is,sh)
call init_kmmod(c,a,m,k,s,is,kmm)
end subroutine init_shot_kmm

subroutine finalize_shot_kmm(sh,kmm)
type(shot2d),         intent(inout) :: sh
type(kmmod2d),        intent(inout) :: kmm
call finalize_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_shot_kmm

!subroutine init_eik_shot_img(c,a,m,o,is,e,sh,i)
!type(coord2d),        intent(in) :: c
!type(aperture),       intent(in) :: a
!type(model2d),        intent(in) :: m
!type(other),          intent(in) :: o
!integer,              intent(in) :: is
!type(eik2d),          intent(out) :: e
!type(shot2d),         intent(out) :: sh
!type(image2d),        intent(out) :: i
!call init_eik(c,a,m,is,e)
!call init_shot(c,a,m,o,is,sh)
!call init_image(c,a,m,o,is,i)
!end subroutine init_eik_shot_img

!subroutine finalize_eik_shot_img(e,sh,i)
!type(eik2d),          intent(inout) :: e
!type(shot2d),         intent(inout) :: sh
!type(image2d),        intent(inout) :: i
!call finalize_eik(e)
!call finalize_shot(sh)
!!call finalize_image(i)
!end subroutine finalize_eik_shot_img

!subroutine init_eik_shot(c,a,m,o,is,e,sh)
!type(coord2d),        intent(in) :: c
!type(aperture),       intent(in) :: a
!type(model2d),        intent(in) :: m
!type(other),          intent(in) :: o
!integer,              intent(in) :: is
!type(eik2d),          intent(out) :: e
!type(shot2d),         intent(out) :: sh
!call init_eik(c,a,m,is,e)
!call init_shot(c,a,m,o,is,sh)
!end subroutine init_eik_shot

!subroutine finalize_eik_shot(e,sh)
!type(eik2d),          intent(inout) :: e
!type(shot2d),         intent(inout) :: sh
!call finalize_eik(e)
!call finalize_shot(sh)
!end subroutine finalize_eik_shot

subroutine init_shot_img_kmm(c,a,m,k,s,o,is,sh,i,kmm)
type(coord2d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(km2d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot2d),         intent(out) :: sh
type(image2d),        intent(out) :: i
type(kmmod2d),        intent(out) :: kmm
call init_shot(c,a,m,o,is,sh)
call init_image(c,a,m,o,is,i)
call init_kmmod(c,a,m,k,s,is,kmm)
end subroutine init_shot_img_kmm

subroutine finalize_shot_img_kmm(sh,i,kmm)
type(shot2d),         intent(inout) :: sh
type(image2d),        intent(inout) :: i
type(kmmod2d),        intent(inout) :: kmm
call finalize_shot(sh)
call finalize_image(i)
call finalize_kmmod(kmm)
end subroutine finalize_shot_img_kmm

subroutine init_pw1_shot_kmm(c,a,m,k,s,o,is,sh,kmm)
type(pw1_coord2d),    intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(km2d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(pw1_shot2d),     intent(out) :: sh
type(kmmod2d),     intent(out) :: kmm
call init_pw1_shot(c,a,m,o,is,sh)
call init_pw1_kmmod(c,a,m,k,s,is,kmm)
end subroutine init_pw1_shot_kmm

subroutine finalize_pw1_shot_kmm(sh,kmm)
type(pw1_shot2d),         intent(inout) :: sh
type(kmmod2d),        intent(inout) :: kmm
call finalize_pw1_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_pw1_shot_kmm

subroutine init_pw1_shot_img_kmm(c,a,m,k,s,o,is,sh,i,kmm)
type(pw1_coord2d),    intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(km2d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(pw1_shot2d),         intent(out) :: sh
type(image2d),        intent(out) :: i
type(kmmod2d),     intent(out) :: kmm
call init_pw1_shot(c,a,m,o,is,sh)
call init_pw1_image(c,a,m,o,is,i)
call init_pw1_kmmod(c,a,m,k,s,is,kmm)
end subroutine init_pw1_shot_img_kmm

subroutine finalize_pw1_shot_img_kmm(sh,i,kmm)
type(pw1_shot2d),     intent(inout) :: sh
type(image2d),        intent(inout) :: i
type(kmmod2d),        intent(inout) :: kmm
call finalize_image(i)
call finalize_pw1_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_pw1_shot_img_kmm

subroutine init_pw2_shot_kmm(c,a,m,k,s,o,is,sh,kmm)
type(pw2_coord2d),    intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(km2d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(pw2_shot2d),     intent(out) :: sh
type(kmmod2d),     intent(out) :: kmm
call init_pw2_shot(c,a,m,o,is,sh)
call init_pw2_kmmod(a,m,k,s,kmm)
end subroutine init_pw2_shot_kmm

subroutine finalize_pw2_shot_kmm(sh,kmm)
type(pw2_shot2d),         intent(inout) :: sh
type(kmmod2d),        intent(inout) :: kmm
call finalize_pw2_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_pw2_shot_kmm

subroutine init_pw2_shot_img_kmm(c,a,m,k,s,o,is,sh,i,kmm)
type(pw2_coord2d),    intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(km2d),           intent(in) :: k
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(pw2_shot2d),     intent(out) :: sh
type(image2d),        intent(out) :: i
type(kmmod2d),        intent(out) :: kmm
call init_pw2_shot(c,a,m,o,is,sh)
call init_pw2_image(a,m,o,i)
call init_pw2_kmmod(a,m,k,s,kmm)
end subroutine init_pw2_shot_img_kmm

subroutine finalize_pw2_shot_img_kmm(sh,i,kmm)
type(pw2_shot2d),     intent(inout) :: sh
type(image2d),        intent(inout) :: i
type(kmmod2d),        intent(inout) :: kmm
call finalize_image(i)
call finalize_pw2_shot(sh)
call finalize_kmmod(kmm)
end subroutine finalize_pw2_shot_img_kmm

subroutine init_shot_gdmod_v1(c,a,m,o,f,g,s,is,sh,gdm)
type(coord2d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(fd2d),           intent(in) :: f
type(gdm2d),          intent(in) :: g
type(source),         intent(in) :: s
type(other),          intent(in) :: o
integer,              intent(in) :: is
type(shot2d),         intent(out) :: sh
type(gdmod2d),        intent(out) :: gdm
call init_shot(c,a,m,o,is,sh)
call init_gdmod_v1(c,a,m,f,g,s,is,gdm)
end subroutine init_shot_gdmod_v1

subroutine init_shot_gdmod_v2(c,a,m,o,g,s,is,sh,gdm)
type(coord2d),        intent(in) :: c
type(aperture),       intent(in) :: a
type(model2d),        intent(in) :: m
type(other),          intent(in) :: o
type(gdm2d),          intent(in) :: g
type(source),         intent(in) :: s
integer,              intent(in) :: is
type(shot2d),         intent(out) :: sh
type(gdmod2d),        intent(out) :: gdm
call init_shot(c,a,m,o,is,sh)
call init_gdmod_v2(c,a,m,g,s,is,gdm)
end subroutine init_shot_gdmod_v2

subroutine finalize_shot_gdmod(sh,gdm)
type(shot2d),         intent(inout) :: sh
type(gdmod2d),        intent(inout) :: gdm
call finalize_shot(sh)
call finalize_gdmod(gdm)
end subroutine finalize_shot_gdmod

!=====================================================================

subroutine coord_min_max(c,a,m,is)
type(coord2d),   intent(in)    :: c
type(aperture),  intent(in)    :: a
type(model2d),   intent(in)    :: m
integer, intent(in) :: is
! get the xmin and xmax for the receiver
xmax = maxval(c%xg(1:c%ng(is),is));  xmin = minval(c%xg(1:c%ng(is),is))
! get the xmin and xmax for source and receiver
xmax = max(xmax, maxval(c%xs(1:c%ng(is),is)))
xmin = min(xmin, minval(c%xs(1:c%ng(is),is)))
! get the zmin and zmax for the receiver
! get the xmin and xmax for source and receiver
zmin = a%z0;   zmax = a%z0+(m%nz-1)*m%dx
call nz_nx(nz,nx,a%aper_z,a%aper_x,a%z0,a%x0,m%dz,m%dx)
end subroutine coord_min_max

!=====================================================================

subroutine ms_coord_min_max(c,a,m,isg)
type(ms_coord2d),intent(in)    :: c
type(aperture),  intent(in)    :: a
type(model2d),   intent(in)    :: m
integer, intent(in) :: isg
! get the xmin and xmax for the receiver
xmax = maxval(c%xg(1:c%ng(isg),isg)); xmin = minval(c%xg(1:c%ng(isg),isg))
! get the xmin and xmax for source and receiver
xmax = max(xmax, maxval(c%xs(1:c%ns(isg),isg)))
xmin = min(xmin, minval(c%xs(1:c%ns(isg),isg)))
! get the zmin and zmax for the receiver
! get the zmin and zmax for source and receiver
zmin = a%z0;   zmax = a%z0+(m%nz-1)*m%dx
call nz_nx(nz,nx,a%aper_z,a%aper_x,a%z0,a%x0,m%dz,m%dx)
end subroutine ms_coord_min_max

!=====================================================================

subroutine pw1_coord_min_max(c,a,m,is)
type(pw1_coord2d),intent(in)    :: c
type(aperture),   intent(in)    :: a
type(model2d),    intent(in)    :: m
integer, intent(in) :: is
! get the xmin and xmax for the traces
xmax = maxval(c%xg(1:c%ng(is),is))
xmin = minval(c%xg(1:c%ng(is),is))
! get the zmin and zmax for the traces
!zmax = maxval(c%zg(1:c%ng(is),is))
!zmin = minval(c%zg(1:c%ng(is),is))
zmin = a%z0
zmax = a%z0+(m%nz-1)*m%dx
call nz_nx(nz,nx,a%aper_z,a%aper_x,a%z0,a%x0,m%dz,m%dx)
end subroutine pw1_coord_min_max

!=====================================================================

subroutine pw2_coord_min_max(a,m)
type(aperture),   intent(in)    :: a
type(model2d),    intent(in)    :: m
! get the xmin and xmax for the traces
xmin = a%x0;   xmax = xmin+(m%nx-1)*m%dx
! get the xmin and xmax for the traces
zmin = a%z0;   zmax = zmin+(m%nz-1)*m%dx
call nz_nx(nz,nx,a%aper_z,a%aper_x,a%z0,a%x0,m%dz,m%dx)
end subroutine pw2_coord_min_max

!=====================================================================

!subroutine init_xmin_xmax_zmin_zmax_pw(c,a,m,nx)
!type(pw_coord2d),  intent(in)    :: c
!type(aperture),    intent(in)    :: a
!type(model2d),     intent(in)    :: m
!integer,           intent(out):: nx
!! get the xmin and xmax for the receiver
!xmax = maxval(c%xg(:))
!xmin = minval(c%xg(:))
!xmax = xmax + a%aper - a%x0
!xmin = xmin - a%aper - a%x0
!ix1 = int(xmin/m%dx) + 1
!ix2 = nint(xmax/m%dx) + 1
!nx = ix2 - ix1 + 1
!end subroutine init_xmin_xmax_zmin_zmax_pw

!=====================================================================

subroutine eik_coord_min_max(a,m)
type(aperture),   intent(in)    :: a
type(model2d),    intent(in)    :: m
xmin = a%x0;   xmax = xmin + (m%nx-1.0)*m%dx
zmin = a%z0;   zmax = zmin+(m%nz-1)*m%dz
call nz_nx(nz,nx,a%aper_z,a%aper_x,a%z0,a%x0,m%dz,m%dx)
end subroutine eik_coord_min_max

!=====================================================================

subroutine init_eik(c,a,m,is,e)
type(eik_coord2d), intent(in) :: c
type(aperture),    intent(in) :: a
type(model2d),     intent(in) :: m
integer,           intent(in) :: is
type(eik2d),       intent(out):: e
call copytype(m,e%m)
if (a%isAper) then
   call eik_coord_min_max(a,m)
   e%m%nz=nz; e%m%nx=nx
   call allocate_and_initial(e%m%v,e%ttt,nz,nx)
   call cut(m%v,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,e%m%v)
else
   xmin=a%x0;   zmin=a%z0;   nz=m%nz;   nx=m%nx
   call allocate_and_initial(e%m%v,e%ttt,nz,nx)
   call copy_array(m%v,e%m%v)
endif
e%xs=c%xs(is)-xmin;   e%zs=c%zs(is)-zmin
end subroutine init_eik

subroutine finalize_eik(e)
type(eik2d),       intent(inout) :: e
call deallocate_and_free(e%m%v)
if (allocated(e%ttt)) call deallocate_and_free(e%ttt)
end subroutine finalize_eik

!=====================================================================

subroutine init_ref(c,a,m,i,fs,is,r)
type(coord2d),     intent(in) :: c
type(aperture),    intent(in) :: a
type(model2d),     intent(in) :: m
type(inv),         intent(in) :: i
integer,           intent(in) :: is,fs(:)
type(ref2d),       intent(out):: r
real                          :: ssmax,diff
integer                       :: i2
ssmax=1.0/i%vmin
call copy_variable(m%dx,r%dx,m%dz,r%dz)
!call copytype(m,r%m)
if (a%isAper) then
   !call init_xmin_xmax_zmin_zmax(c,a,m,is,e%m%nz,e%m%nx)
   !call allocate_and_initial(e%m%v,e%ttt,e%m%nz,e%m%nx)
   !call cut(m%v,m%nz,iz1,iz2,e%m%nz,m%nx,ix1,ix2,e%m%nx,e%m%v)
else
   xmin=a%x0-3*m%dx;   zmin=a%z0-3*m%dz
   r%nx=m%nx+6;        r%nz=m%nz+6
   !call allocate_and_initial(r%m%v,r%ttt,r%m%nz,r%m%nx)
   !call padarray(m%v,m%nz,3,3,m%nx,3,3,r%m%v,"coefficient",1.5)
   call allocate_and_initial(r%v,r%ttt,r%nz,r%nx)
   call padarray(m%v,m%nz,3,3,m%nx,3,3,r%v,"coefficient",1.5)
   call allocate_and_initial(r%fs,r%nx)
   call padarray(fs,m%nx,3,3,r%fs,"replicate")
   r%fs=r%fs+3
   do i2=1,r%nx
      !diff=(r%m%v(r%fs(i2),i2)-ssmax)/3.
      !r%m%v(1:r%fs(i2)-3,i2)=ssmax
      !r%m%v(r%fs(i2)-2,i2)=ssmax+diff
      !r%m%v(r%fs(i2)-1,i2)=ssmax+2.*diff
      diff=(r%v(r%fs(i2),i2)-ssmax)/3.
      r%v(1:r%fs(i2)-3,i2)=ssmax
      r%v(r%fs(i2)-2,i2)=ssmax+diff
      r%v(r%fs(i2)-1,i2)=ssmax+2.*diff
   enddo     
   call allocate_and_initial(r%xs,r%zs,r%xg,r%zg,c%ng(is))
   r%ng=c%ng(is)
   r%xs=c%xs(:,is)-xmin;   r%zs=c%zs(:,is)-zmin
   r%xg=c%xg(:,is)-xmin;   r%zg=c%zg(:,is)-zmin
   r%ix1=4;   r%iz1=4;   r%ix2=3+m%nx;   r%iz2=3+m%nz
endif
end subroutine init_ref

!=====================================================================
subroutine init_shot(c,a,m,o,is,s)
type(coord2d),      intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model2d),      intent(in)    :: m
type(other),        intent(in)    :: o
integer ,           intent(in)    :: is
type(shot2d),       intent(inout) :: s
call copy_variable(a%z0,zmin,o%dt_out,s%dt)
s%nt=o%nt_out
if (a%isAper) then
   call coord_min_max(c,a,m,is)
else
   xmin=a%x0;   zmin=a%z0
endif
s%ng=c%ng(is)
call allocate_and_initial(s%xs,s%zs,s%ng)
s%xs=c%xs(1:s%ng,is)-xmin;   s%zs=c%zs(1:s%ng,is)-zmin
call allocate_and_initial(s%xg,s%zg,s%ng)
s%xg=c%xg(1:s%ng,is)-xmin;   s%zg=c%zg(1:s%ng,is)-zmin
call allocate_and_initial(s%sid,s%gid,s%ng)
s%sid=c%sid(1:s%ng,is);   s%gid=c%gid(1:s%ng,is)
end subroutine init_shot

!=====================================================================

subroutine finalize_shot(s)
type(shot2d), intent(inout) :: s
deallocate(s%xs,s%zs,s%xg,s%zg,s%sid,s%gid)
if (allocated(s%seis)) deallocate(s%seis)
if (allocated(s%seis_ps)) deallocate(s%seis_ps)  ! LSMF
end subroutine finalize_shot

!=====================================================================

subroutine init_ms_shot(c,a,m,o,isg,s)
type(ms_coord2d),   intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model2d),      intent(in)    :: m
type(other),        intent(in)    :: o
integer ,           intent(in)    :: isg
type(ms_shot2d),    intent(inout) :: s
call copy_variable(a%z0,zmin,o%dt_out,s%dt)
s%nt=o%nt_out
if (a%isAper) then
   call ms_coord_min_max(c,a,m,isg)
else
   xmin=a%x0;   zmin=a%z0
endif
s%ns=c%ns(isg);   s%ng=c%ng(isg)
call allocate_and_initial(s%xs,s%zs,s%d,s%p,s%ns)
call allocate_and_initial(s%xg,s%zg,s%ng)
s%d=c%d(1:s%ns,isg);   s%p=c%p(1:s%ns,isg)
s%xs=c%xs(1:s%ns,isg)-xmin;   s%zs=c%zs(1:s%ns,isg)-zmin
s%xg=c%xg(1:s%ng,isg)-xmin;   s%zg=c%zg(1:s%ng,isg)-zmin
end subroutine init_ms_shot

!=====================================================================

subroutine finalize_ms_shot(s)
type(ms_shot2d), intent(inout) :: s
deallocate(s%xs,s%zs,s%d,s%p,s%xg,s%zg)
if (allocated(s%seis)) deallocate(s%seis)
end subroutine finalize_ms_shot

!=====================================================================

subroutine init_pw1_shot(c,a,m,o,is,s)
type(pw1_coord2d),  intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model2d),      intent(in)    :: m
type(other),        intent(in)    :: o
integer ,           intent(in)    :: is
type(pw1_shot2d),   intent(inout) :: s
call copy_variable(a%z0,zmin,o%dt_out,s%dt)
call copy_variable(c%zp1(is),s%zp1,c%theta_x1(is),s%theta_x1,&
                   c%xref_1(is),s%xref_1)
call copy_variable(o%nt_out,s%nt,c%proid(is),s%proid)
if (allocated(c%ng)) then
   call copy_variable(c%ng(is),s%ng)
endif
if (a%isAper) then
   call pw1_coord_min_max(c,a,m,is)
else
   xmin=a%x0;   zmin=a%z0
endif
if (allocated(c%ng)) then
   call allocate_and_initial(s%xg,s%zg,s%ng)
   s%xg=c%xg(1:s%ng,is)-xmin;   s%zg=c%zg(1:s%ng,is)-zmin
   call allocate_and_initial(s%gid,s%ng)
   s%gid=c%gid(1:s%ng,is)
endif
end subroutine init_pw1_shot

!=====================================================================

subroutine finalize_pw1_shot(s)
type(pw1_shot2d), intent(inout) :: s
if (allocated(s%xg)) deallocate(s%xg,s%zg,s%gid)
if (allocated(s%seis)) deallocate (s%seis)
end subroutine finalize_pw1_shot

!=====================================================================

subroutine init_pw2_shot(c,a,m,o,is,s)
type(pw2_coord2d),  intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model2d),      intent(in)    :: m
type(other),        intent(in)    :: o
integer ,           intent(in)    :: is
type(pw2_shot2d),   intent(inout) :: s
call copy_variable(a%z0,zmin,o%dt_out,s%dt)
call copy_variable(c%zp1(is),s%zp1,c%theta_x1(is),s%theta_x1,c%xref_1(is),s%xref_1)
call copy_variable(o%nt_out,s%nt,c%proid(is),s%proid)
if (allocated(c%ng)) then
   call copy_variable(c%ng(is),s%ng)
endif
if (a%isAper) then
   call pw2_coord_min_max(a,m)
else
   xmin=a%x0;   zmin=a%z0
endif
if (allocated(c%ng)) then
   call allocate_and_initial(s%gid,s%ng)
   call allocate_and_initial(s%theta_x2,s%xref_2,s%zp2,s%ng)
   s%gid=c%gid(1:s%ng,is)
   s%theta_x2=c%theta_x2(1:s%ng,is)
   s%xref_2=c%xref_2(1:s%ng,is)
   s%zp2=c%zp2(1:s%ng,is)
endif
end subroutine init_pw2_shot

!=====================================================================

subroutine finalize_pw2_shot(s)
type(pw2_shot2d), intent(inout) :: s
if (allocated(s%gid)) then
   deallocate(s%gid)
   deallocate(s%zp2,s%xref_2,s%theta_x2)
endif
if (allocated(s%seis)) deallocate (s%seis)
end subroutine finalize_pw2_shot

!=====================================================================

subroutine init_image(c,a,m,o,is,i)
type(coord2d),      intent(in)    :: c
type(aperture),     intent(in)    :: a
type(model2d),      intent(in)    :: m
type(other),        intent(in)    :: o
type(image2d),      intent(inout) :: i
integer , intent(in) :: is
call copy_variable(a%z0,zmin,o%illum_order,i%illum_order)
call copy_variable(o%isIllum,i%isIllum,o%isIllumReceiver,i%isIllumReceiver,&
                   o%isPreStackImg,i%isPreStackImg)
if (a%isAper) then
   call coord_min_max(c,a,m,is)
   i%nz=nz; i%nx=nx; i%ix1=ix1; i%ix2=ix2; i%iz1=iz1; i%iz2=iz2
else
   i%nx=m%nx; i%ix1=1; i%ix2=m%nx; i%nz=m%nz; i%iz1=1; i%iz2=m%nz
endif
call allocate_and_initial(i%img,i%nz,i%nx)
call allocate_and_initial(i%img_t,i%nz,i%nx)
call allocate_and_initial(i%img_m,i%nz,i%nx)
if (o%isIllum) call allocate_and_initial(i%illum,i%nz,i%nx)
end subroutine init_image

!=====================================================================

subroutine init_ms_image(mc,a,m,o,isg,i)
type(ms_coord2d),   intent(in)    :: mc
type(aperture),     intent(in)    :: a
type(model2d),      intent(in)    :: m
type(other),        intent(in)    :: o
type(image2d),      intent(inout) :: i
integer , intent(in) :: isg
call copy_variable(a%z0,zmin,o%illum_order,i%illum_order)
call copy_variable(o%isIllum,i%isIllum,o%isIllumReceiver,i%isIllumReceiver,&
                   o%isPreStackImg,i%isPreStackImg)
if (a%isAper) then
   call ms_coord_min_max(mc,a,m,isg)
   i%nz=nz; i%nx=nx; i%ix1=ix1; i%ix2=ix2; i%iz1=iz1; i%iz2=iz2
else
   i%nz=m%nz; i%nx=m%nx; i%ix1=1; i%ix2=m%nx; i%iz1=1; i%iz2=m%nz
endif
call allocate_and_initial(i%img,i%nz,i%nx)
if (o%isIllum) call allocate_and_initial(i%illum,i%nz,i%nx)
end subroutine init_ms_image

!=====================================================================

subroutine init_pw1_image(c,a,m,o,ip,i)
type(pw1_coord2d),  intent(in)    :: c
type(aperture),    intent(in)    :: a
type(model2d),     intent(in)    :: m
type(other),       intent(in)    :: o
integer,           intent(in)    :: ip
type(image2d),     intent(inout) :: i
call copy_variable(a%z0,zmin,o%illum_order,i%illum_order)
call copy_variable(o%isIllum,i%isIllum,o%isIllumReceiver,i%isIllumReceiver,&
                   o%isPreStackImg,i%isPreStackImg)
if (a%isAper) then
   call pw1_coord_min_max(c,a,m,ip)
   i%nz=nz; i%nx=nx; i%ix1=ix1; i%ix2=ix2; i%iz1=iz1; i%iz2=iz2
else
   i%nz=m%nz; i%nx=m%nx; i%ix1=1; i%ix2=m%nx; i%iz1=1; i%iz2=m%nz
endif
call allocate_and_initial(i%img,i%nz,i%nx)
if (o%isIllum) call allocate_and_initial(i%illum,i%nz,i%nx)
end subroutine init_pw1_image

!=====================================================================

subroutine init_pw2_image(a,m,o,i)
type(aperture),    intent(in)    :: a
type(model2d),     intent(in)    :: m
type(other),       intent(in)    :: o
type(image2d),     intent(inout) :: i
call copy_variable(a%z0,zmin,o%illum_order,i%illum_order)
call copy_variable(o%isIllum,i%isIllum,o%isIllumReceiver,i%isIllumReceiver,&
                   o%isPreStackImg,i%isPreStackImg)
if (a%isAper) then
   call pw2_coord_min_max(a,m)
   i%nz=nz; i%nx=nx; i%ix1=ix1; i%ix2=ix2; i%iz1=iz1; i%iz2=iz2
else
   i%nz=m%nz; i%nx=m%nx; i%ix1=1; i%ix2=m%nx; i%iz1=1; i%iz2=m%nz
endif
call allocate_and_initial(i%img,i%nz,i%nx)
if (o%isIllum) call allocate_and_initial(i%illum,i%nz,i%nx)
end subroutine init_pw2_image

!====================================================================

subroutine finalize_image(i)
type(image2d), intent(inout) :: i
if (allocated(i%img)) deallocate(i%img)
if (allocated(i%illum)) deallocate(i%illum)
end subroutine finalize_image

!=====================================================================

subroutine init_wf(c,a,m,f,is,wf)
type(coord2d),       intent(in)    :: c
type(aperture),      intent(in)    :: a
type(model2d),       intent(in)    :: m
type(fd2d),          intent(in)    :: f
integer,             intent(in)    :: is
type(wf2d),          intent(inout) :: wf
! setup normal parameters
wf%nt=f%nt
! setup parameters with aperture
if (a%isAper) then
   call coord_min_max(c,a,m,is)
else
   nx=m%nx;   nz=m%nz
endif
wf%nx_2=nx+2;   wf%nz_2=nz+2
call allocate_and_initial(wf%pwf,wf%nz_2,wf%nx_2,wf%nt)
end subroutine init_wf

!=====================================================================

subroutine finalize_wf(wf)
type(wf2d),          intent(inout) :: wf
if (allocated(wf%pwf)) deallocate(wf%pwf)
end subroutine finalize_wf

!=====================================================================

subroutine init_fdmod(c,a,m,f,s,is,fdm)
type(coord2d),   intent(in) :: c
type(aperture),  intent(in) :: a
type(model2d),   intent(in) :: m
type(fd2d),      intent(in) :: f
type(source),    intent(in) :: s
integer,         intent(in) :: is
type(fdmod2d),   intent(inout) :: fdm
call copytype(f,fdm%f);   call copytype(m,fdm%m); npml=f%npml

! setup source
call copytype(s,fdm%s)
call allocate_and_initial(fdm%s%sou,fdm%f%nt)
call wavelet_expand(s%sou,fdm%s%sou,fdm%f%nt)
! setup parameters with aperture
if (a%isAper) then
   call coord_min_max(c,a,m,is)
   fdm%m%nz=nz; fdm%m%nx=nx; nz_pml=nz + 2*npml;   nx_pml=nx + 2*npml
 ! if (fdm%f%isFS) then       ! modified by Bowen Guo
      call allocate_and_initial(fdm%f%fs,nx_pml)
      call cut_and_pad(f%fs,m%nx,ix1,ix2,npml,nx,fdm%f%fs)
 ! endif                      ! modified by Bowen Guo
   if (allocated(m%v)) then
      call allocate_and_initial(fdm%m%v,nz_pml,nx_pml)
      call cut_and_pad(m%v,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%v)
   endif
   if (allocated(m%vp)) then                                ! modified by Bowen Guo
      call allocate_and_initial(fdm%m%vp,nz_pml,nx_pml)     ! modified by Bowen Guo
      call cut_and_pad(m%vp,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%vp) ! modified by Bowen Guo
   endif                                                    ! modified by Bowen Guo
   if (allocated(m%lambda)) then                            ! modified by Bowen Guo
      call allocate_and_initial(fdm%m%lambda,nz_pml,nx_pml) ! modified by Bowen Guo
      call cut_and_pad(m%lambda,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%lambda) ! modified by Bowen Guo
   endif                                                    ! modified by Bowen Guo
   if (allocated(m%mu)) then                                ! modified by Bowen Guo
      call allocate_and_initial(fdm%m%mu,nz_pml,nx_pml)     ! modified by Bowen Guo
      call cut_and_pad(m%mu,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%mu) ! modified by Bowen Guo
   endif                                                    ! modified by Bowen Guo
   if (allocated(m%refl)) then
      call allocate_and_initial(fdm%m%refl,nz_pml,nx_pml)
      call cut_and_pad(m%refl,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%refl)
      fdm%m%refl(1:npml,:)=0.0
     ! fdm%m%refl(nz+npml+1:nz_pml,:)=0.0
     ! fdm%m%refl(:,1:npml)=0.0
     ! fdm%m%refl(:,nx+npml+1:nx_pml)=0.0
   
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(fdm%m%vs,nz_pml,nx_pml)
      call cut_and_pad(m%vs,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%vs)
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(fdm%m%den,nz_pml,nx_pml)
      call cut_and_pad(m%den,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%den)
   endif
else
   nx=m%nx;   nz=m%nz;   nx_pml=nx+2*npml;   nz_pml=nz+2*npml
  ! if (fdm%f%isFS) then      ! modified by Bowen Guo
      call allocate_and_initial(fdm%f%fs,nx_pml)
      call padarray(f%fs,nx,npml,npml,fdm%f%fs,"replicate")
  ! endif                     ! modified by Bowen Guo

  !  if (allocated(m%v)) write(*,*) "****  allocated ***"

    if (allocated(m%v)) then
      call allocate_and_initial(fdm%m%v,nz_pml,nx_pml)
     ! call padarray(m%v,m%nz,npml,npml,nx,npml,npml,fdm%m%v,"replicate")
   !   write(*,*) "m%nz = " ,m%nz
   !   write(*,*) "m%nx = " ,m%nx
   !   write(*,*) "nx = " ,nx
      call padarray(m%v,m%nz,npml,npml,m%nx,npml,npml,fdm%m%v,"replicate")
   endif
   if (allocated(m%vp)) then                              ! modified by Bowen Guo
      call allocate_and_initial(fdm%m%vp,nz_pml,nx_pml)   ! modified by Bowen Guo
      call padarray(m%vp,m%nz,npml,npml,m%nx,npml,npml,fdm%m%vp,"replicate") ! modified by Bowen Guo
   endif                                                  ! modified by Bowen Guo
   if (allocated(m%lambda)) then
      call allocate_and_initial(fdm%m%lambda,nz_pml,nx_pml) !modified by Bowen Guo
      call padarray(m%lambda,m%nz,npml,npml,nx,npml,npml,fdm%m%lambda,"replicate")! modified by Bowen Guo
   endif   ! modified by Bowen Guo
   if (allocated(m%mu)) then                               ! modified by Bowen Guo
      call allocate_and_initial(fdm%m%mu,nz_pml,nx_pml)    ! modified by Bowen Guo
      call padarray(m%mu,m%nz,npml,npml,nx,npml,npml,fdm%m%mu,"replicate") !modified by Bowen Guo
   endif                                                   ! modified by Bowen Guo
   
   if (allocated(m%refl)) then
      call allocate_and_initial(fdm%m%refl,nz_pml,nx_pml)
      call padarray(m%refl,nz,npml,npml,nx,npml,npml,fdm%m%refl,"replicate")
      fdm%m%refl(1:npml,:)=0.0
   !   fdm%m%refl(nz+npml+1:nz_pml,:)=0.0
   !   fdm%m%refl(:,1:npml)=0.0
   !   fdm%m%refl(:,nx+npml+1:nx_pml)=0.0
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(fdm%m%vs,nz_pml,nx_pml)
      call padarray(m%vs,nz,npml,npml,nx,npml,npml,fdm%m%vs,"replicate")
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(fdm%m%den,nz_pml,nx_pml)
      call padarray(m%den,nz,npml,npml,nx,npml,npml,fdm%m%den,"replicate")
   endif
endif
end subroutine init_fdmod

!=====================================================================

subroutine init_ms_fdmod(c,a,m,f,s,isg,fdm)
type(ms_coord2d),intent(in) :: c
type(aperture),  intent(in) :: a
type(model2d),   intent(in) :: m
type(fd2d),      intent(in) :: f
type(source),    intent(in) :: s
integer,         intent(in) :: isg
type(fdmod2d),   intent(inout) :: fdm
call copytype(f,fdm%f);   call copytype(m,fdm%m); npml=f%npml
! setup source
call copytype(s,fdm%s)
call allocate_and_initial(fdm%s%sou,fdm%f%nt)
call wavelet_expand(s%sou,fdm%s%sou,fdm%f%nt)
! setup parameters with aperture
if (a%isAper) then
   call ms_coord_min_max(c,a,m,isg)
   fdm%m%nz=nz;   fdm%m%nx=nx;   nx_pml=nx + 2*npml;   nz_pml=nz + 2*npml
   call allocate_and_initial(fdm%f%fs,nx_pml)
   call cut_and_pad(f%fs,m%nx,ix1,ix2,npml,nx,fdm%f%fs)
   if (allocated(m%v)) then
      call allocate_and_initial(fdm%m%v,nz_pml,nx_pml)
      call cut_and_pad(m%v,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(fdm%m%refl,nz_pml,nx_pml)
      call cut_and_pad(m%refl,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%refl)
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(fdm%m%vs,nz_pml,nx_pml)
      call cut_and_pad(m%vs,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%vs)
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(fdm%m%den,nz_pml,nx_pml)
      call cut_and_pad(m%den,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,fdm%m%den)
   endif
else
   nx=fdm%m%nx;   nz=fdm%m%nz;   nx_pml=nx+2*npml;   nz_pml=nz+2*npml
   call allocate_and_initial(fdm%f%fs,nx_pml)
   call padarray(f%fs,nx,npml,npml,fdm%f%fs,"replicate")
   if (allocated(m%v)) then
      call allocate_and_initial(fdm%m%v,nz_pml,nx_pml)
      call padarray(m%v,m%nz,npml,npml,nx,npml,npml,fdm%m%v,"replicate")
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(fdm%m%refl,nz_pml,nx_pml)
      call padarray(m%refl,nz,npml,npml,nx,npml,npml,fdm%m%refl,"replicate")
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(fdm%m%vs,nz_pml,nx_pml)
      call padarray(m%vs,nz,npml,npml,nx,npml,npml,fdm%m%vs,"replicate")
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(fdm%m%den,nz_pml,nx_pml)
      call padarray(m%den,nz,npml,npml,nx,npml,npml,fdm%m%den,"replicate")
   endif
endif
end subroutine init_ms_fdmod

!====================================================================

subroutine finalize_fdmod(fdm)
type(fdmod2d), intent(inout) :: fdm

if (allocated(fdm%m%v))    deallocate(fdm%m%v)   ! modified by Bowen Guo
if (allocated(fdm%f%fs))   deallocate(fdm%f%fs)  ! modified by Bowen Guo
if (allocated(fdm%s%sou))  deallocate(fdm%s%sou) ! modified by Bowen Guo
if (allocated(fdm%m%refl)) deallocate(fdm%m%refl)
if (allocated(fdm%m%lambda)) deallocate(fdm%m%lambda)   ! modified by Bowen Guo
if (allocated(fdm%m%mu)) deallocate(fdm%m%mu)           ! modified by Bowen Guo
if (allocated(fdm%m%vp)) deallocate(fdm%m%vp)           ! modified by Bowen Guo
if (allocated(fdm%m%vs)) deallocate(fdm%m%vs)
if (allocated(fdm%m%den)) deallocate(fdm%m%den)
end subroutine finalize_fdmod

!=====================================================================

subroutine init_gdmod_v1(c,a,m,f,g,s,is,gdm)
type(coord2d),   intent(in) :: c
type(aperture),  intent(in) :: a
type(model2d),   intent(in) :: m
type(fd2d),      intent(in) :: f
type(gdm2d),     intent(in) :: g
type(source),    intent(in) :: s
integer,         intent(in) :: is
type(gdmod2d),   intent(inout) :: gdm
call copytype(g,gdm%g);call copytype(m,gdm%m);call copytype(f,gdm%f);npml=f%npml
! setup source
call copytype(s,gdm%s)
call allocate_and_initial(gdm%s%sou,gdm%f%nt)
call wavelet_expand(s%sou,gdm%s%sou,gdm%f%nt)
! setup parameters with aperture
if (a%isAper) then
   call coord_min_max(c,a,m,is)
   gdm%m%nz=nz;   gdm%m%nx=nx;   nz_pml=nz + 2*npml;   nx_pml=nx + 2*npml
   call allocate_and_initial(gdm%f%fs,nx_pml)
   call cut_and_pad(f%fs,m%nx,ix1,ix2,npml,nx,gdm%f%fs)
   if (allocated(m%v)) then
      call allocate_and_initial(gdm%m%v,nz_pml,nx_pml)
      call cut_and_pad(m%v,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,gdm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(gdm%m%refl,nz_pml,nx_pml)
      call cut_and_pad(m%refl,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,gdm%m%refl)
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(gdm%m%vs,nz_pml,nx_pml)
      call cut_and_pad(m%vs,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,gdm%m%vs)
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(gdm%m%den,nz_pml,nx_pml)
      call cut_and_pad(m%den,m%nz,m%nx,iz1,iz2,ix1,ix2,npml,nz,nx,gdm%m%den)
   endif
else
   nz=m%nz;   nx=m%nx;   nx_pml=nx+2*npml;   nz_pml=nz+2*npml
   call allocate_and_initial(gdm%f%fs,nx_pml)
   call padarray(f%fs,nx,npml,npml,gdm%f%fs,"replicate")
   if (allocated(m%v)) then
      call allocate_and_initial(gdm%m%v,nz_pml,nx_pml)
      call padarray(m%v,m%nz,npml,npml,nx,npml,npml,gdm%m%v,"replicate")
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(gdm%m%refl,nz_pml,nx_pml)
      call padarray(m%refl,nz,npml,npml,nx,npml,npml,gdm%m%refl,"replicate")
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(gdm%m%vs,nz_pml,nx_pml)
      call padarray(m%vs,nz,npml,npml,nx,npml,npml,gdm%m%vs,"replicate")
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(gdm%m%den,nz_pml,nx_pml)
      call padarray(m%den,nz,npml,npml,nx,npml,npml,gdm%m%den,"replicate")
   endif
endif
end subroutine init_gdmod_v1

!=====================================================================

subroutine init_gdmod_v2(c,a,m,g,s,is,gdm)
type(coord2d),   intent(in) :: c
type(aperture),  intent(in) :: a
type(model2d),   intent(in) :: m
type(gdm2d),     intent(in) :: g
type(source),    intent(in) :: s
integer,         intent(in) :: is
type(gdmod2d),   intent(inout) :: gdm
call copytype(g,gdm%g);   call copytype(m,gdm%m)
! setup source
call copytype(s,gdm%s)
gdm%g%ntw_conv=gdm%g%nt+s%nw-1;   gdm%g%nt_conv=gdm%g%nt*2-1
call allocate_and_initial(gdm%g%wfft,gdm%g%ntw_conv)
call fft_wavelet(s%sou,gdm%g%wfft,s%nw,gdm%g%ntw_conv)
! setup parameters with aperture
if (a%isAper) then
   call coord_min_max(c,a,m,is)
   gdm%m%nz=nz;   gdm%m%nx=nx
   call allocate_and_initial(gdm%g%z0_mig,nx)
   call cut(g%z0_mig,m%nx,ix1,ix2,nx,gdm%g%z0_mig)
   gdm%g%z0_mig=gdm%g%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(gdm%m%v,nz,nx)
      call cut(m%v,m%nz,m%nx,iz1,iz2,ix1,ix2,nz,nx,gdm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(gdm%m%refl,nz,nx)
      call cut(m%refl,m%nz,m%nx,iz1,iz2,ix1,ix2,nz,nx,gdm%m%refl)
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(gdm%m%vs,nz,nx)
      call cut(m%vs,m%nz,m%nx,iz1,iz2,ix1,ix2,nz,nx,gdm%m%vs)
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(gdm%m%den,nz,nx)
      call cut(m%den,m%nz,m%nx,iz1,iz2,ix1,ix2,nz,nx,gdm%m%den)
   endif
else
   nz=gdm%m%nz;   nx=gdm%m%nx
   call allocate_and_initial(gdm%g%z0_mig,nx)
   gdm%g%z0_mig=g%z0_mig-a%z0
   if (allocated(m%v)) then
      call allocate_and_initial(gdm%m%v,nz,nx)
      call copy_array(m%v,gdm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(gdm%m%refl,nz,nx)
      call copy_array(m%refl,gdm%m%refl)
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(gdm%m%vs,nz,nx)
      call copy_array(m%vs,gdm%m%vs)
   endif
   if (allocated(m%den)) then
      call allocate_and_initial(gdm%m%den,nz,nx)
      call copy_array(m%den,gdm%m%den)
   endif
endif
end subroutine init_gdmod_v2

!=====================================================================

subroutine finalize_gdmod(gdm)
type(gdmod2d), intent(inout) :: gdm
if (allocated(gdm%f%fs)) deallocate(gdm%f%fs,gdm%s%sou,gdm%m%v)
if (allocated(gdm%m%refl)) deallocate(gdm%m%refl)
if (allocated(gdm%m%vs)) deallocate(gdm%m%vs)
if (allocated(gdm%m%den)) deallocate(gdm%m%den)
end subroutine finalize_gdmod

!=====================================================================

!subroutine init_water_fdmod2d(c,a,m,f,s,is,fdm)
!type(coord2d),   intent(in) :: c
!type(aperture),  intent(in) :: a
!type(water_model2d),   intent(in) :: m
!type(fd2d),      intent(in) :: f
!type(source),    intent(in) :: s
!integer,         intent(in) :: is
!type(water_fdmod2d),   intent(inout) :: fdm
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
!   call init_xmin_xmax_zmin_zmax(c,a,m%m,is,fdm%m%nx)
!   nx_pml=fdm%m%nx + 2*fdm%f%npml
!   call allocate_and_initial(fdm%f%fs,nx_pml)
!   call cut_and_pad(f%fs,m%nx,ix1,ix2,fdm%f%npml,fdm%m%nx,fdm%f%fs)
!   call allocate_and_initial(fdm%m%v,nz_pml,nx_pml)
!   call cut_and_pad(m%v,m%nx,fdm%m%nz,ix1,ix2,fdm%f%npml,fdm%m%nx,fdm%m%v)
!   call allocate_and_initial(fdm%m%v_w,nzw_pml,nx_pml)
!   call cut_and_pad(m%v_w,m%nx,fdm%m%nzw,ix1,ix2,fdm%f%npml,fdm%m%nx,fdm%m%v_w)
!   if (allocated(m%v_h)) then
!      call allocate_and_initial(fdm%m%v_h,nzw_pml,nx_pml)
!      call cut_and_pad(m%v_h,m%nx,fdm%m%nzw,ix1,ix2,f%npml,fdm%m%nx,fdm%m%v_h)
!   endif
!   if (allocated(m%refl)) then
!      call allocate_and_initial(fdm%m%refl,nz_pml,nx_pml)
!      call cut_and_pad(m%refl,m%nx,fdm%m%nz,ix1,ix2,f%npml,fdm%m%nx,fdm%m%refl)
!   endif
!   if (allocated(m%refl_fs)) then
!      call allocate_and_initial(fdm%m%refl_fs,nzw_pml,nx_pml)
!      call cut_and_pad(m%refl_fs,m%nx,fdm%m%nzw,ix1,ix2,f%npml,fdm%m%nx,fdm%m%refl_fs)
!   endif
!   if (allocated(m%refl_ob)) then
!      call allocate_and_initial(fdm%m%refl_ob,nzw_pml,nx_pml)
!      call cut_and_pad(m%refl_ob,m%nx,fdm%m%nzw,ix1,ix2,f%npml,fdm%m%nx,fdm%m%refl_ob)
!   endif
!else
!   nx_pml=fdm%m%nx+2*fdm%f%npml
!   call allocate_and_initial(fdm%f%fs,nx_pml)
!   call Padding(f%fs,fdm%m%nx,f%npml,f%npml,fdm%f%fs)
!   call allocate_and_initial(fdm%m%v,nz_pml,nx_pml)
!   call Padding(m%v,m%nz,fdm%f%npml,fdm%f%npml,fdm%m%nx,fdm%f%npml,fdm%f%npml,fdm%m%v)
!   call allocate_and_initial(fdm%m%v_w,nzw_pml,nx_pml)
!   call Padding(m%v_w,m%nzw,fdm%f%npml,fdm%f%npml,fdm%m%nx,fdm%f%npml,&
!                fdm%f%npml,fdm%m%v_w)
!   if (allocated(m%refl)) then
!      call allocate_and_initial(fdm%m%refl,nz_pml,nx_pml)
!      call Padding(m%refl,fdm%m%nz,fdm%f%npml,fdm%f%npml,fdm%m%nx,fdm%f%npml,&
!                   fdm%f%npml,fdm%m%refl)
!   endif
!   if (allocated(m%v_h)) then
!      call allocate_and_initial(fdm%m%v_h,nzw_pml,nx_pml)
!      call Padding(m%v_h,fdm%m%nzw,fdm%f%npml,fdm%f%npml,fdm%m%nx,fdm%f%npml,&
!                   fdm%f%npml,fdm%m%v_h)
!   endif
!   if (allocated(m%refl_fs)) then
!      call allocate_and_initial(fdm%m%refl_fs,nzw_pml,nx_pml)
!      call Padding(m%refl_fs,fdm%m%nzw,fdm%f%npml,fdm%f%npml,fdm%m%nx,fdm%f%npml,&
!                   fdm%f%npml,fdm%m%refl_fs)
!   endif
!   if (allocated(m%refl_ob)) then
!      call allocate_and_initial(fdm%m%refl_ob,nzw_pml,nx_pml)
!      call Padding(m%refl_ob,fdm%m%nzw,fdm%f%npml,fdm%f%npml,fdm%m%nx,fdm%f%npml,&
!                   fdm%f%npml,fdm%m%refl_ob)
!   endif
!endif
!end subroutine init_water_fdmod2d

!====================================================================

!subroutine finalize_water_fdmod2d(fdm)
!type(water_fdmod2d), intent(inout) :: fdm
!deallocate(fdm%f%fs,fdm%s%sou,fdm%m%v,fdm%m%v_w)
!if (allocated(fdm%m%refl)) deallocate(fdm%m%v_h)
!if (allocated(fdm%m%refl)) deallocate(fdm%m%refl)
!if (allocated(fdm%m%refl)) deallocate(fdm%m%refl_fs)
!if (allocated(fdm%m%refl)) deallocate(fdm%m%refl_ob)
!end subroutine finalize_water_fdmod2d

!=====================================================================

subroutine init_kmmod(c,a,m,k,s,is,kmm)
type(coord2d),   intent(in) :: c
type(aperture),  intent(in) :: a
type(model2d),   intent(in) :: m
type(km2d),      intent(in) :: k
type(source),    intent(in) :: s
integer,         intent(in) :: is
type(kmmod2d),   intent(inout) :: kmm
call copytype(m,kmm%m);   call copytype(k,kmm%k)
! setup source
call copytype(s,kmm%s)
kmm%k%nt_conv=kmm%k%nt+s%nw-1
call allocate_and_initial(kmm%k%wfft,kmm%k%nt_conv)
call fft_wavelet(s%sou,kmm%k%wfft,s%nw,kmm%k%nt_conv)

! setup parameters with aperture
if (a%isAper) then
   call coord_min_max(c,a,m,is)
   kmm%m%nz=nz; kmm%m%nx=nx 
   call allocate_and_initial(kmm%k%z0_mig,nx)
   call cut(k%z0_mig,m%nx,ix1,ix2,nx,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,nx)
      call cut(m%v,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,kmm%m%v)
   endif
   if (allocated(m%vs)) then
      call allocate_and_initial(kmm%m%vs,nz,nx)
      call cut(m%vs,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,kmm%m%vs)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,nx)
      call cut(m%refl,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,kmm%m%refl)
   endif
   if (allocated(m%refl_ps)) then
      call allocate_and_initial(kmm%m%refl_ps,nz,nx)
      call cut(m%refl_ps,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,kmm%m%refl_ps)
   endif
else
   nx=m%nx;   nz=m%nz
   call allocate_and_initial(kmm%k%z0_mig,nx)
   call copy_array(k%z0_mig,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%vs)) then
      call allocate_and_initial(kmm%m%vs,nz,nx)
      call copy_array(m%vs,kmm%m%vs)
   endif
   
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,nx)
      call copy_array(m%v,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,nx)
      call copy_array(m%refl,kmm%m%refl)
   endif
   if (allocated(m%refl_ps)) then
      call allocate_and_initial(kmm%m%refl_ps,nz,nx)
      call copy_array(m%refl_ps,kmm%m%refl_ps)
   endif
endif
end subroutine init_kmmod

!=====================================================================

subroutine init_pw1_kmmod(c,a,m,k,s,is,kmm)
type(pw1_coord2d),intent(in) :: c
type(aperture),   intent(in) :: a
type(model2d),    intent(in) :: m
type(km2d),       intent(in) :: k
type(source),     intent(in) :: s
integer,          intent(in) :: is
type(kmmod2d),    intent(inout) :: kmm
call copytype(m,kmm%m);   call copytype(k,kmm%k)
! setup source
call copytype(s,kmm%s)
kmm%k%nt_conv=kmm%k%nt+s%nw-1
call allocate_and_initial(kmm%k%wfft,kmm%k%nt_conv)
call fft_wavelet(s%sou,kmm%k%wfft,s%nw,kmm%k%nt_conv)
! setup parameters with aperture
if (a%isAper) then
   call pw1_coord_min_max(c,a,m,is)
   kmm%m%nz=nz;   kmm%m%nx=nx
   call allocate_and_initial(kmm%k%z0_mig,nx)
   call cut(k%z0_mig,m%nx,ix1,ix2,nx,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,nx)
      call cut(m%v,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,nx)
      call cut(m%refl,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,kmm%m%refl)
   endif
else
   nz=m%nz;   nx=m%nx
   call allocate_and_initial(kmm%k%z0_mig,nx)
   call copy_array(k%z0_mig,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,nx)
      call copy_array(m%v,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,nx)
      call copy_array(m%refl,kmm%m%refl)
   endif
endif
end subroutine init_pw1_kmmod

!=====================================================================

subroutine init_pw2_kmmod(a,m,k,s,kmm)
type(aperture),   intent(in) :: a
type(model2d),    intent(in) :: m
type(km2d),       intent(in) :: k
type(source),     intent(in) :: s
type(kmmod2d),    intent(inout) :: kmm
call copytype(m,kmm%m);   call copytype(k,kmm%k)
! setup source
call copytype(s,kmm%s)
kmm%k%nt_conv=kmm%k%nt+s%nw-1
call allocate_and_initial(kmm%k%wfft,kmm%k%nt_conv)
call fft_wavelet(s%sou,kmm%k%wfft,s%nw,kmm%k%nt_conv)
! setup parameters with aperture
if (a%isAper) then
   call pw2_coord_min_max(a,m)
   kmm%m%nz=nz;   kmm%m%nx=nx
   call allocate_and_initial(kmm%k%z0_mig,nx)
   call cut(k%z0_mig,m%nx,ix1,ix2,nx,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,nx)
      call cut(m%v,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,nx)
      call cut(m%refl,m%nz,iz1,iz2,nz,m%nx,ix1,ix2,nx,kmm%m%refl)
   endif
else
   nx=m%nz;   nx=m%nx
   call allocate_and_initial(kmm%k%z0_mig,nx)
   call copy_array(k%z0_mig,kmm%k%z0_mig)
   kmm%k%z0_mig=kmm%k%z0_mig-zmin
   if (allocated(m%v)) then
      call allocate_and_initial(kmm%m%v,nz,nx)
      call copy_array(m%v,kmm%m%v)
   endif
   if (allocated(m%refl)) then
      call allocate_and_initial(kmm%m%refl,nz,nx)
      call copy_array(m%refl,kmm%m%refl)
   endif
endif
end subroutine init_pw2_kmmod

!====================================================================

subroutine finalize_kmmod(kmm)
type(kmmod2d), intent(inout) :: kmm
deallocate(kmm%k%wfft,kmm%k%z0_mig)
if (allocated(kmm%m%v)) deallocate(kmm%m%v)
if (allocated(kmm%m%vs)) deallocate(kmm%m%vs)
if (allocated(kmm%m%refl)) deallocate(kmm%m%refl)
if (allocated(kmm%m%refl_ps)) deallocate(kmm%m%refl_ps)
end subroutine finalize_kmmod

!=====================================================================

!subroutine init_pwkmmod2d(c,a,m,k,s,kmm)
!type(pw_coord2d),intent(in) :: c
!type(aperture),  intent(in) :: a
!type(model2d),   intent(in) :: m
!type(km2d),      intent(in) :: k
!type(source),    intent(in) :: s
!type(pwkmmod2d), intent(inout) :: kmm
!call copytype(m,kmm%m)
!call copytype(k,kmm%k)
! setup source
!call copytype(s,kmm%s)
!call allocate_and_initial(kmm%k%wfft,kmm%k%nt_conv)
!call fft_wavelet(s%sou,kmm%k%wfft,s%nw,kmm%k%nt_conv)
!!call allocate_and_initial(kmm%taper,kmm%ns_taper)
!!call hanning(kmm%taper,kmm%ns_taper)
!! setup parameters with aperture
!if (a%isAper) then
!   call init_xmin_xmax_zmin_zmax(c,a,m,kmm%m%nx)
!   call allocate_and_initial(kmm%k%z0_mig,kmm%m%nx)
!   call cut(k%z0_mig,m%nx,ix1,ix2,kmm%m%nx,kmm%k%z0_mig)
!   kmm%k%z0_mig=kmm%k%z0_mig-zmin
!   call allocate_and_initial(kmm%m%v,kmm%m%nz,kmm%m%nx)
!   call cut(m%v,m%nz,1,m%nz,kmm%m%nz,m%nx,ix1,ix2,kmm%m%nx,kmm%m%v)
!   if (allocated(m%refl)) then
!      call allocate_and_initial(kmm%m%refl,kmm%m%nz,kmm%m%nx)
!      call cut(m%refl,m%nz,1,m%nz,kmm%m%nz,m%nx,ix1,ix2,kmm%m%nx,kmm%m%refl)
!   endif
!else
!   call allocate_and_initial(kmm%k%z0_mig,kmm%m%nx)
!   call copy_array(k%z0_mig,kmm%k%z0_mig)
!   kmm%k%z0_mig=kmm%k%z0_mig-zmin
!   call allocate_and_initial(kmm%m%v,kmm%m%nz,kmm%m%nx)
!   call copy_array(m%v,kmm%m%v)
!   if (allocated(m%refl)) then
!      call allocate_and_initial(kmm%m%refl,kmm%m%nz,kmm%m%nx)
!      call copy_array(m%refl,kmm%m%refl)
!   endif
!endif
!end subroutine init_pwkmmod2d

!!====================================================================
!
!subroutine finalize_pwkmmod2d(kmm)
!type(pwkmmod2d), intent(inout) :: kmm
!deallocate(kmm%k%wfft,kmm%m%v)
!if (allocated(kmm%m%refl)) deallocate(kmm%m%refl)
!!if (allocated(kmm%taper)) deallocate(kmm%taper)
!end subroutine finalize_pwkmmod2d

!====================================================================

subroutine aperture_image_l2g_real(isAper,img,iz1,iz2,ix1,ix2,img_out,nz,nx)
logical,optional,   intent(in) :: isAper
real,               intent(in) :: img(:,:)
integer,            intent(in) :: iz1,iz2,ix1,ix2,nz, nx
real,               intent(out) :: img_out(:,:)
integer ::ix1_in,ix2_in,iz1_in,iz2_in,ix1_out,ix2_out,iz1_out,iz2_out
if (present(isAper).and.isAper) then
   img_out=0.0
   call l2g_i1i2(ix1,ix2,nx,ix1_in,ix2_in,ix1_out,ix2_out)
   call l2g_i1i2(iz1,iz2,nz,iz1_in,iz2_in,iz1_out,iz2_out)
   img_out(iz1_out:iz2_out,ix1_out:ix2_out)=img(iz1_in:iz2_in,ix1_in:ix2_in)
else
   img_out=img
endif
end subroutine aperture_image_l2g_real

subroutine aperture_image_l2g_int(isAper,img,iz1,iz2,ix1,ix2,img_out,nz,nx)
logical,optional,   intent(in) :: isAper
integer,            intent(in) :: img(:,:)
integer,            intent(in) :: iz1,iz2,ix1,ix2,nz, nx
integer,            intent(out) :: img_out(:,:)
integer ::ix1_in,ix2_in,iz1_in,iz2_in,ix1_out,ix2_out,iz1_out,iz2_out
if (present(isAper).and.isAper) then
   img_out=0.0
   call l2g_i1i2(ix1,ix2,nx,ix1_in,ix2_in,ix1_out,ix2_out)
   call l2g_i1i2(iz1,iz2,nz,iz1_in,iz2_in,iz1_out,iz2_out)
   img_out(iz1_out:iz2_out,ix1_out:ix2_out)=img(iz1_in:iz2_in,ix1_in:ix2_in)
else
   img_out=img
endif
end subroutine aperture_image_l2g_int

!====================================================================

subroutine l2g_i1i2(i1,i2,n,i1_in,i2_in,i1_out,i2_out)
integer, intent(in) :: i1,i2,n
integer, intent(out):: i1_in,i2_in,i1_out,i2_out
if (i1<=0) then
   i1_in=-i1+2;   i1_out=1
   if (i2<=n) then
      i2_in=-i1+1+i2;   i2_out=i2
   else
      i2_in=-i1+1+n;   i2_out=n
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

subroutine nz_nx(nz,nx,aper_z,aper_x,z0,x0,dz,dx)
integer, intent(out) :: nz,nx
real,    intent(in)  :: aper_z,aper_x,z0,x0,dz,dx
xmax = xmax + aper_x - x0; xmin = xmin - aper_x - x0
ix1 = int(xmin/dx) + 1; ix2 = nint(xmax/dx) + 1
zmax = zmax + aper_z - z0;   zmin = zmin - aper_z - z0
iz1 = int(zmin/dz) + 1;   iz2 = nint(zmax/dz) + 1
nz = iz2 - iz1 + 1; nx = ix2 - ix1 + 1
end subroutine nz_nx

end module module_aperture2d
