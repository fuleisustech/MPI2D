module module_boundary3d
use module_global, only : PI
use module_datatype
use module_array
use module_utility
implicit none

type bc3d_general
   real, allocatable :: p_top(:,:,:,:),p_bottom(:,:,:,:),&
     p_left(:,:,:,:),p_right(:,:,:,:),p_front(:,:,:,:),&
     p_back(:,:,:,:),p_nt(:,:,:),p_nt_1(:,:,:)
end type bc3d_general

type bc3d_sm
   real, allocatable :: sw_top(:,:,:,:),sw_bottom(:,:,:,:),&
     sw_left(:,:,:,:),sw_right(:,:,:,:),sw_front(:,:,:,:),&
     sw_back(:,:,:,:),sw_nt(:,:,:),sw_nt_1(:,:,:),&
     sh_top(:,:,:,:),sh_bottom(:,:,:,:),sh_left(:,:,:,:),&
     sh_right(:,:,:,:),sh_front(:,:,:,:),&
     sh_back(:,:,:,:),sh_nt(:,:,:),sh_nt_1(:,:,:),&
     p_top(:,:,:,:),p_bottom(:,:,:,:),p_left(:,:,:,:),&
     p_right(:,:,:,:),p_front(:,:,:,:),p_back(:,:,:,:),&
     p_nt(:,:,:),p_nt_1(:,:,:)
end type bc3d_sm

type bc3d_sm1
   real, allocatable :: sw_top(:,:,:,:),sw_bottom(:,:,:,:),&
     sw_left(:,:,:,:),sw_right(:,:,:,:),sw_front(:,:,:,:),&
     sw_back(:,:,:,:),sw_nt(:,:,:),sw_nt_1(:,:,:),&
     sw2_top(:,:,:,:),sw2_bottom(:,:,:,:),sw2_left(:,:,:,:),&
     sw2_right(:,:,:,:),sw2_front(:,:,:,:),sw2_back(:,:,:,:),&
     sw2_nt(:,:,:),sw2_nt_1(:,:,:),&
     p_top(:,:,:,:),p_bottom(:,:,:,:),p_left(:,:,:,:),&
     p_right(:,:,:,:),p_front(:,:,:,:),p_back(:,:,:,:),&
     p_nt(:,:,:),p_nt_1(:,:,:)
end type bc3d_sm1

type bc3d_sm2
   real, allocatable :: sw_top(:,:,:,:),sw_bottom(:,:,:,:),&
     sw_left(:,:,:,:),sw_right(:,:,:,:),sw_front(:,:,:,:),&
     sw_back(:,:,:,:),sw_nt(:,:,:),sw_nt_1(:,:,:),&
     sw2_top(:,:,:,:),sw2_bottom(:,:,:,:),sw2_left(:,:,:,:),&
     sw2_right(:,:,:,:),sw2_front(:,:,:,:),sw2_back(:,:,:,:),&
     sw2_nt(:,:,:),sw2_nt_1(:,:,:),&
     sw3_top(:,:,:,:),sw3_bottom(:,:,:,:),sw3_left(:,:,:,:),&
     sw3_right(:,:,:,:),sw3_front(:,:,:,:),sw3_back(:,:,:,:),&
     sw3_nt(:,:,:),sw3_nt_1(:,:,:),&
     sw4_top(:,:,:,:),sw4_bottom(:,:,:,:),sw4_left(:,:,:,:),&
     sw4_right(:,:,:,:),sw4_front(:,:,:,:),sw4_back(:,:,:,:),&
     sw4_nt(:,:,:),sw4_nt_1(:,:,:),&
     p_top(:,:,:,:),p_bottom(:,:,:,:),p_left(:,:,:,:),&
     p_right(:,:,:,:),p_front(:,:,:,:),p_back(:,:,:,:),&
     p_nt(:,:,:),p_nt_1(:,:,:)
end type bc3d_sm2

integer, private :: nz_pml, nzw_pml, ny_pml, nx_pml, bc_len, iz, iy, ix, bc_len_1
real,    private :: a, xa, kappa

! 3D coordinate define        _________
!           ^ y              /|       /|
!          /                / | top  / |
!         /                /__|_____/  |
!        /                 |  |back |  | left
!       0---------> x      |  |_____|__|
!       |               rig|ht/     |  /
!       |                  | /bottom| /
!       |                  |/_______|/
!       V z                   front

interface bc_init
   module procedure bc_init_general
   module procedure bc_init_sm
   module procedure bc_init_sm1
   module procedure bc_init_sm2
end interface bc_init

interface bc_final
   module procedure bc_final_general
   module procedure bc_final_sm
   module procedure bc_final_sm1
   module procedure bc_final_sm2
end interface bc_final

interface local_variables
   module procedure local_general
   module procedure local_water
end interface local_variables

private local_variables

contains

!=====================================================================

subroutine local_general(nz,ny,nx,npml)
integer, intent(in) :: nz,nx,ny,npml
nz_pml = nz + 2 * npml
ny_pml = ny + 2 * npml
nx_pml = nx + 2 * npml
end subroutine local_general

!=====================================================================

subroutine local_water(nz,nzw,ny,nx,npml)
integer, intent(in) :: nz,nzw,ny,nx,npml
nz_pml  = nz  + 2 * npml
nzw_pml = nzw + 2 * npml
ny_pml  = ny + 2 * npml
nx_pml  = nx  + 2 * npml
end subroutine local_water

!=====================================================================

subroutine bc_init_general(bc,nz,ny,nx,npml,nt,fd_order,mem)
type(bc3d_general), intent(inout) :: bc
integer,  intent(in)  :: nz,ny,nx,npml,nt,fd_order
integer(8),optional,intent(inout) :: mem
bc_len=(fd_order-20)/2+1
call local_variables(nz,ny,nx,npml)
call allocate_and_initial(bc%p_top,bc%p_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%p_left,bc%p_right,nz,ny,bc_len,nt,mem)
call allocate_and_initial(bc%p_front,bc%p_back,nz,bc_len,nx,nt,mem)
call allocate_and_initial(bc%p_nt,bc%p_nt_1,nz_pml,ny_pml,nx_pml,mem)
end subroutine bc_init_general

!=====================================================================

subroutine bc_final_general(bc,mem)
type(bc3d_general) ,intent(inout) :: bc
integer(8),optional,intent(inout) :: mem
call deallocate_and_free(bc%p_top,bc%p_bottom,bc%p_left,bc%p_right,&
           bc%p_front,bc%p_back,mem)
call deallocate_and_free(bc%p_nt,bc%p_nt_1,mem)
end subroutine bc_final_general

!=====================================================================

subroutine bc_init_sm(bc,nz,nzw,ny,nx,npml,nt,fd_order,mem)
type(bc3d_sm),    intent(inout) :: bc
integer,        intent(in)  :: nz,nzw,ny,nx,npml,nt,fd_order
integer(8),optional,intent(inout) :: mem
bc_len=(fd_order-20)/2+1
call local_variables(nz,nzw,ny,nx,npml)
call allocate_and_initial(bc%p_top,bc%p_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%p_left,bc%p_right,nz,ny,bc_len,nt,mem)
call allocate_and_initial(bc%p_front,bc%p_back,nz,bc_len,nx,nt,mem)
call allocate_and_initial(bc%p_nt,bc%p_nt_1,nz_pml,ny_pml,nx_pml,mem)
call allocate_and_initial(bc%sw_top,bc%sw_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%sw_left,bc%sw_right,nzw,ny,bc_len,nt,mem)
call allocate_and_initial(bc%sw_front,bc%sw_back,nzw,bc_len,nx,nt,mem)
call allocate_and_initial(bc%sw_nt,bc%sw_nt_1,nzw_pml,ny_pml,nx_pml,mem)
call allocate_and_initial(bc%sh_top,bc%sh_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%sh_left,bc%sh_right,nzw,ny,bc_len,nt,mem)
call allocate_and_initial(bc%sh_front,bc%sh_back,nzw,bc_len,nx,nt,mem)
call allocate_and_initial(bc%sh_nt,bc%sh_nt_1,nzw_pml,ny_pml,nx_pml,mem)
end subroutine bc_init_sm

!=====================================================================

subroutine bc_final_sm(bc,mem)
type(bc3d_sm), intent(inout) :: bc
integer(8),optional,intent(inout) :: mem
call deallocate_and_free(bc%p_top,bc%p_bottom,bc%p_left,bc%p_right,&
   bc%p_front,bc%p_back,mem)
call deallocate_and_free(bc%sw_top,bc%sw_bottom,bc%sw_left,bc%sw_right,&
   bc%sw_front,bc%sw_back,mem)
call deallocate_and_free(bc%sh_top,bc%sh_bottom,bc%sh_left,bc%sh_right,&
   bc%sh_front,bc%sh_back,mem)
call deallocate_and_free(bc%p_nt,bc%p_nt_1,bc%sw_nt,bc%sw_nt_1,&
  bc%sh_nt,bc%sh_nt_1,mem)
end subroutine bc_final_sm

!=====================================================================

subroutine bc_init_sm1(bc,nz,nzw,ny,nx,npml,nt,fd_order,mem)
type(bc3d_sm1),  intent(inout) :: bc
integer,         intent(in)    :: nz,nzw,ny,nx,npml,nt,fd_order
integer(8),optional,intent(inout) :: mem
bc_len=(fd_order-20)/2+1
call local_variables(nz,nzw,ny,nx,npml)
call allocate_and_initial(bc%p_top,bc%p_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%p_left,bc%p_right,nz,ny,bc_len,nt,mem)
call allocate_and_initial(bc%p_front,bc%p_back,nz,bc_len,nx,nt,mem)
call allocate_and_initial(bc%p_nt,bc%p_nt_1,nz_pml,ny_pml,nx_pml,mem)
call allocate_and_initial(bc%sw_top,bc%sw_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%sw_left,bc%sw_right,nzw,ny,bc_len,nt,mem)
call allocate_and_initial(bc%sw_front,bc%sw_back,nzw,bc_len,nx,nt,mem)
call allocate_and_initial(bc%sw_nt,bc%sw_nt_1,nzw_pml,ny_pml,nx_pml,mem)
call allocate_and_initial(bc%sw2_top,bc%sw2_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%sw2_left,bc%sw2_right,nzw,ny,bc_len,nt,mem)
call allocate_and_initial(bc%sw2_front,bc%sw2_back,nzw,bc_len,nx,nt,mem)
call allocate_and_initial(bc%sw2_nt,bc%sw2_nt_1,nzw_pml,ny_pml,nx_pml,mem)
end subroutine bc_init_sm1

!=====================================================================

subroutine bc_final_sm1(bc,mem)
type(bc3d_sm1), intent(inout) :: bc
integer(8),optional,intent(inout) :: mem
call deallocate_and_free(bc%p_top,bc%p_bottom,bc%p_left,bc%p_right,&
   bc%p_front,bc%p_back,mem)
call deallocate_and_free(bc%sw_top,bc%sw_bottom,bc%sw_left,bc%sw_right,&
   bc%sw_front,bc%sw_back,mem)
call deallocate_and_free(bc%sw2_top,bc%sw2_bottom,bc%sw2_left,&
   bc%sw2_right,bc%sw2_front,bc%sw2_back,mem)
call deallocate_and_free(bc%p_nt,bc%p_nt_1,bc%sw_nt,bc%sw_nt_1,&
  bc%sw2_nt,bc%sw2_nt_1)
end subroutine bc_final_sm1

!=====================================================================

subroutine bc_init_sm2(bc,nz,nzw,ny,nx,npml,nt,fd_order,mem)
type(bc3d_sm2), intent(inout) :: bc
integer,      intent(in)    :: nz,nzw,ny,nx,npml,nt,fd_order
integer(8),optional,intent(inout) :: mem
bc_len=(fd_order-20)/2+1
call local_variables(nz,nzw,nx,npml)
call allocate_and_initial(bc%p_top,bc%p_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%p_left,bc%p_right,nz,ny,bc_len,nt,mem)
call allocate_and_initial(bc%p_front,bc%p_back,nz,bc_len,nx,nt,mem)
call allocate_and_initial(bc%p_nt,bc%p_nt_1,nz_pml,ny_pml,nx_pml,mem)
call allocate_and_initial(bc%sw_top,bc%sw_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%sw_left,bc%sw_right,nzw,ny,bc_len,nt,mem)
call allocate_and_initial(bc%sw_front,bc%sw_back,nzw,bc_len,nx,nt,mem)
call allocate_and_initial(bc%sw_nt,bc%sw_nt_1,nzw_pml,ny_pml,nx_pml,mem)
call allocate_and_initial(bc%sw2_top,bc%sw2_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%sw2_left,bc%sw2_right,nzw,ny,bc_len,nt,mem)
call allocate_and_initial(bc%sw2_front,bc%sw2_back,nzw,bc_len,nx,nt,mem)
call allocate_and_initial(bc%sw2_nt,bc%sw2_nt_1,nzw_pml,ny_pml,nx_pml,mem)
call allocate_and_initial(bc%sw3_top,bc%sw3_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%sw3_left,bc%sw3_right,nzw,ny,bc_len,nt,mem)
call allocate_and_initial(bc%sw3_front,bc%sw3_back,nzw,bc_len,nx,nt,mem)
call allocate_and_initial(bc%sw3_nt,bc%sw3_nt_1,nzw_pml,ny_pml,nx_pml,mem)
call allocate_and_initial(bc%sw4_top,bc%sw4_bottom,bc_len,ny,nx,nt,mem)
call allocate_and_initial(bc%sw4_left,bc%sw4_right,nzw,ny,bc_len,nt,mem)
call allocate_and_initial(bc%sw4_front,bc%sw4_back,nzw,bc_len,nx,nt,mem)
call allocate_and_initial(bc%sw4_nt,bc%sw4_nt_1,nzw_pml,ny_pml,nx_pml,mem)
end subroutine bc_init_sm2

!=====================================================================

subroutine bc_final_sm2(bc,mem)
type(bc3d_sm2), intent(inout) :: bc
integer(8),optional,intent(inout) :: mem
call deallocate_and_free(bc%p_top,bc%p_bottom,bc%p_left,bc%p_right,&
   bc%p_front,bc%p_back,mem)
call deallocate_and_free(bc%sw_top,bc%sw_bottom,bc%sw_left,bc%sw_right,&
   bc%sw_front,bc%sw_back,mem)
call deallocate_and_free(bc%sw2_top,bc%sw2_bottom,bc%sw2_left,&
   bc%sw2_right,bc%sw2_front,bc%sw2_back,mem)
call deallocate_and_free(bc%sw3_top,bc%sw3_bottom,bc%sw3_left,&
   bc%sw3_right,bc%sw3_front,bc%sw3_back,mem)
call deallocate_and_free(bc%sw4_top,bc%sw4_bottom,bc%sw4_left,&
   bc%sw4_right,bc%sw4_front,bc%sw4_back,mem)
call deallocate_and_free(bc%p_nt,bc%p_nt_1,bc%sw_nt,bc%sw_nt_1,&
  bc%sw2_nt,bc%sw2_nt_1,mem)
call deallocate_and_free(bc%sw3_nt,bc%sw3_nt_1,bc%sw4_nt,bc%sw4_nt_1,mem)
end subroutine bc_final_sm2
!=====================================================================

subroutine saveBC3d(p,nz,ny,nx,nzp,nyp,nxp,npml,top,bottom,left,right,&
   front,back,bc_len)
integer, intent(in)    :: nz,ny,nx,nzp,nyp,nxp,npml,bc_len
real,    intent(in)    :: p(:,:,:)
real,    intent(inout) :: top(:,:,:),bottom(:,:,:),left(:,:,:),&
   right(:,:,:),front(:,:,:),back(:,:,:)
bc_len_1=bc_len-1
!$OMP PARALLEL DO PRIVATE(iz,iy,ix)
do ix=1,nx
   do iy=1,ny
      do iz=1,bc_len
         top(iz,iy,ix)    = p(npml+iz-1,npml+iy,npml+ix)
         bottom(iz,iy,ix) = p(nzp+iz-bc_len_1,npml+iy,npml+ix)
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(iz,iy,ix)
do ix=1,bc_len
   do iy=1,ny
      do iz=1,nz
         left(iz,iy,ix)  = p(npml+iz,npml+iy,npml+ix-1)
         right(iz,iy,ix) = p(npml+iz,npml+iy,nxp+ix-bc_len_1)
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(iz,iy,ix)
do ix=1,nx
   do iy=1,bc_len
      do iz=1,nz
         front(iz,iy,ix) = p(npml+iz,npml+iy-1,npml+ix)
         back(iz,iy,ix)  = p(npml+iz,nyp+iy-bc_len_1,npml+ix)
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
end subroutine saveBC3d

!=====================================================================

subroutine loadBC3d(p,nz,ny,nx,nzp,nyp,nxp,npml,top,bottom,left,right,&
   front,back,bc_len)
integer, intent(in)    :: nz,ny,nx,nzp,nyp,nxp,npml,bc_len
real,    intent(in)    :: top(:,:,:),bottom(:,:,:),left(:,:,:),&
   right(:,:,:),front(:,:,:),back(:,:,:)
real,    intent(inout) :: p(:,:,:)
bc_len_1=bc_len-1
!$OMP PARALLEL DO PRIVATE(iz,iy,ix)
do ix=1,nx
   do iy=1,ny
      do iz=1,bc_len
         p(npml+iz-1,npml+iy,npml+ix)       = top(iz,iy,ix)
         p(nzp+iz-bc_len_1,npml+iy,npml+ix) = bottom(iz,iy,ix)
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(iz,iy,ix)
do ix=1,bc_len
   do iy=1,ny
      do iz=1,nz
         p(npml+iz,npml+iy,npml+ix-1)       = left(iz,iy,ix)
         p(npml+iz,npml+iy,nxp+ix-bc_len_1) = right(iz,iy,ix)
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(iz,iy,ix)
do ix=1,nx
   do iy=1,bc_len
      do iz=1,nz
         p(npml+iz,npml+iy-1,npml+ix)       = front(iz,iy,ix)
         p(npml+iz,nyp+iy-bc_len_1,npml+iy) = back(iz,iy,ix)
      enddo
   enddo
enddo   
!$OMP END PARALLEL DO
end subroutine loadBC3d

!=====================================================================
! Absorbing boundary condition
subroutine abc_get_damp1d(npml, dx, cmin, damp)
integer, intent(in)  :: npml
real,    intent(in)  :: dx, cmin
real,    intent(out) :: damp(:)
a = (npml-1)*dx
kappa = 3.0*cmin*log(10000000.0)/(2.0*a)
do ix=1,npml
   xa = real(ix-1)*dx/a
   damp(ix) = kappa*xa*xa
enddo
end subroutine abc_get_damp1d

!===============================================================================

subroutine abc_get_damp3d(nx, ny, nz, npml, dx, cmin, damp)
integer, intent(in)  :: nx, ny, nz, npml
real,    intent(in)  :: dx, cmin
real,    intent(out) :: damp(:,:,:)
real, allocatable    :: damp1d(:)
allocate(damp1d(npml))
! Setup 1D PML damping term
a = (npml-1)*dx
kappa = 3.0*cmin*log(10000000.0)/(2.0*a)
do ix=1,npml
   xa = real(ix-1)*dx/a
   damp1d(ix) = kappa*xa*xa
enddo
! Setup 3D PML damping array
damp = 0.0
iy = 1
do ix=1,npml
   damp(:,:,npml-ix+1) = damp1d(ix)
   damp(:,:,nx+ix+npml) = damp1d(ix)
enddo
do iz=1,npml
   do ix=npml-(iz-1),nx+npml+(iz-1)
      damp(npml-iz+1,:,ix) = damp1d(iz)
      damp(nz+iz+npml,:,ix) = damp1d(iz)
   enddo
enddo
!$OMP PARALLEL DO PRIVATE(iz,iy,ix)
do iy=1,npml
   do ix=npml-(iy-1),nx+npml+(iy-1)
      do iz=npml-(iy-1),nz+npml+(iy-1)
         damp(iz,npml-iy+1,ix) = damp1d(iy)
         damp(iz,ny+iy+npml,ix) = damp1d(iy)
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
deallocate(damp1d)
end subroutine abc_get_damp3d

!=====================================================================

SUBROUTINE cos_coef( a, n )
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL   , INTENT(OUT) :: a(n)
INTEGER :: ii

DO ii=1,n
    a(ii) = 0.5*( 1.0 - COS(PI*ii/(n+1)) )
ENDDO
END SUBROUTINE cos_coef

!=====================================================================

end module module_boundary3d
