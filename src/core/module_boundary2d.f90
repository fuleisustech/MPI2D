module module_boundary2d
use module_datatype
use module_array
use module_utility
implicit none
integer, private :: nz_pml, nzw_pml, nx_pml, bc_len

type bc2d_general
   real, allocatable :: p_top(:,:,:),p_bottom(:,:,:),&
     p_left(:,:,:),p_right(:,:,:),p_nt(:,:),p_nt_1(:,:),&
     u_top(:,:,:),u_bottom(:,:,:),u_left(:,:,:),u_right(:,:,:),&
     u_nt(:,:),u_nt_1(:,:),w_top(:,:,:),w_bottom(:,:,:),&
     w_left(:,:,:),w_right(:,:,:),w_nt(:,:),w_nt_1(:,:)
end type bc2d_general

type bc2d_sm
   real, allocatable :: sw_top(:,:,:),sw_bottom(:,:,:),&
     sw_left(:,:,:),sw_right(:,:,:),sw_nt(:,:),sw_nt_1(:,:),&
     sh_top(:,:,:),sh_bottom(:,:,:),sh_left(:,:,:),&
     sh_right(:,:,:),sh_nt(:,:),sh_nt_1(:,:),&
     p_top(:,:,:),p_bottom(:,:,:),p_left(:,:,:),&
     p_right(:,:,:),p_nt(:,:),p_nt_1(:,:)
end type bc2d_sm

type bc2d_sm1
   real, allocatable :: sw_top(:,:,:),sw_bottom(:,:,:),&
     sw_left(:,:,:),sw_right(:,:,:),sw_nt(:,:),sw_nt_1(:,:),&
     sw2_top(:,:,:),sw2_bottom(:,:,:),sw2_left(:,:,:),&
     sw2_right(:,:,:),sw2_nt(:,:),sw2_nt_1(:,:),&
     p_top(:,:,:),p_bottom(:,:,:),p_left(:,:,:),&
     p_right(:,:,:),p_nt(:,:),p_nt_1(:,:)
end type bc2d_sm1

type bc2d_sm2
   real, allocatable :: sw_top(:,:,:),sw_bottom(:,:,:),&
     sw_left(:,:,:),sw_right(:,:,:),sw_nt(:,:),sw_nt_1(:,:),&
     sw2_top(:,:,:),sw2_bottom(:,:,:),sw2_left(:,:,:),&
     sw2_right(:,:,:),sw2_nt(:,:),sw2_nt_1(:,:),&
     sw3_top(:,:,:),sw3_bottom(:,:,:),sw3_left(:,:,:),&
     sw3_right(:,:,:),sw3_nt(:,:),sw3_nt_1(:,:),&
     sw4_top(:,:,:),sw4_bottom(:,:,:),sw4_left(:,:,:),&
     sw4_right(:,:,:),sw4_nt(:,:),sw4_nt_1(:,:),&
     p_top(:,:,:),p_bottom(:,:,:),p_left(:,:,:),&
     p_right(:,:,:),p_nt(:,:),p_nt_1(:,:)
end type bc2d_sm2

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

subroutine local_general(nz,nx,npml)
integer, intent(in) :: nz,nx,npml
nz_pml = nz + 2 * npml
nx_pml = nx + 2 * npml
end subroutine local_general

!=====================================================================

subroutine local_water(nz,nzw,nx,npml)
integer, intent(in) :: nz,nzw,nx,npml
nz_pml  = nz  + 2 * npml
nzw_pml = nzw + 2 * npml
nx_pml  = nx  + 2 * npml
end subroutine local_water

!=====================================================================

subroutine bc_init_general(bc,nz,nx,npml,nt,fd_order)
type(bc2d_general), intent(inout) :: bc
integer,  intent(in)  :: nz,nx,npml,nt,fd_order
bc_len=(fd_order-20)/2+1
call local_variables(nz,nx,npml)
call allocate_and_initial(bc%p_top,bc%p_bottom,bc_len,nx,nt)
call allocate_and_initial(bc%p_left,bc%p_right,nz,bc_len,nt)
call allocate_and_initial(bc%p_nt,bc%p_nt_1,nz_pml,nx_pml)
end subroutine bc_init_general

!=====================================================================

subroutine bc_final_general(bc)
type(bc2d_general) ,intent(inout) :: bc
deallocate(bc%p_top,bc%p_bottom,bc%p_left,bc%p_right,&
           bc%p_nt,bc%p_nt_1)
end subroutine bc_final_general

!=====================================================================

subroutine bc_init_sm(bc,nz,nzw,nx,npml,nt,fd_order)
type(bc2d_sm),    intent(inout) :: bc
integer,        intent(in)  :: nz,nzw,nx,npml,nt,fd_order
bc_len=(fd_order-20)/2+1
call local_variables(nz,nzw,nx,npml)
call allocate_and_initial(bc%p_top,bc%p_bottom,bc_len,nx,nt)
call allocate_and_initial(bc%p_left,bc%p_right,nz,bc_len,nt)
call allocate_and_initial(bc%p_nt,bc%p_nt_1,nz_pml,nx_pml)
call allocate_and_initial(bc%sw_top,bc%sw_bottom,bc_len,nx,nt)
call allocate_and_initial(bc%sw_left,bc%sw_right,nzw,bc_len,nt)
call allocate_and_initial(bc%sw_nt,bc%sw_nt_1,nzw_pml,nx_pml)
call allocate_and_initial(bc%sh_top,bc%sh_bottom,bc_len,nx,nt)
call allocate_and_initial(bc%sh_left,bc%sh_right,nzw,bc_len,nt)
call allocate_and_initial(bc%sh_nt,bc%sh_nt_1,nzw_pml,nx_pml)
end subroutine bc_init_sm

!=====================================================================

subroutine bc_final_sm(bc)
type(bc2d_sm), intent(inout) :: bc
deallocate(bc%p_top,bc%p_bottom,bc%p_left,bc%p_right,&
           bc%p_nt,bc%p_nt_1,bc%sw_top,bc%sw_bottom,&
           bc%sw_left,bc%sw_right,bc%sw_nt,bc%sw_nt_1,&
           bc%sh_top,bc%sh_bottom,bc%sh_left,bc%sh_right,&
           bc%sh_nt,bc%sh_nt_1)
end subroutine bc_final_sm

!=====================================================================

subroutine bc_init_sm1(bc,nz,nzw,nx,npml,nt,fd_order)
type(bc2d_sm1),    intent(inout) :: bc
integer,         intent(in)    :: nz,nzw,nx,npml,nt,fd_order
bc_len=(fd_order-20)/2+1
call local_variables(nz,nzw,nx,npml)
call allocate_and_initial(bc%p_top,bc%p_bottom,bc_len,nx,nt)
call allocate_and_initial(bc%p_left,bc%p_right,nz,bc_len,nt)
call allocate_and_initial(bc%p_nt,bc%p_nt_1,nz_pml,nx_pml)
call allocate_and_initial(bc%sw_top,bc%sw_bottom,bc_len,nx,nt)
call allocate_and_initial(bc%sw_left,bc%sw_right,nzw,bc_len,nt)
call allocate_and_initial(bc%sw_nt,bc%sw_nt_1,nzw_pml,nx_pml)
call allocate_and_initial(bc%sw2_top,bc%sw2_bottom,bc_len,nx,nt)
call allocate_and_initial(bc%sw2_left,bc%sw2_right,nzw,bc_len,nt)
call allocate_and_initial(bc%sw2_nt,bc%sw2_nt_1,nzw_pml,nx_pml)
end subroutine bc_init_sm1

!=====================================================================

subroutine bc_final_sm1(bc)
type(bc2d_sm1), intent(inout) :: bc
deallocate(bc%p_top,bc%p_bottom,bc%p_left,bc%p_right,&
           bc%p_nt,bc%p_nt_1,bc%sw_top,bc%sw_bottom,&
           bc%sw_left,bc%sw_right,bc%sw_nt,bc%sw_nt_1,&
           bc%sw2_top,bc%sw2_bottom,bc%sw2_left,bc%sw2_right,&
           bc%sw2_nt,bc%sw2_nt_1)
end subroutine bc_final_sm1

!=====================================================================

subroutine bc_init_sm2(bc,nz,nzw,nx,npml,nt,fd_order)
type(bc2d_sm2), intent(inout) :: bc
integer,      intent(in)    :: nz,nzw,nx,npml,nt,fd_order
bc_len=(fd_order-20)/2+1
call local_variables(nz,nzw,nx,npml)
call allocate_and_initial(bc%sw_top,bc%sw_bottom,2,nx,nt)
call allocate_and_initial(bc%sw_left,bc%sw_right,nzw,2,nt)
call allocate_and_initial(bc%sw_nt,bc%sw_nt_1,nzw_pml,nx_pml)
call allocate_and_initial(bc%sw2_top,bc%sw2_bottom,2,nx,nt)
call allocate_and_initial(bc%sw2_left,bc%sw2_right,nzw,2,nt)
call allocate_and_initial(bc%sw2_nt,bc%sw2_nt_1,nzw_pml,nx_pml)
call allocate_and_initial(bc%sw3_top,bc%sw3_bottom,2,nx,nt)
call allocate_and_initial(bc%sw3_left,bc%sw3_right,nzw,2,nt)
call allocate_and_initial(bc%sw3_nt,bc%sw3_nt_1,nzw_pml,nx_pml)
call allocate_and_initial(bc%sw4_top,bc%sw4_bottom,2,nx,nt)
call allocate_and_initial(bc%sw4_left,bc%sw4_right,nzw,2,nt)
call allocate_and_initial(bc%sw4_nt,bc%sw4_nt_1,nzw_pml,nx_pml)
call allocate_and_initial(bc%p_top,bc%p_bottom,2,nx,nt)
call allocate_and_initial(bc%p_left,bc%p_right,nz,2,nt)
call allocate_and_initial(bc%p_nt,bc%p_nt_1,nz_pml,nx_pml)
end subroutine bc_init_sm2

!=====================================================================

subroutine bc_final_sm2(bc)
type(bc2d_sm2), intent(inout) :: bc
deallocate(bc%p_top,bc%p_bottom,bc%p_left,bc%p_right,&
           bc%p_nt,bc%p_nt_1,bc%sw_top,bc%sw_bottom,&
           bc%sw_left,bc%sw_right,bc%sw_nt,bc%sw_nt_1,&
           bc%sw2_top,bc%sw2_bottom,bc%sw2_left,bc%sw2_right,&
           bc%sw2_nt,bc%sw2_nt_1,bc%sw3_top,bc%sw3_bottom,&
           bc%sw3_left,bc%sw3_right,bc%sw3_nt,bc%sw3_nt_1,&
           bc%sw4_top,bc%sw4_bottom,bc%sw4_left,&
           bc%sw4_right,bc%sw4_nt,bc%sw4_nt_1)
end subroutine bc_final_sm2
!=====================================================================

subroutine saveBC2d(p,nz,nx,nzp,nxp,npml,top,bottom,left,right,bc_len)
integer, intent(in)    :: nz,nx,nzp,nxp,npml,bc_len
real,    intent(in)    :: p(:,:)
real,    intent(inout) :: top(:,:),bottom(:,:),left(:,:),right(:,:)
integer                :: ix, iz, bc_len_1
bc_len_1=bc_len-1
do ix=1,nx
   do iz=1,bc_len
      top(iz,ix)    = p(npml+iz-1,npml+ix)
      bottom(iz,ix) = p(nzp+iz-bc_len_1,npml+ix)
   enddo
enddo
do ix=1,bc_len
   do iz=1,nz
      left(iz,ix)  = p(npml+iz,npml+ix-1)
      right(iz,ix) = p(npml+iz,nxp+ix-bc_len_1)
   enddo
enddo
end subroutine saveBC2d

!=====================================================================

subroutine loadBC2d(p,nz,nx,nzp,nxp,npml,top,bottom,left,right,bc_len)
integer, intent(in)    :: nz,nx,nzp,nxp,npml,bc_len
real,    intent(in)    :: top(:,:),bottom(:,:),left(:,:),right(:,:)
real,    intent(inout) :: p(:,:)
integer                :: ix, iz, bc_len_1
bc_len_1=bc_len-1
do ix=1,nx
   do iz=1,bc_len
      p(npml+iz-1,npml+ix)    = top(iz,ix)
      p(nzp+iz-bc_len_1,npml+ix) = bottom(iz,ix)
   enddo
enddo
do ix=1,bc_len
   do iz=1,nz
      p(npml+iz,npml+ix-1)    = left(iz,ix)
      p(npml+iz,nxp+ix-bc_len_1) = right(iz,ix)
   enddo
enddo
end subroutine loadBC2d

!=====================================================================
! Absorbing boundary condition
subroutine abc_get_damp1d(npml, dx, cmin, damp)
integer, intent(in)  :: npml
real,    intent(in)  :: dx, cmin
real,    intent(out) :: damp(:)
integer              :: ix
real                 :: a, xa, kappa
a = (npml-1)*dx
kappa = 3.0*cmin*log(10000000.0)/(2.0*a)
do ix=1,npml
   xa = real(ix-1)*dx/a
   damp(ix) = kappa*xa*xa
enddo
end subroutine abc_get_damp1d

!===============================================================================

subroutine abc_get_damp2d(nx, nz, npml, dx, cmin, damp)
integer, intent(in)  :: nx, nz, npml
real,    intent(in)  :: dx, cmin
real,    intent(out) :: damp(:,:)
integer              :: ix, iz
real                 :: a, xa, kappa
real, allocatable    :: damp1d(:)
allocate(damp1d(npml))
a = (npml-1)*dx
kappa = 3.0*cmin*log(10000000.0)/(2.0*a)
do ix=1,npml
   xa = real(ix-1)*dx/a
   damp1d(ix) = kappa*xa*xa
enddo
! Setup 2D PML damping array
damp = 0.0
do ix=1,npml
   do iz=1,nz+2*npml
      damp(iz,npml-ix+1) = damp1d(ix)
      damp(iz,nx+npml+ix) = damp1d(ix)
   enddo
enddo
do iz=1,npml
   do ix=1+(npml-(iz-1)),nx+npml+(iz-1)
      damp(npml-iz+1,ix) = damp1d(iz)
      damp(nz+npml+iz,ix) = damp1d(iz)
   enddo
enddo
deallocate(damp1d)
end subroutine abc_get_damp2d

!===============================================================================
! make sponge boundary coefficients
!-------------------------------------------------------------------------------
SUBROUTINE spgcoef(spgx,spgz,spcx1,spcx2,spcz1,spcz2,v,nx_pml,nz_pml,&
                    npml,dx,dz,dt)
IMPLICIT NONE
INTEGER   , INTENT(IN)  :: nx_pml,nz_pml,npml
REAL      , INTENT(IN ) :: v(nz_pml,nx_pml)
REAL      , INTENT(OUT) :: spgx(nx_pml), spgz(nz_pml)
REAL      , INTENT(OUT) :: spcz1(nx_pml), spcz2(nx_pml)
REAL      , INTENT(OUT) :: spcx1(nz_pml), spcx2(nz_pml)
REAL      , ALLOCATABLE :: tmp_spz1(:), tmp_spz2(:), tmp_spx(:)
!REAL      , ALLOCATABLE :: spctmp(:,:)
!REAL     :: scoef, acoef, dist, rr, xr, yr, zratio, rcoef, cfmdlng, z
REAL     :: scoef, acoef, rcoef, cfmdlng
REAL      ,INTENT(IN)   :: dx,dz,dt
INTEGER                 :: ix,iz

ALLOCATE( tmp_spz1(npml), tmp_spz2(npml), tmp_spx(npml) )
CALL cos_coef( tmp_spz1, npml )
CALL cos_coef( tmp_spz2, npml )
CALL cos_coef( tmp_spx , npml  )

! coef. for top ----------------------------------------------------------------
cfmdlng = 1.0
!IF ( lmodeling ) cfmdlng = 0.6
spgz(:) = 0.0
scoef = cfmdlng*1000*dt/dz/npml
DO iz=1,npml
    spgz(npml-iz+1    ) = scoef*tmp_spz1(iz)
ENDDO

! coef. for bottom -------------------------------------------------------------
cfmdlng = 0.6
!IF ( lmodeling ) cfmdlng = 0.4
scoef = cfmdlng*1000*dt/dz/npml
DO iz=1,npml
    spgz(nz_pml-npml+iz) = scoef*tmp_spz2(iz)
ENDDO

! coef. for side walls ---------------------------------------------------------
cfmdlng = 0.8
spgx = 0.0
scoef = cfmdlng*1000*dt/dx/npml
DO ix=1,npml
    spgx(npml-ix+1    ) = scoef*tmp_spx(ix)
    spgx(nx_pml-npml+ix) = spgx(npml-ix+1)
ENDDO

DEALLOCATE( tmp_spz1, tmp_spz2, tmp_spx)

!-- Compute top boundary: spcz1(nx_pml) -----------------------------------
DO ix=1,nx_pml
   scoef = 0.3
   spcz1(ix) = SQRT(v(npml+1,ix))*scoef
ENDDO

!-- Compute bottom boundary: spcz2(nxcom,nycom) --------------------------------
DO ix=1,nx_pml
   scoef = 0.3
   spcz2(ix) = SQRT(v(nz_pml-npml+1,ix))*scoef
ENDDO

!-- Compute left boundary: spxc1(nycom,nzcom) ----------------------------------
rcoef = 1.0
acoef = 0.1
DO iz = 1, nz_pml
   scoef = 0.3
   spcx1(iz) = SQRT(v(iz,npml+1))*scoef
ENDDO

!-- Compute right boundary: spcx2(nycom,nzcom) ---------------------------------
DO iz = 1, nz_pml
   scoef=0.3
   spcx2(iz) = SQRT(v(iz,nx_pml-npml+1))*scoef
ENDDO
END SUBROUTINE spgcoef

!=====================================================================

SUBROUTINE cos_coef( a, n )
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL   , INTENT(OUT) :: a(n)
INTEGER :: ii
REAL    :: pi

pi = 4*ATAN(1.0)
DO ii=1,n
    a(ii) = 0.5*( 1.0 - COS(pi*ii/(n+1)) )
ENDDO
END SUBROUTINE cos_coef

!=====================================================================

SUBROUTINE sponge2D( u2, spgx, spgz, spcx1, spcx2, spcz1,&
                     spcz2,nx_pml,nz_pml,npml )
IMPLICIT NONE
INTEGER, INTENT(IN)    :: nx_pml,nz_pml,npml
REAL   , INTENT(IN)    :: spgx(:), spgz(:)
REAL   , INTENT(IN)    :: spcz1(:), spcz2(:)
REAL   , INTENT(IN)    :: spcx1(:), spcx2(:)
REAL   , INTENT(INOUT) :: u2(:,:)
INTEGER                :: ix,iz

!-- top sponge boundary
DO iz=1,npml
    DO ix=1,nx_pml
        u2(iz,ix) = u2(iz,ix) * ( 1.0 - spgz(iz)*spcz1(ix) )
    ENDDO
ENDDO

!-- bottom sponge boundary
DO iz=nz_pml-npml+1,nz_pml
    DO ix=1,nx_pml
        u2(iz,ix) = u2(iz,ix) * ( 1.0 - spgz(iz)*spcz2(ix) )
    ENDDO
ENDDO

DO iz=1,nz_pml
!-- left sponge boundary
    DO ix=1,npml
        u2(iz,ix) = u2(iz,ix) * ( 1.0 - spgx(ix)*spcx1(iz) )
    ENDDO
!-- right sponge boundary
    DO ix=nx_pml-npml+1,nx_pml
        u2(iz,ix) = u2(iz,ix) * ( 1.0 - spgx(ix)*spcx2(iz) )
    ENDDO
ENDDO
END SUBROUTINE sponge2D

end module module_boundary2d
