! Written by Bowen Guo
! @version 1 2014-05-25
! Based on the report 'Elastic wavefield modeling in 3D by fourth-order staggered grid finite difference technique and on the report '3D fourth-Order Staggered-Grid Finite-Differences Schemes: Stability and Grid Dispersion 
module module_e3d

use module_global
use module_array
use module_utility
use module_boundary3d
use module_fd3d
use module_string
use module_io
use module_datatype
use module_aperture3d

implicit none
character(len=100),private :: fd_type,source_type,data_type
integer,private      :: ix,iy,iz,it,isx,isy,isz,ig,nx,ny,nz,npml,nxp,nyp,nzp,nt,dnt_mod,nx_pml,ny_pml,nz_pml,fd_order,nit_display,fs_thick
real,private     :: dtdx,dx,dy,dz,dt,vmin
integer,private,allocatable  :: igx(:),igy(:),igz(:),fs(:,:)
logical,private              :: isWriteData,isFS
real,private,allocatable     :: tau_xx(:,:,:),tau_yy(:,:,:),tau_zz(:,:,:),tau_xy(:,:,:),tau_xz(:,:,:),tau_yz(:,:,:),total(:,:,:)
real,private,allocatable     :: u(:,:,:),v(:,:,:),w(:,:,:)
real,private,allocatable     :: damp(:,:,:),kappa(:,:,:),temp(:,:,:)

integer(kind=MPI_OFFSET_KIND) :: disp

interface e3d_mod
   module procedure modeling
end interface e3d_mod

interface add_source
   module procedure add_source_p
   module procedure add_source_s
   module procedure add_source_tauzz
   module procedure add_source_w
end interface add_source

interface store_data
   module procedure store_data_w
   module procedure store_data_uvw
end interface store_data

interface add_FS
   module procedure add_FS_vel
   module procedure add_FS_stress
end interface add_FS

private initial_variables,finalize_variables,nx_ny_nz_pml,i_v_kernal,setup_source_receiver

contains
!=========================================================================================
subroutine modeling(fdm,s,output_w,output_v,output_u,mem)
type(fdmod3d),intent(inout)  :: fdm
type(shot3d),intent(inout)   :: s
character(len=*),intent(in)  :: output_w,output_v,output_u
integer                      :: fid,it_out
integer(8),optional,intent(inout)  :: mem
call initial_variables(fdm,s,mem)
if (isWriteData) then
   if (data_type.eq."w") then
      open(10,file=output_w,access="direct",recl=I4)
   elseif (data_type.eq."uvw") then
      open(12,file=output_u,access="direct",recl=I4)
      open(11,file=output_v,access="direct",recl=I4)
      open(10,file=output_w,access="direct",recl=I4)
   endif 
else
   s%seis = 0.0
endif
call setup_source_receiver(s,mem)
! Allocate stress wavefield
call allocate_and_initial(tau_xx,tau_yy,tau_zz,tau_xy,tau_xz,tau_yz,nz_pml,ny_pml,nx_pml)
! Allocate velocity wavefield
call allocate_and_initial(u,v,w,nz_pml,ny_pml,nx_pml)
call allocate_and_initial(total,nz_pml,ny_pml,nx_pml)
do it=1,nt
   if (mod(it,100)==1) then
      write(*,*) " The process has finished", (real(it)/real(nt)*100.0),"%"
   endif
   ! Add source
   if (source_type.eq."p")  then
      call add_source(fdm%s%sou(it),tau_xx,tau_yy,tau_zz)
      !call add_source(u,w,fdm%s%sou(it))
   elseif (source_type.eq."tau_zz") then
      call add_source(tau_zz,fdm%s%sou(it))
   elseif (source_type.eq."s") then
      call add_source(u,v,w,fdm%s%sou(it))
   elseif (source_type.eq."w") then
      call add_source(fdm%s%sou(it),w)
   endif
   
   
   ! Update velocity field 
   call e3d_abc_kernal_vel(u,tau_xx,tau_xy,tau_xz,temp,fdm%m%den,dtdx,nx_pml,ny_pml,nz_pml,fd_order,isFS,npml)
   call e3d_abc_kernal_vel(fd_order,v,tau_xy,tau_yy,tau_yz,temp,fdm%m%den,dtdx,nx_pml,ny_pml,nz_pml,isFS,npml)
   call e3d_abc_kernal_vel(nx_pml,ny_pml,nz_pml,w,tau_xz,tau_yz,tau_zz,temp,fdm%m%den,dtdx,fd_order,isFS,npml)
   ! Free surface boundary condition for velocity
   if (isFS) then 
      call add_FS(u,v,w,fdm%m%lambda,fdm%m%mu,fd_order)
   endif
   ! Update Stress field
   call e3d_abc_kernal_stress(total,u,v,w,fdm%m%lambda,dtdx,nx_pml,ny_pml,nz_pml,fd_order,isFS,npml)
   call e3d_abc_kernal_stress(total,fd_order,u,tau_xx,temp,fdm%m%mu,dtdx,nx_pml,ny_pml,nz_pml,isFS,npml)
   call e3d_abc_kernal_stress(fd_order,total,v,tau_yy,temp,fdm%m%mu,dtdx,nx_pml,ny_pml,nz_pml,isFS,npml)
   call e3d_abc_kernal_stress(nx_pml,ny_pml,nz_pml,total,w,tau_zz,temp,fdm%m%mu,dtdx,fd_order,isFS,npml)

   call e3d_abc_kernal_stress(fd_order,nx_pml,tau_xy,u,v,temp,fdm%m%mu,dtdx,ny_pml,nz_pml,isFS,npml)
   call e3d_abc_kernal_stress(fd_order,nx_pml,ny_pml,nz_pml,tau_xz,u,w,temp,fdm%m%mu,dtdx,isFS,npml)
   call e3d_abc_kernal_stress(fd_order,nx_pml,ny_pml,tau_yz,nz_pml,v,w,temp,fdm%m%mu,dtdx,isFS,npml)
   
   ! Free surface boundary condition for stress
   if (isFS) then
      call add_FS(tau_xz,tau_yz,tau_zz,fd_order)
   endif
   ! Write down the wavefield
   if (mod(it,20)==1) then  
    !  call snapshot("u_",u,it,nz_pml,ny_pml,nx_pml) 
    !  call snapshot("v_",v,it,nz_pml,ny_pml,nx_pml) 
    !  call snapshot("w_",w,it,nz_pml,ny_pml,nx_pml) 
    !  call snapshot("tauxx_",tau_xx,it,nz_pml,ny_pml,nx_pml) 
    !  call snapshot("tauyy_",tau_yy,it,nz_pml,ny_pml,nx_pml) 
    !  call snapshot("tauzz_",tau_zz,it,nz_pml,ny_pml,nx_pml) 
    !  call snapshot("tauxy_",tau_xy,it,nz_pml,ny_pml,nx_pml) 
    !  call snapshot("tauxz_",tau_xz,it,nz_pml,ny_pml,nx_pml) 
    !  call snapshot("tauyz_",tau_yz,it,nz_pml,ny_pml,nx_pml) 
   endif
   if (isWriteData) then
      if (dnt_mod.ne.1) then
         it_out=(it-1)/dnt_mod+1
         do ig=1,s%ng
            if (data_type.eq."w") then
               write(10,rec=(ig-1)*s%nt+it_out)  w(igz(ig),igy(ig),igx(ig))
            elseif (data_type.eq."uvw") then
               write(10,rec=(ig-1)*s%nt+it_out)  w(igz(ig),igy(ig),igx(ig))
               write(11,rec=(ig-1)*s%nt+it_out)  v(igz(ig),igy(ig),igx(ig))
               write(12,rec=(ig-1)*s%nt+it_out)  u(igz(ig),igy(ig),igx(ig))
            endif
         enddo
      endif
   else
      if (data_type.eq."w") then
         call store_data(w,s%seis_w,s%ng,it)
      elseif (data_type.eq."uvw") then
         call store_data(u,v,w,s%seis_u,s%seis_v,s%seis_w,s%ng,it)
      endif
   endif
enddo  ! End of time_marching loop
call finalize_variables()
end subroutine modeling
!=============================================================================
subroutine add_FS_vel(u,v,w,lambda,mu,fd_order)
integer,intent(in)    :: fd_order 
real,intent(inout)    :: u(:,:,:),v(:,:,:),w(:,:,:)
real,intent(in)       :: lambda(:,:,:),mu(:,:,:)
real                  :: alpha
do iz=1,fs_thick
   ! $OMP PARALLEM DO PRIVATE(ix,iy,alpha)
   do ix=2,nx_pml-1
      do iy=2,ny_pml-1
         alpha=lambda(npml,iy,ix)/(lambda(npml,iy,ix)+2*mu(npml,iy,ix))
         w(npml-iz,iy,ix)=w(npml+iz-1,iy,ix)+alpha*(u(npml,iy,ix+iz-1)-u(npml,iy,ix-iz)+v(npml,iy+iz-1,ix)-v(npml,iy-iz,ix))
         u(npml-iz,iy,ix)=u(npml+iz-1,iy,ix)+(u(npml+iz,iy,ix)-u(npml-iz+1,iy,ix)+w(npml,iy,ix+iz)-w(npml,iy,ix-iz+1))+(w(npml-1,iy,ix+iz)-w(npml-1,iy,ix-iz+1))
         v(npml-iz,iy,ix)=v(npml+iz-1,iy,ix)+(w(npml,iy+iz,ix)-w(npml,iy-iz+1,ix)+v(npml+iz,iy,ix)-v(npml-iz+1,iy,ix))+(w(npml-1,iy+iz,ix)-w(npml-1,iy-iz+1,ix))
      enddo
   enddo
   ! $OMP END PARALLEL DO
enddo
end subroutine add_FS_vel
!==============================================================================
subroutine add_FS_stress(tau_xz,tau_yz,tau_zz,fd_order)
integer,intent(in)    :: fd_order
real,intent(inout)    :: tau_xz(:,:,:),tau_yz(:,:,:),tau_zz(:,:,:)
do iz=1,fs_thick
   tau_zz(npml,:,:)=0.0
   tau_zz(npml-iz,:,:)=-tau_zz(npml+iz,:,:)
   tau_xz(npml-iz,:,:)=-tau_xz(npml+iz-1,:,:)
   tau_yz(npml-iz,:,:)=-tau_yz(npml+iz-1,:,:)
enddo

end subroutine add_FS_stress
!==============================================================================================
subroutine finalize_variables()
call deallocate_and_free(igx,igz)
call deallocate_and_free(temp,tau_xx,tau_yy,tau_zz,tau_xy)
call deallocate_and_free(tau_yz,tau_xz,u,v,w)
call deallocate_and_free(total)
end subroutine finalize_variables
!==============================================================================================
subroutine store_data_uvw(u,v,w,seis_u,seis_v,seis_w,ng,it)
integer,intent(in)  :: it,ng
real,intent(in)     :: u(:,:,:),v(:,:,:),w(:,:,:)
real,intent(inout)  :: seis_u(:,:),seis_v(:,:),seis_w(:,:)
integer             :: it_out   
do ig=1,ng
   if (dnt_mod.ne.1) then
      if (mod(it,dnt_mod).eq.1) then
         it_out=(it-1)/dnt_mod+1
         seis_u(it_out,ig)=u(igz(ig),igy(ig),igx(ig))
         seis_v(it_out,ig)=v(igz(ig),igy(ig),igx(ig))
         seis_w(it_out,ig)=w(igz(ig),igy(ig),igx(ig))
      endif
   else 
         seis_u(it,ig)=u(igz(ig),igy(ig),igx(ig))
         seis_v(it,ig)=v(igz(ig),igy(ig),igx(ig))
         seis_w(it,ig)=w(igz(ig),igy(ig),igx(ig))
   endif
enddo
end subroutine store_data_uvw
!===============================================================================================
subroutine  store_data_w(w,seis,ng,it)
integer,intent(in) :: it,ng
real,intent(inout)   :: seis(:,:)
real,intent(in)    :: w(:,:,:)
integer            :: it_out
do ig=1,ng
   if (dnt_mod.ne.1) then
      if (mod(it,dnt_mod).eq.1) then   
         it_out=(it-1)/dnt_mod+1
         seis(it_out,ig)=w(igz(ig),igy(ig),igx(ig))
      endif
   else
      seis(it,ig)=w(igz(ig),igy(ig),igx(ig))
   endif
enddo
end subroutine store_data_w 
!===============================================================================================
subroutine add_source_tauzz(tau_zz,s)
real,intent(in)     :: s
real,intent(inout)  :: tau_zz(:,:,:)

tau_zz(isz,isy,isx)=tau_zz(isz,isy,isx)+s
end subroutine add_source_tauzz
!==============================================================================================
subroutine add_source_w(s,w)
real,intent(in)     :: s
real,intent(inout)  :: w(:,:,:)

w(isz,isy,isx)=w(isz,isy,isx)+s
end subroutine add_source_w
!===============================================================================================
subroutine add_source_s(u,v,w,s)
real,intent(inout)   :: u(:,:,:),v(:,:,:),w(:,:,:)
real,intent(in)      :: s
u(isz,isy,isx)=u(isz,isy,isx)+s
v(isz,isy,isx)=v(isz,isy,isx)+s
w(isz,isy,isx)=w(isz,isy,isx)+s
u(isz-1,isy,isx)=u(isz-1,isy,isx)-s
v(isz,isy-1,isz)=v(isz,isy-1,isz)-s
w(isz,isy,isx+1)=w(isz,isy,isx+1)+s
end subroutine add_source_s
!===============================================================================================
subroutine add_source_p(s,tau_xx,tau_yy,tau_zz)
real,intent(in)      :: s
real,intent(inout)   :: tau_xx(:,:,:),tau_zz(:,:,:),tau_yy(:,:,:)
tau_xx(isz,isy,isx)=tau_xx(isz,isy,isx)+s
tau_zz(isz,isy,isx)=tau_zz(isz,isy,isx)+s
tau_yy(isz,isy,isx)=tau_yy(isz,isy,isx)+s
end subroutine add_source_p

!================================================================================================
subroutine setup_source_receiver(seis,mem)
type(shot3d),    intent(in)         :: seis
integer(8),optional,intent(inout)   :: mem
call allocate_and_initial(igx,igy,igz,seis%ng,mem)
! Setup source position
isx = npml+int(seis%xs(1)/dx)+1
isy = npml+int(seis%ys(1)/dy)+1
isz = npml+int(seis%zs(1)/dz)+1
write(*,*) "isx = ",isx
write(*,*) "isy = ",isy
write(*,*) "isz = ",isz

! Source position is shifted downward 1 node if on the FS
if (isz == fs(isy,isx)) isz = isz + 1
! Setup receiver position
!$OMP PARALLEL DO PRIVATE(ig)
do ig=1,seis%ng
   igx(ig) = npml + int(seis%xg(ig)/dx)+1
   igy(ig) = npml + int(seis%yg(ig)/dy)+1
   igz(ig) = npml + int(seis%zg(ig)/dz)+1
   ! Geophone position is shifted downward 1 node if on the FS.
  ! if (igz(ig) == fs(igy(ig),igx(ig))) igz(ig) = igz(ig) + 1
enddo
!$OMP END PARALLEL DO
end subroutine setup_source_receiver
!=========================================================================================
subroutine initial_variables(fdm,s,mem)
type(fdmod3d),intent(in)   :: fdm
type(shot3d),intent(in)    :: s
integer(8),optional,intent(inout) :: mem
call copy_variable(fdm%f%fd_type,fd_type)
call copy_variable(fdm%f%fd_order,fd_order)
call copy_variable(fdm%f%nt,nt,fdm%m%nz,nz,fdm%m%ny,ny,fdm%m%nx,nx,fdm%f%npml,npml)
call copy_variable(fdm%f%dt,dt,fdm%m%dx,dx) 
call copy_variable(fdm%m%dy,dy,fdm%m%dz,dz) 
call copy_variable(fdm%f%isWriteData,isWriteData)
call copy_variable(fdm%f%isFS,isFS)

source_type=fdm%f%source_type;data_type=fdm%f%data_type
vmin=fdm%m%vmin
nit_display=int(nt/100)
dnt_mod=nint(s%dt/fdm%f%dt)
fs_thick=(fd_order-20)/2
call i_v_kernal(vmin,fdm%f%fs,mem)
end subroutine initial_variables
!=========================================================================================
subroutine nx_ny_nz_pml
nx_pml=nx+2*npml
ny_pml=ny+2*npml
nz_pml=nz+2*npml
nxp=nx+npml
nyp=ny+npml
nzp=nz+npml
end subroutine nx_ny_nz_pml
!=========================================================================================
subroutine i_v_kernal(v_min,fs_in,mem)
real,intent(in)     :: v_min
integer,intent(in)  :: fs_in(:,:)
integer(8),optional,intent(inout) :: mem
call nx_ny_nz_pml
dtdx=dt/dx
! Allocate free surface boundary
call allocate_and_initial(fs,ny_pml,nx_pml)
fs = fs_in
! Absorbing boundary condition is used
call allocate_and_initial(damp,kappa,temp,nz_pml,ny_pml,nx_pml,mem)
call abc_get_damp3d(nx,ny,nz,npml,dx,vmin,damp)
kappa=damp*dt
temp=1-kappa


call deallocate_and_free(damp,kappa,mem)
end subroutine i_v_kernal
!=============================================================================================

end module module_e3d 


