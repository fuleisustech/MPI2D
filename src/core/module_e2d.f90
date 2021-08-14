! written by Bowen Guo according to Xing Wang 's structure
module module_e2d

use module_global
use module_datatype
use module_array
use module_utility
use module_boundary2d
use module_fd2d
use module_string
use module_io

implicit none
character(len=100),private  :: fd_type,source_type,data_type
integer,private             :: ix,iz,it,isx,isz,ig,nx,nz,npml,nxp,nzp,&
nt,dnt_mod,nx_pml,nz_pml,ic,fd_order,bc_type,bc_len,fs_thick
real, private   :: dtdx,dx,dz,dt,vmin,c1t2
integer,private,allocatable :: igx(:),igz(:),fs(:)
logical,private :: isFS, isDipole, isSaveBC, isSaveWF, isWriteData, isReadData,sourcestress
real,private,allocatable    :: u(:,:),w(:,:),tau_xx(:,:),tau_xz(:,:),tau_zz(:,:)
real,private,allocatable    :: beta_dt(:,:)
real,private,allocatable    :: alpha_vel(:,:),alpha_lambda(:,:),alpha_mu(:,:),temp(:,:)
real,private,allocatable    :: damp(:,:),kappa(:,:)
integer(kind=MPI_OFFSET_KIND) :: disp
interface e2d_mod
   module procedure modeling
   module procedure modeling_fp
  ! module procedure modeling_bp
end interface e2d_mod


interface add_source
   module procedure add_source_p
   module procedure add_source_s
   module procedure add_source_tauzz
   module procedure add_source_w
end interface add_source


interface store_data
   module procedure store_data_w
   module procedure store_data_uw
   module procedure store_data_p
end interface store_data



private initial_variable,finalize_variable,nx_nz_pml,i_v_kernal,&
        setup_source_receiver

contains
!==========================================================================================
subroutine modeling(fdm,s,output_u,output_w,output_p)
type(fdmod2d),intent(inout)      :: fdm
type(shot2d) ,intent(inout)      :: s
!type(bc2d_general),intent(inout) :: bc
!type(wf2d),        intent(inout) :: wf
character(len=*), intent(in)     :: output_u,output_w,output_p
integer                          :: fid,it_out
call initial_variable(fdm,s)
if (isWriteData) then
   if (data_type.eq."w") then
      open(10,file=output_w,access="direct",recl=I4)
   elseif (data_type.eq."uw") then
      open(11,file=output_u,access="direct",recl=I4)
      open(10,file=output_w,access="direct",recl=I4)
   elseif (data_type.eq."p" ) then
      open(12,file=output_p,access="direct",recl=I4)
   endif    
else 
   s%seis=0.0
endif 
call setup_source_receiver(s)
call allocate_and_initial(tau_xx,tau_xz,tau_zz,u,w,nz_pml,nx_pml)
do it=1,nt
   
   ! Add source
   if (source_type.eq."P") then
      call add_source(tau_xx,tau_zz,fdm%s%sou(it))
   elseif (source_type.eq."tau_zz") then
      call add_source(tau_zz,fdm%s%sou(it))
   elseif (source_type.eq."S") then
      call add_source(fdm%s%sou(it),u,w)
   elseif (source_type.eq."w") then
      call add_source(fdm%s%sou(it),w)
   endif
   ! update velocity field 
   call e2d_abc_kernal_vel(u,w,tau_xx,tau_xz,tau_zz,temp,alpha_vel,nz_pml,nx_pml,fd_order,isFS,npml)
      
   ! Free surface Boundary Condition
   if (isFS) then 
      call add_FS_vel()
   endif
   ! update stress field
   !write(*,*) size(tau_xx,2)

   call e2d_abc_kernal_stress(u,w,tau_xx,tau_xz,tau_zz,temp,alpha_lambda,alpha_mu,nz_pml,nx_pml,fd_order,isFS,npml)
   
   if (isFS) then 
      call add_FS_stress()
   endif
   !if (mod(it,10)==1) call snapshot("w_",w,it,nz_pml,nx_pml)
   
   
   ! Record data
   if (isWriteData) then 
      if (dnt_mod.ne.1) then
         do ig=1,s%ng
            if (data_type.eq."w") then 
               write(10,rec=(ig-1)*s%nt+it_out)  w(igz(ig),igx(ig))
            elseif (data_type.eq."uw") then
               write(10,rec=(ig-1)*s%nt+it_out)  w(igz(ig),igx(ig))
               write(11,rec=(ig-1)*s%nt+it_out)  u(igz(ig),igx(ig))
            !elseif (data_type.eq."p") then
            !   write(12,rec=(ig-1)*s%nt+it_out) -0.5*(tau_xx(igz(ig),igx(ig))+tau_zz(igz(ig),igx(ig))) 
            endif
         enddo
       endif
   else
      if (data_type.eq."w") then
         call store_data(w,s%seis_w,s%ng,it)
      elseif (data_type.eq."uw") then
         call store_data(u,w,s%seis_u,s%seis_w,s%ng,it)
      elseif (data_type.eq."p") then
         call store_data(tau_xx,tau_zz,s%seis_p,s%ng,it)
      endif
   endif
enddo    ! End of time_marching loop
call finalize_variable()
call deallocate_and_free(u,w,tau_xx,tau_xz,tau_zz)
end subroutine modeling

!==========================================================================================
subroutine modeling_fp(fdm,s,wf)
type(fdmod2d),intent(inout)      :: fdm
type(shot2d) ,intent(inout)      :: s
!type(bc2d_general),intent(inout) :: bc
type(wf2d),        intent(inout) :: wf

call initial_variable(fdm,s)
call setup_source_receiver(s)
call allocate_and_initial(tau_xx,tau_xz,tau_zz,u,w,nz_pml,nx_pml)
do it=1,nt

   ! Add source
   if (source_type.eq."P") then
      call add_source(tau_xx,tau_zz,fdm%s%sou(it))
   elseif (source_type.eq."tau_zz") then
      call add_source(tau_zz,fdm%s%sou(it))
   elseif (source_type.eq."S") then
      call add_source(fdm%s%sou(it),u,w)
   elseif (source_type.eq."w") then
      call add_source(fdm%s%sou(it),w)
   endif

  ! update velocity field 
   call e2d_abc_kernal_vel(u,w,tau_xx,tau_xz,tau_zz,temp,alpha_vel,nz_pml,nx_pml,fd_order,isFS,npml)

   ! Free surface Boundary Condition
   if (isFS) then
      call add_FS_vel()
   endif
   ! update stress field
   call e2d_abc_kernal_stress(u,w,tau_xx,tau_xz,tau_zz,temp,alpha_lambda,alpha_mu,nz_pml,nx_pml,fd_order,isFS,npml)

   if (isFS) then
      call add_FS_stress()
   endif
   !if (mod(it,10)==1) call snapshot("w_",w,it,nz_pml,nx_pml)

  ! Record w component wavefield
    call store_wf_w(wf,w,it,nz,nx,npml)
enddo    ! End of time_marching loop
call finalize_variable()
call deallocate_and_free(u,w,tau_xx,tau_xz,tau_zz)
end subroutine modeling_fp
 
!==========================================================================================
subroutine add_FS_vel()
iz=npml
do ix=1,nx_pml-1
   w(iz-1,ix)=w(iz,ix)+alpha_lambda(iz,ix)/(alpha_lambda(iz,ix)+2*alpha_mu(iz,ix))*(u(iz,ix+1)-u(iz,ix))
enddo
if (fd_order.eq.24) then
   do ix=3,nx_pml-1
      u(iz-1,ix)= u(iz+2,ix)+S41/S42*(u(iz+1,ix)-u(iz,ix)+w(iz,ix)-w(iz,ix-1))+w(iz,ix+1)-w(iz,ix-2)
      w(iz-2,ix)= w(iz+1,ix)+S41/S42*(w(iz,ix)-w(iz-1,ix)+alpha_lambda(iz,ix)/(alpha_lambda(iz,ix)+2*alpha_mu(iz,ix))*(u(iz,ix)-u(iz,ix-1)))+&
                 alpha_lambda(iz,ix)/(alpha_lambda(iz,ix)+2*alpha_mu(iz,ix))*(u(iz,ix+1)-u(iz,ix-2))
   enddo
endif
end subroutine add_FS_vel
!==========================================================================================
subroutine add_FS_stress()
iz=npml
tau_zz(iz,:)=0.0;tau_xz(iz,:)=0.0;tau_xz(iz-1,:)=-tau_xz(iz+1,:)
if (fd_order.eq.24) then
   tau_zz(iz-1,:)=-tau_zz(iz+1,:)
   tau_xz(iz-1,:)=-tau_xz(iz+1,:)
endif
end subroutine add_FS_stress
!============================================================================================

subroutine store_data_uw(u,w,seis_u,seis_w,ng,it)
integer,intent(in)  :: it,ng
real,intent(in)     :: u(:,:),w(:,:)
real,intent(inout)  :: seis_u(:,:),seis_w(:,:)
integer             :: it_out
do ig=1,ng
   if (dnt_mod.ne.1) then
      if (mod(it,dnt_mod).eq.1) then
         it_out=(it-1)/dnt_mod+1
         seis_u(it_out,ig)=u(igz(ig),igx(ig))
         seis_w(it_out,ig)=w(igz(ig),igx(ig))
      endif
   else
         seis_u(it,ig)=u(igz(ig),igx(ig))
         seis_w(it,ig)=w(igz(ig),igx(ig))
   endif
enddo
end subroutine store_data_uw
!===============================================================================================
subroutine  store_data_w(w,seis,ng,it)
integer,intent(in) :: it,ng
real,intent(inout)   :: seis(:,:)
real,intent(in)    :: w(:,:)
integer            :: it_out
do ig=1,ng
   if (dnt_mod.ne.1) then
      if (mod(it,dnt_mod).eq.1) then
         it_out=(it-1)/dnt_mod+1
         seis(it_out,ig)=w(igz(ig),igx(ig))
      endif
   else
      seis(it,ig)=w(igz(ig),igx(ig))
   endif
enddo
end subroutine store_data_w

!===============================================================================================
subroutine  store_data_p(tau_xx,tau_zz,seis,ng,it)
integer,intent(in) :: it,ng
real,intent(inout)   :: seis(:,:)
real,intent(in)    :: tau_xx(:,:),tau_zz(:,:)
integer            :: it_out
do ig=1,ng
   seis(it,ig)=0.5*(tau_xx(igz(ig),igx(ig))+tau_zz(igz(ig),igx(ig)))
enddo
end subroutine store_data_p

!================================================
subroutine  store_wf_w(wf,w,it,nz,nx,npml)
integer,intent(in) :: it,nz,nx,npml
real,intent(in)    :: w(:,:)
type(wf2d), intent(inout) :: wf

integer  :: ix,iz
do ix=1,nx
 do iz=1,nz
  wf%wwf(it,iz,ix)=w(iz+npml,ix+npml);
 enddo
enddo
end subroutine store_wf_w

!========================================================================================
subroutine add_source_tauzz(tau_zz,s)

real, intent(in)    :: s
real, intent(inout) :: tau_zz(:,:)

tau_zz(isz,isx)=tau_zz(isz,isx)+s

end subroutine add_source_tauzz
!==================================================================================
subroutine add_source_w(s,w)
real,    intent(in)     :: s
real,    intent(inout)  :: w(:,:)

w(isz,isx)=w(isz,isx)+s
end subroutine add_source_w

!====================================================================================
subroutine add_source_p(tau_xx,tau_zz,s)
real,    intent(in)     :: s
real,    intent(inout)  :: tau_xx(:,:),tau_zz(:,:)

tau_xx(isz,isx)=tau_xx(isz,isx)+s
tau_zz(isz,isx)=tau_zz(isz,isx)+s
end subroutine add_source_p

!====================================================================================
subroutine add_source_n(tau_xx,tau_zz,nx,nz,nps,s) ! Add noise source by Lei Fu 2015-04-07
! nps : number source per shot
integer, intent(in)     :: nx,nz,nps
real,    intent(in)     :: s(:)
real,    intent(inout)  :: tau_xx(:,:),tau_zz(:,:)
integer                 :: k
real                    :: rand1, rand2
do k=1,nps
   call random_number(rand1)
   call random_number(rand2)
   isz = floor(rand1*nz)
   isx = floor(rand2*nx)
   tau_xx(isz,isx)=tau_xx(isz,isx)+s(k)
   tau_zz(isz,isx)=tau_zz(isz,isx)+s(k)
enddo
end subroutine add_source_n

!====================================================================================
subroutine add_source_s(s,u,w)
real,    intent(in)     :: s
real,    intent(inout)  :: u(:,:),w(:,:)

u(isz,isx)=u(isz,isx)+s
u(isz-1,isx)=u(isz-1,isx)-s
w(isz,isx)=w(isz,isx)+s
w(isz,isx+1)=w(isz,isx+1)+s
end subroutine add_source_s

!==============================================================================================
subroutine initial_variable(fdm,s)
type(fdmod2d),intent(inout)        :: fdm
type(shot2d),intent(in)         :: s
fd_type=fdm%f%fd_type;fd_order=fdm%f%fd_order; bc_type=fdm%f%bc_type
nt=fdm%f%nt; nz=fdm%m%nz; nx=fdm%m%nx; npml=fdm%f%npml
dt=fdm%f%dt; dx=fdm%m%dx
isFS=fdm%f%isFS; isDipole=fdm%f%isDipole; isWriteData=fdm%f%isWriteData;
isReadData=fdm%f%isReadData; isSaveBC=fdm%f%isSaveBC; isSaveWF=fdm%f%isSaveWF
source_type=fdm%f%source_type;data_type=fdm%f%data_type
dnt_mod=nint(s%dt/fdm%f%dt)
call i_v_kernal(fdm%m%vp,fdm%m%vs,fdm%m%lambda,fdm%m%mu,fdm%m%den,fdm%f%fs)
end subroutine initial_variable
!===========================================================================================
subroutine finalize_variable()
call deallocate_and_free(igx,igz,fs)
call deallocate_and_free(alpha_vel,alpha_lambda,alpha_mu,temp,beta_dt)
end subroutine finalize_variable
!===========================================================================================
subroutine i_v_kernal(vp,vs,lambda,mu,den,fs_in)
real,intent(inout)    :: vp(:,:),vs(:,:),lambda(:,:),mu(:,:),den(:,:)
integer,intent(in)    :: fs_in(:)
call nx_nz_pml
bc_len=(fd_order-20)/2    
fs_thick=bc_len-1
dtdx=dt/dx
call allocate_and_initial(beta_dt,nz_pml,nx_pml)
beta_dt=1   ! perhaps further construction is required
call allocate_and_initial(fs,nx_pml)
fs=fs_in
call allocate_and_initial(damp,alpha_vel,alpha_lambda,alpha_mu,kappa,temp,nz_pml,nx_pml) 
vmin=min(minval(vp),minval(vs))                                                 
call abc_get_damp2d(nx,nz,npml,dx,vmin,damp)
alpha_vel=dtdx/den                             
alpha_lambda=lambda*dtdx
alpha_mu=mu*dtdx                      
kappa=damp*dt                                
temp=1.0-kappa
!call snapshot("alpha_vel",alpha_vel,1,nz_pml,nx_pml)
!call snapshot("alpha_lambda",alpha_lambda,1,nz_pml,nx_pml)
!call snapshot("alpha_mu",alpha_mu,1,nz_pml,nx_pml)
call deallocate_and_free(damp,kappa)
end subroutine i_v_kernal         
!==========================================================================================
subroutine nx_nz_pml
nx_pml = nx + 2 * npml
nz_pml = nz + 2 * npml
nxp = nx + npml
nzp = nz + npml
end subroutine nx_nz_pml
!==========================================================================================
subroutine setup_source_receiver(s)
type(shot2d),       intent(in) :: s
call allocate_and_initial(igx,igz,s%ng)
! Setup source position
isx = npml+int(s%xs(1)/dx)+1
isz = npml+int(s%zs(1)/dx)+1
write(*,*) "isx=",isx,"isz=",isz
write(*,*) "s%xs(1)=",s%xs(1)/dx,"s%zs(1)=",s%zs(1)/dx
write(*,*) "npml=",npml
write(*,*) "dx=",dx
! Source position is shifted downward 1 node if on the FS.
if (isz == fs(isx)) isz = isz + 1
!write(*,*) "fs(isx)===",fs(isx)
! Setup receiver position
do ig=1,s%ng
   igx(ig) = npml + int(s%xg(ig)/dx)+1
   igz(ig) = npml + int(s%zg(ig)/dx)+1
   ! Geophone position is shifted downward 1 node if on the FS.
   if (igz(ig) == fs(igx(ig))) igz(ig) = igz(ig) + 1
enddo
end subroutine setup_source_receiver
!=========================================================================================
subroutine cal_lambda_mu(vp,vs,den,lambda,mu)
real,intent(in)             :: vp(:,:)
real,intent(in)             :: vs(:,:)
real,intent(in)             :: den(:,:)
real,intent(out)            :: lambda(:,:)
real,intent(out)            :: mu(:,:)

mu=(vs**2)*den
lambda=(vp**2)*den-2*mu

if (minval(lambda).lt.0.0) then
   write(*,*) "Vp/Vs is phyiscally impossible"
endif
end subroutine cal_lambda_mu
!=====================================================================================






end module module_e2d
