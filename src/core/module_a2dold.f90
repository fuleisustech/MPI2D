module module2d_mp_i_v_a2d

use module_global
use module_datatype
use module_array
use module_utility
use module_boundary2d
use module_fd2d
use module_string
use module_io

!use module_parser  !Added by Bowen Guo

implicit none
character(len=100),private :: fd_type
integer,private :: ix,iz,it,isx,isz,ig,nx,nz,npml,nxp,nzp,&
   nt,dnt_mod,nx_pml,nz_pml,ic,fd_order,bc_type,bc_len,fs_thick
real, private   :: dtdx,dx,dz,dt,vmin,c1t2
integer,private,allocatable :: igx(:),igz(:),fs(:)
logical,private :: isFS, isDipole, isSaveBC, isSaveWF, isWriteData, isReadData
real,private,target,allocatable :: p0(:,:),p1(:,:),p(:,:),&
   q0(:,:),q1(:,:),q(:,:)
real,private,pointer :: pp0(:,:),pp1(:,:),pp(:,:),pq0(:,:),pq1(:,:),pq(:,:)
real,private,allocatable :: u(:,:),w(:,:)                             ! modified by Bowen Guo ! used in staggered grid
real,private,allocatable :: uq(:,:),wq(:,:),ppq(:,:)                  ! modified by Bowen Guo ! used in staggered grid
real,private, allocatable :: beta_dt(:,:)
! only use in illumination
real,private,allocatable :: illum_source(:,:),illum_receiver(:,:)
! only use in Absorbing Boundary Condition 2nd Order
real,private,allocatable :: alpha(:,:),temp1(:,:),temp2(:,:),temp3(:,:),temp4(:,:)! modified by Bowen Guo  
! only use in Absorbing Boundary Condition SG
real,private,allocatable :: alpha1(:,:),alpha2(:,:),temp(:,:) 
! only use in Sponge BC
real,private,allocatable :: spgx(:),spgz(:),spcx1(:),spcx2(:),&
   spcz1(:),spcz2(:),beta_dtdx(:,:),u1(:,:)
real,   private, allocatable :: damp(:,:), kappa(:,:) ! only use in abc initialize
real,   private, allocatable :: p_top(:,:,:),p_bottom(:,:,:),p_left(:,:,:),&
   p_right(:,:,:),p_nt(:,:),p_nt_1(:,:), pwf(:,:,:) ! only use in rtm code
real,   private, allocatable :: p_1(:,:)
real,   private, allocatable :: u_top(:,:,:),u_bottom(:,:,:),u_left(:,:,:),& ! modified by Bowen Guo ! used in staggered rtm
   u_right(:,:,:)                                                            ! modified by Bowen Guo ! used in staggered rtm
real,   private, allocatable :: w_top(:,:,:),w_bottom(:,:,:),w_left(:,:,:),& ! modified by Bowen Guo ! used in staggered rtm
   w_right(:,:,:)                                                            ! modified by Bowen Guo ! used in staggered rtm      
real,   private, allocatable :: u_nt(:,:), w_nt(:,:),uwf(:,:,:),wwf(:,:,:)   ! modified by Bowen Guo ! used in staggered rtm
integer(kind=MPI_OFFSET_KIND) :: disp
integer :: icheck

interface initial_variables
   module procedure initial_variables_mod
   !module procedure initial_variables_adcig   ! Added by Bowen Guo
end interface initial_variables

interface a2d_mod
   module procedure modeling
end interface a2d_mod

interface a2d_mod_SG             ! modified by Bowen Guo
  module procedure modeling_SG   ! modified by Bowen Guo
end interface a2d_mod_SG         ! modified by Bowen Guo

interface a2d_bmod
   module procedure linear_modeling
end interface a2d_bmod

interface a2d_bmod_SG                  ! modified by Bowen Guo
   module procedure linear_modeling_SG ! modified by Bowen GUo
end interface a2d_bmod_SG              ! modified by Bowen Guo

interface a2d_mig
   module procedure rtm
   !module procedure rtm_adcig    ! rtm to form angle domain common image gathers
end interface a2d_mig

interface a2d_mig_SG                ! modified by Bowen Guo
   module procedure rtm_SG          ! modified by Bowen Guo
end interface a2d_mig_SG            ! modified by Bowen Guo

interface add_source
   module procedure add_source_2nd_mod
   module procedure add_source_2nd_mig
!   module procedure add_source_sg_mod   
!   module procedure add_source_sg_mig
end interface add_source

private initial_variables,finalize_variables,nx_nz_pml,i_v_kernal,&
   a2d_abc_kernal,add_source,subtract_source,store_data,&
   add_FS,pertubation,setup_source_receiver,saveBC2d,loadBC2d

contains

!=====================================================================

subroutine add_source_2nd_mod(p,beta_dt,s,isFS,isDipole)
logical, intent(in) :: isFS, isDipole
real,    intent(in) :: beta_dt(:,:),s
real,    intent(inout) :: p(:,:)
p(isz,isx)=p(isz,isx)+beta_dt(isz,isx)*s
if (.not.isFS) then
   if (isDipole) then
      p(isz-2,isx) = p(isz-2,isx) - beta_dt(isz-2,isx)*s 
   endif
endif
end subroutine add_source_2nd_mod

!=====================================================================
subroutine add_source_2nd_mig(p,beta_dt,s,isDipole)
logical, intent(in) :: isDipole
real,    intent(in) :: beta_dt(:,:),s
real,    intent(inout) :: p(:,:)
p(isz,isx) = p(isz,isx) + beta_dt(isz,isx)*s
if (isDipole) then
   p(isz-2,isx) = p(isz-2,isx) - beta_dt(isz-2,isx)*s
endif
end subroutine add_source_2nd_mig

!=====================================================================

subroutine subtract_source(p0,beta_dt,s,it,isDipole)
integer, intent(in) :: it
logical, intent(in) :: isDipole
real,    intent(in) :: beta_dt(:,:),s(:)
real,    intent(inout) :: p0(:,:)
if (fd_type.eq."2ND")  then                                    ! modified by Bowen Guo
   p0(isz,isx) = p0(isz,isx) - s(it+2) * beta_dt(isz,isx)
elseif (fd_type.eq."SG") then                              ! modified by Bowen Guo
   p0(isz,isx) = p0(isz,isx) - s(it+1) * beta_dt(isz,isx)  ! modified by Bowen Guo
endif
if (fd_type.eq."2ND")   then    ! modified by Bowen Guo
   if (isDipole) then
      p0(isz-2,isx) = p0(isz-2,isx) + s(it+2) * beta_dt(isz-2,isx)
   endif                        !  modified by Bowen Guo
elseif (fd_type.eq."SG") then   !  modified by Bowen Guo
   if (isDipole) then
      p0(isz-2,isx) = p0(isz-2,isx) + s(it+1) * beta_dt(isz-2,isx)
   endif
endif                          ! modified by Bowen Guo
end subroutine subtract_source

!=====================================================================

subroutine add_seis(q,beta_dt,ng,seis,it)
integer, intent(in) :: ng,it
real,    intent(in) :: beta_dt(:,:),seis(:,:)
real,    intent(inout) :: q(:,:)
do ig=1,ng
   q(igz(ig),igx(ig)) = q(igz(ig),igx(ig)) &
                        + beta_dt(igz(ig),igx(ig)) * seis(it,ig)
enddo
end subroutine add_seis

!=====================================================================

subroutine add_FS(p,npml,fs_thick)
integer, intent(in) :: npml,fs_thick
real,    intent(inout) :: p(:,:)
integer :: iz
p(npml+1,:)=0.0
do iz=1,fs_thick
   p(npml+1-iz,:)=-p(npml+1+iz,:)
enddo
end subroutine add_FS

!=====================================================================

subroutine store_data(p,ng,it,seis)
integer, intent(in) :: ng,it
real,    intent(in) :: p(:,:)
real,    intent(inout) :: seis(:,:)
integer             :: it_out
do ig=1,ng
   if (dnt_mod.ne.1) then
      if (mod(it,dnt_mod).eq.1) then
         it_out=(it-1)/dnt_mod+1 
         seis(it_out,ig) = p(igz(ig),igx(ig))
      endif
   else
      seis(it,ig) = p(igz(ig),igx(ig))
   endif
enddo
end subroutine store_data

!=====================================================================

subroutine pertubation(ic,refl,beta_dt,p,p1,q,nz_pml,nx_pml)
integer, intent(in) :: ic, nz_pml,nx_pml
real,    intent(in) :: refl(:,:),beta_dt(:,:),p(:,:),p1(:,:)
real,    intent(inout) :: q(:,:)
if (ic.eq.0) then
   forall (iz=1:nz_pml,ix=1:nx_pml)
      q(iz,ix)=q(iz,ix)+p1(iz,ix)*refl(iz,ix)*beta_dt(iz,ix)
   endforall
elseif (ic.eq.1) then
   forall (iz=1:nz_pml,ix=1:nx_pml)
      q(iz,ix)=q(iz,ix)+(p(iz,ix)-p1(iz,ix))*refl(iz,ix)*beta_dt(iz,ix)
   endforall
elseif (ic.eq.2) then
   forall (iz=2:nz_pml-1,ix=2:nx_pml-1)
      q(iz,ix)=q(iz,ix)+(p1(iz-1,ix)+p1(iz+1,ix)+p1(iz,ix-1)+p1(iz,ix+1)&
                        -4.0*p1(iz,ix))*refl(iz,ix)*beta_dt(iz,ix)
   endforall
endif
end subroutine pertubation
!=====================================================================

subroutine image_condition(ic,img,nzp,nxp,npml,p0,p1,q0)
integer, intent(in) :: ic,nzp,nxp,npml
real,    intent(in) :: p0(:,:),p1(:,:),q0(:,:)
real,    intent(inout) :: img(:,:)
if (ic.eq.0) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      img(iz-npml,ix-npml) = img(iz-npml,ix-npml) &
                                 + p1(iz,ix) * q0(iz,ix)
   endforall
elseif (ic.eq.1) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      img(iz-npml,ix-npml) = img(iz-npml,ix-npml) + &
            (p0(iz,ix)- p1(iz,ix)) * q0(iz,ix)
   endforall
elseif (ic.eq.2) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      img(iz-npml,ix-npml) = img(iz-npml,ix-npml) &
            + (p1(iz-1,ix) + p1(iz+1,ix) + p1(iz,ix-1) &
             + p1(iz,ix+1)-4.0 * p1(iz,ix)) * q0(iz,ix)
   endforall
endif
end subroutine image_condition

!=====================================================================
!subroutine modeling(fdm,bc,wf,s,output)
subroutine modeling(fdm,bc,wf,s,output)
type(fdmod2d),      intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(bc2d_general), intent(inout) :: bc
type(wf2d),         intent(inout) :: wf
character(len=*),   intent(in)    :: output
integer                           :: fid,it_out,mod_type

write(*,*) "check point 1111"


call initial_variables(fdm,s)

write(*,*) "check point 3333"
if (allocated(bc%p_top).or.allocated(wf%pwf)) then
   if (allocated(bc%p_top)) then
      mod_type=2
   else
      mod_type=3
   endif
else
   mod_type=1
endif
write(*,*) "check point 2222"

if (isWriteData) then
   call MPI_FILE_OPEN(MPI_COMM_SELF, output, &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, &
         MPI_INFO_NULL, fid, ierr)
else
   s%seis = 0.0
endif

call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
do it=1,nt
   if (bc_type==1) then
      if (fd_type.eq."2ND") then 
         call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
      elseif (fd_type.eq."Lax_wen") then 
         call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,temp3,temp4,fdm%m%v,dx,dz,nz_pml,nx_pml,fd_order)
      endif 
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif 
   call add_source(pp,beta_dt,fdm%s%sou(it),isFS,isDipole)
   if (isFS) then
      call add_FS(pp,npml,fs_thick)
   endif
   if (mod_type.eq.2) then
      call saveBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
               bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   elseif (mod_type.eq.3) then
      wf%pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   endif
  ! if (mod(it,50)==1) call snapshot("p_",p,it,nz_pml,nx_pml)
   if (isWriteData) then
      if (dnt_mod.ne.1) then
         if (mod(it,dnt_mod).eq.1) then
            it_out=(it-1)/dnt_mod+1
            do ig=1,s%ng
               disp=((ig-1)*s%nt+(it_out-1))*4
               call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
               call MPI_FILE_WRITE(fid, pp(igz(ig),igx(ig)),1, MPI_REAL, &
                         MPI_STATUS_IGNORE, ierr)
            enddo
         endif
      else
         do ig=1,s%ng
            disp=((ig-1)*s%nt+(it-1))*4
            call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
            call MPI_FILE_WRITE(fid, pp(igz(ig),igx(ig)),1, MPI_REAL, &
                      MPI_STATUS_IGNORE, ierr)
         enddo
      endif
   else
      call store_data(pp,s%ng,it,s%seis)
   endif
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
if (mod_type.eq.2) then
   call copy_array(pp0,bc%p_nt_1,p1,bc%p_nt)
endif
if (isWriteData) then
   call MPI_FILE_CLOSE(fid,ierr)
endif
call finalize_variables()
call deallocate_and_free(p,p0,p1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine modeling

!=====================================================================
! modified by Bowen Guo
subroutine modeling_SG(fdm,bc,wf,s,output)
type(fdmod2d),          intent(inout) :: fdm
type(shot2d),           intent(inout) :: s
type(bc2d_general), intent(inout)     :: bc
type(wf2d),             intent(inout) :: wf
character(len=*),       intent(in)    :: output
integer                               :: fid,it_out,mod_type

call initial_variables(fdm,s)
if (allocated(bc%p_top).or.allocated(wf%pwf)) then
   if (allocated(bc%p_top)) then
      mod_type=2
   else
      mod_type=3
   endif
else
    mod_type=1
endif
if (isWriteData) then
   call MPI_FILE_OPEN(MPI_COMM_SELF,output,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                      MPI_INFO_NULL,fid,ierr)
else
  s%seis=0.0
endif

call setup_source_receiver(s)
call allocate_and_initial(p,u,w,nz_pml,nx_pml)
do it=1,nt
   ! update pressure field
   call a2d_abc_kernal_SG_p(u,w,p,temp,alpha1,nz_pml,nx_pml,fd_order)
   call add_source(p,beta_dt,fdm%s%sou(it),isFS,isDipole)
   if (isFS) then
      call add_FS(p,npml,fs_thick)
   endif
   ! update velocity field
   call a2d_abc_kernal_SG_vel(u,w,p,temp,alpha2,nz_pml,nx_pml,fd_order)
   if (mod_type.eq.2) then
      ! save pressure boundary
      call saveBC2d(p,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
           bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
      ! save particle velocity boundary
      call saveBC2d(u,nz,nx,nzp,nxp,npml,bc%u_top(:,:,it),bc%u_bottom(:,:,it),&
               bc%u_left(:,:,it),bc%u_right(:,:,it),bc_len)
      call saveBC2d(w,nz,nx,nzp,nxp,npml,bc%w_top(:,:,it),bc%w_bottom(:,:,it),&
               bc%w_left(:,:,it),bc%w_right(:,:,it),bc_len)
   elseif (mod_type.eq.3) then
      wf%pwf(:,:,it)=p(npml:nzp+1,npml:nxp+1)
      wf%uwf(:,:,it)=u(npml:nzp+1,npml:nxp+1)
      wf%wwf(:,:,it)=w(npml:nzp+1,npml:nxp+1)
   endif
   ! if (mod(it,1000)==1) call_snapshot("p_:",p,it,nz_pml,nx_pml)
   if (isWriteData) then  
      if (dnt_mod.ne.1) then
         it_out=(it-1)/dnt_mod+1
         do ig=1,s%ng
            disp=((ig-1)*s%nt+(it_out-1))*4
            call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
            call MPI_FILE_WRITE(fid,p(igz(ig),igx(ig)),1,MPI_REAL, &
                               MPI_STATUS_IGNORE,ierr)
         enddo
       endif
   else
       call store_data(p,s%ng,it,s%seis)
   endif
enddo       ! End of time_marching loop
if (mod_type.eq.2) then
   call copy_array(p,bc%p_nt,u,bc%u_nt,w,bc%w_nt)
endif
if (isWriteData) then
   call MPI_FILE_CLOSE(fid,ierr)
endif
call finalize_variables()
call deallocate_and_free(p,u,w)
end subroutine modeling_SG
!====================================================================

subroutine check_time_now()
call date_and_time(date,time_now)
write(*,*)icheck,time_now
!call flush(6)
icheck=icheck+1
end subroutine check_time_now

!=====================================================================

subroutine linear_modeling(fdm,bc,wf,s,output)
type(fdmod2d),       intent(inout) :: fdm
type(shot2d),        intent(inout) :: s
type(bc2d_general),  intent(inout) :: bc
type(wf2d),          intent(inout) :: wf
character(len=*),    intent(in)    :: output
integer                            :: fid,it_out,mod_type
call initial_variables(fdm,s)
if (allocated(bc%p_top).or.allocated(wf%pwf)) then
   if (allocated(bc%p_top)) then
      mod_type=2
   else
      mod_type=3
   endif
else
   mod_type=1
endif
if (isWriteData) then
   call MPI_FILE_OPEN(MPI_COMM_SELF, output, &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, &
         MPI_INFO_NULL, fid, ierr)
else
   s%seis = 0.0
endif
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
do it=1,nt
   !--------- Pressure field --------------------------------------------
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_source(pp,beta_dt,fdm%s%sou(it),isFS,isDipole)
   if (mod_type.eq.2) then
      call saveBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
               bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   elseif (mod_type.eq.3) then
      wf%pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   endif
   !--------- Pertubation field -----------------------------------------
   if (bc_type==1) then
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pq,pq0,pq1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   if (isWriteData) then
      if (dnt_mod.ne.1) then
         if (mod(it,dnt_mod).eq.1) then
            it_out=(it-1)/dnt_mod+1
            do ig=1,s%ng
               disp=((ig-1)*s%nt+(it_out-1))*4
               call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
               call MPI_FILE_WRITE(fid, pq(igz(ig),igx(ig)),1, MPI_REAL, &
                         MPI_STATUS_IGNORE, ierr)
            enddo
         endif
      else
         do ig=1,s%ng
            disp=((ig-1)*s%nt+(it-1))*4
            call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
            call MPI_FILE_WRITE(fid, pq(igz(ig),igx(ig)),1, MPI_REAL, &
                      MPI_STATUS_IGNORE, ierr)
         enddo
      endif
   else
      call store_data(pq,s%ng,it,s%seis)
   endif
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
if (mod_type.eq.2) then
   call copy_array(pp0,bc%p_nt_1,p1,bc%p_nt)
endif
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
if (bc_type==2) call deallocate_and_free(u1)
if (isWriteData) then
   call MPI_FILE_CLOSE(fid,ierr)
endif
end subroutine linear_modeling

!=====================================================================
! modified by Bowen Guo
subroutine linear_modeling_SG(fdm,bc,wf,s,output)
type(fdmod2d),     intent(inout)   :: fdm
type(shot2d),      intent(inout)   :: s
type(bc2d_general),intent(inout)   :: bc
type(wf2d),        intent(inout)   :: wf
character(len=*),  intent(inout)   :: output
integer                            :: fid,it_out,mod_type
call initial_variables(fdm,s)
if (allocated(bc%p_top).or.allocated(wf%pwf)) then
   if (allocated(bc%p_top)) then
      mod_type=2
   else
      mod_type=3
   endif
else 
   mod_type=1
endif
if (isWriteData) then
   call MPI_FILE_OPEN(MPI_COMM_SELF,output,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                      MPI_INFO_NULL,fid,ierr)
else 
   s%seis=0.0
endif
call setup_source_receiver(s)
call allocate_and_initial(u,w,p,uq,wq,ppq,p_1,nz_pml,nx_pml)
do it=1,nt
   ! update pressure field
   p_1=p     !  save the wavefield before
   call a2d_abc_kernal_SG_p(u,w,p,temp,alpha1,nz_pml,nx_pml,fd_order)
   call add_source(p,beta_dt,fdm%s%sou(it),isFS,isDipole)
   if (isFS) then
      call add_FS(p,npml,fs_thick)
   endif
   ! update velocity field
   call a2d_abc_kernal_SG_vel(u,w,p,temp,alpha2,nz_pml,nx_pml,fd_order)
   if (mod_type.eq.2) then
      ! save pressure boundary
      call saveBC2d(p,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it), &                   bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
      ! save partile velocity field
      call saveBC2d(u,nz,nx,nzp,nxp,npml,bc%u_top(:,:,it),bc%u_bottom(:,:,it), &
                   bc%u_left(:,:,it),bc%u_right(:,:,it),bc_len)
      call saveBC2d(w,nz,nx,nzp,nxp,npml,bc%w_top(:,:,it),bc%w_bottom(:,:,it), &                   bc%w_left(:,:,it),bc%w_right(:,:,it),bc_len)
   elseif (mod_type.eq.3) then
      wf%pwf(:,:,it)=p(npml:nzp+1,npml:nxp+1)
   endif
   !------------------Pertubation field --------------------------------
   call a2d_abc_kernal_SG_p(uq,wq,ppq,temp,alpha1,nz_pml,nx_pml,fd_order)
   call a2d_abc_kernal_SG_vel(uq,wq,ppq,temp,alpha2,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,p_1,p,ppq,nz_pml,nx_pml)
   if (isWriteData) then
      if (dnt_mod.ne.1)  then
         if (mod(it,dnt_mod).eq.1) then
            it_out=(it-1)/dnt_mod+1
            do ig=1,s%ng
               disp=((ig-1)*s%nt+(it_out-1))*4
               call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
               call MPI_FILE_WRITE(fid,ppq(igz(ig),igx(ig)),1,MPI_REAL, &
                                  MPI_STATUS_IGNORE,ierr)
            enddo
         endif
      else 
         do ig=1,s%ng
            disp=((ig-1)*s%nt+(it-1))*4
            call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
            call MPI_FILE_WRITE(fid,ppq(igz(ig),igx(ig)),1,MPI_REAL,& 
                               MPI_STATUS_IGNORE,ierr)
         enddo
      endif
   else
      call store_data(ppq,s%ng,it,s%seis)
   endif
enddo ! End of time-marching loop
if (mod_type.eq.2) then
   c2d_mp_i_vall copy_array(p_1,bc%p_nt_1,p1,bc%p_nt)
endif
call finalize_variables()
call deallocate_and_free(u,w,p,uq,wq,ppq)
if (isWriteData) then
   call MPI_FILE_CLOSE(fid,ierr)
endif
end subroutine linear_modeling_SG 










   













 














!=====================================================================
subroutine rtm(s,bc,wf,fdm,img,output)
type(fdmod2d),      intent(inout)  :: fdm
type(shot2d),       intent(inout)  :: s
type(bc2d_general), intent(inout)  :: bc
type(wf2d),         intent(inout)  :: wf
type(image2d),      intent(inout)  :: img
character(len=*),   intent(in)     :: output
integer                            :: fid,rtm_type
call initial_variables(fdm,s)
call setup_source_receiver(s)
if (isSaveBC) then
   if (allocated(bc%p_top)) then
      rtm_type=1 ! read BC from outside
   else
      rtm_type=2 ! generate BC inside 
   endif
elseif (isSaveWF) then
   if (allocated(wf%pwf)) then
      rtm_type=3 ! read WF from outside
   else
      rtm_type=4 ! generate WF inside
   endif
endif
img%img=0.0
if (img%isIllum) then
   img%img=0.0
endif
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
if (isReadData) then
   call MPI_FILE_OPEN(MPI_COMM_SELF, output, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
endif
if (rtm_type==2) then
   call allocate_and_initial(p_top,p_bottom,bc_len,nx,nt)
   call allocate_and_initial(p_left,p_right,nz,bc_len,nt)
   call allocate_and_initial(p_nt,p_nt_1,nz_pml,nx_pml)
elseif (rtm_type==4) then
   call allocate_and_initial(pwf,nz+2,nx+2,nt)
endif

! Forward propagate to save Forward field
if ((rtm_type.eq.2) .or. (rtm_type.eq.4)) then
   do it=1,nt
      !if (rank.eq.0) write(*,*)"it1",it; call flush(6)
      if (bc_type==1) then
         call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
      else
         call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
                 spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
      endif
      call add_source(pp,beta_dt,fdm%s%sou(it),isDipole)
      if (rtm_type.eq.2) then
         call saveBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
                  p_left(:,:,it),p_right(:,:,it),bc_len)
      else 
         pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
      endif
      call wf_refresh(pp0,pp1,pp)
   enddo  ! End of time-marching loop
   if (rtm_type.eq.2) then
      call copy_array(pp1,p_nt,pp0,p_nt_1)
   else
      pp0(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt)
      pp1(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt-1)
   endif
endif
   
if (rtm_type.eq.1) then
   call copy_array(bc%p_nt,pp0,bc%p_nt_1,pp1)
elseif (rtm_type.eq.2) then
   call copy_array(p_nt,pp0,p_nt_1,pp1)
elseif (rtm_type.eq.3) then
   pp1(npml:nzp+1,npml:nxp+1)=wf%pwf(:,:,nt)
   pp0(npml:nzp+1,npml:nxp+1)=wf%pwf(:,:,nt-1)
else
   pp1(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt)
   pp0(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt-1)
endif

! Back propagate 
do it=nt-2,1,-1
   !if (rank.eq.0) write(*,*)"it2",it; call flush(6)
   if ((rtm_type.eq.3).or.(rtm_type.eq.4)) then
      if (rtm_type.eq.4) then
         pp(npml:nzp+1,npml:nxp+1)=pwf(:,:,it)
      else
         pp(npml:nzp+1,npml:nxp+1)=wf%pwf(:,:,it)
      endif
   else
      !-.-.-.-.-.-. Reconstruct Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_source(pp0,beta_dt,fdm%s%sou,it,isDipole)
      if (bc_type==1) then
         call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
      else
         call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,npml,nz_pml,nx_pml,fd_order)
      endif
      if (rtm_type.eq.2) then
         call loadBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
                  p_left(:,:,it),p_right(:,:,it),bc_len)
      else
         call loadBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
                  bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
      endif
   endif
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   if (bc_type==1) then
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pq,pq0,pq1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   if (isReadData) then
      do ig=1,s%ng
         disp=((ig-1)*s%nt+(it-1))*4
         call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
         call MPI_FILE_READ(fid,pq(igz(ig),igx(ig)),1, MPI_REAL, &
                  MPI_STATUS_IGNORE, ierr)
      enddo
   else
      call add_seis(pq,beta_dt,s%ng,s%seis,it)
   endif
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   !if (mod(it,1000)==1) call snapshot("p_",pp,it,nz_pml,nx_pml,npml)
   !if (mod(it,1000)==1) call snapshot("q_",pq,it,nz_pml,nx_pml,npml)
   if (img%isIllum) then
      forall (iz=npml+1:nzp,ix=npml+1:nxp)
         img%illum(iz-npml,ix-npml) = img%illum(iz-npml,ix-npml) &
                                       + (abs(pp(iz,ix)))**img%illum_order
      end forall
   endif
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
if (bc_type==2) call deallocate_and_free(u1)
if (rtm_type.eq.2) then
   call deallocate_and_free(p_top,p_bottom,p_left,p_right)
   call deallocate_and_free(p_nt,p_nt_1)
elseif (rtm_type.eq.4) then
   call deallocate_and_free(pwf)
endif
if (isReadData) then
   call MPI_FILE_CLOSE(fid,ierr)
endif
end subroutine rtm

!=====================================================================
! modified by Bowen Guo
subroutine rtm_SG(s,bc,wf,fdm,img,output)
type(fdmod2d), intent(inout)       :: fdm
type(shot2d),  intent(inout)       :: s
type(bc2d_general), intent(inout)  :: bc
type(wf2d),         intent(inout)  :: wf
type(image2d),      intent(inout)  :: img
character(len=*),   intent(in)     :: output
integer                            :: fid, rtm_type
call initial_variables(fdm,s)
call setup_source_receiver(s)
if (isSaveBC) then 
   if ((allocated(bc%p_top)).or.(allocated(bc%u_top))) then
      rtm_type=1   ! read BC from outside
   else 
      rtm_type=2   ! generate BC inside
   endif 
elseif (isSaveWF) then
   if (allocated(wf%pwf)) then
      rtm_type=3   ! read WF from outside
   else 
      rtm_type=4   ! generate WF inside
   endif
endif  
img%img=0.0
if (img%isIllum) then
   img%illum=0.0
endif
call allocate_and_initial(u,w,p,uq,wq,ppq,p_1,nz_pml,nx_pml)
!uq,wq,pq,p0 are used in the backprogagated wave field
if (isReadData) then
   call MPI_FILE_OPEN(MPI_COMM_SELF,output,MPI_MODE_RDONLY,MPI_INFO_NULL,fid,ierr)
endif
if (rtm_type==2) then
   call allocate_and_initial(p_top,p_bottom,bc_len,nx,nt)
   call allocate_and_initial(p_left,p_right,nz,bc_len,nt)
   call allocate_and_initial(u_top,u_bottom,bc_len,nx,nt)
   call allocate_and_initial(u_left,u_right,nz,bc_len,nt)
   call allocate_and_initial(w_top,w_bottom,bc_len,nx,nt)
   call allocate_and_initial(w_left,w_right,nz,bc_len,nt)

   call allocate_and_initial(p_nt,u_nt,w_nt,nz_pml,nx_pml)
elseif (rtm_type==4) then 
   call allocate_and_initial(pwf,nz+2,nx+2,nt)
   call allocate_and_initial(uwf,nz+2,nx+2,nt)
   call allocate_and_initial(wwf,nz+2,nx+2,nt)
endif

! Forward propagate to save forward field
if (( rtm_type.eq.2) .or. (rtm_type.eq.4)) then ! no modeling has been done
   do it=1,nt
      ! if (rank.eq.0) write(*,*) "it1",it; call flush(6)
      if (bc_type.eq.1) then 
         ! update pressure field 
         call a2d_abc_kernal_SG_p(u,w,p,temp,alpha1,nz_pml,nx_pml,fd_order)
      endif
      call add_source(p,beta_dt,fdm%s%sou(it),isDipole)
      if (bc_type.eq.1) then 
         ! update velocity field
         call a2d_abc_kernal_SG_vel(u,w,p,temp,alpha2,nz_pml,nx_pml,fd_order)
      endif 
      if (rtm_type.eq.2) then 
         call saveBC2d(p,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it), &
                      p_left(:,:,it),p_right(:,:,it),bc_len)
         call saveBC2d(u,nz,nx,nzp,nxp,npml,u_top(:,:,it),u_bottom(:,:,it), &
                      u_left(:,:,it),u_right(:,:,it),bc_len)
         call saveBC2d(w,nz,nx,nzp,nxp,npml,w_top(:,:,it),w_bottom(:,:,it), &
                      w_left(:,:,it),w_right(:,:,it),bc_len)
      ! only u,w  boundary is saved and needed to be loaded in reconstruction 
      else 
         pwf(:,:,it)=p(npml:nzp+1,npml:nxp+1)
         uwf(:,:,it)=u(npml:nzp+1,npml:nxp+1)
         wwf(:,:,it)=w(npml:nzp+1,npml:nxp+1)
      endif
   enddo  ! End of time_marching loop
   if (rtm_type.eq.2) then
      call copy_array(p,p_nt,u,u_nt,w,w_nt) ! save last wavefields
   endif
endif
 
! load last wavefield as the first wavefield during reconstruction
if (rtm_type.eq.1) then
   call copy_array(bc%p_nt,p,bc%u_nt,u,bc%w_nt,w)
elseif (rtm_type.eq.2) then
   call copy_array(p_nt,p,u_nt,u,w_nt,w)
elseif (rtm_type.eq.3) then
   p(npml:nzp+1,npml:nxp+1)=wf%pwf(:,:,nt)
   u(npml:nzp+1,npml:nxp+1)=wf%uwf(:,:,nt)
   w(npml:nzp+1,npml:nxp+1)=wf%wwf(:,:,nt)
else 
   p(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt)
   u(npml:nzp+1,npml:nxp+1)=uwf(:,:,nt)
   w(npml:nzp+1,npml:nxp+1)=wwf(:,:,nt)
endif

! back progagate
do it=nt-1,1,-1
   ! if (rank .eq.0) write(*,*) "it2",it; call flush(6) 
   if ((rtm_type.eq.3) .or. (rtm_type.eq.4)) then ! the wavefield already exists, no reconstruction
      if (rtm_type.eq.4) then
         p(npml:nzp+1,npml:nxp+1)=pwf(:,:,it)
      else 
         p(npml:nzp+1,npml:nxp+1)=wf%pwf(:,:,it)
      endif
   else 
   !=========================Reconstruct pressure field=========================
      if (ic.eq.1)   p_1=p   ! note the p wavefield before 
   
   ! update pressure field
      if (bc_type.eq.1) then
         call a2d_abc_kernal_SG_p(u,w,p,temp,alpha1,nz_pml,nx_pml,fd_order)
      endif
   
   ! load pressure field
      if (rtm_type.eq.2) then
         call loadBC2d(p,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),p_left(:,:,it)&
                      ,p_right(:,:,it),bc_len)
      else
         call loadBC2d(p,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),bc%p_left&
                   (:,:,it),bc%p_right(:,:,it),bc_len)
      endif
   ! substract source from the p
      call subtract_source(p,beta_dt,fdm%s%sou,it,isDipole)
   ! update velocity field
      if (bc_type.eq.1) then
         call a2d_abc_kernal_SG_vel(u,w,p,temp,alpha2,nz_pml,nx_pml,fd_order)
      endif  
   ! load velocity boundary
   !   if (rtm_type.eq.2) then
   !      call loadBC2d(u,nz,nx,nzp,nxp,npml,u_top(:,:,it),u_bottom(:,:,it),u_left(:,:,it)&
   !                   ,u_right(:,:,it),bc_len)
   !      call loadBC2d(w,nz,nx,nzp,nxp,npml,w_top(:,:,it),w_bottom(:,:,it),w_left(:,:,it)&
   !                   ,w_right(:,:,it),bc_len)
   !   else
   !      call loadBC2d(u,nz,nx,nzp,nxp,npml,bc%u_top(:,:,it),bc%u_bottom(:,:,it),bc%u_left&
   !                (:,:,it),bc%u_right(:,:,it),bc_len)
   !      call loadBC2d(w,nz,nx,nzp,nxp,npml,bc%w_top(:,:,it),bc%w_bottom(:,:,it),bc%w_left&
   !           (:,:,it),bc%w_right(:,:,it),bc_len)
   !   endif

   endif  
   !========================Back propagate Receiver Pressure field=================
   
   if (bc_type.eq.1) then
      call a2d_abc_kernal_SG_p(uq,wq,ppq,temp,alpha1,nz_pml,nx_pml,fd_order)      
   endif
   if (isReadData) then
      do ig=1,s%ng
         disp=((ig-1)*s%nt+(it-1))*4
         call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
         call MPI_FILE_READ(fid,ppq(igz(ig),igx(ig)),1,MPI_REAL,MPI_STATUS_IGNORE,ierr)
      enddo
   else
      call add_seis(ppq,beta_dt,s%ng,s%seis,it)
   endif
   if (bc_type.eq.1) then
      call a2d_abc_kernal_SG_vel(uq,wq,ppq,temp,alpha2,nz_pml,nx_pml,fd_order)
   endif
   call image_condition(ic,img%img,nzp,nxp,npml,p_1,p,ppq)
  ! if (mod(it,50).eq.1) call snapshot("p_",p,it,nz_pml,nx_pml,npml)
  ! if (mod(it,50).eq.1) call snapshot("uq_",uq,it,nz_pml,nx_pml,npml)
   if (img%isILLum) then
      forall (iz=npml+1:nzp,ix=npml+1:nxp)
         img%illum(iz-npml,ix-npml)=img%illum(iz-npml,ix-npml)+(abs(pp(iz,ix)))**&
                                    img%illum_order
      endforall
   endif
enddo
call finalize_variables()
call deallocate_and_free(u,w,p,uq,wq,ppq)
if (rtm_type.eq.2) then
   call deallocate_and_free(p_top,p_bottom,p_left,p_right)
   call deallocate_and_free(u_top,u_bottom,u_left,u_right)
   call deallocate_and_free(w_top,w_bottom,w_left,w_right)
   
   call deallocate_and_free(p_nt,u_nt,w_nt)
elseif (rtm_type.eq.4) then
   call deallocate_and_free(pwf,uwf,wwf)
endif
if (isReadData) then
   call MPI_FILE_CLOSE(fid,ierr)
endif
end subroutine rtm_SG   
!==================================================================================
! Use rtm to calculate angle domain common image gathers
! subroutine rtm_adcig(s,bc,wf,fdm,img,output)
! type(fdmod2d),intent(inout)     :: fdm
! type(image2d),intent(inout)     :: img
! type(shot2d),intent(inout)      :: s
! type(bc_general),intent(inout)  :: bc
! type(wf2d), intent(inout)       :: wf
! character(len=*), intent(inout) :: output
! integer                         :: fid,it_out,mod_type
 
! call initial_variables(fdm,s,img)















!=====================================================================================
subroutine setup_source_receiver(s)
type(shot2d),       intent(in) :: s
call allocate_and_initial(igx,igz,s%ng)
! Setup source position
isx = npml+int(s%xs(1)/dx)+1
isz = npml+int(s%zs(1)/dx)+1
!write(*,*) "isx=",isx,"isz=",isz
!write(*,*) "s%xs(1)=",s%xs(1)/dx,"s%zs(1)=",s%zs(1)/dx
!write(*,*) "npml=",npml
!write(*,*) "dx=",dx
! Source position is shifted downward 1 node if on the FS.
if (isz == fs(isx)) isz = isz + 1
! Setup receiver position
do ig=1,s%ng
   igx(ig) = npml + int(s%xg(ig)/dx)+1
   igz(ig) = npml + int(s%zg(ig)/dx)+1
   ! Geophone position is shifted downward 1 node if on the FS.
   if (igz(ig) == fs(igx(ig))) igz(ig) = igz(ig) + 1
enddo
end subroutine setup_source_receiver

!=====================================================================

subroutine initial_variables_mod(fdm,s)
type(fdmod2d),intent(in) :: fdm
type(shot2d), intent(in) :: s
write(*,*) "101010101010"
fd_type=fdm%f%fd_type; fd_order=fdm%f%fd_order; bc_type=fdm%f%bc_type
nt=fdm%f%nt; nz=fdm%m%nz; nx=fdm%m%nx; npml=fdm%f%npml
dt=fdm%f%dt; dx=fdm%m%dx; dz=fdm%m%dz;
isFS=fdm%f%isFS; isDipole=fdm%f%isDipole; isWriteData=fdm%f%isWriteData; 
isReadData=fdm%f%isReadData; isSaveBC=fdm%f%isSaveBC; isSaveWF=fdm%f%isSaveWF
dnt_mod=nint(s%dt/fdm%f%dt)
write(*,*) "202020202020"
!den=fdm%m%den       ! modified by Bowen Guo
call i_v_kernal(fdm%m%v,fdm%m%den,fdm%f%fs)


end subroutine initial_variables_mod

!=====================================================================

subroutine finalize_variables()
if (fd_type.eq."2ND") then               ! modified by Bowen Guo
   if (bc_type==1) then
      call deallocate_and_free(igx,igz,fs)
      call deallocate_and_free(alpha,temp1,temp2,beta_dt)
   endif
   if (bc_type==2) then 
      call deallocate_and_free(igx,igz,fs)
      call deallocate_and_free(beta_dtdx,beta_dt,u1)
   endif
elseif (fd_type.eq."SG") then              ! modified by Bowen Guo
   call deallocate_and_free(igx,igz,fs)    ! modified by Bowen Guo
   call deallocate_and_free(alpha1,alpha2,temp,beta_dt) ! modified by Bowen Guo
endif  
end subroutine finalize_variables

!=====================================================================

subroutine nx_nz_pml
nx_pml = nx + 2 * npml
nz_pml = nz + 2 * npml
nxp = nx + npml
nzp = nz + npml
end subroutine nx_nz_pml

!=====================================================================

subroutine i_v_kernal(v,den,fs_in)
real,intent(in) :: v(:,:),den(:,:)
integer, intent(in) :: fs_in(:)
call nx_nz_pml
if ((fd_type.eq."2ND").or.(fd_type.eq."Lax_wen")) then             ! modified by Bowen Guo
   bc_len=(fd_order-20)/2+1
elseif (fd_type.eq."SG") then         ! modified by Bowen Guo
   bc_len=(fd_order-20)/2              ! modified by Bowen Guo
endif                                  ! modified by Bowen Guo

write(*,*)"66666666666"
fs_thick=bc_len-1
dtdx = dt/dx

write(*,*) "777777777777777777"
call allocate_and_initial(beta_dt,nz_pml,nx_pml)
if ((fd_type.eq."2ND").or.(fd_type.eq."Lax_wen"))  then
   beta_dt=(v*dt)**2.0
!write(*,*) 'dt = ',dt*v(20,60)
elseif (fd_type.eq."SG") then
   beta_dt=(v**2.0)*dt
endif
write(*,*) "888888888888888888"
call allocate_and_initial(fs,nx_pml)
fs = fs_in
if (fd_type.eq."2ND") then ! 2nd Order FD
   if (fd_order==22) then
      c1t2=C21*2.0
   elseif (fd_order==24) then
      c1t2=C41*2.0
   elseif (fd_order==26) then
      c1t2=C61*2.0
   elseif (fd_order==28) then
      c1t2=C81*2.0
   endif
write(*,*) "9999999999999999999999"
   if (bc_type==1) then ! Absorbing Boundary Condition
      call allocate_and_initial(damp,alpha,kappa,temp1,temp2,nz_pml,nx_pml)
      vmin=minval(v)
      call abc_get_damp2d(nx,nz,npml,dx,vmin,damp)
      alpha=(v*dtdx)**2.0
      kappa=damp*dt
!      call write_binfile('kappa.bin',kappa,nz+2*npml,nx+2*npml)
      temp1=2.0+c1t2*alpha-kappa
      temp2=1.0-kappa
      call deallocate_and_free(damp,kappa)
   elseif (bc_type==2) then ! Sponge Boundary Condition
      call allocate_and_initial(beta_dtdx,u1,nz_pml,nx_pml)
      call allocate_and_initial(spgx,spcz1,spcz2,nx_pml)
      call allocate_and_initial(spgz,spcx1,spcx2,nz_pml)
      call spgcoef(spgx,spgz,spcx1,spcx2,spcz1,spcz2,v,&
                   nx_pml,nz_pml,npml,dx,dx,dt)
      beta_dtdx=(v*dtdx)**2.0
!      call write_binfile("beta_dtdx.bin",beta_dtdx,nz_pml,nx_pml)
!      call write_binfile("vel.bin",v,nz_pml,nx_pml)
   endif
elseif (fd_type.eq."SG")  then! Staggerred FD 
   if (bc_type==1) then ! Absorbing Boundary Condition                       ! modified by Bowen Guo
      call allocate_and_initial(damp,alpha1,alpha2,kappa,temp,nz_pml,nx_pml) ! modified by Bowen Guo
      vmin=minval(v)                                                         ! modified by Biwen Guo
      call abc_get_damp2d(nx,nz,npml,dx,vmin,damp) ! modified by Bowen Guo
      alpha1=dtdx*(v**2)*den                       ! modified by Bowen Guo
      alpha2=dtdx/den                              ! modified by Bowen Guo
      kappa=damp*dt                               ! modified by Bowen Guo
      temp=1.0-kappa                              ! modified by Bowen Guo
!      if (rank.eq.0) call write_binfile('temp_sg.bin',temp,nz_pml,nx_pml)
      call deallocate_and_free(damp,kappa)         ! modified by Bowen Guo
   elseif (bc_type==2) then ! Sponge Boundary Condition
      write(*,*)"Under construction!"
   endif
elseif (fd_type.eq."Lax_wen") then                                        ! modified by Bowen Guo
   call allocate_and_initial(damp,kappa,temp,nz_pml,nx_pml) ! modified by Bowen Guo
   call allocate_and_initial(temp1,temp2,temp3,temp4,nz_pml,nx_pml) ! modified by Bowen Guo
   vmin=minval(v)
   call abc_get_damp2d(nx,nz,npml,dx,vmin,damp)
   kappa=damp*dt
   !temp=1.0
   !temp1=temp*(dt**2) 
   temp1=dt**2
   temp2=kappa-1
   temp3=2-kappa
   temp4=-(v**2)*(dt**4)/12
endif 
end subroutine i_v_kernal

end module module_a2d
