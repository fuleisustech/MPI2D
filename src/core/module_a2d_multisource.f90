module module_a2d_multisource

use module_global
use module_datatype
use module_array
use module_utility
use module_math
use module_boundary2d
use module_fd2d
use module_string
use module_io

implicit none
character(len=100) :: fd_type
integer,private :: ix,iz,it,nx,nz,npml,nxp,nzp,ns,ng,&
   nt,dnt_mod,nx_pml,nz_pml,ic,fd_order,bc_type,bc_len,fs_thick
real, private   :: dtdx,dx,dz,dt,vmin,c1t2
integer,private,allocatable :: isx(:),isz(:),igx(:),igz(:),fs(:),idelay(:)
logical,private :: isFS, isDipole, isSaveBC, isSaveWF
real,private,target,allocatable :: p0(:,:),p1(:,:),p(:,:),&
   q0(:,:),q1(:,:),q(:,:)
real,private,pointer :: pp0(:,:),pp1(:,:),pp(:,:),pq0(:,:),pq1(:,:),pq(:,:)
real,private, allocatable :: beta_dt(:,:), polarity(:), taper(:)
! only use in illumination
real,private,allocatable :: illum_source(:,:),illum_receiver(:,:)
! only use in Absorbing Boundary Condition 
real,private,allocatable :: alpha(:,:),temp1(:,:),temp2(:,:) 
! only use in Sponge BC
real,private,allocatable :: spgx(:),spgz(:),spcx1(:),spcx2(:),&
   spcz1(:),spcz2(:),beta_dtdx(:,:),u1(:,:)
real,   private, allocatable :: damp(:,:), kappa(:,:) ! only use in abc initialize
real,   private, allocatable :: p_top(:,:,:),p_bottom(:,:,:),p_left(:,:,:),&
   p_right(:,:,:),p_nt(:,:),p_nt_1(:,:), pwf(:,:,:) ! only use in rtm code
integer :: icheck

interface initial_variables
   module procedure initial_variables_mod
end interface initial_variables

interface a2d_ms_mod
   module procedure modeling
   module procedure modeling_saveBC
   module procedure modeling_saveWF
end interface a2d_ms_mod

interface a2d_ms_bmod
   module procedure linear_modeling
   module procedure linear_modeling_saveBC
   module procedure linear_modeling_saveWF
end interface a2d_ms_bmod

interface a2d_ms_mig
   module procedure rtm
   module procedure rtm_saveBC
   module procedure rtm_saveWF
end interface a2d_ms_mig

interface a2d_ms_illum
   module procedure illum
end interface a2d_ms_illum

interface add_source
   module procedure add_source_mod
   module procedure add_source_mig
end interface add_source

private initial_variables,finalize_variables,nx_nz_pml,i_v_kernal,&
   a2d_abc_kernal,add_source,subtract_source,store_data,&
   add_FS,pertubation,setup_source_receiver,saveBC2d,loadBC2d

contains

!=====================================================================

subroutine add_source_mod(p,beta_dt,s,it,isFS,isDipole)
logical, intent(in) :: isFS, isDipole
real,    intent(in) :: beta_dt(:,:),s(:)
integer, intent(in) :: it
real,    intent(inout) :: p(:,:)
integer :: is, it_delay
do is=1,ns
   it_delay=it-idelay(is)
   if (it_delay > 0) then 
      p(isz(is),isx(is)) = p(isz(is),isx(is)) &
         + beta_dt(isz(is),isx(is))*s(it_delay)*polarity(is)*taper(is)
      if (.not.isFS) then
         if (isDipole) then
            p(isz(is)-2,isx(is)) = p(isz(is)-2,isx(is)) &
               - beta_dt(isz(is)-2,isx(is))*s(it_delay)*polarity(is)*taper(is)
         endif
      endif
   endif
enddo
end subroutine add_source_mod

!=====================================================================

subroutine add_source_mig(p,beta_dt,s,it,isDipole)
logical, intent(in) :: isDipole
real,    intent(in) :: beta_dt(:,:),s(:)
integer, intent(in) :: it
real,    intent(inout) :: p(:,:)
integer :: is, it_delay
do is=1,ns
   it_delay=it-idelay(is)
   if (it_delay > 0) then 
      p(isz(is),isx(is)) = p(isz(is),isx(is)) &
         + beta_dt(isz(is),isx(is))*s(it_delay)*polarity(is)*taper(is)
      if (isDipole) then
         p(isz(is)-2,isx(is)) = p(isz(is)-2,isx(is)) &
            - beta_dt(isz(is)-2,isx(is))*s(it_delay)*polarity(is)*taper(is)
      endif
   endif
enddo
end subroutine add_source_mig

!=====================================================================

subroutine subtract_source(p0,beta_dt,s,it,isDipole)
integer, intent(in) :: it
logical, intent(in) :: isDipole
real,    intent(in) :: beta_dt(:,:),s(:)
real,    intent(inout) :: p0(:,:)
integer :: is, it_delay
do is=1,ns
   it_delay=it-idelay(is)
   if (it_delay > -2) then 
      p0(isz(is),isx(is)) = p0(isz(is),isx(is)) - &
         s(it_delay+2) * beta_dt(isz(is),isx(is))*polarity(is)*taper(is)
      if (isDipole) then
         p0(isz(is)-2,isx(is)) = p0(isz(is)-2,isx(is)) &
            + s(it_delay+2) * beta_dt(isz(is)-2,isx(is))*polarity(is)*taper(is)
      endif
   endif
enddo
end subroutine subtract_source

!=====================================================================

subroutine add_seis(q,beta_dt,ng,seis,it)
integer, intent(in) :: ng,it
real,    intent(in) :: beta_dt(:,:),seis(:,:)
real,    intent(inout) :: q(:,:)
integer :: ig
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
integer             :: it_out, ig
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

subroutine modeling(fdm,xs_taper,s)
type(fdmod2d),      intent(inout) :: fdm
type(ms_shot2d),    intent(inout) :: s
real,               intent(in)    :: xs_taper
s%seis = 0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
do it=1,nt
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif 
   call add_source(pp,beta_dt,fdm%s%sou,it,isFS,isDipole)
   if (isFS) call add_FS(pp,npml,fs_thick)
   !if (mod(it,100)==1) call snapshot("p_",p,it,nz_pml,nx_pml)
   call store_data(pp,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine modeling

!=====================================================================

subroutine check_time_now()
call date_and_time(date,time_now)
write(*,*)icheck,time_now
!call flush(6)
icheck=icheck+1
end subroutine check_time_now

!=====================================================================

subroutine modeling_saveBC(fdm,bc,xs_taper,s)
type(fdmod2d),      intent(inout) :: fdm
type(ms_shot2d),    intent(inout) :: s
type(bc2d_general), intent(inout) :: bc
real,               intent(in)    :: xs_taper
s%seis = 0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
do it=1,nt
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_source(pp,beta_dt,fdm%s%sou,it,isFS,isDipole)
   ! Free Surface
   if (isFS) call add_FS(pp,npml,fs_thick)
   call store_data(pp,s%ng,it,s%seis)
   call saveBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
               bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   !call wf_refresh(p0,p1,p)
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
call copy_array(pp0,bc%p_nt_1,p1,bc%p_nt)
call finalize_variables()
call deallocate_and_free(p,p0,p1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine modeling_saveBC

!=====================================================================

subroutine modeling_saveWF(fdm,xs_taper,wf,s)
type(fdmod2d),    intent(inout) :: fdm
type(ms_shot2d),  intent(inout) :: s
type(wf2d),       intent(inout) :: wf
real,             intent(in)    :: xs_taper
s%seis = 0.0
wf%pwf=0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
do it=1,nt
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_source(pp,beta_dt,fdm%s%sou,it,isFS,isDipole)
   ! Free Surface
   if (isFS) call add_FS(pp,npml,fs_thick)
   call store_data(pp,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp)
   !call wf_refresh(p0,p1,p)
   wf%pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine modeling_saveWF

!=====================================================================

subroutine linear_modeling(fdm,xs_taper,s)
type(fdmod2d),       intent(inout) :: fdm
type(ms_shot2d),     intent(inout) :: s
real,                intent(in)    :: xs_taper
s%seis = 0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
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
   call add_source(pp,beta_dt,fdm%s%sou,it,isFS,isDipole)
   !--------- Pertubation field -----------------------------------------
   if (bc_type==1) then
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pq,pq0,pq1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine linear_modeling

!=====================================================================

subroutine linear_modeling_saveBC(fdm,bc,xs_taper,s)
type(fdmod2d),      intent(inout) :: fdm
type(ms_shot2d),    intent(inout) :: s
type(bc2d_general), intent(inout) :: bc
real,               intent(in)    :: xs_taper
s%seis = 0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
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
   call add_source(pp,beta_dt,fdm%s%sou,it,isFS,isDipole)
   call saveBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
                  bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   !--------- Pertubation field -----------------------------------------
   if (bc_type==1) then
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pq,pq0,pq1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call copy_array(p0,bc%p_nt_1,p1,bc%p_nt)
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine linear_modeling_saveBC

!=====================================================================

subroutine linear_modeling_saveWF(fdm,wf,xs_taper,s)
type(fdmod2d),      intent(inout) :: fdm
type(ms_shot2d),    intent(inout) :: s
type(wf2d),         intent(inout) :: wf
real,               intent(in)    :: xs_taper
s%seis = 0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
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
   call add_source(pp,beta_dt,fdm%s%sou,it,isFS,isDipole)
   wf%pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   !--------- Pertubation field -----------------------------------------
   if (bc_type==1) then
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pq,pq0,pq1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine linear_modeling_saveWF

!=====================================================================

subroutine rtm(s,fdm,xs_taper,img)
type(fdmod2d),      intent(inout) :: fdm
type(ms_shot2d),    intent(inout) :: s
type(image2d),      intent(inout) :: img
real,               intent(in)    :: xs_taper
img%img=0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
if (isSaveBC) then
   call allocate_and_initial(p_top,p_bottom,bc_len,nx,nt)
   call allocate_and_initial(p_left,p_right,nz,bc_len,nt)
   call allocate_and_initial(p_nt,p_nt_1,nz_pml,nx_pml)
endif
if (isSaveWF) then
   call allocate_and_initial(pwf,nz+2,nx+2,nt)
endif
! Forward propagate to save BC
do it=1,nt
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_source(pp,beta_dt,fdm%s%sou,it,isDipole)
   if (isSaveBC) then
      call saveBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
                  p_left(:,:,it),p_right(:,:,it),bc_len)
   endif
   if (isSaveWF) then
      pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   endif
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
if (isSaveBC) then
   call copy_array(pp1,p_nt,pp0,p_nt_1)
   call copy_array(p_nt,pp0,p_nt_1,pp1)
endif
if (isSaveWF) then
   call initial(pp0,pp1,pp)
   pp0(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt)
   pp1(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt-1)
endif   
! Back propagate 
do it=nt-2,1,-1
   if (isSaveBC) then
      !-.-.-.-.-.-. Reconstruct Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_source(pp0,beta_dt,fdm%s%sou,it,isDipole)
      if (bc_type==1) then
         call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
      else
         call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,npml,nz_pml,nx_pml,fd_order)
      endif
      call loadBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
                  p_left(:,:,it),p_right(:,:,it),bc_len)
   endif
   if (isSaveWF) then
      pp(npml:nzp+1,npml:nxp+1)=pwf(:,:,it)
   endif
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   if (bc_type==1) then
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pq,pq0,pq1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_seis(pq,beta_dt,s%ng,s%seis,it)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
if (bc_type==2) call deallocate_and_free(u1)
if (isSaveBC) then
   call deallocate_and_free(p_top,p_bottom,p_left,p_right)
   call deallocate_and_free(p_nt,p_nt_1)
endif
if (isSaveWF) call deallocate_and_free(pwf)
end subroutine rtm

!=====================================================================

subroutine rtm_saveBC(s,bc,fdm,xs_taper,img)
type(fdmod2d),      intent(inout) :: fdm
type(ms_shot2d),    intent(inout) :: s
type(image2d),      intent(inout) :: img
type(bc2d_general), intent(in)    :: bc
real,               intent(in)    :: xs_taper
img%img=0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
call copy_array(bc%p_nt,p0,bc%p_nt_1,p1)
do it=nt-2,1,-1
   !-.-.-.-.-.-. Reconstruct Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_source(pp0,beta_dt,fdm%s%sou,it,isDipole)
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,npml,nz_pml,nx_pml,fd_order)
   endif
   call loadBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
               bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   if (bc_type==1) then
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pq,pq0,pq1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_seis(pq,beta_dt,s%ng,s%seis,it)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine rtm_saveBC

!=====================================================================
subroutine rtm_saveWF(s,wf,fdm,xs_taper,img)
type(fdmod2d),      intent(inout) :: fdm
type(ms_shot2d),    intent(inout) :: s
type(image2d),      intent(inout) :: img
type(wf2d),         intent(in)    :: wf
real,               intent(in)    :: xs_taper
img%img=0.0
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
pp1(npml:nzp+1,npml:nxp+1)=wf%pwf(:,:,nt)
pp0(npml:nzp+1,npml:nxp+1)=wf%pwf(:,:,nt-1)
do it=nt-2,1,-1
   ! ------ Load P field ------
   pp(npml:nzp+1,npml:nxp+1)=wf%pwf(:,:,it)
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   if (bc_type==1) then
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pq,pq0,pq1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_seis(pq,beta_dt,s%ng,s%seis,it)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine rtm_saveWF

!=====================================================================

subroutine illum(fdm,xs_taper,s,img)
type(fdmod2d),       intent(in)    :: fdm
type(ms_shot2d),     intent(in)    :: s
type(image2d),       intent(inout) :: img
real,                intent(in)    :: xs_taper
integer :: ig
if (.not.allocated(img%img)) then
   call allocate_and_initial(img%illum,img%nz,img%nx)
else
   img%illum = 0.0
endif
call initial_variables(fdm,s)
call setup_source_receiver(s,xs_taper)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
call allocate_and_initial(illum_source,nz,nx)
do it=1,nt
   !--------- Pressure field --------------------------------------------
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_source(pp,beta_dt,fdm%s%sou,it,isFS,isDipole)
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      illum_source(iz-npml,ix-npml) = illum_source(iz-npml,ix-npml) &
                                    + (abs(pp(iz,ix)))**img%illum_order
   end forall
   call wf_refresh(pp0,pp1,pp)
   !call wf_refresh(p0,p1,p)
enddo  ! End of time-marching loop
if (img%isIllumReceiver) then
   call initial(pp0,pp,pp1)
   call allocate_and_initial(illum_receiver,nz,nx)
   do it=1,nt
      if (bc_type==1) then
         call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
      else
         call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
                 spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
      endif
      do ig=1,s%ng
         pp(igz(ig),igx(ig)) = pp(igz(ig),igx(ig)) + beta_dt(igz(ig),igx(ig))*fdm%s%sou(it)
      enddo
      forall (iz=npml+1:nzp,ix=npml+1:nxp)
         illum_receiver(iz-npml,ix-npml) = illum_receiver(iz-npml,ix-npml) &
                                         + (abs(pp(iz,ix)))**img%illum_order
      end forall
      call wf_refresh(pp0,pp1,pp)
   enddo  ! End of time-marching loop
endif
if (img%isIllumReceiver) then
   img%illum=illum_source*illum_receiver
   call deallocate_and_free(illum_receiver)
else
   img%illum=illum_source
endif
call finalize_variables()
call deallocate_and_free(p,p0,p1,illum_source)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine illum

!=====================================================================

subroutine setup_source_receiver(s,xs_taper)
type(ms_shot2d),       intent(in) :: s
real,                  intent(in) :: xs_taper
integer :: is, ig, n_taper, itaper
real,allocatable :: taper_tmp(:),taper_tmp_i(:)
real :: xs_taper_max,xs_taper_min
call allocate_and_initial(isx,isz,idelay,ns)
call allocate_and_initial(polarity,ns)
! Setup source position
do is=1,ns
   isx(is) = npml + int(s%xs(is)/dx)+1
   isz(is) = npml + int(s%zs(is)/dx)+1
   ! Source position is shifted downward 1 node if on the FS.
   if (isz(is) == fs(isx(is))) isz(is) = isz(is) + 1
   idelay(is) = nint(s%d(is)/dt)
   polarity(is) = s%p(is)
enddo
call allocate_and_initial(igx,igz,s%ng)
! Setup receiver position
do ig=1,s%ng
   igx(ig) = npml + int(s%xg(ig)/dx)+1
   igz(ig) = npml + int(s%zg(ig)/dx)+1
   ! Geophone position is shifted downward 1 node if on the FS.
   if (igz(ig) == fs(igx(ig))) igz(ig) = igz(ig) + 1
enddo
! add the tapering for the source on the eadge
call allocate_and_initial(taper,ns)
taper=1.0
if (xs_taper>0.01) then
   n_taper=nint(xs_taper/dx)+1
   call allocate_and_initial(taper_tmp,taper_tmp_i,n_taper)
   taper_tmp=hanning(n_taper)
   taper_tmp_i=taper_tmp(n_taper:1:-1)
   xs_taper_max=maxval(s%xs(1:ns)) - xs_taper
   xs_taper_min=minval(s%xs(1:ns)) + xs_taper
   do is=1,ns
      if (s%xs(is) < xs_taper_min) then
         itaper=nint((xs_taper_min-s%xs(is))/dx)+1
         taper(is)=taper_tmp_i(itaper)
      elseif (s%xs(is) > xs_taper_max) then
         itaper=nint((s%xs(is)-xs_taper_max)/dx)+1
         taper(is)=taper_tmp_i(itaper)
      endif
   enddo
   deallocate(taper_tmp,taper_tmp_i)
endif 
end subroutine setup_source_receiver

!=====================================================================

subroutine initial_variables_mod(fdm,s)
type(fdmod2d),   intent(in) :: fdm
type(ms_shot2d), intent(in) :: s
call copy_variable(fdm%f%fd_type,fd_type)
call copy_variable(fdm%f%fd_order,fd_order,fdm%f%bc_type,bc_type)
call copy_variable(fdm%f%nt,nt,fdm%m%nz,nz,fdm%m%nx,nx,fdm%f%npml,npml)
call copy_variable(fdm%f%dt,dt,fdm%m%dx,dx)
call copy_variable(fdm%f%isFS,isFS,fdm%f%isDipole,isDipole)
call copy_variable(fdm%f%isSaveBC,isSaveBC,fdm%f%isSaveWF,isSaveWF)
dnt_mod=nint(s%dt/fdm%f%dt)
call i_v_kernal(fdm%m%v,fdm%f%fs)
call copy_variable(s%ns,ns,s%ng,ng)
end subroutine initial_variables_mod

!=====================================================================

subroutine finalize_variables()
if (bc_type==1) then
   call deallocate_and_free(isx,isz,igx,igz,fs,idelay)
   call deallocate_and_free(alpha,temp1,temp2,beta_dt)
   call deallocate_and_free(polarity)
endif
if (bc_type==2) then 
   call deallocate_and_free(isx,isz,igx,igz,fs,idelay)
   call deallocate_and_free(beta_dtdx,beta_dt,u1)
   call deallocate_and_free(polarity)
endif
if (allocated(taper)) call deallocate_and_free(taper)
if (allocated(polarity)) call deallocate_and_free(polarity)
if (allocated(idelay)) call deallocate_and_free(idelay)
end subroutine finalize_variables

!=====================================================================

subroutine nx_nz_pml
nx_pml = nx + 2 * npml
nz_pml = nz + 2 * npml
nxp = nx + npml
nzp = nz + npml
end subroutine nx_nz_pml

!=====================================================================

subroutine i_v_kernal(v,fs_in)
real,intent(in) :: v(:,:)
integer, intent(in) :: fs_in(:)
call nx_nz_pml
bc_len=(fd_order-20)/2+1
fs_thick=bc_len-1
dtdx = dt/dx
call allocate_and_initial(beta_dt,nz_pml,nx_pml)
beta_dt=(v*dt)**2.0
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
   if (bc_type==1) then ! Absorbing Boundary Condition
      call allocate_and_initial(damp,alpha,kappa,temp1,temp2,nz_pml,nx_pml)
      vmin=minval(v)
      call abc_get_damp2d(nx,nz,npml,dx,vmin,damp)
      alpha=(v*dtdx)**2.0
      kappa=damp*dt
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
else ! Staggerred FD
   write(*,*)"Under construction!"
endif
end subroutine i_v_kernal

end module module_a2d_multisource
