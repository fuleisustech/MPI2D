module module_a3d

use module_global
use module_datatype
use module_array
use module_utility
use module_boundary3d
use module_fd3d
use module_string
use module_io
use module_aperture3d     !modified by Bowen Guo
implicit none
character(len=100), private :: fd_type
integer,private :: ix,iy,iz,it,isx,isy,isz,ig,nx,ny,nz,npml,nxp,nyp,nzp,&
   nt,dnt_mod,nx_pml,ny_pml,nz_pml,ic,fd_order,bc_type,bc_len,fs_thick,&
   nit_display,it_display
real, private   :: dtdx,dx,dy,dz,dt,vmin,c1t3,beta_dt
integer,private,allocatable :: igx(:),igy(:),igz(:),fs(:,:)
logical,private :: isFS, isDipole, isSaveBC, isSaveWF, isWriteData
real,private,target,allocatable :: p0(:,:,:),p1(:,:,:),p(:,:,:),&
   q0(:,:,:),q1(:,:,:),q(:,:,:)
real,private,pointer :: pp0(:,:,:),pp1(:,:,:),pp(:,:,:),&
   pq0(:,:,:),pq1(:,:,:),pq(:,:,:)
! only use in illumination
real,private,allocatable :: illum_source(:,:,:),illum_receiver(:,:,:)
! only use in Absorbing Boundary Condition 
real,private,allocatable :: kappa(:,:,:)
!real,private,allocatable :: spgx(:),spgy(:),spgz(:),spcx1(:),spcx2(:),&
!   spcy1(:),spcy2(:),spcz1(:),spcz2(:),beta_dtdx(:,:,:),u1(:,:,:)
real,   private, allocatable :: damp(:,:,:) ! only use in abc initialize
real,   private, allocatable :: temp1(:,:,:),temp2(:,:,:),temp3(:,:,:),temp4(:,:,:)! Bowen Guo
real,   private, allocatable :: p_top(:,:,:,:),p_bottom(:,:,:,:),&
   p_left(:,:,:,:),p_right(:,:,:,:),p_front(:,:,:,:),p_back(:,:,:,:),&
   p_nt(:,:,:),p_nt_1(:,:,:),pwf(:,:,:,:) ! only use in rtm code
integer(kind=MPI_OFFSET_KIND) :: disp
integer :: icheck

interface a3d_mod
   module procedure modeling
   module procedure modeling_saveBC
   module procedure modeling_saveWF
end interface a3d_mod

interface a3d_bmod
   module procedure linear_modeling
   module procedure linear_modeling_saveBC
   module procedure linear_modeling_saveWF
end interface a3d_bmod

interface a3d_mig
   module procedure rtm
   module procedure rtm_saveBC
   module procedure rtm_saveWF
end interface a3d_mig

interface a3d_illum
   module procedure illum
end interface a3d_illum

interface add_source
   module procedure add_source_mod
   module procedure add_source_mig
end interface add_source

private initial_variables,finalize_variables,nx_ny_nz_pml,i_v_kernal,&
   a3d_abc_kernal,add_source,subtract_source,store_data,&
   add_FS,pertubation,setup_source_receiver,saveBC3d,loadBC3d

contains

!=====================================================================

subroutine add_source_mod(p,v,dt,s,isFS,isDipole)
logical, intent(in) :: isFS, isDipole
real,    intent(in) :: v(:,:,:),dt,s
real,    intent(inout) :: p(:,:,:)
beta_dt=(v(isz,isy,isx)*dt)**2.0
p(isz,isy,isx) = p(isz,isy,isx) + beta_dt*s
if (.not.isFS) then
   if (isDipole) then
      beta_dt=(v(isz-2,isy,isx)*dt)**2.0
      p(isz-2,isy,isx) = p(isz-2,isy,isx) - beta_dt*s
   endif
endif
end subroutine add_source_mod

!=====================================================================

subroutine add_source_mig(p,v,dt,s,isDipole)
logical, intent(in) :: isDipole
real,    intent(in) :: v(:,:,:),dt,s
real,    intent(inout) :: p(:,:,:)
beta_dt=(v(isz,isy,isx)*dt)**2.0
p(isz,isy,isx) = p(isz,isy,isx) + beta_dt*s
if (isDipole) then
   beta_dt=(v(isz-2,isy,isx)*dt)**2.0
   p(isz-2,isy,isx) = p(isz-2,isy,isx) - beta_dt*s
endif
end subroutine add_source_mig

!=====================================================================

subroutine subtract_source(p0,v,dt,s,it,isDipole)
integer, intent(in) :: it
logical, intent(in) :: isDipole
real,    intent(in) :: v(:,:,:),dt,s(:)
real,    intent(inout) :: p0(:,:,:)
beta_dt=(v(isz,isy,isx)*dt)**2.0
p0(isz,isy,isx) = p0(isz,isy,isx) - s(it+2) * beta_dt
if (isDipole) then
   beta_dt=(v(isz-2,isy,isx)*dt)**2.0
   p0(isz-2,isy,isx) = p0(isz-2,isy,isx) + s(it+2) * beta_dt
endif
end subroutine subtract_source

!=====================================================================

subroutine add_seis(q,v,dt,ng,seis,it)
integer, intent(in) :: ng,it
real,    intent(in) :: v(:,:,:),dt,seis(:,:)
real,    intent(inout) :: q(:,:,:)
!$OMP PARALLEL DO PRIVATE(ig,beta_dt)
do ig=1,ng
   beta_dt=(v(igz(ig),igy(ig),igx(ig))*dt)**2.0
   q(igz(ig),igy(ig),igx(ig)) = q(igz(ig),igy(ig),igx(ig)) &
                        + beta_dt * seis(it,ig)
enddo
!$OMP END PARALLEL DO
end subroutine add_seis

!=====================================================================

subroutine add_FS(p,npml,fs_thick)
integer, intent(in) :: npml,fs_thick
real,    intent(inout) :: p(:,:,:)
integer :: iz
p(npml+1,:,:)=0.0
!$OMP PARALLEL DO PRIVATE(iz)
do iz=1,fs_thick
   p(npml+1-iz,:,:)=-p(npml+1+iz,:,:)
enddo
!$OMP END PARALLEL DO
end subroutine add_FS

!=====================================================================

subroutine store_data(p,ng,it,seis)
integer, intent(in) :: ng,it
real,    intent(in) :: p(:,:,:)
real,    intent(inout) :: seis(:,:)
integer             :: it_out
do ig=1,ng
   if (dnt_mod.ne.1) then
      if (mod(it,dnt_mod).eq.1) then
         it_out=(it-1)/dnt_mod+1
         seis(it_out,ig) = p(igz(ig),igy(ig),igx(ig))
      endif
   else
      seis(it,ig) = p(igz(ig),igy(ig),igx(ig))
   endif    
enddo
end subroutine store_data

!=====================================================================

subroutine pertubation(ic,refl,v,dt,p,p1,q,nz_pml,ny_pml,nx_pml)
integer, intent(in) :: ic, nz_pml,ny_pml,nx_pml
real,    intent(in) :: refl(:,:,:),v(:,:,:),p(:,:,:),p1(:,:,:),dt
real,    intent(inout) :: q(:,:,:)
if (ic.eq.0) then
   !forall (iz=1:nz_pml,iy=1:ny_pml,ix=1:nx_pml)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,beta_dt)
   do ix=1,nx_pml
      do iy=1,ny_pml
         do iz=1,nz_pml
            beta_dt=(v(iz,iy,ix)*dt)**2.0
            q(iz,iy,ix)=q(iz,iy,ix)+p1(iz,iy,ix)*refl(iz,iy,ix)*beta_dt
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !endforall
elseif (ic.eq.1) then
   !forall (iz=1:nz_pml,iy=1:ny_pml,ix=1:nx_pml)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,beta_dt)
   do ix=1,nx_pml
      do iy=1,ny_pml
         do iz=1,nz_pml
            beta_dt=(v(iz,iy,ix)*dt)**2.0
            q(iz,iy,ix)=q(iz,iy,ix)+(p(iz,iy,ix)-p1(iz,iy,ix))&
                        *refl(iz,iy,ix)*beta_dt
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !endforall
elseif (ic.eq.2) then
   !forall (iz=2:nz_pml-1,iy=2:ny_pml-1,ix=2:nx_pml-1)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,beta_dt)
   do ix=1,nx_pml
      do iy=1,ny_pml
         do iz=1,nz_pml
            beta_dt=(v(iz,iy,ix)*dt)**2.0
            q(iz,iy,ix)=q(iz,iy,ix)&
                       +(p1(iz-1,iy,ix)+p1(iz+1,iy,ix)&
                       +p1(iz,iy-1,ix)+p1(iz,iy+1,ix)&
                       +p1(iz,iy,ix-1)+p1(iz,iy,ix+1)&
                       -6.0*p1(iz,iy,ix))*refl(iz,iy,ix)*beta_dt
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !endforall
endif
end subroutine pertubation

!=====================================================================

subroutine image_condition(ic,img,nzp,nyp,nxp,npml,p0,p1,q0)
integer, intent(in) :: ic,nzp,nyp,nxp,npml
real,    intent(in) :: p0(:,:,:),p1(:,:,:),q0(:,:,:)
real,    intent(inout) :: img(:,:,:)
if (ic.eq.0) then
   !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix)
   do ix=npml+1,nxp
      do iy=npml+1,nyp
         do iz=npml+1,nzp
            img(iz-npml,iy-npml,ix-npml) = img(iz-npml,iy-npml,ix-npml) &
                                         + p1(iz,iy,ix) * q0(iz,iy,ix)
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !endforall
elseif (ic.eq.1) then
   !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix)
   do ix=npml+1,nxp
      do iy=npml+1,nyp
         do iz=npml+1,nzp
            img(iz-npml,iy-npml,ix-npml) = img(iz-npml,iy-npml,ix-npml) + &
                   (p0(iz,iy,ix)- p1(iz,iy,ix)) * q0(iz,iy,ix)
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !endforall
elseif (ic.eq.2) then
   !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix)
   do ix=npml+1,nxp
      do iy=npml+1,nyp
         do iz=npml+1,nzp
            img(iz-npml,iy-npml,ix-npml) = img(iz-npml,iy-npml,ix-npml) &
                  + (p1(iz-1,iy,ix) + p1(iz+1,iy,ix) + p1(iz,iy-1,ix) &
                    +p1(iz,iy+1,ix) + p1(iz,iy,ix-1) + p1(iz,iy,ix+1)&
                    -6.0 * p1(iz,iy,ix)) * q0(iz,iy,ix)
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !endforall
endif
end subroutine image_condition

!=====================================================================

subroutine modeling(fdm,s,output,mem)
type(fdmod3d),      intent(inout) :: fdm
type(shot3d),       intent(inout) :: s
character(len=*),   intent(in)    :: output
integer(8),optional,intent(inout) :: mem
integer                           :: it_out
call initial_variables(fdm,s,mem)
!if (isWriteData) then
!   call MPI_FILE_OPEN(MPI_COMM_SELF, output, &
!         MPI_MODE_WRONLY + MPI_MODE_CREATE, &
!         MPI_INFO_NULL, fid, ierr)
if (isWriteData) then
   open(13,file=output,access="direct",recl=I4)
else
   s%seis = 0.0
endif
call setup_source_receiver(s,mem)
call allocate_and_initial(p,p0,p1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml,mem)
do it=1,nt
   if (bc_type==1) then
      if (fd_type.eq."2ND") then 
         call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
      elseif (fd_type.eq."Lax_wen") then 
         call a3d_abc_kernal(pp,pp0,pp1,temp1,temp2,temp3,temp4,fdm%m%v,dx,dy,dz,nz_pml,&
              ny_pml,nx_pml,fd_order)
      endif
      !call a3d_abc_kernal(p,p0,p1,temp1,temp2,temp3,alpha,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif 
   call add_source(pp,fdm%m%v,dt,fdm%s%sou(it),isFS,isDipole)
   !call add_source(p,beta_dt,model%s(it),isFS,isDipole)
   if (isFS) call add_FS(pp,npml,fs_thick)
   if (mod(it,nit_display)==1) then
      call date_and_time(date,time_now);it_display=int(it/nit_display)+1
      call message_all("rank=",rank,"time=",time_now,it_display,"/100 done")
   endif
   !if (mod(it,100)==1) call snapshot("p_",p,it,nz_pml,ny_pml,nx_pml)
   if (isWriteData) then
      if (dnt_mod.ne.1) then
         if (mod(it,dnt_mod).eq.1) then
            it_out=(it-1)/dnt_mod+1
            do ig=1,s%ng
               write(13,rec=(ig-1)*s%nt+it_out)pp(igz(ig),igy(ig),igx(ig))
               !disp=((ig-1)*s%nt+(it_out-1))*4
               !call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
               !call MPI_FILE_WRITE(fid, pp(igz(ig),igy(ig),igx(ig)),1, MPI_REAL, &
               !          MPI_STATUS_IGNORE, ierr)
            enddo
         endif
      else
         do ig=1,s%ng
            write(13,rec=(ig-1)*s%nt+it)pp(igz(ig),igy(ig),igx(ig))
            !disp=((ig-1)*s%nt+(it-1))*4
            !call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
            !call MPI_FILE_WRITE(fid, pp(igz(ig),igy(ig),igx(ig)),1, MPI_REAL, &
            !          MPI_STATUS_IGNORE, ierr)
         enddo
      endif
   else
      call store_data(pp,s%ng,it,s%seis)
   endif
   !call store_data(p,seis%ng,it,seis%seis)
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
if (isWriteData) then
   close(13)
endif
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,mem)
!if (bc_type==2) deallocate(u1)
end subroutine modeling

!=====================================================================

subroutine check_time_now()
call date_and_time(date,time_now)
write(*,*)icheck,time_now
!call flush(6)
icheck=icheck+1
end subroutine check_time_now

!=====================================================================

subroutine modeling_saveBC(fdm,bc,s,mem)
type(fdmod3d),      intent(inout) :: fdm
type(shot3d),       intent(inout) :: s
type(bc3d_general), intent(inout) :: bc
integer(8),optional,intent(inout) :: mem
s%seis = 0.0
call initial_variables(fdm,s,mem)
call setup_source_receiver(s,mem)
call allocate_and_initial(p,p0,p1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml,mem)
do it=1,nt
   if (bc_type==1) then
      call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a2d_sponge_kernal(p,p0,p1,u1,c1t2,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_source(pp,fdm%m%v,dt,fdm%s%sou(it),isFS,isDipole)
   ! Free Surface
   if (isFS) call add_FS(pp,npml,fs_thick)
   call store_data(pp,s%ng,it,s%seis)
   call saveBC3d(pp,nz,ny,nx,nzp,nyp,nxp,npml,bc%p_top(:,:,:,it),&
      bc%p_bottom(:,:,:,it),bc%p_left(:,:,:,it),bc%p_right(:,:,:,it),&
      bc%p_front(:,:,:,it),bc%p_back(:,:,:,it),bc_len)
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
call copy_array(pp0,bc%p_nt_1,pp1,bc%p_nt)
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,mem)
!if (bc_type==2) deallocate(u1)
end subroutine modeling_saveBC

!=====================================================================

subroutine modeling_saveWF(fdm,wf,s,mem)
type(fdmod3d),    intent(inout) :: fdm
type(shot3d),     intent(inout) :: s
type(wf3d),       intent(inout) :: wf
integer(8),optional,intent(inout) :: mem
s%seis = 0.0
wf%pwf=0.0
call initial_variables(fdm,s,mem)
call setup_source_receiver(s,mem)
call allocate_and_initial(p,p0,p1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml,mem)
do it=1,nt
   if (bc_type==1) then
      call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_source(pp,fdm%m%v,dt,fdm%s%sou(it),isFS,isDipole)
   ! Free Surface
   if (isFS) call add_FS(pp,npml,fs_thick)
   call store_data(pp,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp)
   wf%pwf(:,:,:,it)=pp(npml:nzp+1,npml:nyp+1,npml:nxp+1)
enddo  ! End of time-marching loop
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,mem)
!if (bc_type==2) deallocate(u1)
end subroutine modeling_saveWF

!=====================================================================

subroutine linear_modeling(fdm,s,mem)
type(fdmod3d),       intent(inout) :: fdm
type(shot3d),        intent(inout) :: s
integer(8),optional,intent(inout) :: mem
s%seis = 0.0
call initial_variables(fdm,s,mem)
call setup_source_receiver(s,mem)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml,mem)
do it=1,nt
   !--------- Pressure field --------------------------------------------
   if (bc_type==1) then
      call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_source(pp,fdm%m%v,dt,fdm%s%sou(it),isFS,isDipole)
   !--------- Pertubation field -----------------------------------------
   if (bc_type==1) then
      call a3d_abc_kernal(pq,pq0,pq1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(q,q0,q1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call pertubation(ic,fdm%m%refl,fdm%m%v,dt,pp,pp1,pq,nz_pml,ny_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,q,q0,q1,mem)
!if (bc_type==2) deallocate(u1)
end subroutine linear_modeling

!=====================================================================

subroutine linear_modeling_saveBC(fdm,bc,s,mem)
type(fdmod3d),      intent(inout) :: fdm
type(shot3d),       intent(inout) :: s
type(bc3d_general), intent(inout) :: bc
integer(8),optional,intent(inout) :: mem
s%seis = 0.0
call initial_variables(fdm,s,mem)
call setup_source_receiver(s,mem)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml,mem)
do it=1,nt
   !--------- Pressure field --------------------------------------------
   if (bc_type==1) then
      call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_source(pp,fdm%m%v,dt,fdm%s%sou(it),isFS,isDipole)
   call saveBC3d(pp,nz,ny,nx,nzp,nyp,nxp,npml,bc%p_top(:,:,:,it),&
      bc%p_bottom(:,:,:,it),bc%p_left(:,:,:,it),bc%p_right(:,:,:,it),&
      bc%p_front(:,:,:,it),bc%p_back(:,:,:,it),bc_len)
   !--------- Pertubation field -----------------------------------------
   if (bc_type==1) then
      call a3d_abc_kernal(pq,pq0,pq1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(q,q0,q1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call pertubation(ic,fdm%m%refl,fdm%m%v,dt,pp,pp1,pq,nz_pml,ny_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call copy_array(pp0,bc%p_nt_1,pp1,bc%p_nt)
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,q,q0,q1,mem)
!if (bc_type==2) deallocate(u1)
end subroutine linear_modeling_saveBC

!=====================================================================

subroutine linear_modeling_saveWF(fdm,wf,s,mem)
type(fdmod3d),      intent(inout) :: fdm
type(shot3d),       intent(inout) :: s
type(wf3d),         intent(inout) :: wf
integer(8),optional,intent(inout) :: mem
s%seis = 0.0
call initial_variables(fdm,s,mem)
call setup_source_receiver(s,mem)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml,mem)
do it=1,nt
   !--------- Pressure field --------------------------------------------
   if (bc_type==1) then
      call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_source(pp,fdm%m%v,dt,fdm%s%sou(it),isFS,isDipole)
   wf%pwf(:,:,:,it)=pp(npml:nzp+1,npml:nyp+1,npml:nxp+1)
   !--------- Pertubation field -----------------------------------------
   if (bc_type==1) then
      call a3d_abc_kernal(pq,pq0,pq1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(q,q0,q1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call pertubation(ic,fdm%m%refl,fdm%m%v,dt,pp,pp1,pq,nz_pml,ny_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,q,q0,q1,mem)
!if (bc_type==2) deallocate(u1)
end subroutine linear_modeling_saveWF

!=====================================================================

subroutine rtm(s,fdm,img,mem)
type(fdmod3d),      intent(inout) :: fdm
type(shot3d),       intent(inout) :: s
type(image3d),      intent(inout) :: img
integer(8),optional,intent(inout) :: mem
img%img=0.0
call initial_variables(fdm,s,mem)
call setup_source_receiver(s,mem)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml,mem)
if (isSaveBC) then
   call allocate_and_initial(p_top,p_bottom,bc_len,ny,nx,nt,mem)
   call allocate_and_initial(p_left,p_right,nz,ny,bc_len,nt,mem)
   call allocate_and_initial(p_front,p_back,nz,bc_len,nx,nt,mem)
   call allocate_and_initial(p_nt,p_nt_1,nz_pml,ny_pml,nx_pml,mem)
endif
if (isSaveWF) then
   call allocate_and_initial(pwf,nz+2,ny+2,nx+2,nt,mem)
endif
! Forward propagate to save BC
do it=1,nt
   if (bc_type==1) then
      call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(pp,p0,p1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_source(pp,fdm%m%v,dt,fdm%s%sou(it),isDipole)
   if (isSaveBC) then
      call saveBC3d(pp,nz,ny,nx,nzp,nyp,nxp,npml,p_top(:,:,:,it),&
         p_bottom(:,:,:,it),p_left(:,:,:,it),p_right(:,:,:,it),&
         p_front(:,:,:,it),p_back(:,:,:,it),bc_len)
   endif
   if (isSaveWF) then
      pwf(:,:,:,it)=pp(npml:nzp+1,npml:nyp+1,npml:nxp+1)
   endif
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
if (isSaveBC) then
   call copy_array(pp1,p_nt,pp0,p_nt_1)
   call copy_array(p_nt,pp0,p_nt_1,pp1)
endif
if (isSaveWF) then
   call initial(pp0,pp1,pp)
   pp0(npml:nzp+1,npml:nyp+1,npml:nxp+1)=pwf(:,:,:,nt)
   pp1(npml:nzp+1,npml:nyp+1,npml:nxp+1)=pwf(:,:,:,nt-1)
endif   
! Back propagate 
do it=nt-2,1,-1
   if (isSaveBC) then
      !-.-.-.-.-.-. Reconstruct Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_source(pp0,fdm%m%v,dt,fdm%s%sou,it,isDipole)
      if (bc_type==1) then
         call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,npml,nzp,nyp,nxp,fd_order)
      else
         !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,npml,nz_pml,ny_pml,nx_pml,fd_order)
      endif
      call loadBC3d(pp,nz,ny,nx,nzp,nyp,nxp,npml,p_top(:,:,:,it),p_bottom(:,:,:,it),&
         p_left(:,:,:,it),p_right(:,:,:,it),p_front(:,:,:,it),p_back(:,:,:,it),bc_len)
   endif
   if (isSaveWF) then
      pp(npml:nzp+1,npml:nyp+1,npml:nxp+1)=pwf(:,:,:,it)
   endif
!   call snapshot("p_",p,it,nz_pml,nx_pml)
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   if (bc_type==1) then
      call a3d_abc_kernal(pq,pq0,pq1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(q,q0,q1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_seis(pq,fdm%m%v,dt,s%ng,s%seis,it)
!   call snapshot("q_",q,it,nz_pml,nx_pml)
   call image_condition(ic,img%img,nzp,nyp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,q,q0,q1,mem)
!if (bc_type==2) deallocate(u1)
if (isSaveBC) then
   call deallocate_and_free(p_top,p_bottom,p_left,p_right,p_front,p_back,mem)
   call deallocate_and_free(p_nt,p_nt_1,mem)
endif
if (isSaveWF) call deallocate_and_free(pwf,mem)
end subroutine rtm

!=====================================================================

subroutine rtm_saveBC(s,bc,fdm,img,mem)
type(fdmod3d),        intent(inout) :: fdm
type(shot3d),         intent(inout) :: s
type(image3d),        intent(inout) :: img
type(bc3d_general),   intent(in)    :: bc
integer(8),optional,intent(inout) :: mem
img%img=0.0
call initial_variables(fdm,s,mem)
call setup_source_receiver(s,mem)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml,mem)
call copy_array(bc%p_nt,pp0,bc%p_nt_1,pp1)
do it=nt-2,1,-1
   !-.-.-.-.-.-. Reconstruct Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_source(pp0,fdm%m%v,dt,fdm%s%sou,it,isDipole)
   if (bc_type==1) then
      call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,npml,nzp,nyp,nxp,fd_order)
   else
      !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call loadBC3d(pp,nz,ny,nx,nzp,nyp,nxp,npml,bc%p_top(:,:,:,it),bc%p_bottom(:,:,:,it),&
      bc%p_left(:,:,:,it),bc%p_right(:,:,:,it),bc%p_front(:,:,:,it),&
      bc%p_back(:,:,:,it),bc_len)
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   if (bc_type==1) then
      call a3d_abc_kernal(pq,pq0,pq1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(q,q0,q1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif
   call add_seis(pq,fdm%m%v,dt,s%ng,s%seis,it)
   call image_condition(ic,img%img,nzp,nyp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,q,q0,q1,mem)
!if (bc_type==2) deallocate(u1)
end subroutine rtm_saveBC

!=====================================================================
subroutine rtm_saveWF(s,wf,fdm,img,mem)
type(fdmod3d),      intent(inout) :: fdm
type(shot3d),       intent(inout) :: s
type(image3d),      intent(inout) :: img
type(wf3d),         intent(in)    :: wf
integer(8),optional,intent(inout) :: mem
img%img=0.0
call initial_variables(fdm,s,mem)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml)
pp1(npml:nzp+1,npml:nyp+1,npml:nxp+1)=wf%pwf(:,:,:,nt)
pp0(npml:nzp+1,npml:nyp+1,npml:nxp+1)=wf%pwf(:,:,:,nt-1)
do it=nt-2,1,-1
   ! ------ Load P field ------
   pp(npml:nzp+1,npml:nyp+1,npml:nxp+1)=wf%pwf(:,:,:,it)
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   if (bc_type==1) then
      call a3d_abc_kernal(pq,pq0,pq1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(q,q0,q1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_seis(pq,fdm%m%v,dt,s%ng,s%seis,it)
   call image_condition(ic,img%img,nzp,nyp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1,q,q0,q1)
!if (bc_type==2) deallocate(u1)
end subroutine rtm_saveWF

!=====================================================================

subroutine illum(fdm,s,img,mem)
type(fdmod3d),       intent(in)    :: fdm
type(shot3d),        intent(in)    :: s
type(image3d),       intent(inout) :: img
integer(8),optional,intent(inout) :: mem
if (.not.allocated(img%img)) then
   call allocate_and_initial(img%illum,img%nz,img%ny,img%nx,mem)
else
   img%illum = 0.0
endif
call initial_variables(fdm,s)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,nz_pml,ny_pml,nx_pml,mem)
call wf_assign(p,pp,p0,pp0,p1,pp1)
!if (bc_type==2) call allocate_and_initial(u1,nz_pml,ny_pml,nx_pml)
call allocate_and_initial(illum_source,nz,ny,nx,mem)
do it=1,nt
   !--------- Pressure field --------------------------------------------
   if (bc_type==1) then
      call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,nz_pml,ny_pml,nx_pml,fd_order)
   else
      !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
      !        spcy1,spcy2,spcz1,spcz2,npml,nz_pml,ny_pml,nx_pml,fd_order)
   endif
   call add_source(pp,fdm%m%v,dt,fdm%s%sou(it),isFS,isDipole)
   !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix)
   do ix=npml+1,nxp
      do iy=npml+1,nyp
         do iz=npml+1,nzp
            illum_source(iz-npml,iy-npml,ix-npml) &
          = illum_source(iz-npml,iy-npml,ix-npml) &
          + (abs(p(iz,iy,ix)))**img%illum_order
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
if (img%isIllumReceiver) then
   call initial(pp0,pp,pp1)
   call allocate_and_initial(illum_receiver,nz,ny,nx,mem)
   do it=1,nt
      if (bc_type==1) then
         call a3d_abc_kernal(pp,pp0,pp1,fdm%m%v,kappa,dtdx,c1t3,&
                 nz_pml,ny_pml,nx_pml,fd_order)
      else
         !call a3d_sponge_kernal(p,p0,p1,u1,c1t3,beta_dtdx,spgx,spgy,spgz,spcx1,spcx2,&
         !        spcx1,spcx2,spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
      endif
      !$OMP PARALLEL DO PRIVATE(iz,iy,ix,beta_dt)
      do ig=1,s%ng
         beta_dt=(fdm%m%v(igz(ig),igy(ig),igx(ig))*dt)**2.0
         pp(igz(ig),igy(ig),igx(ig)) = pp(igz(ig),igy(ig),igx(ig)) &
                                    + beta_dt*fdm%s%sou(it)
      enddo
      !$OMP END PARALLEL DO
      !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
      !$OMP PARALLEL DO PRIVATE(iz,iy,ix)
      do ix=npml+1,nxp
         do iy=npml+1,nyp
            do iz=npml+1,nzp
               illum_receiver(iz-npml,iy-npml,ix-npml) &
             = illum_receiver(iz-npml,iy-npml,ix-npml) &
             + (abs(pp(iz,iy,ix)))**img%illum_order
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !end forall
      call wf_refresh(pp0,pp1,pp)
   enddo  ! End of time-marching loop
endif
if (img%isIllumReceiver) then
   img%illum=illum_source*illum_receiver
   call deallocate_and_free(illum_receiver,mem)
else
   img%illum=illum_source
endif
call finalize_variables(mem)
call deallocate_and_free(p,p0,p1)
call deallocate_and_free(illum_source)
!if (bc_type==2) deallocate(u1)
end subroutine illum

!=====================================================================

subroutine setup_source_receiver(seis,mem)
type(shot3d),       intent(in) :: seis
integer(8),optional,intent(inout) :: mem
call allocate_and_initial(igx,igy,igz,seis%ng,mem)
! Setup source position
isx = npml+int(seis%xs(1)/dx)+1
isy = npml+int(seis%ys(1)/dx)+1
isz = npml+int(seis%zs(1)/dx)+1
! Source position is shifted downward 1 node if on the FS.
if (isz == fs(isy,isx)) isz = isz + 1
! Setup receiver position
!$OMP PARALLEL DO PRIVATE(ig)
do ig=1,seis%ng
   igx(ig) = npml + int(seis%xg(ig)/dx)+1
   igy(ig) = npml + int(seis%yg(ig)/dx)+1
   igz(ig) = npml + int(seis%zg(ig)/dx)+1
   ! Geophone position is shifted downward 1 node if on the FS.
   if (igz(ig) == fs(igy(ig),igx(ig))) igz(ig) = igz(ig) + 1
enddo
!$OMP END PARALLEL DO
end subroutine setup_source_receiver

!=====================================================================

subroutine initial_variables(fdm,s,mem)
type(fdmod3d),intent(in) :: fdm
type(shot3d), intent(in) :: s
integer(8),optional,intent(inout) :: mem
call copy_variable(fdm%f%fd_type,fd_type)
call copy_variable(fdm%f%fd_order,fd_order,fdm%f%bc_type,bc_type)
call copy_variable(fdm%f%nt,nt,fdm%m%nz,nz,fdm%m%ny,ny,fdm%m%nx,nx,fdm%f%npml,npml)
call copy_variable(fdm%f%dt,dt,fdm%m%dx,dx)
call copy_variable(fdm%m%dy,dy,fdm%m%dz,dz) ! modified by Bowen Guo
call copy_variable(fdm%f%isFS,isFS,fdm%f%isDipole,isDipole,fdm%f%isWriteData,isWriteData)
call copy_variable(fdm%f%isSaveBC,isSaveBC,fdm%f%isSaveWF,isSaveWF)
nit_display=int(nt/100)
dnt_mod=nint(s%dt/fdm%f%dt)
call i_v_kernal(fdm%m%v,fdm%f%fs,mem)
end subroutine initial_variables

!=====================================================================

subroutine finalize_variables(mem)
integer(8),optional,intent(inout) :: mem
if (bc_type==1) then
   call deallocate_and_free(igx,igy,igz,mem)
   call deallocate_and_free(kappa,mem)
   call deallocate_and_free(fs,mem)
endif
!if (bc_type==2) deallocate(igx,igz,beta_dtdx,beta_dt,u1,fs)
end subroutine finalize_variables

!=====================================================================

subroutine nx_ny_nz_pml
nx_pml = nx + 2 * npml
ny_pml = ny + 2 * npml
nz_pml = nz + 2 * npml
nxp = nx + npml
nyp = ny + npml
nzp = nz + npml
end subroutine nx_ny_nz_pml

!=====================================================================

subroutine i_v_kernal(v,fs_in,mem)
real,intent(in) :: v(:,:,:)
integer, intent(in) :: fs_in(:,:)
integer(8),optional,intent(inout) :: mem
call nx_ny_nz_pml
bc_len=(fd_order-20)/2+1
fs_thick=bc_len-1
dtdx = dt/dx
call allocate_and_initial(fs,ny_pml,nx_pml,mem)
fs = fs_in
if (fd_type.eq."2ND") then ! 2nd Order FD
   if (fd_order==22) then
      c1t3=C21*3.0
   elseif (fd_order==24) then
      c1t3=C41*3.0
   elseif (fd_order==26) then
      c1t3=C61*3.0
   elseif (fd_order==28) then
      c1t3=C81*3.0
   endif
   if (bc_type==1) then ! Absorbing Boundary Condition
      call allocate_and_initial(damp,kappa,nz_pml,ny_pml,nx_pml,mem)
      vmin=minval(v)
      call abc_get_damp3d(nx,ny,nz,npml,dx,vmin,damp)
      kappa=damp*dt
      !write(*,*)size(damp,1),size(damp,2),size(damp,3)
      !write(*,*)vmin
      !call write_binfile("vel.bin",v)
      !call write_binfile("damp.bin",damp)
      call deallocate_and_free(damp,mem)
   elseif (bc_type==2) then ! Sponge Boundary Condition
   !   call allocate_and_initial(beta_dtdx,u1,nz_pml,ny_pml,nx_pml)
   !   call allocate_and_initial(spgx,spcz1,spcz2,nx_pml)
   !   call allocate_and_initial(spgz,spcx1,spcx2,nz_pml)
   !   call spgcoef(spgx,spgy,spgz,spcx1,spcx2,spcy1,spcy2,spcz1,spcz2,v,&
   !                nx_pml,ny_pml,nz_pml,npml,dx,dx,dx,dt)
   !   beta_dtdx=(v*dtdx)**2.0
!      call write_binfile("beta_dtdx.bin",beta_dtdx,nz_pml,nx_pml)
!      call write_binfile("vel.bin",v,nz_pml,nx_pml)
   endif
elseif (fd_type.eq."SG") then 
   write(*,*)"Under construction!"

elseif (fd_type.eq."Lax_wen") then 
   call allocate_and_initial(damp,kappa,nz_pml,ny_pml,nx_pml,mem)
   call allocate_and_initial(temp1,temp2,temp3,temp4,nz_pml,ny_pml,nx_pml,mem)
   vmin=minval(v)
   call abc_get_damp3d(nx,ny,nz,npml,dx,vmin,damp)
   kappa=damp*dt
   temp1=dt**2
   temp2=kappa-1
   temp3=2-kappa
   temp4=-(v**2)*(dt**4)/12
   call deallocate_and_free(damp,mem)
endif

end subroutine i_v_kernal

end module module_a3d
