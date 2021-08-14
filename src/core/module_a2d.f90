!  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
!  King Abdullah University of Science and Technology, All rights reserved.
!
!  Sponsors of CSIM are granted a non-exclusive, irrevocable royalty free
!  world-wide license to use this software and associated documentation files
!  (the "Software"), in their business, including the rights to modify and
!  distribute the Software to their affiliates, partners, clients or consultants
!  as necessary in the conduct of the sponsors business, all without accounting
!  to the King Abdullah University of Science and Technology, subject to the
!  following conditions:
!
!  The above copyright notice and this permission notice shall be included
!  in all copies or substantial portions of the Software.
!
!  Warranty Disclaimer:
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!  DEALINGS IN THE SOFTWARE.
!
!  file:        module_a2d.f90
!  @author:     Xin Wang
!  email:       netmacross@gmail.com, Gerard.Schuster@kaust.edu.sa
!  last update: 2014.03.20 (Xin)   
!  purpose:     module for 2D acoustic wave equation fd solution

module module_a2d

use module_global
use module_datatype
use module_array
use module_utility
use module_boundary2d
use module_fd2d
use module_string
use module_io


implicit none
character(len=100),private :: fd_type
integer,private :: ix,iz,it,isx,isz,ig,nx,nz,npml,nx1,nz1,nxp,nzp,&
   nt,dnt_mod,nx_pml,nz_pml,ic,fd_order,bc_type,bc_len,fs_thick
real, private   :: dtdx,dx,dz,dt,vmin,c1t2
integer,private,allocatable :: igx(:),igz(:),fs(:)
logical,private :: isFS, isDipole, isSaveBC, isSaveWF, isWriteData, isReadData
real,private,target,allocatable :: p0(:,:),p1(:,:),p(:,:),&
   q0(:,:),q1(:,:),q(:,:)
real,private,pointer :: pp0(:,:),pp1(:,:),pp(:,:),pq0(:,:),pq1(:,:),pq(:,:)
real,private, allocatable :: beta_dt(:,:)
! only use in illumination
real,private,allocatable :: illum_source(:,:),illum_receiver(:,:)
! only use in Absorbing Boundary Condition 
real,private,allocatable :: alpha(:,:),temp1(:,:),temp2(:,:) 
! only use in Sponge BC
real,private,allocatable :: spgx(:),spgz(:),spcx1(:),spcx2(:),&
   spcz1(:),spcz2(:),beta_dtdx(:,:),u1(:,:)
real,   private, allocatable :: damp(:,:), kappa(:,:) ! only use in abc initialize
real,   private, allocatable :: p_top(:,:,:),p_bottom(:,:,:),p_left(:,:,:),&
   p_right(:,:,:),p_nt(:,:),p_nt_1(:,:), pwf(:,:,:),bwf(:,:,:) ! only use in rtm code
integer(kind=MPI_OFFSET_KIND) :: disp

interface initial_variables
   module procedure initial_variables_mod
end interface initial_variables

interface a2d_mod
   module procedure modeling
end interface a2d_mod

interface a2d_bmod
   module procedure linear_modeling
end interface a2d_bmod

interface a2d_mig
   module procedure rtm
   module procedure rtm_h
   module procedure rtm_hil
end interface a2d_mig

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

!if (rank.eq.0) then
!write(*,*),isz,isx
!write(*,*),s
!endif
!write(*,*), size(beta_dt)

p(isz,isx) = p(isz,isx) + beta_dt(isz,isx)*s
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
p0(isz,isx) = p0(isz,isx) - s(it+2) * beta_dt(isz,isx)
if (isDipole) then
   p0(isz-2,isx) = p0(isz-2,isx) + s(it+2) * beta_dt(isz-2,isx)
endif
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
! image condition for separated wavefield by Hilbert transform
! here the wavefield is divided into up- and down- going wavefield
subroutine image_condition_hil(ic,img0,img1,img2,nzp,nxp,npml,p0,p1,q0)
integer, intent(in)    :: ic,nzp,nxp,npml
real,    intent(in)    :: p0(:,:),p1(:,:),q0(:,:)
real,    intent(inout) :: img0(:,:), img1(:,:), img2(:,:)
real,    allocatable   :: Hzp(:,:),Hzq(:,:)
integer                :: nz,nx
nz = nzp - npml;
nx = nxp - npml;


allocate(Hzp(nz,nx))
allocate(Hzq(nz,nx))


Hzp = p1(npml+1:nzp,npml+1:nxp)
Hzq = q0(npml+1:nzp,npml+1:nxp)

   !write(*,*) "c0"
call Matrix_hil(Hzp,nz,nx)
call Matrix_hil(Hzq,nz,nx)

   !write(*,*) "c1"

if (ic.eq.0) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      img0(iz-npml,ix-npml) = img0(iz-npml,ix-npml) &
                                 + p1(iz,ix) * q0(iz,ix)
      img2(iz-npml,ix-npml) = img2(iz-npml,ix-npml) &
                                 + p1(iz,ix) * q0(iz,ix)-Hzp(iz-npml,ix-npml)*Hzq(iz-npml,ix-npml)
      img1(iz-npml,ix-npml) = img1(iz-npml,ix-npml) &
                                 + p1(iz,ix) * q0(iz,ix)+Hzp(iz-npml,ix-npml)*Hzq(iz-npml,ix-npml)
   endforall
elseif (ic.eq.1) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      img0(iz-npml,ix-npml) = img0(iz-npml,ix-npml) + &
            (p0(iz,ix)- p1(iz,ix)) * q0(iz,ix)
   endforall
elseif (ic.eq.2) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      img0(iz-npml,ix-npml) = img0(iz-npml,ix-npml) &
            + (p1(iz-1,ix) + p1(iz+1,ix) + p1(iz,ix-1) &
             + p1(iz,ix+1)-4.0 * p1(iz,ix)) * q0(iz,ix)
   endforall
endif

end subroutine image_condition_hil


!==========================================================
subroutine image_condition_h(ic,img,nzp,nxp,npml,p0,p1,q0,nh)
integer, intent(in) :: ic,nzp,nxp,npml
real,    intent(in) :: p0(:,:),p1(:,:),q0(:,:)
real,    intent(inout) :: img(:,:,:)
integer :: nh,ih
integer,allocatable :: h(:)

! Image condition
!nh = size(img%img_h,3)
!call allocate_and_initial(h,nh)
allocate(h(nh))
do ih = 1,nh
   h(ih) = -(nh-1)/2 +ih -1
enddo


if (ic.eq.0) then
do ih=1,nh
   forall (iz=npml+1:nzp,ix=npml+1+(nh-1)/2:nxp-(nh-1)/2)
      img(iz-npml,ix-npml,ih) = img(iz-npml,ix-npml,ih) &
                                 + p1(iz,ix-h(ih)) * q0(iz,ix+h(ih))
   endforall
enddo

elseif (ic.eq.1) then
do ih=1,nh
   forall (iz=npml+1:nzp,ix=npml+1+(nh-1)/2:nxp-(nh-1)/2)
      img(iz-npml,ix-npml,ih) = img(iz-npml,ix-npml,ih) + &
            (p0(iz,ix-h(ih))- p1(iz,ix-h(ih))) * q0(iz,ix+h(ih))
   endforall
enddo

elseif (ic.eq.2) then
do ih=1,nh
   forall (iz=npml+1:nzp,ix=npml+1+(nh-1)/2:nxp-(nh-1)/2)
      img(iz-npml,ix-npml,ih) = img(iz-npml,ix-npml,ih) &
            + (p1(iz-1,ix-h(ih)) + p1(iz+1,ix-h(ih)) + p1(iz,ix-1-h(ih)) &
             + p1(iz,ix+1-h(ih))-4.0 * p1(iz,ix-h(ih))) * q0(iz,ix+h(ih))
   endforall
enddo

endif
end subroutine image_condition_h

!=====================================================================

subroutine modeling(fdm,bc,wf,s,output)
type(fdmod2d),      intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(bc2d_general), intent(inout) :: bc
type(wf2d),         intent(inout) :: wf
character(len=*),   intent(in)    :: output
integer                           :: fid,it_out,mod_type
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
   !write(*,*)size(fdm%s%sou); call flush(6)
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
   img%illum=0.0
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
subroutine rtm_hil(s,bc,wf,fdm,img,output,is_Hilbert)
type(fdmod2d),      intent(inout)  :: fdm
type(shot2d),       intent(inout)  :: s
type(bc2d_general), intent(inout)  :: bc
type(wf2d),         intent(inout)  :: wf
type(image2d),      intent(inout)  :: img
character(len=*),   intent(in)     :: output
integer                            :: fid,rtm_type
logical                            :: is_Hilbert

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
img%img_t=0.0
img%img_m=0.0

if (img%isIllum) then
   img%illum=0.0
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
   !write(*,*) "b0"
   !call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   if (mod(it,8)==1) call image_condition_hil(ic,img%img,img%img_t,img%img_m,nzp,nxp,npml,pp0,pp1,pq0)
   !write(*,*) "b1"
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
end subroutine rtm_hil
!=====================================================================


subroutine rtm_h(s,bc,wf,fdm,img,output,nh)
type(fdmod2d),      intent(inout)  :: fdm
type(shot2d),       intent(inout)  :: s
type(bc2d_general), intent(inout)  :: bc
type(wf2d),         intent(inout)  :: wf
type(image2d),      intent(inout)  :: img
character(len=*),   intent(in)     :: output
integer                            :: fid,rtm_type,nh
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
   img%illum=0.0
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
   !bwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   !call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
  ! if (rank.eq.0) write(*,), size(img%img_h,3)

   call image_condition_h(ic,img%img_h,nzp,nxp,npml,pp0,pp1,pq0,nh)
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
end subroutine rtm_h

!=====================================================================

subroutine geometric(fdm,s,geo)
type(fdmod2d),      intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
real,               intent(out)   :: geo(:,:)
integer                           :: mod_type
call initial_variables(fdm,s)
mod_type=1; geo=-100000.0
call setup_source_receiver(s)
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
   !if (mod(it,100)==1) call snapshot("p_",p,it,nz_pml,nx_pml)
   !write(*,*)size(fdm%s%sou); call flush(6)
   call add_source(pp,beta_dt,fdm%s%sou(it),isFS,isDipole)
   if (isFS) then
      call add_FS(pp,npml,fs_thick)
   endif
   !do ix=1,nx;   do iz=1,nz
   !   if (geo(iz,ix)<pp(iz+npml,ix+npml)) then
   !      geo(iz,ix)=pp(iz+npml,ix+npml)
   !   endif
   !enddo;   enddo;
   geo=max(geo,p(nz1:nzp,nx1:nxp))
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine geometric

!=====================================================================

subroutine setup_source_receiver(s)
type(shot2d),       intent(in) :: s
call allocate_and_initial(igx,igz,s%ng)
! Setup source position
isx = npml+int(s%xs(1)/dx)+1
isz = npml+int(s%zs(1)/dx)+1
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
fd_type=fdm%f%fd_type; fd_order=fdm%f%fd_order; bc_type=fdm%f%bc_type
nt=fdm%f%nt; nz=fdm%m%nz; nx=fdm%m%nx; npml=fdm%f%npml
dt=fdm%f%dt; dx=fdm%m%dx
isFS=fdm%f%isFS; isDipole=fdm%f%isDipole; isWriteData=fdm%f%isWriteData; 
isReadData=fdm%f%isReadData; isSaveBC=fdm%f%isSaveBC; isSaveWF=fdm%f%isSaveWF
!dnt_mod=nint(s%dt/fdm%f%dt)
dnt_mod=1
!write(*,*) "v=",fdm%m%v(50,120)
!write(*,*) "v=",fdm%m%v(50,121)
!write(*,*) "s-dt=",s%dt
!write(*,*) "fdmdt=",fdm%f%dt
call i_v_kernal(fdm%m%v,fdm%f%fs)
end subroutine initial_variables_mod

!=====================================================================

subroutine finalize_variables()
if (bc_type==1) then
   call deallocate_and_free(igx,igz,fs)
   call deallocate_and_free(alpha,temp1,temp2,beta_dt)
endif
if (bc_type==2) then 
   call deallocate_and_free(igx,igz,fs)
   call deallocate_and_free(beta_dtdx,beta_dt,u1)
endif
end subroutine finalize_variables

!=====================================================================

subroutine nx_nz_pml
nx_pml = nx + 2 * npml;   nz_pml = nz + 2 * npml
nxp = nx + npml;   nzp = nz + npml
nx1 = npml + 1;   nz1 = npml + 1
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

end module module_a2d
