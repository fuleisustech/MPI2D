module module_a2d_mul

use module_global
use module_datatype
use module_array
use module_utility
use module_boundary2d
use module_fd2d
use module_string
use module_io

implicit none
integer,private :: ix,iz,it,isx,isz,ig,nx,nz,nzw,npml,nxp,nzp,nzwp,&
  nt,nx_pml,nz_pml,nzw_pml,ic,fd_type,fd_order,bc_type,bc_len,fs_thick
real, private   :: dtdx,dx,dz,dt,vmin,c1t2
integer,private, allocatable :: igx(:),igz(:),fs(:)
logical,private :: isFS, isDipole, isSaveBC, isSaveWF
real,private,target,allocatable :: p0(:,:),p1(:,:),p(:,:),&
   sh0(:,:),sh1(:,:),sh(:,:),sw0(:,:),sw1(:,:),sw(:,:),&
   swh0(:,:),swh1(:,:),swh(:,:),q0(:,:),q1(:,:),q(:,:),&
   rh0(:,:),rh1(:,:),rh(:,:),rw0(:,:),rw1(:,:),rw(:,:),&
   rwh0(:,:),rwh1(:,:),rwh(:,:)
real, private,pointer :: pp0(:,:),pp1(:,:),pp(:,:),&
   psh0(:,:),psh1(:,:),psh(:,:),psw0(:,:),psw1(:,:),psw(:,:),&
   pswh0(:,:),pswh1(:,:),pswh(:,:),pq0(:,:),pq1(:,:),pq(:,:),&
   prh0(:,:),prh1(:,:),prh(:,:),prw0(:,:),prw1(:,:),prw(:,:),&
   prwh0(:,:),prwh1(:,:),prwh(:,:)
real, private, allocatable :: beta_dt(:,:),beta_dt_w(:,:),beta_dt_h(:,:)
real,private,target,allocatable :: sw2_0(:,:),sw2_1(:,:),sw2(:,:),&
   sw3_0(:,:),sw3_1(:,:),sw3(:,:),sw4_0(:,:),sw4_1(:,:),sw4(:,:),&
   sw5_0(:,:),sw5_1(:,:),sw5(:,:),sw6_0(:,:),sw6_1(:,:),sw6(:,:)
real,private,pointer :: psw2_0(:,:),psw2_1(:,:),psw2(:,:),&
   psw3_0(:,:),psw3_1(:,:),psw3(:,:),psw4_0(:,:),psw4_1(:,:),psw4(:,:),&
   psw5_0(:,:),psw5_1(:,:),psw5(:,:),psw6_0(:,:),psw6_1(:,:),psw6(:,:)
! only use in abc kind
real,   private, allocatable ::alpha(:,:),temp1(:,:),&
   temp2(:,:),alpha_w(:,:),temp1_w(:,:),temp2_w(:,:),&
   alpha_h(:,:),temp1_h(:,:),temp2_h(:,:)
! only use in sponge kind
real,   private, allocatable ::beta_dtdx(:,:),&
   beta_dtdx_w(:,:),beta_dtdx_h(:,:),u1(:,:),u2(:,:),&
   spgx(:),spgz(:),spcx1(:),spcx2(:),spcz1(:),spcz2(:),&
   spgxw(:),spgzw(:),spcx1w(:),spcx2w(:),spcz1w(:),spcz2w(:),&
   spgxh(:),spgzh(:),spcx1h(:),spcx2h(:),spcz1h(:),spcz2h(:)
! only in initialize
real,   private, allocatable :: damp(:,:), kappa(:,:),&
   damp_w(:,:),kappa_w(:,:),damp_h(:,:),kappa_h(:,:) 
! only use in rtm withour saveBC or WF
real,   private, allocatable :: p_top(:,:,:),p_bottom(:,:,:),&
   p_left(:,:,:),p_right(:,:,:),p_nt(:,:),p_nt_1(:,:),&
   sw_top(:,:,:),sw_bottom(:,:,:),sw_left(:,:,:),sw_right(:,:,:),&
   sw_nt(:,:),sw_nt_1(:,:),sh_top(:,:,:),sh_bottom(:,:,:),&
   sh_left(:,:,:),sh_right(:,:,:),sh_nt(:,:),sh_nt_1(:,:),pwf(:,:,:)
! only use in illumination
real,   private, allocatable :: illum_source(:,:),illum_receiver(:,:)
integer :: icheck

interface initial_variables
   module procedure initial_variables_mod
end interface initial_variables

!=====================================================================
! source side multiples
! ====================================== free surface
!  s----->/\             g
!   \ sh /  \           /                ic: p*q
!    \  /    \         /
!     \/sw    \       /
! -------------\-----/-----------------  ocean bottom
!               \   /                    bc: bc2d_sm -- sw,sh,p
!               p\ /q                    wf: p
!                 x                      x - image point        
! In saving WF aproach, same as conventional RTM, only save p
interface a2d_sm_bmod
   module procedure sm_modeling
   module procedure sm_modeling_saveBC
   module procedure sm_modeling_saveWF
end interface a2d_sm_bmod

interface a2d_sm_mig
   module procedure sm_rtm
   module procedure sm_rtm_saveBC
   module procedure sm_rtm_saveWF
end interface a2d_sm_mig

interface a2d_sm_illum
   module procedure sm_illum
end interface a2d_sm_illum

! receiver side multiples
! ====================================== free surface
!  s             /\----->g
!   \           /  \ rh /                ic: p*q
!    \         /    \  /
!     \       /      \/rw
! -----\-----/-------------------------- ocean bottom
!       \   /                            bc: bc2d_general -- p
!       p\ /q                            wf: p    
!         x                              x - image point        
interface a2d_rm_bmod
   module procedure rm_modeling
   module procedure rm_modeling_saveBC
   module procedure rm_modeling_saveWF
end interface a2d_rm_bmod

interface a2d_rm_mig
   module procedure rm_rtm
   module procedure rm_rtm_saveBC
   module procedure rm_rtm_saveWF
end interface a2d_rm_mig

interface a2d_rm_illum
   module procedure rm_illum
end interface a2d_rm_illum

!=============================================================
! both side multiples
! ====================================== free surface
!  s----->/\             /\----->g
!   \ sh /  \           /  \ rh /        ic: p*q
!    \  /    \         /    \  /
!     \/sw    \       /    rw\/
! -------------\-----/-----------------  ocean bottom
!               \   /                    bc: bc2d_sm -- sw,sh,p
!               p\ /q                    wf: p
!                 x                      x - image point        
! In saving WF aproach, same as conventional RTM, only save p
interface a2d_srm_bmod
   module procedure srm_modeling
   module procedure srm_modeling_saveBC
   module procedure srm_modeling_saveWF
end interface a2d_srm_bmod

interface a2d_srm_mig
   module procedure srm_rtm
   module procedure srm_rtm_saveBC
   module procedure rm_rtm_saveWF ! same as rm rtm
end interface a2d_srm_mig

interface a2d_srm_illum
   module procedure srm_illum
end interface a2d_srm_illum

!interface a2d_sm1_bmod
!   module procedure sm1_modeling
!   module procedure sm1_modeling_saveBC
!   module procedure sm1_modeling_saveWF
!end interface a2d_sm1_bmod

!interface a2d_sm1_mig
!   module procedure sm1_rtm
!   module procedure sm1_rtm_saveBC
!   module procedure sm1_rtm_saveWF
!end interface a2d_sm1_mig

!interface a2d_sm1_illum
!   module procedure sm1_illum
!end interface a2d_sm1_illum

!interface a2d_sm1_bmod
!   module procedure sm1_modeling
!   module procedure sm1_modeling_saveBC
!   module procedure sm1_modeling_saveWF
!end interface a2d_sm1_bmod

!interface a2d_sm1_rtm
!   module procedure sm1_rtm
!   module procedure sm1_rtm_saveBC
!   module procedure sm1_rtm_saveWF
!end interface a2d_sm1_rtm

interface add_source
   module procedure add_source_FS
   module procedure add_source_noFS
end interface add_source

interface subtract_source
   module procedure subtract_source_FS
   module procedure subtract_source_noFS
end interface subtract_source

private initial_variables,finalize_variables,nx_nz_pml,i_v_kernal,&
   a2d_abc_kernal,a2d_sponge_kernal,add_source,subtract_source,&
   pertubation,subtract_wavefield,store_data,add_FS,&
   setup_source_receiver,saveBC2d,loadBC2d

contains

!=====================================================================

subroutine add_source_FS(p,beta_dt,s)
real,    intent(in) :: beta_dt(:,:),s
real,    intent(inout) :: p(:,:)
p(isz,isx) = p(isz,isx) + beta_dt(isz,isx)*s
end subroutine add_source_FS

!=====================================================================

subroutine add_source_noFS(p,beta_dt,s,isDipole)
logical, intent(in) :: isDipole
real,    intent(in) :: beta_dt(:,:),s
real,    intent(inout) :: p(:,:)
p(isz,isx) = p(isz,isx) + beta_dt(isz,isx)*s
if (isDipole) then
   p(isz-2,isx) = p(isz-2,isx) - beta_dt(isz-2,isx)*s
endif
end subroutine add_source_noFS

!=====================================================================

subroutine subtract_source_FS(p0,beta_dt,s,it)
integer, intent(in) :: it
real,    intent(in) :: beta_dt(:,:),s(:)
real,    intent(inout) :: p0(:,:)
p0(isz,isx) = p0(isz,isx) - s(it+2) * beta_dt(isz,isx)
end subroutine subtract_source_FS

!=====================================================================

subroutine subtract_source_noFS(p0,beta_dt,s,it,isDipole)
integer, intent(in) :: it
logical, intent(in) :: isDipole
real,    intent(in) :: beta_dt(:,:),s(:)
real,    intent(inout) :: p0(:,:)
p0(isz,isx) = p0(isz,isx) - s(it+2) * beta_dt(isz,isx)
if (isDipole) then
   p0(isz-2,isx) = p0(isz-2,isx) + s(it+2)*beta_dt(isz-2,isx)
endif
end subroutine subtract_source_noFS

!=====================================================================

subroutine subtract_wavefield(p,b,b1,beta_dt,refl,nz_pml,nx_pml,ic)
real, intent(inout) :: p(:,:),b(:,:),b1(:,:),beta_dt(:,:),refl(:,:)
integer, intent(in) :: nz_pml,nx_pml,ic
if (ic.eq.0) then
   forall (iz=1:nz_pml,ix=1:nx_pml) 
      p(iz,ix)=p(iz,ix) - b1(iz,ix)*refl(iz,ix)*beta_dt(iz,ix)
   end forall
elseif (ic.eq.1) then
   forall (iz=1:nz_pml,ix=1:nx_pml) 
      p(iz,ix)=p(iz,ix)-(b(iz,ix)-b1(iz,ix))*refl(iz,ix)*beta_dt(iz,ix)
   end forall
elseif (ic.eq.2) then
   forall (iz=1:nz_pml,ix=1:nx_pml) 
      p(iz,ix)=p(iz,ix)-(b1(iz-1,ix)+b1(iz+1,ix)+b1(iz,ix-1)+b1(iz,ix+1)&
                    -4.0*b1(iz,ix))*refl(iz,ix)*beta_dt(iz,ix)
   end forall
endif
end subroutine subtract_wavefield
      
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
do ig=1,ng
   seis(it,ig) = p(igz(ig),igx(ig))
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
   forall (iz=2:nz_pml-1,ix=2:nx_pml)
      q(iz,ix)=q(iz,ix)+(p1(iz-1,ix)+p1(iz+1,ix)+p1(iz,ix-1)+p1(iz,ix+1)&
                    -4.0*p1(iz,ix))*refl(iz,ix)*beta_dt(iz,ix)
   endforall
endif
end subroutine pertubation

!=====================================================================

subroutine image_condition(ic,img,nzp,nxp,npml,p0,p1,q0)
integer, intent(in) :: ic,nxp,nzp,npml
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

subroutine check_time_now()
call date_and_time(date,time_now)
write(*,*)icheck,time_now
call flush(6)
icheck=icheck+1
end subroutine check_time_now

!=====================================================================
! source side multiples
subroutine sm_modeling(fdm,s)
type(water_fdmod2d), intent(inout) :: fdm
type(shot2d),        intent(inout) :: s
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,&
               sw,psw,sw0,psw0,sw1,psw1,sh,psh,sh0,psh0,sh1,psh1,&
               swh,pswh,swh0,pswh0,swh1,pswh1)
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   !call snapshot("sw_",psw,it,nzw_pml,nx_pml)
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   !call snapshot("sh_",psh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   pswh = psw - psh
   !call snapshot("swh_",pswh,it,nzw_pml,nx_pml)
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pswh,pswh1,pp,nzw_pml,nx_pml)
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   !call snapshot("q_",pq,it,nz_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,&
                   pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1,sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
end subroutine sm_modeling

!=====================================================================

subroutine sm_modeling_saveBC(fdm,bc,s)
type(water_fdmod2d), intent(inout) :: fdm
type(shot2d),      intent(inout) :: s
type(bc2d_sm),     intent(inout) :: bc
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,sw,psw,sw0,psw0,&
        sw1,psw1,sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   call saveBC2d(psw,nzw,nx,nzwp,nxp,npml,bc%sw_top(:,:,it),&
      bc%sw_bottom(:,:,it),bc%sw_left(:,:,it),bc%sw_right(:,:,it),bc_len)
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   call saveBC2d(psh,nzw,nx,nzwp,nxp,npml,bc%sh_top(:,:,it),&
      bc%sh_bottom(:,:,it),bc%sh_left(:,:,it),bc%sh_right(:,:,it),bc_len)
   ! ------ mute the direct wave
   pswh = psw - psh
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pswh,pswh1,pp,nzw_pml,nx_pml)
   call saveBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
      bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nzw_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,&
                   pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call copy_array(pp1,bc%p_nt,psw1,bc%sw_nt,psh1,bc%sh_nt)
call copy_array(pp0,bc%p_nt_1,psw0,bc%sw_nt_1,psh0,bc%sh_nt_1)
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1,sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
end subroutine sm_modeling_saveBC

!=====================================================================

subroutine sm_modeling_saveWF(fdm,wf,s)
type(water_fdmod2d),intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(wf2d), intent(inout) :: wf
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,sw,psw,sw0,psw0,&
        sw1,psw1,sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
     nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
     nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   ! ------ mute the direct wave
   pswh = psw - psh
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pswh,pswh1,pp,nzw_pml,nx_pml)
   wf%p(:,:,it) = pp
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nzw_pml,nx_pml)
   call store_data(pq,s%ng,it,s%seis)
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1,sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
end subroutine sm_modeling_saveWF

!=====================================================================

subroutine sm_rtm(s,fdm,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
img%img=0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,sw,psw,sw0,psw0,&
        sw1,psw1,sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
if (isSaveBC) then
   call allocate_and_initial(p_top,p_bottom,sw_top,sw_bottom,sh_top,sh_bottom,&
                             bc_len,nx,nt)
   call allocate_and_initial(p_left,p_right,nz,bc_len,nt)
   call allocate_and_initial(p_nt,p_nt_1,nz_pml,nx_pml)
   call allocate_and_initial(sw_left,sw_right,sh_left,sh_right,nz,bc_len,nt)
   call allocate_and_initial(sw_nt,sw_nt_1,sh_nt,sh_nt_1,nzw_pml,nx_pml)
endif
if (isSaveWF) then
   call allocate_and_initial(pwf,nz+2,nx+2,nt)
endif
! Forward propagate to save BC or WF
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   if (isSaveBC) then
      call saveBC2d(psw,nzw,nx,nzwp,nxp,npml,sw_top(:,:,it),sw_bottom(:,:,it),&
              sw_left(:,:,it),sw_right(:,:,it),bc_len)
   endif
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   !call snapshot("sh_f_",sh,it,nzw_pml,nx_pml)
   if (isSaveBC) then
      call saveBC2d(psh,nzw,nx,nzwp,nxp,npml,sh_top(:,:,it),sh_bottom(:,:,it),&
              sh_left(:,:,it),sh_right(:,:,it),bc_len)
   endif
   ! ------ mute the direct wave
   pswh = psw - psh
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt,pswh,pswh1,pp,nzw_pml,nx_pml)
   if (isSaveBC) then 
     call saveBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
             p_left(:,:,it),p_right(:,:,it),bc_len)
   endif
   if (isSaveWF) then
      pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   endif
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp)
enddo

if (isSaveBC) then
   call copy_array(pp1,p_nt,psw1,sw_nt,psh1,sh_nt)
   call copy_array(pp0,p_nt_1,psw0,sw_nt_1,psh0,sh_nt_1)
   call copy_array(p_nt,pp0,sw_nt,psw0,sh_nt,psh0)
   call copy_array(p_nt_1,pp1,sw_nt_1,psw1,sh_nt_1,psh1)
   call initial(psw,psh,pswh)
   pswh0=psw0-psh0
   pswh1=psw1-psh1
endif
if (isSaveWF) then
   call deallocate_and_free(sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
   call initial(pp0,pp1,pp)
   pp0(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt)
   pp1(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt-1)
endif

! Back propagate
do it=nt-2,1,-1
   if (isSaveBC) then
      !-.-.-.-.-.-. Reconstruct P Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_wavefield(pp0,pswh,pswh1,beta_dt_w,fdm%m%refl_fs,nzw_pml,nx_pml,ic)
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
      call loadBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
              p_left(:,:,it),p_right(:,:,it),bc_len)
      !-.-.-.-.-.-. Reconstruct SW Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_source(psw0,beta_dt_w,fdm%s%sou,it)  
      call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,npml,nzwp,nxp,fd_order)
      call loadBC2d(psw,nzw,nx,nzwp,nxp,npml,sw_top(:,:,it),sw_bottom(:,:,it),&
              sw_left(:,:,it),sw_right(:,:,it),bc_len)
      !call snapshot("sw_",psw,it,nzw_pml,nx_pml)
      !-.-.-.-.-.-. Reconstruct SH Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_source(psh0,beta_dt_h,fdm%s%sou,it)  
      call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,npml,nzwp,nxp,fd_order)
      call loadBC2d(psh,nzw,nx,nzwp,nxp,npml,sh_top(:,:,it),sh_bottom(:,:,it),&
                  sh_left(:,:,it),sh_right(:,:,it),bc_len)
      !call snapshot("sh_",psh,it,nzw_pml,nx_pml)
      pswh = psw - psh
      !call snapshot("swh_",pswh,it,nzw_pml,nx_pml)
      call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh)
   endif
   if (isSaveWF) then
      pp(npml:nzp+1,npml:nxp+1)=pwf(:,:,it)
   endif
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-..-.-.-
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call add_seis(pq,beta_dt,s%ng,s%seis,it)
   !call snapshot("q_",pq,it,nz_pml,nx_pml)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
if (isSaveBC) then
   call deallocate_and_free(sw_top,sw_bottom,sw_left,sw_right)
   call deallocate_and_free(sh_top,sh_bottom,sh_left,sh_right)
   call deallocate_and_free(p_top,p_bottom,p_left,p_right)
   call deallocate_and_free(sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
   call deallocate_and_free(p_nt,p_nt_1,sw_nt,sw_nt_1,sh_nt,sh_nt_1)
else
   call deallocate_and_free(pwf)
endif
call deallocate_and_free(p,p0,p1,q,q0,q1)
end subroutine sm_rtm

!=====================================================================

subroutine sm_rtm_saveBC(s,bc,fdm,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
type(bc2d_sm),        intent(inout) :: bc
img%img=0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,sw,psw,sw0,psw0,&
        sw1,psw1,sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
! the last two time step
call copy_array(bc%p_nt,pp0,bc%p_nt_1,pp1,bc%sw_nt,psw0)
call copy_array(bc%sw_nt_1,psw1,bc%sh_nt,psh0,bc%sh_nt_1,psh1)
pswh0=psw0-psh0
pswh1=psw1-psh1
do it=nt-2,1,-1
   !-.-.-.-.-.-. Reconstruct P Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_wavefield(pp0,pswh,pswh1,beta_dt_w,fdm%m%refl_fs,nzw_pml,nx_pml,ic)
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
   call loadBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
           bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   !-.-.-.-.-.-. Reconstruct SW Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_source(psw0,beta_dt_w,fdm%s%sou,it)  
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,npml,nzwp,nxp,fd_order)
   call loadBC2d(psw,nzw,nx,nzwp,nxp,npml,bc%sw_top(:,:,it),bc%sw_bottom(:,:,it),&
           bc%sw_left(:,:,it),bc%sw_right(:,:,it),bc_len)
   !-.-.-.-.-.-. Reconstruct SH Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_source(psh0,beta_dt_h,fdm%s%sou,it)  
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,npml,nzwp,nxp,fd_order)
   call loadBC2d(psh,nzw,nx,nzp,nxp,npml,bc%sh_top(:,:,it),bc%sh_bottom(:,:,it),&
           bc%sh_left(:,:,it),bc%sh_right(:,:,it),bc_len)
   pswh = psw - psh
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-.-.-.-.-.-.-.-.-
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call add_seis(pq,beta_dt,s%ng,s%seis,it)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1,p,p0,p1,q,q0,q1)
end subroutine sm_rtm_saveBC

!=====================================================================

subroutine sm_rtm_saveWF(s,wf,fdm,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
type(wf2d),   intent(inout) :: wf
img%img=0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1)
call copy_array(wf%p(:,:,nt),pp0,wf%p(:,:,nt-1),pp1)
do it=nt-2,1,-1
   !-.-.-.-.-.-. Reconstruct P Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
   pp(npml:nzp+1,npml:nxp+1)=wf%p(:,:,it)
   !-.-.-.-.-.-. Back Propagate Receiver Pressure field .-.-.-.-.-
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,npml,nz_pml,nx_pml,fd_order)
   call add_seis(pq,beta_dt,s%ng,s%seis,it)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1)
end subroutine sm_rtm_saveWF

!=====================================================================

subroutine sm_illum(fdm,s,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
if (.not.allocated(img%img)) then
   call allocate_and_initial(img%illum,img%nz,img%nx)
else
   img%illum = 0.0
endif
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,sw,psw,sw0,psw0,sw1,psw1,&
        sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
call allocate_and_initial(illum_source,nz,nx)
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   ! ------ mute the direct wave
   pswh = psw - psh
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pswh,pswh1,pp,nzw_pml,nx_pml)
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      illum_source(iz-npml,ix-npml) = illum_source(iz-npml,ix-npml) &
                                    + (abs(pp(iz,ix)))**img%illum_order
   end forall
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp)
enddo
call deallocate_and_free(sw0,sw1,sw,sh0,sh1,sh,swh0,swh1,swh)
if (img%isIllumReceiver) then
   call initial(pp,pp0,pp1)
   call allocate_and_initial(illum_receiver,nz,nx)
   do it=1,nt
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
      do ig=1,s%ng
         pp(igz(ig),igx(ig))=pp(igz(ig),igx(ig))+beta_dt(igz(ig),igx(ig))*fdm%s%sou(it)
      enddo
      !call snapshot("q_",q,it,nz_pml,nx_pml)
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
end subroutine sm_illum

!=====================================================================
! receiver side multiples
subroutine rm_modeling(fdm,s)
type(water_fdmod2d), intent(inout) :: fdm
type(shot2d),        intent(inout) :: s
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,&
        rw1,prw1,rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
do it=1,nt
   ! ------ P , add source from WAVELET
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call add_source(pp,beta_dt,fdm%s%sou(it),isDipole)
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   !call snapshot("q_",pq,it,nz_pml,nx_pml)
   ! ------ RW in water, add source from Q ------
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pq,pq1,prw,nzw_pml,nx_pml)
   call add_FS(prw,npml,fs_thick)
   !call snapshot("rw_",prw,it,nzw_pml,nx_pml)
   ! ------ RH in water, add source from Q ------
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_h,pq,pq1,prh,nzw_pml,nx_pml)
   call add_FS(prh,npml,fs_thick)
   !call snapshot("rh_",prh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   !call snapshot("rwh_",prwh,it,nzw_pml,nx_pml)
   call store_data(prwh,s%ng,it,s%seis)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1,rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1)
end subroutine rm_modeling

!=====================================================================
subroutine rm_modeling_saveBC(fdm,bc,s)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),         intent(inout) :: s
type(bc2d_general),   intent(inout) :: bc
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,&
        rw1,prw1,rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
do it=1,nt
   ! ------ P , add source from WAVELET
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call add_source(pp,beta_dt,fdm%s%sou(it),isDipole)
   call saveBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),&
      bc%p_bottom(:,:,it),bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   !call snapshot("q_",pq,it,nz_pml,nx_pml)
   ! ------ RW in water, add source from Q ------
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pq,pq1,prw,nzw_pml,nx_pml)
   call add_FS(rw,npml,fs_thick)
   !call snapshot("rw_",prw,it,nzw_pml,nx_pml)
   ! ------ RH in water, add source from Q ------
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_h,pq,pq1,prh,nzw_pml,nx_pml)
   call add_FS(prh,npml,fs_thick)
   !call snapshot("rh_",prh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   !call snapshot("rwh_",rwh,it,nzw_pml,nx_pml)
   call store_data(prwh,s%ng,it,s%seis)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call copy_array(pp1,bc%p_nt,pp0,bc%p_nt_1)
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1,rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1)
end subroutine rm_modeling_saveBC

!=====================================================================
subroutine rm_modeling_saveWF(fdm,wf,s)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),         intent(inout) :: s
type(wf2d),           intent(inout) :: wf
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,&
        rw1,prw1,rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
do it=1,nt
   ! ------ P , add source from WAVELET
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call add_source(pp,beta_dt,fdm%s%sou(it),isDipole)
   wf%p(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   !call snapshot("q_",pq,it,nz_pml,nx_pml)
   ! ------ RW in water, add source from Q ------
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pq,pq1,prw,nzw_pml,nx_pml)
   call add_FS(prw,npml,fs_thick)
   !call snapshot("rw_",prw,it,nzw_pml,nx_pml)
   ! ------ RH in water, add source from Q ------
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_h,pq,pq1,prh,nzw_pml,nx_pml)
   call add_FS(prh,npml,fs_thick)
   !call snapshot("rh_",prh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   !call snapshot("rwh_",prwh,it,nzw_pml,nx_pml)
   call store_data(prwh,s%ng,it,s%seis)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,q,q0,q1,rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1)
end subroutine rm_modeling_saveWF

!=====================================================================

subroutine rm_rtm(s,fdm,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
img%img=0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,&
        rw1,prw1,rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
if (isSaveBC) then
   call allocate_and_initial(p_top,p_bottom,bc_len,nx,nt)
   call allocate_and_initial(p_left,p_right,nz,bc_len,nt)
   call allocate_and_initial(p_nt,p_nt_1,nz_pml,nx_pml)
endif
if (isSaveWF) then
   call allocate_and_initial(pwf,nz+2,nx+2,nt)
endif
! Forward propagate to save BC or WF
do it=1,nt
   ! ------ P , add source from WAVELET
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call add_source(pp,beta_dt,fdm%s%sou(it),isDipole)
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   if (isSaveBC) then
     call saveBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
             p_left(:,:,it),p_right(:,:,it),bc_len)
   endif
   if (isSaveWF) then
      pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   endif
   call wf_refresh(pp0,pp1,pp)
enddo

if (isSaveBC) then
   call copy_array(pp1,p_nt,pp0,p_nt_1)
   call copy_array(p_nt,pp0,p_nt_1,pp1)
endif
if (isSaveWF) then
   call initial(pp,pp0,pp1)
   pp0(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt)
   pp1(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt-1)
endif
! Back propagate
do it=nt-2,1,-1
   if (isSaveBC) then
      !-.-.-.-.-.-. Reconstruct P Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_source(pp0,beta_dt,fdm%s%sou,it,isDipole)
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
      call loadBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
              p_left(:,:,it),p_right(:,:,it),bc_len)
   endif
   if (isSaveWF) then
      pp(npml:nzp+1,npml:nxp+1)=pwf(:,:,it)
   endif
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   !-.-.-.-.-.-. Back propagate RW field in water, add source from SEISMOGRAM
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
   call add_seis(prw,beta_dt_w,s%ng,s%seis,it)
   call add_FS(prw,npml,fs_thick)
   !call snapshot("rw_",prw,it,nzw_pml,nx_pml)
   !-.-.-.-.-.-. Back propagate RH field in water, add source from SEISMOGRAM
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
   call add_seis(prh,beta_dt_h,s%ng,s%seis,it)
   call add_FS(prh,npml,fs_thick)
   !call snapshot("rh_",prh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   ! ------ Q r, add source from rwh ------
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,prwh,prwh1,pq,nzw_pml,nx_pml)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
if (isSaveBC) then
   call deallocate_and_free(p_top,p_bottom,p_left,p_right)
   call deallocate_and_free(p_nt,p_nt_1)
else
   call deallocate_and_free(pwf)
endif
call deallocate_and_free(rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1,p,p0,p1,q,q0,q1)
end subroutine rm_rtm

!=====================================================================

subroutine rm_rtm_saveBC(s,bc,fdm,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
type(bc2d_general),   intent(inout) :: bc
img%img=0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,&
        rw1,prw1,rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
call copy_array(bc%p_nt,pp0,bc%p_nt_1,pp1)
! Back propagate
do it=nt-2,1,-1
   !-.-.-.-.-.-. Reconstruct P Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_source(pp0,beta_dt,fdm%s%sou,it,isDipole)
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
   call loadBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
           bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   if (isSaveWF) then
      pp(npml:nzp+1,npml:nxp+1)=pwf(:,:,it)
   endif
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   !-.-.-.-.-.-. Back propagate RW field in water, add source from SEISMOGRAM
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
   call add_seis(prw,beta_dt_w,s%ng,s%seis,it)
   call add_FS(prw,npml,fs_thick)
   !call snapshot("rw_",prw,it,nzw_pml,nx_pml)
   !-.-.-.-.-.-. Back propagate RH field in water, add source from SEISMOGRAM
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
   call add_seis(prh,beta_dt_h,s%ng,s%seis,it)
   call add_FS(prh,npml,fs_thick)
   !call snapshot("rh_",prh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   ! ------ Q r, add source from rwh ------
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,prwh,prwh1,pq,nzw_pml,nx_pml)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1,p,p0,p1,q,q0,q1)
end subroutine rm_rtm_saveBC

!=====================================================================

subroutine rm_rtm_saveWF(s,wf,fdm,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
type(wf2d), intent(inout) :: wf
img%img=0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,&
        rw1,prw1,rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
pp0(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt)
pp1(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt-1)
! Back propagate
do it=nt-2,1,-1
   pp(npml:nzp+1,npml:nxp+1)=pwf(:,:,it)
   !call snapshot("p_",p,it,nz_pml,nx_pml)
   !-.-.-.-.-.-. Back propagate RW field in water, add source from SEISMOGRAM
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
   call add_seis(prw,beta_dt_w,s%ng,s%seis,it)
   call add_FS(prw,npml,fs_thick)
   !call snapshot("rw_",prw,it,nzw_pml,nx_pml)
   !-.-.-.-.-.-. Back propagate RH field in water, add source from SEISMOGRAM
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
   call add_seis(prh,beta_dt_h,s%ng,s%seis,it)
   call add_FS(prh,npml,fs_thick)
   !call snapshot("rh_",prh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   ! ------ Q r, add source from rwh ------
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,prwh,prwh1,pq,nzw_pml,nx_pml)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1,p,p0,p1,q,q0,q1)
end subroutine rm_rtm_saveWF

!=====================================================================
! receiver side multiples
subroutine rm_illum(fdm,s,img)
type(water_fdmod2d), intent(inout) :: fdm
type(shot2d),      intent(inout) :: s
type(image2d),     intent(inout) :: img
if (.not.allocated(img%img)) then
   call allocate_and_initial(img%illum,img%nz,img%nx)
else
   img%illum = 0.0
endif
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call allocate_and_initial(illum_source,nz,nx)
call wf_assign(p,pp,p0,pp0,p1,pp1)
do it=1,nt
   ! ------ P , add source from WAVELET
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call add_source(pp,beta_dt,fdm%s%sou(it),isDipole)
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      illum_source(iz-npml,ix-npml) = illum_source(iz-npml,ix-npml) &
                                    + (abs(pp(iz,ix)))**img%illum_order
   end forall
   call wf_refresh(pp0,pp1,pp)
enddo
call initial(pp,pp0,pp1)
if (img%isIllumReceiver) then
   call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
   call allocate_and_initial(rwh,rwh0,rwh1,nzw_pml,nx_pml)
   call wf_assign(rw,prw,rw0,prw0,rw1,prw1,rh,prh,rh0,prh0,rh1,prh1,&
           rwh,prwh,rwh0,prwh0,rwh1,prwh1)
   call allocate_and_initial(illum_receiver,nz,nx)
   do it=1,nt
      !-.-.-.-.-.-. Back propagate RW field in water, add source from SEISMOGRAM
      call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
      do ig=1,s%ng
         prw(igz(ig),igx(ig))=prw(igz(ig),igx(ig))+beta_dt_w(igz(ig),igx(ig))*fdm%s%sou(it)
      enddo
      call add_FS(prw,npml,fs_thick)
      !-.-.-.-.-.-. Back propagate RH field in water, add source from SEISMOGRAM
      call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
      do ig=1,s%ng
         prw(igz(ig),igx(ig)) = prw(igz(ig),igx(ig)) + beta_dt_w(igz(ig),igx(ig))*fdm%s%sou(it)
      enddo
      call add_FS(prh,npml,fs_thick)
      ! ------ mute the direct wave
      prwh = prw - prh
      ! ------ P, add source from rwh ------
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
      call pertubation(ic,fdm%m%refl_fs,beta_dt_w,prwh,prwh1,pp,nzw_pml,nx_pml)
      forall (iz=npml+1:nzp,ix=npml+1:nxp)
         illum_receiver(iz-npml,ix-npml) = illum_receiver(iz-npml,ix-npml) &
                                         + (abs(pp(iz,ix)))**img%illum_order
      end forall
      call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pp0,pp1,pp)
   enddo  ! End of time-marching loop
endif
if (img%isIllumReceiver) then
   img%illum=illum_source*illum_receiver
   call deallocate_and_free(rw,rw0,rw1,rh,rh0,rh1)
   call deallocate_and_free(illum_receiver)
else
   img%illum=illum_source
endif
call finalize_variables()
call deallocate_and_free(p,p0,p1,illum_source)
end subroutine rm_illum

!=====================================================================
! both side multiples
subroutine srm_modeling(fdm,s)
type(water_fdmod2d), intent(inout) :: fdm
type(shot2d),        intent(inout) :: s
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,sw,psw,sw0,psw0,sw1,psw1,&
        sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
call wf_assign(q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,rw1,prw1,&
        rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   !call snapshot("sw_",psw,it,nzw_pml,nx_pml)
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   !call snapshot("sh_",psh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   pswh = psw - psh
   !call snapshot("swh_",pswh,it,nzw_pml,nx_pml)
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pswh,pswh1,pp,nzw_pml,nx_pml)
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   !call snapshot("q_",pq,it,nz_pml,nx_pml)
   ! ------ RW reflected by FS, add source from Q
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pq,pq1,prw,nzw_pml,nx_pml)
   ! ------ RH reflected by FS, add source from Q
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_h,pq,pq1,prh,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   call store_data(prwh,s%ng,it,s%seis)
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
call deallocate_and_free(q,q0,q1,rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1)
end subroutine srm_modeling

!=====================================================================
! both side multiples
subroutine srm_modeling_saveBC(fdm,bc,s)
type(water_fdmod2d), intent(inout) :: fdm
type(shot2d),        intent(inout) :: s
type(bc2d_sm),       intent(inout) :: bc
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,sw,psw,sw0,psw0,sw1,psw1,&
        sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
call wf_assign(q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,rw1,prw1,&
        rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   call saveBC2d(psw,nzw,nx,nzwp,nxp,npml,bc%sw_top(:,:,it),&
      bc%sw_bottom(:,:,it),bc%sw_left(:,:,it),bc%sw_right(:,:,it),bc_len)
   !call snapshot("sw_",psw,it,nzw_pml,nx_pml)
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   call saveBC2d(psh,nzw,nx,nzwp,nxp,npml,bc%sh_top(:,:,it),&
      bc%sh_bottom(:,:,it),bc%sh_left(:,:,it),bc%sh_right(:,:,it),bc_len)
   !call snapshot("sh_",psh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   pswh = psw - psh
   !call snapshot("swh_",pswh,it,nzw_pml,nx_pml)
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pswh,pswh1,pp,nzw_pml,nx_pml)
   call saveBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),&
      bc%p_bottom(:,:,it),bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,p,p1,q,nz_pml,nx_pml)
   !call snapshot("q_",pq,it,nz_pml,nx_pml)
   ! ------ RW reflected by FS, add source from Q
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pq,pq1,prw,nzw_pml,nx_pml)
   ! ------ RH reflected by FS, add source from Q
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_h,pq,pq1,prh,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   call store_data(prwh,s%ng,it,s%seis)
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pq0,pq1,pq)
enddo  ! End of time-marching loop
call copy_array(pp1,bc%p_nt,psw1,bc%sw_nt,psh1,bc%sh_nt)
call copy_array(pp0,bc%p_nt_1,psw0,bc%sw_nt_1,psh0,bc%sh_nt_1)
call finalize_variables()
call deallocate_and_free(p,p0,p1,sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
call deallocate_and_free(q,q0,q1,rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1)
end subroutine srm_modeling_saveBC

!=====================================================================
! both side multiples
subroutine srm_modeling_saveWF(fdm,wf,s)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),         intent(inout) :: s
type(wf2d),   intent(inout) :: wf
s%seis = 0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,sw,psw,sw0,psw0,sw1,psw1,&
        sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
call wf_assign(q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,rw1,prw1,&
        rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   !call snapshot("sw_",psw,it,nzw_pml,nx_pml)
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   !call snapshot("sh_",psh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   pswh = psw - psh
   !call snapshot("swh_",pswh,it,nzw_pml,nx_pml)
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pswh,pswh1,pp,nzw_pml,nx_pml)
   wf%p(:,:,it) = pp
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   ! ------ Q reflected by model reflector, add source from P
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl,beta_dt,pp,pp1,pq,nz_pml,nx_pml)
   !call snapshot("q_",pq,it,nz_pml,nx_pml)
   ! ------ RW reflected by FS, add source from Q
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pq,pq1,prw,nzw_pml,nx_pml)
   ! ------ RH reflected by FS, add source from Q
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_h,pq,pq1,prh,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   call store_data(prwh,s%ng,it,s%seis)
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(p,p0,p1,sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
call deallocate_and_free(q,q0,q1,rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1)
end subroutine srm_modeling_saveWF

!=====================================================================

subroutine srm_rtm(s,fdm,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
img%img=0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,sw,psw,sw0,psw0,sw1,psw1,&
        sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
call wf_assign(q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,rw1,prw1,&
        rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
if (isSaveBC) then
   call allocate_and_initial(p_top,p_bottom,sw_top,sw_bottom,sh_top,sh_bottom,bc_len,nx,nt)
   call allocate_and_initial(p_left,p_right,nz,bc_len,nt)
   call allocate_and_initial(p_nt,p_nt_1,nz_pml,nx_pml)
   call allocate_and_initial(sw_left,sw_right,sh_left,sh_right,nz,bc_len,nt)
   call allocate_and_initial(sw_nt,sw_nt_1,sh_nt,sh_nt_1,nzw_pml,nx_pml)
endif
if (isSaveWF) then
   call allocate_and_initial(pwf,nz+2,nx+2,nt)
endif
! Forward propagate to save BC or WF
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   if (isSaveBC) then
      call saveBC2d(psw,nzw,nx,nzwp,nxp,npml,sw_top(:,:,it),sw_bottom(:,:,it),&
              sw_left(:,:,it),sw_right(:,:,it),bc_len)
   endif
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   !call snapshot("sh_f_",psh,it,nzw_pml,nx_pml)
   if (isSaveBC) then
      call saveBC2d(psh,nzw,nx,nzwp,nxp,npml,sh_top(:,:,it),sh_bottom(:,:,it),&
              sh_left(:,:,it),sh_right(:,:,it),bc_len)
   endif
   ! ------ mute the direct wave
   pswh = psw - psh
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt,pswh,pswh1,pp,nzw_pml,nx_pml)
   if (isSaveBC) then 
     call saveBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
             p_left(:,:,it),p_right(:,:,it),bc_len)
   endif
   if (isSaveWF) then
      pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   endif
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp)
enddo

if (isSaveBC) then
   call copy_array(pp1,p_nt,psw1,sw_nt,psh1,sh_nt)
   call copy_array(pp0,p_nt_1,psw0,sw_nt_1,psh0,sh_nt_1)
   call copy_array(p_nt,pp0,sw_nt,psw0,sh_nt,psh0)
   call copy_array(p_nt_1,pp1,sw_nt_1,psw1,sh_nt_1,psh1)
   call initial(psw,psh,pswh)
   pswh0 = psw0 - psh0
   pswh1 = psw1 - psh1
endif
if (isSaveWF) then
   call deallocate_and_free(sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
   call initial(pp0,pp1)
   pp0(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt)
   pp1(npml:nzp+1,npml:nxp+1)=pwf(:,:,nt-1)
endif

! Back propagate
pp=0.0
do it=nt-2,1,-1
   if (isSaveBC) then
      !-.-.-.-.-.-. Reconstruct P Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_wavefield(pp0,pswh,pswh1,beta_dt_w,fdm%m%refl_fs,nzw_pml,nx_pml,ic)
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
      call loadBC2d(pp,nz,nx,nzp,nxp,npml,p_top(:,:,it),p_bottom(:,:,it),&
              p_left(:,:,it),p_right(:,:,it),bc_len)
      !-.-.-.-.-.-. Reconstruct SW Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_source(psw0,beta_dt_w,fdm%s%sou,it)  
      call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,npml,nzwp,nxp,fd_order)
      call loadBC2d(psw,nzw,nx,nzwp,nxp,npml,sw_top(:,:,it),sw_bottom(:,:,it),&
              sw_left(:,:,it),sw_right(:,:,it),bc_len)
      !call snapshot("sw_",psw,it,nzw_pml,nx_pml)
      !-.-.-.-.-.-. Reconstruct SH Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
      call subtract_source(psh0,beta_dt_h,fdm%s%sou,it)  
      call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,npml,nzwp,nxp,fd_order)
      call loadBC2d(psh,nzw,nx,nzwp,nxp,npml,sh_top(:,:,it),sh_bottom(:,:,it),&
                  sh_left(:,:,it),sh_right(:,:,it),bc_len)
      !call snapshot("sh_",psh,it,nzw_pml,nx_pml)
      pswh = psw - psh
      !call snapshot("swh_",pswh,it,nzw_pml,nx_pml)
      call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh)
   endif
   if (isSaveWF) then
      pp(npml:nzp+1,npml:nxp+1)=pwf(:,:,it)
   endif
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   !-.-.-.-.-.-. Back Propagate RW field .-.-..-.-.-
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
   call add_seis(prw,beta_dt_w,s%ng,s%seis,it)
   !call snapshot("rw_",prw,it,nzw_pml,nx_pml)
   !-.-.-.-.-.-. Back Propagate RH field .-.-..-.-.-
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
   call add_seis(prh,beta_dt_h,s%ng,s%seis,it)
   !call snapshot("rh_",prh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   ! ------ Q, add source from rwh ------
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,prwh,prwh1,pq,nzw_pml,nx_pml)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pp0,pp1,pp,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
if (isSaveBC) then
   call deallocate_and_free(sw_top,sw_bottom,sw_left,sw_right)
   call deallocate_and_free(sh_top,sh_bottom,sh_left,sh_right)
   call deallocate_and_free(p_top,p_bottom,p_left,p_right)
   call deallocate_and_free(sw_nt,sw_nt_1,sh_nt,sh_nt_1,p_nt,p_nt_1)
   call deallocate_and_free(sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1)
else
   call deallocate_and_free(pwf)
endif
call deallocate_and_free(rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1,p,p0,p1,q,q0,q1)
end subroutine srm_rtm

!=====================================================================

subroutine srm_rtm_saveBC(s,bc,fdm,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
type(bc2d_sm),      intent(inout) :: bc
img%img=0.0
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,q,q0,q1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,rwh,rwh0,rwh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,sw,psw,sw0,psw0,sw1,psw1,&
        sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
call wf_assign(q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,rw1,prw1,&
        rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
call copy_array(bc%p_nt,pp0,bc%p_nt_1,pp1,bc%sw_nt,psw0)
call copy_array(bc%sw_nt_1,psw1,bc%sh_nt,psh0,bc%sh_nt_1,psh1)
pswh0 = psw0 - psh0
pswh1 = psw1 - psh1
do it=nt-2,1,-1
   !-.-.-.-.-.-. Reconstruct P Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_wavefield(pp0,pswh,pswh1,beta_dt_w,fdm%m%refl_fs,nzw_pml,nx_pml,ic)
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,npml,nzp,nxp,fd_order)
   call loadBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
           bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   !-.-.-.-.-.-. Reconstruct SW Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_source(psw0,beta_dt_w,fdm%s%sou,it)  
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,npml,nzwp,nxp,fd_order)
   call loadBC2d(psw,nzw,nx,nzwp,nxp,npml,bc%sw_top(:,:,it),bc%sw_bottom(:,:,it),&
           bc%sw_left(:,:,it),bc%sw_right(:,:,it),bc_len)
   !-.-.-.-.-.-. Reconstruct SH Pressure .-.-.-.-.-.-.-.-.-.-.-.-.-
   call subtract_source(psh0,beta_dt_h,fdm%s%sou,it)  
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,npml,nzwp,nxp,fd_order)
   call loadBC2d(psh,nzw,nx,nzp,nxp,npml,bc%sh_top(:,:,it),bc%sh_bottom(:,:,it),&
           bc%sh_left(:,:,it),bc%sh_right(:,:,it),bc_len)
   pswh = psw - psh
   !-.-.-.-.-.-. Back Propagate RW field .-.-..-.-.-
   call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
   call add_seis(prw,beta_dt_w,s%ng,s%seis,it)
   !call snapshot("rw_",prw,it,nzw_pml,nx_pml)
   !-.-.-.-.-.-. Back Propagate RH field .-.-..-.-.-
   call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
   call add_seis(prh,beta_dt_h,s%ng,s%seis,it)
   !call snapshot("rh_",prh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   prwh = prw - prh
   ! ------ Q, add source from rwh ------
   call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,prwh,prwh1,pq,nzw_pml,nx_pml)
   call image_condition(ic,img%img,nzp,nxp,npml,pp0,pp1,pq0)
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp)
   call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pq0,pq1,pq)
enddo  ! End of time-marching loop
call finalize_variables()
call deallocate_and_free(sw,sw0,sw1,sh,sh0,sh1,swh,swh0,swh1,p,p0,p1)
call deallocate_and_free(rw,rw0,rw1,rh,rh0,rh1,rwh,rwh0,rwh1,q,q0,q1)
end subroutine srm_rtm_saveBC

!=====================================================================

subroutine srm_illum(fdm,s,img)
type(water_fdmod2d),  intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(image2d),      intent(inout) :: img
if (.not.allocated(img%img)) then
   call allocate_and_initial(img%illum,img%nz,img%nx)
else
   img%illum = 0.0
endif
call initial_variables(fdm)
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call allocate_and_initial(sw,sw0,sw1,sh,sh0,sh1,nzw_pml,nx_pml)
call allocate_and_initial(swh,swh0,swh1,nzw_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1,sw,psw,sw0,psw0,sw1,psw1,&
        sh,psh,sh0,psh0,sh1,psh1,swh,pswh,swh0,pswh0,swh1,pswh1)
call allocate_and_initial(illum_source,nz,nx)
do it=1,nt
   ! ------ SW in water, add source from WAVELET ------
   call a2d_abc_kernal(psw,psw0,psw1,temp1_w,temp2_w,alpha_w,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psw,beta_dt_w,fdm%s%sou(it))
   call add_FS(psw,npml,fs_thick)
   !call snapshot("sw_",psw,it,nzw_pml,nx_pml)
   ! ------ SH in water, add source from WAVELET ------
   call a2d_abc_kernal(psh,psh0,psh1,temp1_h,temp2_h,alpha_h,nzw_pml,&
      nx_pml,fd_order)
   call add_source(psh,beta_dt_h,fdm%s%sou(it))
   call add_FS(psh,npml,fs_thick)
   !call snapshot("sh_",psh,it,nzw_pml,nx_pml)
   ! ------ mute the direct wave
   pswh = psw - psh
   !call snapshot("swh_",pswh,it,nzw_pml,nx_pml)
   ! ------ P reflected by FS, add source from SWH
   call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   call pertubation(ic,fdm%m%refl_fs,beta_dt_w,pswh,pswh1,pp,nzw_pml,nx_pml)
   !call snapshot("p_",pp,it,nz_pml,nx_pml)
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      illum_source(iz-npml,ix-npml) = illum_source(iz-npml,ix-npml) &
                                    + (abs(pp(iz,ix)))**img%illum_order
   end forall
   call wf_refresh(psw0,psw1,psw,psh0,psh1,psh,pswh0,pswh1,pswh,pp0,pp1,pp)
enddo
deallocate(sw0,sw1,sw,sh0,sh1,sh,swh0,swh1,swh,p0,p1,p)
if (img%isIllumReceiver) then
   call allocate_and_initial(q,q0,q1,nz_pml,nx_pml)
   call allocate_and_initial(rw,rw0,rw1,rh,rh0,rh1,nzw_pml,nx_pml)
   call allocate_and_initial(rwh,rwh0,rwh1,nzw_pml,nx_pml)
   call wf_assign(q,pq,q0,pq0,q1,pq1,rw,prw,rw0,prw0,rw1,prw1,&
           rh,prh,rh0,prh0,rh1,prh1,rwh,prwh,rwh0,prwh0,rwh1,prwh1)
   call allocate_and_initial(illum_receiver,nz,nx)
   do it=1,nt
      !-.-.-.-.-.-. Back propagate RW field in water, add source from SEISMOGRAM
      call a2d_abc_kernal(prw,prw0,prw1,temp1_w,temp2_w,alpha_w,nzw_pml,nx_pml,fd_order)
      do ig=1,s%ng
         prw(igz(ig),igx(ig))=prw(igz(ig),igx(ig))+beta_dt_w(igz(ig),igx(ig))*fdm%s%sou(it)
      enddo
      call add_FS(prw,npml,fs_thick)
      !-.-.-.-.-.-. Back propagate RH field in water, add source from SEISMOGRAM
      call a2d_abc_kernal(prh,prh0,prh1,temp1_h,temp2_h,alpha_h,nzw_pml,nx_pml,fd_order)
      do ig=1,s%ng
         prw(igz(ig),igx(ig)) = prw(igz(ig),igx(ig)) + beta_dt_w(igz(ig),igx(ig))*fdm%s%sou(it)
      enddo
      call add_FS(prh,npml,fs_thick)
      ! ------ mute the direct wave
      prwh = prw - prh
      ! ------ Q, add source from rwh ------
      call a2d_abc_kernal(pq,pq0,pq1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
      call pertubation(ic,fdm%m%refl_fs,beta_dt_w,prwh,prwh1,q,nzw_pml,nx_pml)
      forall (iz=npml+1:nzp,ix=npml+1:nxp)
         illum_receiver(iz-npml,ix-npml) = illum_receiver(iz-npml,ix-npml) &
                                         + (abs(q(iz,ix)))**img%illum_order
      end forall
      call wf_refresh(prw0,prw1,prw,prh0,prh1,prh,prwh0,prwh1,prwh,pq0,pq1,pq)
   enddo ! End of time-marching loop
endif
if (img%isIllumReceiver) then
   img%illum=illum_source*illum_receiver
   call deallocate_and_free(rw0,rw1,rw,rh0,rh1,rh,rwh0,rwh1,rwh,q0,q1,q)
   call deallocate_and_free(illum_receiver)
else
   img%illum=illum_source
endif
call finalize_variables()
call deallocate_and_free(illum_source)
end subroutine srm_illum

!=====================================================================

subroutine setup_source_receiver(s)
type(shot2d),       intent(in) :: s
allocate(igx(s%ng),igz(s%ng))
! Setup source position
isx = npml+int(s%xs/dx)+1
isz = npml+int(s%zs/dx)+1
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

subroutine initial_variables_mod(fdm)
type(water_fdmod2d),intent(in) :: fdm
call copy_variable(fdm%f%fd_type,fd_type,fdm%f%fd_order,fd_order)
call copy_variable(fdm%f%ic,ic,fdm%f%bc_type,bc_type)
call copy_variable(fdm%f%nt,nt,fdm%m%nz,nz,fdm%m%nzw,nzw,fdm%m%nx,nx,fdm%f%npml,npml)
call copy_variable(fdm%f%dt,dt,fdm%m%dx,dx)
call copy_variable(fdm%f%isDipole,isDipole)
call i_v_kernal(fdm%m%v,fdm%m%v_w,fdm%m%v_h,fdm%f%fs)
end subroutine initial_variables_mod

!=====================================================================

subroutine finalize_variables()
call deallocate_and_free(igx,igz,fs) 
call deallocate_and_free(beta_dt,beta_dt_w,beta_dt_h)
if (bc_type==1) call deallocate_and_free(alpha,temp1,temp2,&
   alpha_w,temp1_w,temp2_w,alpha_h,temp1_h,temp2_h)
if (bc_type==2) then 
   call deallocate_and_free(beta_dtdx,beta_dtdx_w,beta_dtdx_h)
   call deallocate_and_free(u1,u2)
endif
end subroutine finalize_variables

!=====================================================================

subroutine nx_nz_pml(nz,nzw,nx,npml,nz_pml,nzw_pml,nx_pml,nzp,nzwp,nxp)
integer, intent(in)  :: nz,nzw,nx,npml
integer, intent(out) :: nz_pml,nzw_pml,nx_pml,nzp,nzwp,nxp
nx_pml = nx + 2 * npml
nz_pml = nz + 2 * npml
nzw_pml = nzw + 2 * npml
nxp = nx + npml
nzp = nz + npml
nzwp = nzw + npml
end subroutine nx_nz_pml

!=====================================================================

subroutine i_v_kernal(v,v_w,v_h,fs_in)
real,intent(in) :: v(:,:),v_w(:,:),v_h(:,:)
integer, intent(in) :: fs_in(:)
call nx_nz_pml(nz,nzw,nx,npml,nz_pml,nzw_pml,nx_pml,nzp,nzwp,nxp)
bc_len=(fd_order-20)/2+1
fs_thick=bc_len-1
dtdx = dt/dx
call allocate_and_initial(beta_dt,nz_pml,nx_pml)
call allocate_and_initial(beta_dt_w,beta_dt_h,nzw_pml,nx_pml)
beta_dt=(v*dt)**2.0
beta_dt_w=(v_w*dt)**2.0
beta_dt_h=(v_h*dt)**2.0
call allocate_and_initial(fs,nx_pml)
fs = fs_in
if (fd_type==1) then ! Second Order FD
   if (bc_type==1) then ! Absorbing Boundary Condition
     ! for v
      call allocate_and_initial(damp,alpha,kappa,temp1,temp2,nz_pml,nx_pml)
      vmin=minval(v)
      call abc_get_damp2d(nx,nz,npml,dx,vmin,damp)
      alpha=(v*dtdx)**2.0
      kappa=damp*dt
      ! for v_w
      call allocate_and_initial(damp_w,alpha_w,kappa_w,temp1_w,temp2_w,nzw_pml,nx_pml)
      vmin=minval(v_w)
      call abc_get_damp2d(nx,nzw,npml,dx,vmin,damp_w)
      alpha_w=(v_w*dtdx)**2.0
      kappa_w=damp_w*dt
      ! for v_h
      call allocate_and_initial(damp_h,alpha_h,kappa_h,temp1_h,temp2_h,nzw_pml,nx_pml)
      vmin=minval(v_h)
      call abc_get_damp2d(nx,nzw,npml,dx,vmin,damp_h)
      alpha_h=(v_h*dtdx)**2.0
      kappa_h=damp_h*dt
      if (fd_order==22) then
         temp1  =2.0+2.0*C21*alpha  -kappa
         temp1_w=2.0+2.0*C21*alpha_w-kappa_w
         temp1_h=2.0+2.0*C21*alpha_h-kappa_h
      elseif (fd_order==24) then
         temp1  =2.0+2.0*C41*alpha  -kappa 
         temp1_w=2.0+2.0*C41*alpha_w-kappa_w
         temp1_h=2.0+2.0*C41*alpha_h-kappa_h
      elseif (fd_order==26) then
         temp1  =2.0+2.0*C61*alpha  -kappa  
         temp1_w=2.0+2.0*C61*alpha_w-kappa_w
         temp1_h=2.0+2.0*C61*alpha_h-kappa_h
      elseif (fd_order==28) then
         temp1  =2.0+2.0*C81*alpha  -kappa
         temp1_w=2.0+2.0*C81*alpha_w-kappa_w
         temp1_h=2.0+2.0*C81*alpha_h-kappa_h
      endif
      temp2  =1.0-kappa
      temp2_w=1.0-kappa_w
      temp2_h=1.0-kappa_h
      call deallocate_and_free(damp,kappa,damp_w,kappa_w,damp_h,kappa_h)
   elseif (bc_type==2) then ! Sponge Boundary Condition
      call allocate_and_initial(beta_dtdx,u1,nz_pml,nx_pml)
      call allocate_and_initial(beta_dtdx_w,beta_dtdx_h,u2,nzw_pml,nx_pml)
      call allocate_and_initial(spgx,spcz1,spcz2,nx_pml)
      call allocate_and_initial(spgz,spcx1,spcx2,nz_pml)
      call allocate_and_initial(spgxw,spcz1w,spcz2w,spgxh,spcz1h,spcz2h,nx_pml)
      call allocate_and_initial(spgzw,spcx1w,spcx2w,spgzh,spcx1h,spcx2h,nzw_pml)
      call spgcoef(spgx,spgz,spcx1,spcx2,spcz1,spcz2,v,&
                   nx_pml,nz_pml,npml,dx,dx,dt)
      call spgcoef(spgxw,spgzw,spcx1w,spcx2w,spcz1w,spcz2w,v_w,&
                   nx_pml,nzw_pml,npml,dx,dx,dt)
      call spgcoef(spgxh,spgzh,spcx1h,spcx2h,spcz1h,spcz2h,v_h,&
                   nx_pml,nzw_pml,npml,dx,dx,dt)
      beta_dtdx  =(v*dtdx)**2.0
      beta_dtdx_w=(v_w*dtdx)**2.0
      beta_dtdx_h=(v_h*dtdx)**2.0
!      call write_binfile("beta_dtdx.bin",beta_dtdx,nz_pml,nx_pml)
!      call write_binfile("vel.bin",v,nz_pml,nx_pml)
      if (fd_order==22) then
         c1t2=C21*2.0
      elseif (fd_order==24) then
         c1t2=C41*2.0
      elseif (fd_order==26) then
         c1t2=C61*2.0
      elseif (fd_order==28) then
         c1t2=C81*2.0
      endif
   endif
elseif (fd_type==2) then ! Staggerred FD
   write(*,*)"Under construction!"
endif
end subroutine i_v_kernal

end module module_a2d_mul
