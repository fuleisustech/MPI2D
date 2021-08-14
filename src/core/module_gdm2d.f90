module module_gdm2d

use module_global
use module_datatype
use module_array
use module_utility
use module_boundary2d
use module_fd2d
use module_string
use module_io
use module_fftw

implicit none
character(len=100),private :: fd_type,sort_in
integer,private :: ix,iz,it,isx,isz,ig,nx,nz,nx1,nz1,nx2,nz2,npml,nxp,nzp,nt,ntw_conv,&
   nt_conv,dnt_mod,nx_pml,nz_pml,ic,fd_order,bc_type,bc_len,fs_thick,nxp1,nzp1
real, private   :: dtdx,dx,dz,dt,vmin,c1t2
integer,private,allocatable :: igx(:),igz(:),fs(:),iz0(:)
logical,private :: isFS, isDipole, isSaveBC, isSaveWF, isWriteData, isReadData
real,private,target,allocatable :: p0(:,:),p1(:,:),p(:,:),&
   q0(:,:),q1(:,:),q(:,:)
real,private,pointer :: pp0(:,:),pp1(:,:),pp(:,:),pq0(:,:),pq1(:,:),pq(:,:)
real,private,allocatable :: beta_dt(:,:),gf_s(:),gf_g(:),gf_sg(:),gf_sg_xcorr(:)
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
integer(kind=MPI_OFFSET_KIND) :: disp
integer :: icheck
character(len=slen),private :: input_temp2, input, output_temp2, output

interface gdm2d_gfs
   module procedure gfs_cal
end interface gdm2d_gfs

interface gdm2d_mod
   module procedure linear_modeling
end interface gdm2d_mod

interface gdm2d_mig
   module procedure migration
end interface gdm2d_mig

interface add_source
   module procedure add_source_mod
   module procedure add_source_mig
end interface add_source

private initial_variables_gfs,finalize_variables_gfs,&
   initial_variables_mod,finalize_variables_mod,nx_nz_pml,i_v_kernal,&
   a2d_abc_kernal,add_source,subtract_source,store_data,&
   add_FS,pertubation,setup_source_receiver,saveBC2d,loadBC2d

contains

!=====================================================================

subroutine add_source_mod(p,beta_dt,s,isFS,isDipole)
logical, intent(in) :: isFS, isDipole
real,    intent(in) :: beta_dt(:,:),s
real,    intent(inout) :: p(:,:)
p(isz,isx) = p(isz,isx) + beta_dt(isz,isx)*s
if (.not.isFS) then
   if (isDipole) then
      p(isz-2,isx) = p(isz-2,isx) - beta_dt(isz-2,isx)*s
   endif
endif
end subroutine add_source_mod

!=====================================================================

subroutine add_source_mig(p,beta_dt,s,isDipole)
logical, intent(in) :: isDipole
real,    intent(in) :: beta_dt(:,:),s
real,    intent(inout) :: p(:,:)
p(isz,isx) = p(isz,isx) + beta_dt(isz,isx)*s
if (isDipole) then
   p(isz-2,isx) = p(isz-2,isx) - beta_dt(isz-2,isx)*s
endif
end subroutine add_source_mig

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

subroutine add_FS(p,fs_thick)
integer, intent(in) :: fs_thick
real,    intent(inout) :: p(:,:)
integer :: iz
p(nz1,:)=0.0
do iz=1,fs_thick
   p(nz1-iz,:)=-p(nz1+iz,:)
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
   forall (iz=nz1:nzp,ix=nx1:nxp)
      img(iz-npml,ix-npml) = img(iz-npml,ix-npml) &
                                 + p1(iz,ix) * q0(iz,ix)
   endforall
elseif (ic.eq.1) then
   forall (iz=nz1:nzp,ix=nx1:nxp)
      img(iz-npml,ix-npml) = img(iz-npml,ix-npml) + &
            (p0(iz,ix)- p1(iz,ix)) * q0(iz,ix)
   endforall
elseif (ic.eq.2) then
   forall (iz=nz1:nzp,ix=nx1:nxp)
      img(iz-npml,ix-npml) = img(iz-npml,ix-npml) &
            + (p1(iz-1,ix) + p1(iz+1,ix) + p1(iz,ix-1) &
             + p1(iz,ix+1)-4.0 * p1(iz,ix)) * q0(iz,ix)
   endforall
endif
end subroutine image_condition

!=====================================================================

subroutine gfs_cal(gdm,s,output_temp1)
type(gdmod2d),      intent(inout) :: gdm
type(shot2d),       intent(inout) :: s
character(len=*),   intent(in)    :: output_temp1
integer :: it_out,it,ix,iz,ntmp1,ntmp2,ntmp
!call initial_variables_gfs(gdm,s)
call initial_variables_gfs(gdm)
!if (isWriteData) then
!   call MPI_FILE_OPEN(MPI_COMM_SELF, output, &
!         MPI_MODE_WRONLY + MPI_MODE_CREATE, &
!         MPI_INFO_NULL, fid, ierr)
!else
!   s%seis = 0.0
!endif
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1)
call filename(output, output_temp1,s%sid(1),".bin")
open(11,file=output,access="direct",recl=nz*I4)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
do it=1,nt
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
   endif 
   call add_source(pp,beta_dt,gdm%s%sou(it),isFS,isDipole)
   if (isFS) then
      call add_FS(pp,fs_thick)
   endif
   !if (mod(it,1000)==1) call snapshot("p_",p,it,nz_pml,nx_pml)
   !call write_binfile(output,pp(nz1:nzp,nx1:nxp),nz,nx)  
   if (dnt_mod/=1) then
      if (mod(it,dnt_mod)==1) then
         it_out=(it-1)/dnt_mod+1
         call message(0,"it_out",it_out)
         write(11,rec=it_out)pp(nzp1:nzp,nx1:nxp)
      endif
   else
      ntmp1=nx*(it-1)
      do ix=nx1,nxp
         ntmp2=ix-nx1
         ntmp=ntmp1+ntmp2
         write(11,rec=ntmp+1)(pp(iz,ix),iz=nz1,nzp)
      enddo
   endif
   !call write_binfile(output,pp(npml:nz1,npml:nx1),nz2,nx2)  
   !if (isWriteData) then
   !   if (dnt_mod.ne.1) then
   !      if (mod(it,dnt_mod).eq.1) then
   !         it_out=(it-1)/dnt_mod+1
   !         do ig=1,s%ng
   !            disp=((ig-1)*s%nt+(it_out-1))*4
   !            call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
   !            call MPI_FILE_WRITE(fid, pp(igz(ig),igx(ig)),1, MPI_REAL, &
   !                      MPI_STATUS_IGNORE, ierr)
   !         enddo
   !      endif
   !   else
   !      do ig=1,s%ng
   !         disp=((ig-1)*s%nt+(it-1))*4
   !         call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
   !         call MPI_FILE_WRITE(fid, pp(igz(ig),igx(ig)),1, MPI_REAL, &
   !                   MPI_STATUS_IGNORE, ierr)
   !      enddo
   !   endif
   !else
   !   call store_data(pp,s%ng,it,s%seis)
   !endif
   call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
close(11)
call finalize_variables_gfs()
call deallocate_and_free(p,p0,p1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine gfs_cal

!=====================================================================

subroutine check_time_now()
call date_and_time(date,time_now)
write(*,*)icheck,time_now
!call flush(6)
icheck=icheck+1
end subroutine check_time_now

!=====================================================================

subroutine linear_modeling(gdm,s,input_temp1)
type(gdmod2d),       intent(in)    :: gdm
type(shot2d),        intent(inout) :: s
character(len=*),    intent(in)    :: input_temp1
integer                            :: fid, itrace, ix, iz
call initial_variables_mod(gdm,s)
s%seis = 0.0
do itrace=1,s%ng
   do ix=1,nx
      do iz=iz0(ix),nz
         disp=((ix-1)*nz+iz-1)*nt*4
         call filename(input,input_temp1,s%sid(itrace),".bin")
         call MPI_FILE_OPEN(MPI_COMM_SELF, input, MPI_MODE_RDONLY, &
                            MPI_INFO_NULL, fid, ierr)
         call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
         call MPI_FILE_READ(fid,gf_s,nt,MPI_REAL,MPI_STATUS_IGNORE, ierr)
         call MPI_FILE_CLOSE(fid,ierr)
         call filename(input,input_temp1,s%gid(itrace),".bin")
         call MPI_FILE_OPEN(MPI_COMM_SELF, input, MPI_MODE_RDONLY, &
                            MPI_INFO_NULL, fid, ierr)
         call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
         call MPI_FILE_READ(fid,gf_g,nt,MPI_REAL,MPI_STATUS_IGNORE, ierr)
         call MPI_FILE_CLOSE(fid,ierr)
         call conv_gf(gf_s,gf_g,gf_sg,nt,nt_conv)
         call xcorr_wavelet(gf_sg,gf_sg_xcorr,gdm%g%wfft,nt,ntw_conv)
         s%seis(:,itrace)=s%seis(:,itrace)+gdm%m%refl(iz,ix)*gf_sg_xcorr
      enddo
   enddo
enddo
call finalize_variables_mod()
end subroutine linear_modeling

!=====================================================================

!subroutine migration(gdm,s,img,input_temp1,output)
subroutine migration(gdm,s,img,input_temp1)
type(gdmod2d),       intent(in)    :: gdm
type(shot2d),        intent(in)    :: s
type(image2d),       intent(inout) :: img
character(len=*),    intent(in)    :: input_temp1!,output
integer                            :: fid, itrace, ix, iz
call initial_variables_mod(gdm,s)
img%img = 0.0
if (img%isIllum) then
   img%illum=0.0
endif
do itrace=1,s%ng
   do ix=1,nx
      do iz=iz0(ix),nz
         disp=((ix-1)*nz+iz-1)*nt*4
         call filename(input,input_temp1,s%sid(itrace),".bin")
         call MPI_FILE_OPEN(MPI_COMM_SELF, input, MPI_MODE_RDONLY, &
                            MPI_INFO_NULL, fid, ierr)
         call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
         call MPI_FILE_READ(fid,gf_s,nt,MPI_REAL,MPI_STATUS_IGNORE, ierr)
         call MPI_FILE_CLOSE(fid,ierr)
         call filename(input,input_temp1,s%gid(itrace),".bin")
         call MPI_FILE_OPEN(MPI_COMM_SELF, input, MPI_MODE_RDONLY, &
                            MPI_INFO_NULL, fid, ierr)
         call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
         call MPI_FILE_READ(fid,gf_g,nt,MPI_REAL,MPI_STATUS_IGNORE, ierr)
         call MPI_FILE_CLOSE(fid,ierr)
         call conv_gf(gf_s,gf_g,gf_sg,nt,nt_conv)
         call xcorr_wavelet(gf_sg,gf_sg_xcorr,gdm%g%wfft,nt,ntw_conv)
         img%img(iz,ix)=img%img(iz,ix)+sum(gf_sg_xcorr*s%seis(:,itrace))
         if (img%isIllum) then
            img%illum(iz,ix)=img%illum(iz,ix)+sqrt(sum(gf_sg_xcorr*gf_sg_xcorr))
         endif
      enddo
   enddo
enddo
if (img%isIllum) then
   img%illum=merge(1.0,img%illum,img%illum<0.000000001)
endif
call finalize_variables_mod()
end subroutine migration

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

!subroutine initial_variables_gfs(gdm,s)
subroutine initial_variables_gfs(gdm)
type(gdmod2d),intent(in) :: gdm
!type(shot2d), intent(in) :: s
call copy_variable(gdm%f%fd_type,fd_type,gdm%g%sort_in,sort_in)
call copy_variable(gdm%f%fd_order,fd_order,gdm%f%bc_type,bc_type)
call copy_variable(gdm%f%nt,nt,gdm%m%nz,nz,gdm%m%nx,nx,gdm%f%npml,npml)
call copy_variable(gdm%f%dt,dt,gdm%m%dx,dx)
call copy_variable(gdm%f%isFS,isFS,gdm%f%isDipole,isDipole,&
                   gdm%f%isWriteData,isWriteData)
call copy_variable(gdm%f%isSaveBC,isSaveBC,gdm%f%isSaveWF,isSaveWF)
!write(*,*) 'g%dt=====================',gdm%g%dt
!write(*,*) 'f%dt=====================',gdm%f%dt
dnt_mod=nint(gdm%g%dt/gdm%f%dt)
call i_v_kernal(gdm%m%v,gdm%f%fs)
end subroutine initial_variables_gfs

!=====================================================================

subroutine initial_variables_mod(gdm,s)
type(gdmod2d),intent(in) :: gdm
type(shot2d), intent(in) :: s
call copy_variable(s%nt,nt,gdm%m%nz,nz,gdm%m%nx,nx)
call copy_variable(gdm%g%ntw_conv,ntw_conv,gdm%g%nt_conv,nt_conv)
call copy_variable(s%dt,dt,gdm%m%dx,dx,gdm%m%dz,dz)
call copy_variable(gdm%f%isWriteData,isWriteData)
call allocate_and_initial(gf_s,gf_g,gf_sg,gf_sg_xcorr,nt)
call allocate_and_initial(iz0,nx)
iz0=nint(gdm%g%z0_mig/dx)+1
end subroutine initial_variables_mod

subroutine finalize_variables_gfs()
if (bc_type==1) then
      call deallocate_and_free(igx,igz,fs)
      call deallocate_and_free(alpha,temp1,temp2,beta_dt)
endif
if (bc_type==2) then 
      call deallocate_and_free(igx,igz,fs)
      call deallocate_and_free(beta_dtdx,beta_dt,u1)
endif
end subroutine finalize_variables_gfs

subroutine finalize_variables_mod
deallocate(iz0,gf_s,gf_g,gf_sg,gf_sg_xcorr)
end subroutine finalize_variables_mod

!=====================================================================

subroutine nx_nz_pml
nx_pml = nx + 2 * npml
nz_pml = nz + 2 * npml
nx1 = npml + 1
nz1 = npml + 1
nx2 = nx + 2
nz2 = nz + 2
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

end module module_gdm2d
