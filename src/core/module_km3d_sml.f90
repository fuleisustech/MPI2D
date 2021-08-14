module module_km3d
use module_global, only : date,time_now
use module_datatype
use module_utility
use module_eikfd3d_dir
use module_eikfd3d_fir
use module_fftw
use module_io
use module_sf90_mpi

implicit none
real, allocatable, private :: ttt_s(:,:,:),ttt_g(:,:,:),trace(:),trace_temp(:)
integer, allocatable, private :: iz0(:,:)
character(len=100),private :: eik_type, sort_in
logical, private :: isTheta, isDip, isSaveTTT, isSaveBothTTT
real, private :: dx, dy,dz,dt, theta, dip, theta2, dip2
integer, private :: nx, ny, nz, nt, nt_conv, ng

private :: initial_variable, finalize_variable,&
   read_ttt, calculate_ttt, mod_kernal, mig_kernal, &
   mod_kernal_control, mig_kernal_control      

contains

!------------------------------------------------------------------------

subroutine km3d_mod(kmm,s,sttt_in,gttt_in)
type(kmmod3d), intent(in)  :: kmm
type(shot3d),  intent(inout) :: s
character(len=100), intent(in)  :: sttt_in,gttt_in
integer                    :: itrace
s%seis=0.0
call initial_variable(kmm,s)
if (sort_in.eq."CSG") then
   if (isSaveTTT) then
      if(isSaveBothTTT) then
         call read_ttt(sttt_in,s%sid(1),ttt_s)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%ys(1),s%xs(1))
      endif
   else
      call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%ys(1),s%xs(1))
   endif
   do itrace=1,s%ng
      if (isSaveTTT) then
         call read_ttt(gttt_in,s%gid(itrace),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(itrace),s%yg(itrace),s%xg(itrace))
      endif
      if (isTheta) then
         call mod_kernal_control(trace,kmm%m%refl,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      else
         !call mod_kernal(trace,kmm%m%refl,kmm%k%wfft)
         call mod_kernal(trace,kmm%m%refl,kmm%k%wfft,s%xs(1),s%ys(1),s%zs(1),&
                  s%xg(itrace),s%yg(itrace),s%zg(itrace))
      endif
      s%seis(:,itrace)=trace
   enddo
elseif (sort_in.eq."CRG") then
   if (isSaveTTT) then
      if(isSaveBothTTT) then
         call read_ttt(gttt_in,s%gid(1),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(1),s%yg(1),s%xg(1))
      endif
   else
      call calculate_ttt(kmm%m%v,ttt_g,s%zg(1),s%yg(1),s%xg(1))
   endif
   do itrace=1,s%ng
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(itrace),ttt_s)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(itrace),s%ys(itrace),s%xs(itrace))
      endif
      if (isTheta) then
         call mod_kernal_control(trace,kmm%m%refl,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      else
         !call mod_kernal(trace,kmm%m%refl,kmm%k%wfft)
         call mod_kernal(trace,kmm%m%refl,kmm%k%wfft,s%xs(itrace),s%ys(itrace),&
                  s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      endif
      s%seis(:,itrace)=trace
   enddo
else
   do itrace=1,s%ng
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(itrace),ttt_s)
         call read_ttt(gttt_in,s%gid(itrace),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(itrace),s%ys(itrace),s%xs(itrace))
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(itrace),s%yg(itrace),s%xg(itrace))
      endif
      if (isTheta) then
         call mod_kernal_control(trace,kmm%m%refl,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      else
         !call mod_kernal(trace,kmm%m%refl,kmm%k%wfft)
         call mod_kernal(trace,kmm%m%refl,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      endif
      s%seis(:,itrace)=trace
   enddo
endif
call finalize_variable
end subroutine km3d_mod

!------------------------------------------------------------------------

subroutine km3d_mig(kmm,s,img,sttt_in,gttt_in)
type(kmmod3d), intent(in)    :: kmm
type(shot3d),  intent(in)    :: s
type(image3d), intent(inout) :: img 
character(len=100), intent(in)  :: sttt_in,gttt_in
integer                      :: itrace
call initial_variable(kmm,s)
img%img=0.0
if (img%isIllum) then
   img%illum=0.0
endif
if (sort_in.eq."CSG") then
   if (isSaveTTT) then
      if(isSaveBothTTT) then
         call read_ttt(sttt_in,s%sid(1),ttt_s)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%ys(1),s%xs(1))
      endif
   else
      call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%ys(1),s%xs(1))
   endif
   do itrace=1,s%ng
      trace=s%seis(:,itrace)
      if (isSaveTTT) then
         call read_ttt(gttt_in,s%gid(itrace),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(itrace),s%yg(itrace),s%xg(itrace))
      endif
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,trace,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      else
         call mig_kernal(img%img,img%isIllum,img%illum,trace,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      endif
   enddo
elseif (sort_in.eq."CRG") then
   if (isSaveTTT) then
      if(isSaveBothTTT) then
         call read_ttt(gttt_in,s%gid(1),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(1),s%yg(1),s%xg(1))
      endif
   else
      call calculate_ttt(kmm%m%v,ttt_g,s%zg(1),s%yg(1),s%xg(1))
   endif
   do itrace=1,s%ng
      trace=s%seis(:,itrace)
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(itrace),ttt_s)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(itrace),s%ys(itrace),s%xs(itrace))
      endif
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,trace,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      else
         call mig_kernal(img%img,img%isIllum,img%illum,trace,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      endif
   enddo
else
   do itrace=1,s%ng
      trace=s%seis(:,itrace)
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(itrace),ttt_s)
         call read_ttt(gttt_in,s%gid(itrace),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(itrace),s%ys(itrace),s%xs(itrace))
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(itrace),s%yg(itrace),s%xg(itrace))
      endif
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,trace,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      else
         call mig_kernal(img%img,img%isIllum,img%illum,trace,kmm%k%wfft,s%xs(itrace),&
            s%ys(itrace),s%zs(itrace),s%xg(itrace),s%yg(itrace),s%zg(itrace))
      endif
   enddo
endif
if (img%isIllum) then
   img%illum=merge(1.0,img%illum,img%illum<0.000000001)
endif
call finalize_variable
end subroutine km3d_mig

!------------------------------------------------------------------------

subroutine initial_variable(kmm,s)
type(kmmod3d),  intent(in) :: kmm
type(shot3d),   intent(in) :: s
call copy_variable(kmm%k%nt,nt,kmm%k%nt_conv,nt_conv) 
call copy_variable(kmm%k%isTheta,isTheta,kmm%k%isDip,isDip)
call copy_variable(kmm%k%isSaveTTT,isSaveTTT,kmm%k%isSaveBothTTT,isSaveBothTTT)
call copy_variable(kmm%k%eik_type,eik_type,kmm%k%sort_in,sort_in)
call copy_variable(kmm%k%dt,dt,kmm%k%theta,theta,kmm%k%dip,dip)
call copy_variable(kmm%m%dx,dx,kmm%m%dy,dy,kmm%m%dz,dz)
call copy_variable(kmm%m%nx,nx,kmm%m%ny,ny,kmm%m%nz,nz,s%ng,ng)
call allocate_and_initial(ttt_s,ttt_g,nz,ny,nx)
call allocate_and_initial(trace,nt)
call allocate_and_initial(trace_temp,nt)
call allocate_and_initial(iz0,ny,nx)
iz0=nint(kmm%k%z0_mig/dx)+1
theta2=(sin(theta/180.0*PI))**2.0
dip2=(sin(dip/180.0*PI))**2.0
!call copy_variable(s%xs,xs,s%zs,zs)
end subroutine initial_variable

subroutine finalize_variable
call deallocate_and_free(ttt_s,ttt_g)
call deallocate_and_free(trace)
call deallocate_and_free(trace_temp)
call deallocate_and_free(iz0)
end subroutine finalize_variable

subroutine read_ttt(ttt_in,id,ttt)
integer, intent(in) :: id
character(len=100),intent(in) :: ttt_in
real, intent(out) :: ttt(:,:,:)
character(len=200) :: ttt_name
integer :: fid
call filename(ttt_name, ttt_in, id, ".bin")
call MPI_FILE_OPEN(MPI_COMM_SELF, ttt_name, &
                   MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,ttt,nz*ny*nx,MPI_REAL,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine read_ttt

subroutine calculate_ttt(s,ttt,zs,ys,xs)
real, intent(in) :: s(:,:,:),zs,ys,xs
real, intent(out):: ttt(:,:,:)
if (eik_type.eq."DIR") then
   call eikfd3d_dir(s,ttt,zs,ys,xs,nz,ny,nx,dx)
else
   call eikfd3d_fir(s,ttt,zs,ys,xs,nz,ny,nx,dx,dx,dx)
endif
end subroutine calculate_ttt

subroutine mod_kernal(trace,refl,wfft,sx,sy,sz,gx,gy,gz)
real, intent(in)    :: refl(:,:,:),sx,sy,sz,gx,gy,gz
real, intent(out)   :: trace(:)
!real,allocatable :: trace_temp(:)
complex, intent(in) :: wfft(:)
integer :: ix, iy, iz, it0, it1
real    :: x,y,z,xsx2,xgx2,ysy2,ygy2,zsz2,zgz2,xsx2ysy2,xgx2ygy2,&
   dist_s,dist_g,oblq_s,oblq_g
trace_temp=0.0
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,x,xsx2,xgx2,y,ysy2,ygy2,xsx2ysy2,&
!$OMP    it0,it1,z,zsz2,zgz2,dist_s,dist_g,oblq_s,oblq_g)
do ix=1,nx
   x=(ix-1)*dx
   xsx2=(x-sx)**2.0
   xgx2=(x-gx)**2.0
   do iy=1,ny
      y=(iy-1)*dx
      ysy2=(y-sy)**2.0
      ygy2=(y-gy)**2.0
      xsx2ysy2=xsx2+ysy2
      xgx2ygy2=xgx2+ygy2
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            z=(iz-1)*dx
            zsz2=(z-sz)**2.0
            zgz2=(z-gz)**2.0
            dist_s=max(sqrt(xsx2ysy2+zsz2),dx)
            dist_g=max(sqrt(xgx2ygy2+zgz2),dx)
            if (dist_s < dx+0.0005) then
               !oblq(iz,iy,ix)=1.0
               oblq_s=0.0
            else
               oblq_s=abs(z)/dist_s
            endif
            if (dist_g < dx+0.0005) then
               !oblq(iz,iy,ix)=1.0
               oblq_g=0.0
            else
               oblq_g=abs(z)/dist_g
            endif          
            trace_temp(it1)=refl(iz,iy,ix)*oblq_s*oblq_g/(dist_s*dist_g)+trace_temp(it1)
                     !/(dist_s(iz,iy,ix)+dist_g(iz,iy,ix))+trace_temp(it1)
         endif
!         trace_temp(it)=refl(iz,iy,ix)*oblq_s(iz,iy,ix)*oblq_g(iz,iy,ix)&
!                     /(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))+trace_temp(it)
         !endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
call conv_wavelet(trace_temp,trace,wfft,nt)
end subroutine mod_kernal

subroutine mig_kernal(img,isIllum,illum,trace,wfft,sx,sy,sz,gx,gy,gz)
real,    intent(inout) :: img(:,:,:),illum(:,:,:)
logical, intent(in)    :: isIllum
real,    intent(in)    :: trace(:),sx,sy,sz,gx,gy,gz
complex, intent(in)    :: wfft(:)
real    :: x,y,z,xsx2,xgx2,ysy2,ygy2,zsz2,zgz2,xsx2ysy2,xgx2ygy2,&
   dist_s,dist_g,oblq_s,oblq_g
integer                :: ix, iy, iz, it0, it1
trace_temp=0.0
call xcorr_wavelet(trace,trace_temp,wfft,nt)
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,x,xsx2,xgx2,y,ysy2,ygy2,xsx2ysy2,&
!$OMP    it0,it1,z,zsz2,zgz2,dist_s,dist_g,oblq_s,oblq_g)
do ix=1,nx
   x=(ix-1)*dx
   xsx2=(x-sx)**2.0
   xgx2=(x-gx)**2.0
   do iy=1,ny
      y=(iy-1)*dx
      ysy2=(y-sy)**2.0
      ygy2=(y-gy)**2.0
      xsx2ysy2=xsx2+ysy2
      xgx2ygy2=xgx2+ygy2
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            z=(iz-1)*dx
            zsz2=(z-sz)**2.0
            zgz2=(z-gz)**2.0
            dist_s=max((xsx2ysy2+zsz2)**0.5,dx)
            dist_g=max((xgx2ygy2+zgz2)**0.5,dx)
            if (dist_s < dx+0.0005) then
               !oblq(iz,iy,ix)=1.0
               oblq_s=0.0
            else
               oblq_s=abs(z)/dist_s
            endif
            if (dist_g < dx+0.0005) then
               !oblq(iz,iy,ix)=1.0
               oblq_g=0.0
            else
               oblq_g=abs(z)/dist_g
            endif
            img(iz,iy,ix)=trace_temp(it1)*oblq_s*oblq_g/(dist_s*dist_g)+img(iz,iy,ix)
                  !/(dist_s(iz,iy,ix)+dist_g(iz,iy,ix))+img(iz,iy,ix)
            if (isIllum) then
               illum(iz,iy,ix)=illum(iz,iy,ix)+1.0/(dist_s*dist_g)
            endif
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
end subroutine mig_kernal

subroutine mod_kernal_control(trace,refl,wfft,xs,ys,zs,xg,yg,zg)
real, intent(in)    :: refl(:,:,:),xs,ys,zs,xg,yg,zg
real, intent(inout) :: trace(:)
complex, intent(in) :: wfft(:)
integer             :: ix, iy, iz, it0, it1
real    :: x,y,z,xsx2,xgx2,ysy2,ygy2,zsz2,zgz2,xsx2ysy2,xgx2ygy2,&
   dist_s,dist_g,oblq_s,oblq_g,xm,ym,zm,xmx2,ymy2,zmz2,xymxy2,theta_s,&
   theta_g,dip_sg

trace_temp=0.0
xm=0.5*(xs+xg)  
ym=0.5*(ys+yg)  
zm=0.5*(zs+zg)  
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,x,y,z,xsx2,ysy2,zsz2,xgx2,ygy2,zgz2,&
!$OMP    xmx2,ymy2,zmz2,xsx2ysy2,xgx2ygy2,xymxy2,theta_s,theta_g,dip_sg,it0,it1)
do ix=1,nx
   x=(ix-1)*dx
   xsx2=(xs-x)**2.0
   xgx2=(xg-x)**2.0
   xmx2=(xm-x)**2.0
   do iy=1,ny
      y=(iy-1)*dx
      ysy2=(ys-y)**2.0
      ygy2=(yg-y)**2.0
      ymy2=(ym-y)**2.0
      xsx2ysy2=xsx2+ysy2
      xgx2ygy2=xgx2+ygy2
      xymxy2=xmx2+ymy2
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            z=(iz-1)*dx
            zsz2=(zs-z)**2.0
            theta_s=xsx2ysy2/(xsx2ysy2+zsz2)
            if (theta_s<=theta2) then
               zgz2=(zg-z)**2.0
               theta_g=xgx2ygy2/(xgx2ygy2+zgz2)
               if (theta_g<=theta2) then
                  zmz2=(zm-z)**2.0
                  dip_sg=xymxy2/(xymxy2+zmz2)
                  if (dip_sg<dip2) then
                     dist_s=max(sqrt(xsx2ysy2+zsz2),dx)
                     dist_g=max(sqrt(xgx2ygy2+zgz2),dx)
                     if (dist_s < dx+0.0005) then
                        !oblq(iz,iy,ix)=1.0
                        oblq_s=0.0
                     else
                        oblq_s=abs(z)/dist_s
                     endif
                     if (dist_g < dx+0.0005) then
                        !oblq(iz,iy,ix)=1.0
                        oblq_g=0.0
                     else
                        oblq_g=abs(z)/dist_g
                     endif
                     trace_temp(it1)=refl(iz,iy,ix)*oblq_s*oblq_g/(dist_s*dist_g)+trace_temp(it1)
                     !/(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))+trace_temp(it1)
                  endif
               endif
            endif
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
call conv_wavelet(trace_temp,trace,wfft,nt)
end subroutine mod_kernal_control

subroutine mig_kernal_control(img,isIllum,illum,trace,wfft,xs,ys,zs,xg,yg,zg)
real, intent(in)    :: xs,ys,zs,xg,yg,zg,trace(:)
real, intent(inout) :: img(:,:,:),illum(:,:,:)
logical, intent(in) :: isIllum
complex, intent(in) :: wfft(:)
integer             :: ix,iy,iz,it0,it1
real    :: x,y,z,xsx2,xgx2,ysy2,ygy2,zsz2,zgz2,xsx2ysy2,xgx2ygy2,&
   dist_s,dist_g,oblq_s,oblq_g,xm,ym,zm,xmx2,ymy2,zmz2,xymxy2,theta_s,&
   theta_g,dip_sg
call xcorr_wavelet(trace,trace_temp,wfft,nt)
xm=0.5*(xs+xg)  
ym=0.5*(ys+yg)  
zm=0.5*(zs+zg)  
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,x,y,z,xsx2,ysy2,zsz2,xgx2,ygy2,zgz2,&
!$OMP    xmx2,ymy2,zmz2,xsx2ysy2,xgx2ygy2,xymxy2,theta_s,theta_g,dip_sg,it0,it1)
do ix=1,nx
   x=(ix-1)*dx
   xsx2=(xs-x)**2.0
   xgx2=(xg-x)**2.0
   xmx2=(xm-x)**2.0
   do iy=1,ny
      y=(iy-1)*dx
      ysy2=(ys-y)**2.0
      xsx2ysy2=xsx2+ysy2
      ygy2=(yg-y)**2.0
      xgx2ygy2=xgx2+ygy2
      ymy2=(ym-y)**2.0
      xymxy2=xmx2+ymy2
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            z=(iz-1)*dx
            zsz2=(zs-z)**2.0
            zgz2=(zg-z)**2.0
            theta_s=xsx2ysy2/(xsx2ysy2+zsz2)
            theta_g=xgx2ygy2/(xgx2ygy2+zgz2)
            if (theta_s<=theta2 .and. theta_g<=theta2) then
               zmz2=(zm-z)**2.0
               dip_sg=xymxy2/(xymxy2+zmz2)
               if (dip_sg<dip2) then
                  dist_s=max(sqrt(xsx2ysy2+zsz2),dx)
                  dist_g=max(sqrt(xgx2ygy2+zgz2),dx)
                  img(iz,iy,ix)=trace_temp(it1)*oblq_s*oblq_g&
                                /sqrt(dist_s*dist_g)+img(iz,iy,ix)
                      !/(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))+img(iz,iy,ix)
                  if (isIllum) then
                     illum(iz,iy,ix)=illum(iz,iy,ix)+1.0/(dist_s*dist_g)
                  endif
               endif
            endif
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
end subroutine mig_kernal_control

end module module_km3d
