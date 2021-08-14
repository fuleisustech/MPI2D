module module_pwkm3d
use module_global, only : slen,date,time_now
use module_datatype
use module_eikfd3d_dir
use module_eikfd3d_fir
use module_pweikfd3d_dir
use module_fftw
use module_io
use module_sf90_mpi
use module_string
use module_utility

implicit none
real, allocatable, private :: ttt_s(:,:,:),ttt_g(:,:,:),dist_g(:,:,:),&
   oblq_g(:,:,:),trace(:),trace_temp(:)
integer, allocatable, private :: iz0(:,:)
character(len=100),private :: eik_type
character(len=slen),private :: output
logical, private :: isTheta, isDip, isSaveTTT,isSaveBothTTT
real, private :: dx, dy, dz, dt, theta, theta2, &
   x, y, z, sxx2, syy2, szz2, gxx2, gyy2, gzz2, sxy2, gxy2, &
   xm, ym, zm, xmx2, ymy2, zmz2, xymxy2, theta_s, theta_g, zs, ys, xs
integer, private :: nx, ny, nz, nt, nt_conv, it, ng

private :: initial_variable, finalize_variable, read_ttt, calculate_ttt,&
   calculate_pwttt,pw1mod_kernal, pw1mig_kernal, pw1mod_kernal_control,&
   pw1mig_kernal_control, pw2mod_kernal, pw2mig_kernal  

interface initial_variable
   module procedure pw1_initial_variable
   module procedure pw2_initial_variable
end interface initial_variable

contains

!------------------------------------------------------------------------

subroutine pw1km3d_mod(kmm,s,pw1ttt_in,gttt_in)
type(kmmod3d), intent(in)  :: kmm
character(len=100), intent(in)  :: pw1ttt_in,gttt_in
type(pw1_shot3d),  intent(inout) :: s
integer                    :: itrace
call initial_variable(kmm,s)
s%seis=0.0
if (isSaveTTT.and.isSaveBothTTT) then
   call read_ttt(pw1ttt_in,s%proid,ttt_s)
else
   call calculate_pwttt(kmm%m%v,ttt_s,s%theta_x1,s%xref_1,&
                        s%theta_y1,s%yref_1,s%zp1)
endif
do itrace=1,s%ng
   if (isSaveTTT) then
      call read_ttt(gttt_in,s%gid(itrace),ttt_g)
   else
      call calculate_ttt(kmm%m%v,ttt_g,s%zg(itrace),s%yg(itrace),s%xg(itrace))
   endif
   call dist_cal(dist_g,nz,ny,nx,dz,dy,dx,s%xg(itrace),s%yg(itrace),s%zg(itrace))
   call oblq_cal(dist_g,oblq_g,nz,ny,nx,dz)
   if (isTheta) then
      !call pw1mod_kernal_control(trace,kmm%m%refl,kmm%k%wfft,s%zg(itrace),&
      !        s%yg(itrace),s%xg(itrace))
      call pw1mod_kernal_control(s%seis(:,itrace),kmm%m%refl,kmm%k%wfft,&
              s%zg(itrace),s%yg(itrace),s%xg(itrace))
   else
      !call pw1mod_kernal(trace,kmm%m%refl,kmm%k%wfft)
      call pw1mod_kernal(s%seis(:,itrace),kmm%m%refl,kmm%k%wfft)
   endif
   !s%seis(:,itrace)=trace
enddo
call finalize_variable
end subroutine pw1km3d_mod

!------------------------------------------------------------------------

subroutine pw1km3d_mig(kmm,s,img,pw1ttt_in,gttt_in)
type(kmmod3d), intent(in)    :: kmm
type(pw1_shot3d),  intent(in):: s
character(len=100),intent(in):: pw1ttt_in,gttt_in
type(image3d), intent(inout) :: img 
integer                      :: itrace
call initial_variable(kmm,s)
img%img=0.0
if (img%isIllum) then
   img%illum=0.0
endif
if (isSaveTTT.and.isSaveBothTTT) then
   call read_ttt(pw1ttt_in,s%proid,ttt_s)
else
   call calculate_pwttt(kmm%m%v,ttt_s,s%theta_x1,s%xref_1,&
                        s%theta_y1,s%yref_1,s%zp1)
endif
do itrace=1,s%ng
   !trace=s%seis(:,itrace)
   if (isSaveTTT) then
      call read_ttt(gttt_in,s%gid(itrace),ttt_g)
   else
      call calculate_ttt(kmm%m%v,ttt_g,s%zg(itrace),s%yg(itrace),s%xg(itrace))
   endif
   call dist_cal(dist_g,nz,ny,nx,dz,dy,dx,s%xg(itrace),s%yg(itrace),s%zg(itrace))
   call oblq_cal(dist_g,oblq_g,nz,ny,nx,dz)
   if (isTheta) then
      !call pw1mig_kernal_control(img%img,img%isIllum,img%illum,trace,kmm%k%wfft,&
      !        s%zg(itrace),s%yg(itrace),s%xg(itrace))
      call pw1mig_kernal_control(img%img,img%isIllum,img%illum,s%seis(:,itrace),&
              kmm%k%wfft,s%zg(itrace),s%yg(itrace),s%xg(itrace))
   else
      !call pw1mig_kernal(img%img,img%isIllum,img%illum,trace,kmm%k%wfft)
      call pw1mig_kernal(img%img,img%isIllum,img%illum,s%seis(:,itrace),kmm%k%wfft)
   endif
enddo
if (img%isIllum) then
   img%illum=merge(1.0,img%illum,img%illum<0.000000001)
endif
call finalize_variable
end subroutine pw1km3d_mig

!------------------------------------------------------------------------

subroutine pw2km3d_mod(kmm,s,pw1ttt_in,pw2ttt_in)
type(kmmod3d), intent(in)  :: kmm
character(len=100), intent(in)  :: pw1ttt_in, pw2ttt_in
type(pw2_shot3d),  intent(inout) :: s
integer                    :: itrace
call initial_variable(kmm,s)
s%seis=0.0
if (isSaveTTT.and.isSaveBothTTT) then
   call read_ttt(pw1ttt_in,s%proid,ttt_s)
else
   call calculate_pwttt(kmm%m%v,ttt_s,s%theta_x1,s%xref_1,&
                        s%theta_y1,s%yref_1,s%zp1)
endif
do itrace=1,s%ng
   !trace=0.0
   if (isSaveTTT) then
      call read_ttt(pw2ttt_in,s%gid(itrace),ttt_g)
   else
      call calculate_pwttt(kmm%m%v,ttt_g,s%theta_x2,s%xref_2,&
                           s%theta_y2,s%yref_2,s%zp2)
   endif
   !call pw2mod_kernal(trace,kmm%m%refl,kmm%k%wfft)
   call pw2mod_kernal(s%seis(:,itrace),kmm%m%refl,kmm%k%wfft)
   !s%seis(:,itrace)=trace
enddo
call finalize_variable
end subroutine pw2km3d_mod

!------------------------------------------------------------------------

subroutine pw2km3d_mig(kmm,s,img,pw1ttt_in,pw2ttt_in)
type(kmmod3d), intent(in)    :: kmm
type(pw2_shot3d),  intent(in):: s
character(len=100), intent(in)  :: pw1ttt_in,pw2ttt_in
type(image3d), intent(inout) :: img 
integer                      :: itrace
call initial_variable(kmm,s)
img%img=0.0
if (isSaveTTT.and.isSaveBothTTT) then
   call read_ttt(pw1ttt_in,s%proid,ttt_s)
else
   call calculate_pwttt(kmm%m%v,ttt_s,s%theta_x1,s%xref_1,&
                        s%theta_y1,s%yref_1,s%zp1)
endif
do itrace=1,s%ng
   !trace=s%seis(:,itrace)
   if (isSaveTTT) then
      call read_ttt(pw2ttt_in,s%gid(itrace),ttt_g)
   else
      call calculate_pwttt(kmm%m%v,ttt_g,s%theta_x2,s%xref_2,&
                           s%theta_y2,s%yref_2,s%zp2)
   endif
   !call pw2mig_kernal(img%img,trace,kmm%k%wfft)
   call pw2mig_kernal(img%img,s%seis(:,itrace),kmm%k%wfft)
enddo
call finalize_variable
end subroutine pw2km3d_mig

!------------------------------------------------------------------------

subroutine pw1_initial_variable(kmm,s)
type(kmmod3d),  intent(in) :: kmm
type(pw1_shot3d),   intent(in) :: s
call copy_variable(kmm%k%nt,nt,kmm%k%nt_conv,nt_conv) 
call copy_variable(kmm%k%isTheta,isTheta)
call copy_variable(kmm%k%eik_type,eik_type)
call copy_variable(kmm%k%isSaveTTT,isSaveTTT)
call copy_variable(kmm%k%isSaveBothTTT,isSaveBothTTT)
call copy_variable(kmm%k%dt,dt,kmm%k%theta,theta)
call copy_variable(kmm%m%dx,dx,kmm%m%dy,dy,kmm%m%dz,dz)
call copy_variable(kmm%m%nx,nx,kmm%m%ny,ny,kmm%m%nz,nz,s%ng,ng)
call allocate_and_initial(ttt_s,ttt_g,dist_g,oblq_g,nz,ny,nx)
call allocate_and_initial(trace,nt)
call allocate_and_initial(trace_temp,nt)
call allocate_and_initial(iz0,ny,nx)
iz0=nint(kmm%k%z0_mig/dx)+1
theta2=(sin(theta/180.0*PI))**2.0
end subroutine pw1_initial_variable

subroutine pw2_initial_variable(kmm,s)
type(kmmod3d),  intent(in) :: kmm
type(pw2_shot3d),   intent(in) :: s
call copy_variable(kmm%k%nt,nt,kmm%k%nt_conv,nt_conv) 
call copy_variable(kmm%k%eik_type,eik_type)
call copy_variable(kmm%k%isSaveTTT,isSaveTTT)
call copy_variable(kmm%k%isSaveBothTTT,isSaveBothTTT)
call copy_variable(kmm%m%dx,dx,kmm%m%dy,dy,kmm%m%dz,dz)
call copy_variable(kmm%k%dt,dt)
call copy_variable(kmm%m%nx,nx,kmm%m%ny,ny,kmm%m%nz,nz,s%ng,ng)
call allocate_and_initial(ttt_s,ttt_g,nz,ny,nx)
call allocate_and_initial(trace,nt)
call allocate_and_initial(trace_temp,nt)
call allocate_and_initial(iz0,ny,nx)
iz0=nint(kmm%k%z0_mig/dx)+1
end subroutine pw2_initial_variable

subroutine finalize_variable
call deallocate_and_free(ttt_s,ttt_g,dist_g,oblq_g)
call deallocate_and_free(trace)
call deallocate_and_free(trace_temp)
call deallocate_and_free(iz0)
end subroutine finalize_variable

subroutine read_ttt(ttt_in,id,ttt)
character(len=*),intent(in)  :: ttt_in
integer,         intent(in)  :: id
real,            intent(out) :: ttt(:,:,:)
integer :: fid
character(len=slen) :: output
call filename(output, ttt_in, id, ".bin")
call sf90_readself(output,ttt,nz,ny,nx)
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

subroutine calculate_pwttt(s,ttt,theta_x,xref,theta_y,yref,zp)
real, intent(in) :: s(:,:,:),theta_x,xref,theta_y,yref,zp
real, intent(out):: ttt(:,:,:)
!if (eik_type.eq."DIR") then
call pweikfd3d_dir(s,ttt,theta_x,xref,theta_y,yref,zp,nz,ny,nx,dx)
!else
!   call eikfd2d_fir(s,ttt,zs,xs,nz,nx,dx,dx)
!endif
end subroutine calculate_pwttt

subroutine pw1mod_kernal(trace,refl,wfft)
real, intent(in)    :: refl(:,:,:)
real, intent(out)   :: trace(:)
complex, intent(in) :: wfft(:)
integer             :: ix, iy, iz, it0, it1
trace_temp=0.0
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,it0,it1)
do ix=1,nx
   do iy=1,ny
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            trace_temp(it1)=refl(iz,iy,ix)*oblq_g(iz,iy,ix)&
                           /dist_g(iz,iy,ix)+trace_temp(it1)
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine pw1mod_kernal

subroutine pw1mig_kernal(img,isIllum,illum,trace,wfft)
real,    intent(inout) :: img(:,:,:),illum(:,:,:)
logical, intent(in)    :: isIllum
real,    intent(in)    :: trace(:)
complex, intent(in)    :: wfft(:)
integer                :: ix, iy, iz, it0, it1
real :: x,y,z,xgx2,ygy2,zgz2,xgx2ygy2,dist_g,oblq_g
trace_temp=0.0
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,it0,it1)
do ix=1,nx
   do iy=1,ny
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            img(iz,iy,ix)=trace_temp(it)*oblq_g(iz,iy,ix)&
                                        /dist_g(iz,iy,ix)+img(iz,iy,ix)
            if (isIllum) then
               illum(iz,iy,ix)=illum(iz,iy,ix)+1.0/dist_g(iz,iy,ix)
            endif
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
end subroutine pw1mig_kernal

subroutine pw1mod_kernal_control(trace,refl,wfft,xg,yg,zg)
real, intent(in)    :: refl(:,:,:),xg,yg,zg
real, intent(out)   :: trace(:)
complex, intent(in) :: wfft(:)
integer             :: ix, iy, iz, it0, it1
real :: x,y,xgx2,ygy2,xgx2ygy2
trace_temp=0.0
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,x,xgx2,y,ygy2,xgx2ygy2,&
!$OMP    it0,it1,theta_g)
do ix=1,nx
   x=(ix-1)*dx
   xgx2=(x-xg)**2.0
   do iy=1,ny
      y=(iy-1)*dx
      ygy2=(y-yg)**2.0
      xgx2ygy2=xgx2+ygy2
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            theta_g=xgx2ygy2/dist_g(iz,iy,ix)
            if (theta_g<=theta2) then
               trace_temp(it1)=refl(iz,iy,ix)*oblq_g(iz,iy,ix)&
                              /dist_g(iz,iy,ix)+trace_temp(it1)
            endif
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine pw1mod_kernal_control

subroutine pw1mig_kernal_control(img,isIllum,illum,trace,wfft,xg,yg,zg)
real, intent(in)    :: trace(:),xg,yg,zg
real, intent(inout) :: img(:,:,:),illum(:,:,:)
logical,intent(in)  :: isIllum
complex, intent(in) :: wfft(:)
integer             :: ix,iy,iz,it0,it1
real :: x,y,xgx2,ygy2,xgx2ygy2,theta_g
trace_temp=0.0
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,x,xgx2,y,ygy2,xgx2ygy2,&
!$OMP    it0,it1,theta_g)
do ix=1,nx
   x=(ix-1)*dx
   xgx2=(x-xg)**2.0
   do iy=1,ny
      y=(iy-1)*dx
      ygy2=(y-yg)**2.0
      xgx2ygy2=xgx2+ygy2
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            theta_g=xgx2ygy2/dist_g(iz,iy,ix)
            if (theta_g<=theta2) then
               img(iz,iy,ix)=trace_temp(it)*oblq_g(iz,iy,ix)&
                            /dist_g(iz,iy,ix)+img(iz,iy,ix)
            endif
            if (isIllum) then
              illum(iz,iy,ix)=illum(iz,iy,ix)+1.0/dist_g(iz,iy,ix)
            endif
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
end subroutine pw1mig_kernal_control

subroutine pw2mod_kernal(trace,refl,wfft)
real, intent(in)    :: refl(:,:,:)
real, intent(out)   :: trace(:)
complex, intent(in) :: wfft(:)
integer             :: ix, iy, iz, it0, it1
trace_temp=0.0
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,it0,it1)
do ix=1,nx
   do iy=1,ny
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            trace_temp(it)=refl(iz,iy,ix)+trace_temp(it)
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine pw2mod_kernal

subroutine pw2mig_kernal(img,trace,wfft)
real,    intent(inout) :: img(:,:,:)
real,    intent(in)    :: trace(:)
complex, intent(in)    :: wfft(:)
integer                :: ix, iy, iz, it0, it1
trace_temp=0.0
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,it0,it1)
do ix=1,nx
   do iy=1,ny
      do iz=iz0(iy,ix),nz
         it0=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)
         it1=it0+1
         if (it0>0 .and. it0<=nt-1) then
            img(iz,iy,ix)=trace_temp(it)+img(iz,iy,ix)
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
end subroutine pw2mig_kernal

end module module_pwkm3d
