! Big memory version
module module_km3d
use module_global, only : PI!date,time_now
use module_datatype
use module_eikfd3d_dir
use module_eikfd3d_fir
use module_fftw
use module_utility
use module_io
use module_sf90_io

implicit none
real,   allocatable,private :: ttt_s(:,:,:),ttt_g(:,:,:),&
   dist_s(:,:,:),dist_g(:,:,:),oblq_s(:,:,:),oblq_g(:,:,:),trace(:)
integer,allocatable,private :: iz0(:,:)
character(len=100), private :: eik_type, sort_in,dfor
logical, private :: isTheta,isDip,isSaveTTT,isSaveBothTTT,&
   isWriteData,isReadData
real,private :: dz,dy,dx,dt,theta,theta2,dip,dip2,dip_sg, &
   x,y,z,sxx2,syy2,szz2,gxx2,gyy2,gzz2,sxy2,gxy2,xm,ym,zm,&
   xmx2,ymy2,zmz2,xymxy2,theta_s,theta_g,vmin,vmax,dfor_coe
integer,private :: nx,ny,nz,nt,nt_conv,ng

private :: initial_variable,finalize_variable,read_ttt,calculate_ttt,&
   mod_kernal,mig_kernal,mod_kernal_control,mig_kernal_control      

contains

!------------------------------------------------------------------------

subroutine km3d_mod(kmm,s,output,sttt_in,gttt_in)
type(kmmod3d),    intent(in)    :: kmm
type(shot3d),     intent(inout) :: s
character(len=*), intent(in)    :: output,sttt_in,gttt_in
integer                         :: ig
call initial_variable(kmm,s)
if (isWriteData) open(11,file=output,access="direct",recl=I4*nt,action="write")
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
   call dist_cal(dist_s,nz,ny,nx,dz,dy,dx,s%zs(1),s%ys(1),s%xs(1))
   call oblq_cal(dist_s,oblq_s,nz,ny,nx,dz)
   do ig=1,ng
      if (isSaveTTT) then
         call read_ttt(gttt_in,s%gid(ig),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%yg(ig),s%xg(ig))
      endif
      call dist_cal(dist_g,nz,ny,nx,dz,dy,dx,s%zg(ig),s%yg(ig),s%xg(ig))
      call oblq_cal(dist_g,oblq_g,nz,ny,nx,dz)
      if (isTheta) then
         call mod_kernal_control(kmm%m%refl,kmm%k%wfft,s%zs(ig),&
            s%ys(ig),s%xs(ig),s%zg(ig),s%yg(ig),s%xg(ig))        
      else
         call mod_kernal(kmm%m%refl,kmm%k%wfft)
      endif        
      if (isWriteData) then
         write(11,rec=ig)trace
      else
         s%seis(:,ig)=trace
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
   call dist_cal(dist_g,nz,ny,nx,dz,dy,dx,s%zg(1),s%yg(1),s%xg(1))
   call oblq_cal(dist_g,oblq_g,nz,ny,nx,dz)
   do ig=1,ng
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(ig),ttt_s)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(ig),s%ys(ig),s%xs(ig))
      endif
      call dist_cal(dist_s,nz,ny,nx,dz,dy,dx,s%zs(ig),s%ys(ig),s%xs(ig))
      call oblq_cal(dist_s,oblq_s,nz,ny,nx,dz)
      if (isTheta) then
         call mod_kernal_control(kmm%m%refl,kmm%k%wfft,s%zs(ig),&
            s%ys(ig),s%xs(ig),s%zg(ig),s%yg(ig),s%xg(ig))
      else
         call mod_kernal(kmm%m%refl,kmm%k%wfft)
      endif
      if (isWriteData) then
         write(11,rec=ig)trace
      else
         s%seis(:,ig)=trace
      endif
   enddo
else
   do ig=1,ng
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(ig),ttt_s)
         call read_ttt(gttt_in,s%gid(ig),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(ig),s%ys(ig),s%xs(ig))
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%yg(ig),s%xg(ig))
      endif
      call dist_cal(dist_s,nz,ny,nx,dz,dy,dx,s%zs(ig),s%ys(ig),s%xs(ig))
      call oblq_cal(dist_s,oblq_s,nz,ny,nx,dz)
      call dist_cal(dist_g,nz,ny,nx,dz,dy,dx,s%zg(ig),s%yg(ig),s%xg(ig))
      call oblq_cal(dist_g,oblq_g,nz,ny,nx,dz)
      if (isTheta) then
         call mod_kernal_control(kmm%m%refl,kmm%k%wfft,s%zs(ig),&
            s%ys(ig),s%xs(ig),s%zg(ig),s%yg(ig),s%xg(ig))
      else
         call mod_kernal(kmm%m%refl,kmm%k%wfft)
      endif
      if (isWriteData) then
         write(11,rec=ig)trace
      else
         s%seis(:,ig)=trace
      endif
   enddo
endif
if (isWriteData) close(11)
call finalize_variable
end subroutine km3d_mod

!------------------------------------------------------------------------

subroutine km3d_mig(kmm,s,img,output,sttt_in,gttt_in)
type(kmmod3d),    intent(in)    :: kmm
type(shot3d),     intent(in)    :: s
type(image3d),    intent(inout) :: img 
character(len=*), intent(in)    :: output,sttt_in,gttt_in
integer                         :: ig
img%img=0.0
call initial_variable(kmm,s)
if (isReadData) open(13,file=output,access="direct",recl=I4*nt,action="read")
if (img%isIllum) img%illum=0.0
if (sort_in.eq."CSG") then
   if (isSaveTTT) then
      if (isSaveBothTTT) then
         call read_ttt(sttt_in,s%sid(1),ttt_s)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%ys(1),s%xs(1))
      endif
   else
      call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%ys(1),s%xs(1))
   endif
   call dist_cal(dist_s,nz,ny,nx,dz,dy,dx,s%zs(1),s%ys(1),s%xs(1))
   call oblq_cal(dist_s,oblq_s,nz,ny,nx,dz)
   do ig=1,ng
      if (isReadData) then
         read(13,rec=ig)trace
      else
         trace=s%seis(:,ig)
      endif
      if (isSaveTTT) then
         call read_ttt(gttt_in,s%gid(ig),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%yg(ig),s%xg(ig))
      endif
      call dist_cal(dist_g,nz,ny,nx,dz,dy,dx,s%zg(ig),s%yg(ig),s%xg(ig))
      call oblq_cal(dist_g,oblq_g,nz,ny,nx,dz)
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,kmm%k%wfft,&
            s%zs(ig),s%ys(ig),s%xs(ig),s%zg(ig),s%yg(ig),s%xg(ig))
      else
         call mig_kernal(img%img,img%isIllum,img%illum,kmm%k%wfft)
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
   call dist_cal(dist_g,nz,ny,nx,dz,dy,dx,s%zg(1),s%yg(1),s%xg(1))
   call oblq_cal(dist_g,oblq_g,nz,ny,nx,dz)
   do ig=1,ng
      if (isReadData) then
         read(13,rec=ig)trace
      else
         trace=s%seis(:,ig)
      endif
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(ig),ttt_s)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(ig),s%ys(ig),s%xs(ig))
      endif
      call dist_cal(dist_s,nz,ny,nx,dz,dy,dx,s%zs(ig),s%ys(ig),s%xs(ig))
      call oblq_cal(dist_s,oblq_s,nz,ny,nx,dz)
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,kmm%k%wfft,&
            s%zs(ig),s%ys(ig),s%xs(ig),s%zg(ig),s%yg(ig),s%xg(ig))
      else
         call mig_kernal(img%img,img%isIllum,img%illum,kmm%k%wfft)
      endif
   enddo
else
   do ig=1,ng
      if (isReadData) then
         read(13,rec=ig)trace
      else
         trace=s%seis(:,ig)
      endif
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(ig),ttt_s)
         call read_ttt(gttt_in,s%gid(ig),ttt_g)
      else
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(ig),s%ys(ig),s%xs(ig))
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%yg(ig),s%xg(ig))
      endif
      call dist_cal(dist_s,nz,ny,nx,dz,dy,dx,s%zs(ig),s%ys(ig),s%xs(ig))
      call oblq_cal(dist_s,oblq_s,nz,ny,nx,dz)
      call dist_cal(dist_g,nz,ny,nx,dz,dy,dx,s%zg(ig),s%yg(ig),s%xg(ig))
      call oblq_cal(dist_g,oblq_g,nz,ny,nx,dz)
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,kmm%k%wfft,&
            s%zs(ig),s%ys(ig),s%xs(ig),s%zg(ig),s%yg(ig),s%xg(ig))
      else
         call mig_kernal(img%img,img%isIllum,img%illum,kmm%k%wfft)
      endif
   enddo
endif
if (isReadData) close(13)
if (img%isIllum) img%illum=merge(1.0,img%illum,img%illum<0.000000001)
call finalize_variable
end subroutine km3d_mig

!------------------------------------------------------------------------

subroutine initial_variable(kmm,s)
type(kmmod3d),  intent(in) :: kmm
type(shot3d),   intent(in) :: s
call copy_variable(kmm%k%nt,nt,kmm%k%nt_conv,nt_conv) 
call copy_variable(kmm%k%isTheta,isTheta,kmm%k%isDip,isDip)
call copy_variable(kmm%k%isWriteData,isWriteData,kmm%k%isReadData,isReadData)
call copy_variable(kmm%k%isSaveTTT,isSaveTTT,kmm%k%isSaveBothTTT,isSaveBothTTT)
call copy_variable(kmm%k%eik_type,eik_type,kmm%k%sort_in,sort_in)
call copy_variable(kmm%k%dt,dt,kmm%m%dx,dx,kmm%m%dy,dy,kmm%m%dz,dz)
call copy_variable(kmm%m%nx,nx,kmm%m%ny,ny,kmm%m%nz,nz,s%ng,ng)
call copy_variable(kmm%k%dfor,dfor)
if (dfor.eq."SHORT") call copy_variable(kmm%k%dfor_coe,dfor_coe)
call allocate_and_initial(ttt_s,ttt_g,dist_s,nz,ny,nx)
call allocate_and_initial(dist_g,oblq_s,oblq_g,nz,ny,nx)
call allocate_and_initial(trace,nt)
call allocate_and_initial(iz0,ny,nx)
iz0=nint(kmm%k%z0_mig/dz)+1
if (isTheta) theta=kmm%k%theta; theta2=(sin(theta/180.0*PI))**2.0
if (isDip) dip=kmm%k%dip; dip2=(sin(dip/180.0*PI))**2.0
if (allocated(kmm%m%v)) then
   if (eik_type.eq."FIR") then
      vmin=1.0/maxval(kmm%m%v); vmax=1.0/minval(kmm%m%v)
   endif
endif
end subroutine initial_variable

subroutine finalize_variable
call deallocate_and_free(ttt_s,ttt_g,dist_s,dist_g,oblq_s,oblq_g)
call deallocate_and_free(trace)
call deallocate_and_free(iz0)
end subroutine finalize_variable

subroutine read_ttt(ttt_in,id,ttt)
integer, intent(in) :: id
character(len=*),intent(in) :: ttt_in
real, intent(out) :: ttt(:,:,:)
character(len=slen) :: ttt_name
call filename(ttt_name, ttt_in, id, ".bin")
if (dfor.eq."SHORT") then
   call sf90_readself(ttt_name,ttt,nz,ny,nx,dfor_coe)
else
   call sf90_readself(ttt_name,ttt,nz,ny,nx)
endif
end subroutine read_ttt

subroutine calculate_ttt(s,ttt,zs,ys,xs)
real, intent(in) :: s(:,:,:),zs,ys,xs
real, intent(out):: ttt(:,:,:)
if (eik_type.eq."DIR") then
   call eikfd3d_dir(s,ttt,zs,ys,xs,nz,ny,nx,dx)
else
   call eikfd3d_fir(s,ttt,zs,ys,xs,nz,ny,nx,dx,vmax,vmin)
endif
end subroutine calculate_ttt

subroutine mod_kernal(refl,wfft)
real, intent(in)    :: refl(:,:,:)
complex, intent(in) :: wfft(:)
real :: trace_temp(nt) 
integer             :: ix, iy, iz, it
trace_temp=0.0
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,it)
do ix=1,nx
   do iy=1,ny
      do iz=iz0(iy,ix),nz
         it=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)+1
         if (it<=nt) then
            trace_temp(it)=refl(iz,iy,ix)*oblq_s(iz,iy,ix)*oblq_g(iz,iy,ix)&
                     /sqrt(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))+trace_temp(it)
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine mod_kernal

subroutine mig_kernal(img,isIllum,illum,wfft)
real,    intent(inout) :: img(:,:,:),illum(:,:,:)
logical, intent(in)    :: isIllum
complex, intent(in)    :: wfft(:)
integer                :: ix, iy, iz, it
real :: trace_temp(nt)
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,it)
do ix=1,nx
   do iy=1,ny
      do iz=iz0(iy,ix),nz
         it=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)+1
         if (it<=nt) then
            img(iz,iy,ix)=trace_temp(it)*oblq_s(iz,iy,ix)*oblq_g(iz,iy,ix)&
                  /sqrt(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))+img(iz,iy,ix)
            if (isIllum) then
               illum(iz,iy,ix)=illum(iz,iy,ix)&
                              +1.0/sqrt(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))
            endif
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
end subroutine mig_kernal

subroutine mod_kernal_control(refl,wfft,zs,ys,xs,zg,yg,xg)
real, intent(in)    :: refl(:,:,:),zs,ys,xs,zg,yg,xg
complex, intent(in) :: wfft(:)
integer             :: ix, iy, iz, it
real :: trace_temp(nt)
trace_temp=0.0
xm=0.5*(xs+xg);   ym=0.5*(ys+yg);  zm=0.5*(zs+zg)  
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,x,y,z,sxx2,syy2,szz2,gxx2,gyy2,gzz2,&
!$OMP xmx2,ymy2,zmz2,sxy2,gxy2,xymxy2,theta_s,theta_g,dip_sg,it)
do ix=1,nx
   x=(ix-1)*dx; sxx2=(xs-x)**2.0; gxx2=(xg-x)**2.0; xmx2=(xm-x)**2.0
   do iy=1,ny
      y=(iy-1)*dy; syy2=(ys-y)**2.0; gyy2=(yg-y)**2.0; ymy2=(ym-y)**2.0
      sxy2=sxx2+syy2; gxy2=gxx2+gyy2; xymxy2=xmx2+ymy2
      do iz=iz0(iy,ix),nz
         z=(iz-1)*dz; szz2=(zs-z)**2.0; gzz2=(zg-z)**2.0
         theta_s=sxy2/(sxy2+szz2); theta_g=gxy2/(gxy2+gzz2)
         if (theta_s<=theta2 .and. theta_g<=theta2) then
            zmz2=(zm-z)**2.0; dip_sg=xymxy2/(xymxy2+zmz2)
            if (dip_sg<dip2) then
               it=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)+1
               if (it<=nt) then
                  trace_temp(it)=refl(iz,iy,ix)*oblq_s(iz,iy,ix)*oblq_g(iz,iy,ix)&
                     /sqrt(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))+trace_temp(it)
               endif
            endif
         endif
      enddo
   enddo
enddo
!$OMP END PARALLEL DO
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine mod_kernal_control

subroutine mig_kernal_control(img,isIllum,illum,wfft,zs,ys,xs,zg,yg,xg)
real, intent(in)    :: zs,ys,xs,zg,yg,xg
real, intent(inout) :: img(:,:,:),illum(:,:,:)
logical, intent(in) :: isIllum
complex, intent(in) :: wfft(:)
real :: trace_temp(nt)
integer             :: ix,iy,iz,it
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
xm=0.5*(xs+xg);   ym=0.5*(ys+yg);   zm=0.5*(zs+zg)  
!$OMP PARALLEL DO PRIVATE(iz,iy,ix,x,y,z,sxx2,syy2,szz2,gxx2,gyy2,gzz2,&
!$OMP xmx2,ymy2,zmz2,sxy2,gxy2,xymxy2,theta_s,theta_g,dip_sg,it)
do ix=1,nx
   x=(ix-1)*dx; sxx2=(xs-x)**2.0; gxx2=(xg-x)**2.0; xmx2=(xm-x)**2.0
   do iy=1,ny
      y=(iy-1)*dy; syy2=(ys-y)**2.0; sxy2=sxx2+syy2; gyy2=(yg-y)**2.0
      gxy2=gxx2+gyy2; ymy2=(ym-y)**2.0; xymxy2=xmx2+ymy2
      do iz=iz0(iy,ix),nz
         z=(iz-1)*dz; szz2=(zs-z)**2.0; gzz2=(zg-z)**2.0
         theta_s=sxy2/(sxy2+szz2); theta_g=gxy2/(gxy2+gzz2)
         if (theta_s<=theta2 .and. theta_g<=theta2) then
            zmz2=(zm-z)**2.0; dip_sg=xymxy2/(xymxy2+zmz2)
            if (dip_sg<dip2) then
               it=nint((ttt_s(iz,iy,ix)+ttt_g(iz,iy,ix))/dt)+1
               if (it<=nt) then
                  img(iz,iy,ix)=trace_temp(it)*oblq_s(iz,iy,ix)*oblq_g(iz,iy,ix)&
                      /sqrt(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))+img(iz,iy,ix)
                  if (isIllum) then
                     illum(iz,iy,ix)=illum(iz,iy,ix)&
                                    +1.0/sqrt(dist_s(iz,iy,ix)*dist_g(iz,iy,ix))
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
