module module_pwkm2d
use module_global, only : PI
use module_datatype
use module_eikfd2d_dir
use module_eikfd2d_fir
use module_pweikfd2d_dir
use module_fftw
use module_io
use module_sf90_io
use module_string

implicit none
real, allocatable, private :: ttt_s(:,:),ttt_g(:,:),dist_g(:,:),&
   oblq_g(:,:),trace(:)
character(len=100) ,private :: eik_type,dfor
integer, allocatable, private :: iz0(:)
logical, private :: isTheta,isSaveTTT,isSaveBothTTT,isReadData,isWriteData
real, private :: dx,dz,dt,theta,theta2,x,z,gxx2,gzz2,theta_g,dfor_coe
integer, private :: nx, nz, nt, nt_conv, ng

private :: initial_variable, finalize_variable,read_ttt, calculate_ttt, &
   calculate_pwttt,pw1_mod_kernal, pw1_mig_kernal, pw1_mod_kernal_control,&
   pw1_mig_kernal_control, pw2_mod_kernal, pw2_mig_kernal        

interface initial_variable
   module procedure pw1_initial_variable
   module procedure pw2_initial_variable
end interface initial_variable

contains

!------------------------------------------------------------------------

subroutine pw1km2d_mod(kmm,s,output,pw1ttt_in,gttt_in)
type(kmmod2d), intent(in)  :: kmm
character(len=*), intent(in)    :: output,pw1ttt_in,gttt_in
type(pw1_shot2d), intent(inout) :: s
integer                    :: ig
call initial_variable(kmm,s)
if (isWriteData) open(11,file=output,access="direct",recl=I4*nt,action="write")
if (isSaveTTT.and.isSaveBothTTT) then
   call read_ttt(pw1ttt_in,s%proid,ttt_s)
else
   call calculate_pwttt(kmm%m%v,ttt_s,s%theta_x1,s%xref_1,s%zp1)
endif
do ig=1,ng
   if (isSaveTTT) then
      call read_ttt(gttt_in,s%gid(ig),ttt_g)
   else ! calculate the source side traveltime table
      call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%xg(ig))
   endif
   call dist_cal(dist_g,nz,nx,dz,dx,s%zg(ig),s%xg(ig))
   call oblq_cal(dist_g,oblq_g,nz,nx,dz)   
   if (isTheta) then
      call pw1_mod_kernal_control(kmm%m%refl,kmm%k%wfft,s%zg(ig),s%xg(ig))
   else
      call pw1_mod_kernal(kmm%m%refl,kmm%k%wfft)
   endif
   if (isWriteData) then
      write(11,rec=ig)trace
   else
      s%seis(:,ig)=trace
   endif
enddo
if (isWriteData) close(11)
call finalize_variable
end subroutine pw1km2d_mod

!------------------------------------------------------------------------

subroutine pw1km2d_mig(kmm,s,img,output,pw1ttt_in,gttt_in)
type(kmmod2d),     intent(in)    :: kmm
type(pw1_shot2d),  intent(in)    :: s
character(len=*),  intent(in)    :: output,pw1ttt_in,gttt_in
type(image2d),  intent(inout)    :: img 
integer                          :: ig
img%img=0.0
call initial_variable(kmm,s)
if (isReadData) open(13,file=output,access="direct",recl=I4*nt,action="read")
if (img%isIllum) img%illum=0.0
if (isSaveTTT.and.isSaveBothTTT) then
   call read_ttt(pw1ttt_in,s%proid,ttt_s)
else
   call calculate_pwttt(kmm%m%v,ttt_s,s%theta_x1,s%xref_1,s%zp1)
endif
do ig=1,ng
   if (isReadData) then
      read(13,rec=ig)trace
   else
      trace=s%seis(:,ig)
   endif
   if (isSaveTTT) then
      call read_ttt(gttt_in,s%gid(ig),ttt_g)
   else ! calculate the source side traveltime table
      call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%xg(ig))
   endif
   call dist_cal(dist_g,nz,nx,dz,dx,s%zg(ig),s%xg(ig))   
   call oblq_cal(dist_g,oblq_g,nz,nx,dx)   
   if (isTheta) then
      call pw1_mig_kernal_control(img%img,img%isIllum,img%illum,&
                                  kmm%k%wfft,s%zg(ig),s%xg(ig))
   else
      call pw1_mig_kernal(img%img,img%isIllum,img%illum,kmm%k%wfft)
   endif
enddo
if (isReadData) close(13)
if (img%isIllum) img%illum=merge(1.0,img%illum,img%illum<0.0000000001)
call finalize_variable
end subroutine pw1km2d_mig

!------------------------------------------------------------------------

subroutine pw2km2d_mod(kmm,s,output,pw1ttt_in,pw2ttt_in)
type(kmmod2d), intent(in)  :: kmm
character(len=*) :: output,pw1ttt_in,pw2ttt_in
type(pw2_shot2d),  intent(inout) :: s
integer :: ig
call initial_variable(kmm,s)
if (isWriteData) open(11,file=output,access="direct",recl=I4*nt,action="write")
if (isSaveTTT.and.isSaveBothTTT) then
   call read_ttt(pw1ttt_in,s%proid,ttt_s)
else
   call calculate_pwttt(kmm%m%v,ttt_s,s%theta_x1,s%xref_1,s%zp1)
endif
do ig=1,ng
   if (isSaveTTT) then
      call read_ttt(pw2ttt_in,s%gid(ig),ttt_g)
   else
      call calculate_pwttt(kmm%m%v,ttt_g,s%theta_x2(ig),s%xref_2(ig),s%zp2(ig))
   endif
   call pw2_mod_kernal(kmm%m%refl,kmm%k%wfft)
   if (isWriteData) then
      write(11,rec=ig)trace
   else
      s%seis(:,ig)=trace
   endif
enddo
if (isWriteData) close(11)
call finalize_variable
end subroutine pw2km2d_mod

!------------------------------------------------------------------------

subroutine pw2km2d_mig(kmm,s,img,output,pw1ttt_in,pw2ttt_in)
type(kmmod2d),      intent(in)    :: kmm
type(pw2_shot2d),   intent(in)    :: s
character(len=*),   intent(in)    :: output,pw1ttt_in,pw2ttt_in
type(image2d),      intent(inout) :: img 
integer                           :: ig
img%img=0.0
call initial_variable(kmm,s)
if (isReadData) open(13,file=output,access="direct",recl=I4*nt,action="read")
if (isSaveTTT.and.isSaveBothTTT) then
   call read_ttt(pw1ttt_in,s%proid,ttt_s)
else
   call calculate_pwttt(kmm%m%v,ttt_s,s%theta_x1,s%xref_1,s%zp1)
endif
do ig=1,ng
   if (isReadData) then
      read(13,rec=ig)trace
   else
      trace=s%seis(:,ig)
   endif
   if (isSaveTTT) then
      call read_ttt(pw2ttt_in,s%gid(ig),ttt_g)
   else
      call calculate_pwttt(kmm%m%v,ttt_g,s%theta_x2(ig),s%xref_2(ig),s%zp2(ig))
   endif
   call pw2_mig_kernal(img%img,kmm%k%wfft)
enddo
if (isReadData) close(13)
if (img%isIllum) img%illum=merge(1.0,img%illum,img%illum<0.0000000001)
call finalize_variable
end subroutine pw2km2d_mig

!------------------------------------------------------------------------

subroutine pw1_initial_variable(kmm,s)
type(kmmod2d),  intent(in) :: kmm
type(pw1_shot2d),   intent(in) :: s
call copy_variable(kmm%k%nt,nt,kmm%k%nt_conv,nt_conv) 
call copy_variable(kmm%k%isTheta,isTheta,kmm%k%isSaveTTT,isSaveTTT,&
                   kmm%k%isSaveBothTTT,isSaveBothTTT)
call copy_variable(kmm%k%isWriteData,isWriteData,kmm%k%isReadData,isReadData)
call copy_variable(kmm%k%eik_type,eik_type) 
call copy_variable(kmm%k%dt,dt,kmm%m%dx,dx,kmm%m%dz,dz)
call copy_variable(kmm%m%nx,nx,kmm%m%nz,nz,s%ng,ng)
call copy_variable(kmm%k%dfor,dfor)
if (dfor.eq."SHORT") call copy_variable(kmm%k%dfor_coe,dfor_coe)
call allocate_and_initial(ttt_s,ttt_g,dist_g,oblq_g,nz,nx)
call allocate_and_initial(trace,nt)
call allocate_and_initial(iz0,nx)
iz0=nint(kmm%k%z0_mig/dz)+1
if (isTheta) then
   theta=kmm%k%theta; theta2=(sin(theta/180.0*PI))**2.0
endif
end subroutine pw1_initial_variable

!------------------------------------------------------------------------

subroutine pw2_initial_variable(kmm,s)
type(kmmod2d),  intent(in) :: kmm
type(pw2_shot2d),   intent(in) :: s
call copy_variable(kmm%k%nt,nt,kmm%k%nt_conv,nt_conv) 
call copy_variable(kmm%k%isSaveTTT,isSaveTTT,&
                   kmm%k%isSaveBothTTT,isSaveBothTTT)
call copy_variable(kmm%k%isWriteData,isWriteData,kmm%k%isReadData,isReadData)
call copy_variable(kmm%k%eik_type,eik_type)
call copy_variable(kmm%k%dt,dt,kmm%m%dx,dx,kmm%m%dz,dz)
call copy_variable(kmm%m%nx,nx,kmm%m%nz,nz,s%ng,ng)
call copy_variable(kmm%k%dfor,dfor)
if (dfor.eq."SHORT") call copy_variable(kmm%k%dfor_coe,dfor_coe)
call allocate_and_initial(ttt_s,ttt_g,nz,nx)
call allocate_and_initial(trace,nt)
call allocate_and_initial(iz0,nx)
iz0=nint(kmm%k%z0_mig/dx)+1
end subroutine pw2_initial_variable

subroutine finalize_variable
call deallocate_and_free(ttt_s,ttt_g)
if (allocated(dist_g)) then
   call deallocate_and_free(dist_g,oblq_g)
endif
call deallocate_and_free(trace)
call deallocate_and_free(iz0)
end subroutine finalize_variable

subroutine read_ttt(ttt_in,id,ttt)
integer, intent(in) :: id
character(len=100),intent(in) :: ttt_in
real, intent(out) :: ttt(:,:)
character(len=slen) :: output
call filename(output, ttt_in, id, ".bin")
if (dfor.eq."SHORT") then
   call sf90_readself(output,ttt,nz,nx,dfor_coe)
else
   call sf90_readself(output,ttt,nz,nx)
endif
end subroutine read_ttt

subroutine calculate_ttt(s,ttt,zs,xs)
real, intent(in) :: s(:,:),zs,xs
real, intent(out):: ttt(:,:)
if (eik_type.eq."DIR") then
   call eikfd2d_dir(s,ttt,zs,xs,nz,nx,dx)
else
   call eikfd2d_fir(s,ttt,zs,xs,nz,nx,dx,dx)
endif
end subroutine calculate_ttt

subroutine calculate_pwttt(s,ttt,theta_x,xref,zp)
real, intent(in) :: s(:,:),theta_x,xref,zp
real, intent(out):: ttt(:,:)
!if (eik_type.eq."DIR") then
call pweikfd2d_dir(s,ttt,theta_x,xref,zp,nz,nx,dx)
!else
!   call eikfd2d_fir(s,ttt,zs,xs,nz,nx,dx,dx)
!endif
end subroutine calculate_pwttt

subroutine pw1_mod_kernal(refl,wfft)
real, intent(in)    :: refl(:,:)
complex, intent(in) :: wfft(:)
integer             :: ix, iz, it
real :: trace_temp(nt)
trace_temp=0.0
do ix=1,nx
   do iz=iz0(ix),nz
      it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
      if (it<=nt) then
         trace_temp(it)=refl(iz,ix)*oblq_g(iz,ix) &
                        /sqrt(dist_g(iz,ix))+trace_temp(it)
      endif
   enddo
enddo
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine pw1_mod_kernal

subroutine pw1_mig_kernal(img,isIllum,illum,wfft)
real,    intent(inout) :: img(:,:),illum(:,:)
logical, intent(in)    :: isIllum
complex, intent(in)    :: wfft(:)
integer                :: ix, iz, it
real :: trace_temp(nt)
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
do ix=1,nx
   do iz=iz0(ix),nz
      it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
      if (it<=nt) then
         img(iz,ix)=trace_temp(it)*oblq_g(iz,ix)/sqrt(dist_g(iz,ix))+img(iz,ix)
         if (isIllum) then
            illum(iz,ix)=illum(iz,ix)+1.0/sqrt(dist_g(iz,ix))
         endif
      endif
   enddo
enddo
end subroutine pw1_mig_kernal

subroutine pw1_mod_kernal_control(refl,wfft,zg,xg)
real, intent(in)    :: refl(:,:),zg,xg
complex, intent(in) :: wfft(:)
integer             :: ix, iz, it
real :: trace_temp(nt)
trace_temp=0.0
do ix=1,nx
   x=(ix-1)*dx;   gxx2=(xg-x)**2.0
   do iz=iz0(ix),nz
      z=(iz-1)*dz;   gzz2=(zg-z)**2.0;   theta_g=gxx2/(gxx2+gzz2)
      if (theta_g<=theta2) then
         it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
         if (it<=nt) then
            trace_temp(it)=refl(iz,ix)*oblq_g(iz,ix)  &
                           /sqrt(dist_g(iz,ix))+trace_temp(it)
         endif
      endif
   enddo
enddo
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine pw1_mod_kernal_control

subroutine pw1_mig_kernal_control(img,isIllum,illum,wfft,zg,xg)
real, intent(in)    :: zg,xg
real, intent(inout) :: img(:,:),illum(:,:)
logical,intent(in)  :: isIllum
complex, intent(in) :: wfft(:)
integer             :: ix,iz,it
real :: trace_temp(nt)
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
do ix=1,nx
   x=(ix-1)*dx;   gxx2=(xg-x)**2.0
   do iz=iz0(ix),nz
      z=(iz-1)*dx;   gzz2=(zg-z)**2.0;   theta_g=gxx2/(gxx2+gzz2)
      if (theta_g<=theta2) then
         it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
         if (it<=nt) then
            img(iz,ix)=trace_temp(it)*oblq_g(iz,ix) &
                       /sqrt(dist_g(iz,ix))+img(iz,ix)
            if (isIllum) then
               illum(iz,ix)=illum(iz,ix)+1.0/sqrt(dist_g(iz,ix))
            endif
         endif
      endif
   enddo
enddo
end subroutine pw1_mig_kernal_control

subroutine pw2_mod_kernal(refl,wfft)
real, intent(in)    :: refl(:,:)
complex, intent(in) :: wfft(:)
integer             :: ix, iz, it
real :: trace_temp(nt)
trace_temp=0.0
do ix=1,nx
   do iz=iz0(ix),nz
      it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
      if (it<=nt) then
         trace_temp(it)=refl(iz,ix)+trace_temp(it)
      endif
   enddo
enddo
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine pw2_mod_kernal

subroutine pw2_mig_kernal(img,wfft)
real,    intent(inout) :: img(:,:)
complex, intent(in)    :: wfft(:)
integer                :: ix, iz, it
real :: trace_temp(nt)
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
do ix=1,nx
   do iz=iz0(ix),nz
      it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
      if (it<=nt) then
         img(iz,ix)=trace_temp(it)+img(iz,ix)
      endif
   enddo
enddo
end subroutine pw2_mig_kernal

end module module_pwkm2d
