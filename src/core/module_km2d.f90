module module_km2d
use module_global, only : PI
use module_utility
use module_datatype
use module_eikfd2d_dir
use module_eikfd2d_fir
use module_fftw
use module_io
use module_sf90_io
use module_string

implicit none
real, allocatable, private :: ttt_s(:,:),ttt_g(:,:),dist_s(:,:),dist_g(:,:),&
   oblq_s(:,:),oblq_g(:,:),trace(:)
integer, allocatable, private :: iz0(:)
character(len=100),private :: eik_type,sort_in,dfor
logical, private :: isTheta,isDip,isSaveTTT,isSaveBothTTT,isWriteData,isReadData,issvs,isgvs
real, private :: dx,dz,dt,theta,dip,theta2,dip2,dip_sg,x,z,xm,zm,sxx2,&
   gxx2,szz2,gzz2,xmx2,zmz2,theta_s,theta_g,dfor_coe
integer, private :: nx, nz, nt, nt_conv, ng

private :: initial_variable, finalize_variable, read_ttt, calculate_ttt,&
   mod_kernal, mig_kernal, mod_kernal_control, mig_kernal_control       

contains

!------------------------------------------------------------------------

subroutine km2d_mod(kmm,s,output,sttt_in,gttt_in)
type(kmmod2d), intent(in)      :: kmm
character(len=*),   intent(in) :: output,sttt_in,gttt_in
type(shot2d),  intent(inout)   :: s
integer                        :: ig 
call initial_variable(kmm,s)
if (isWriteData) open(11,file=output,access="direct",recl=I4*nt,action="write")
if (sort_in.eq."CSG") then
   if (isSaveTTT) then
      if(isSaveBothTTT) then
         call read_ttt(sttt_in,s%sid(1),ttt_s)
      else 
         if (issvs) then                                          ! LSMF
            call calculate_ttt(kmm%m%vs,ttt_s,s%zs(1),s%xs(1))    ! LSMF
         else                                                     ! LSMF
            call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%xs(1))
         endif                                                    ! LSMF
      endif
   else
         if (issvs) then                                          ! LSMF
            call calculate_ttt(kmm%m%vs,ttt_s,s%zs(1),s%xs(1))    ! LSMF
         else                                                     ! LSMF
            call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%xs(1))
         endif                                                    ! LSMF
   endif
   call dist_cal(dist_s,nz,nx,dz,dx,s%zs(1),s%xs(1))   
   call oblq_cal(dist_s,oblq_s,nz,nx,dz)   
   do ig=1,ng
      if (isSaveTTT) then
         call read_ttt(gttt_in,s%gid(ig),ttt_g)
      else
         if (isgvs) then                                          ! LSMF
            call calculate_ttt(kmm%m%vs,ttt_g,s%zg(ig),s%xg(ig))  ! LSMF
         else                                                     ! LSMF
            call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%xg(ig))
         endif                                                    ! LSMF
      endif
      call dist_cal(dist_g,nz,nx,dz,dx,s%zg(ig),s%xg(ig))   
      call oblq_cal(dist_g,oblq_g,nz,nx,dz)   
      if (isTheta) then
         if ((.not.(issvs)).and.((.not.isgvs))) then               ! LSMF
            call mod_kernal_control(kmm%m%refl,kmm%k%wfft,&
            s%zs(ig),s%xs(ig),s%zg(ig),s%xg(ig))
         elseif ((.not.(issvs)).and.((isgvs))) then                ! LSMF
            call mod_kernal_control(kmm%m%refl_ps,kmm%k%wfft,&     ! LSMF
            s%zs(ig),s%xs(ig),s%zg(ig),s%xg(ig))                   ! LSMF
         endif                                                     ! LSMF
      else
         if ((.not.(issvs)).and.((.not.isgvs))) then               ! LSMF
            call mod_kernal(kmm%m%refl,kmm%k%wfft)
         elseif ((.not.(issvs)).and.((isgvs))) then                     ! LSMF
            call mod_kernal(kmm%m%refl_ps,kmm%k%wfft)              ! LSMF
         endif
      endif
      if (isWriteData) then
         write(11,rec=ig)trace
      else
         if ((.not.(issvs)).and.(isgvs)) then     ! LSMF
            s%seis_ps(:,ig)=trace                  ! LSMF
         else                                    ! LSMF
            s%seis(:,ig)=trace
         endif
      endif
   enddo
elseif (sort_in.eq."CRG") then      
   if (isSaveTTT) then
      if(isSaveBothTTT) then
         call read_ttt(gttt_in,s%gid(1),ttt_g)
      else 
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(1),s%xg(1))
      endif
   else
      call calculate_ttt(kmm%m%v,ttt_g,s%zg(1),s%xg(1))
   endif
   call dist_cal(dist_g,nz,nx,dz,dx,s%zg(1),s%xg(1))   
   call oblq_cal(dist_g,oblq_g,nz,nx,dx)   
   do ig=1,ng
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(ig),ttt_s)
      else 
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(ig),s%xs(ig))
      endif
      call dist_cal(dist_s,nz,nx,dz,dx,s%zs(ig),s%xs(ig))   
      call oblq_cal(dist_s,oblq_s,nz,nx,dz)   
      if (isTheta) then
         call mod_kernal_control(kmm%m%refl,kmm%k%wfft,&
            s%zs(ig),s%xs(ig),s%zg(ig),s%xg(ig))
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
         if (issvs) then
            call calculate_ttt(kmm%m%vs,ttt_s,s%zs(ig),s%xs(ig))
         else
            call calculate_ttt(kmm%m%v,ttt_s,s%zs(ig),s%xs(ig))
         endif
         if (isgvs) then
            call calculate_ttt(kmm%m%vs,ttt_g,s%zg(ig),s%xg(ig))
         else 
            call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%xg(ig))
         endif
      endif
      call dist_cal(dist_s,nz,nx,dz,dx,s%zs(ig),s%xs(ig))
      call oblq_cal(dist_s,oblq_s,nz,nx,dz)   
      call dist_cal(dist_g,nz,nx,dz,dx,s%zg(ig),s%xg(ig))
      call oblq_cal(dist_g,oblq_g,nz,nx,dz)   
      if (isTheta) then
         if ((.not.issvs).and.(.not.isgvs)) then            !LSMF
            call mod_kernal_control(kmm%m%refl,kmm%k%wfft,&
                 s%zs(ig),s%xs(ig),s%zg(ig),s%xg(ig))
         elseif ((.not.issvs).and.(isgvs)) then             !LSMF
            call mod_kernal_control(kmm%m%refl_ps,kmm%k%wfft,& !LMSF
                 s%zs(ig),s%xs(ig),s%zg(ig),s%xg(ig))       !LSMF
         endif
      else
         if ((.not.issvs).and.(.not.isgvs)) then       !LSMF
           ! call mod_kernal(kmm%m%refl,kmm%k%wfft)
            call mod_kernal_zo(kmm%m%refl,kmm%k%wfft,s%xs(ig),s%xg(ig))
         elseif ((.not.issvs).and.(isgvs)) then        !LSMF
            call mod_kernal(kmm%m%refl_ps,kmm%k%wfft)  !LSMF
         endif
      endif
      if (isWriteData) then
         write(11,rec=ig)trace
      else
         if ((.not.issvs).and.(isgvs))then
            s%seis_ps(:,ig)=trace
         else
            s%seis(:,ig)=trace           ! LSMF
         endif
      endif
   enddo
endif
if (isWriteData) close(11)
call finalize_variable
end subroutine km2d_mod

!------------------------------------------------------------------------

subroutine km2d_mig(kmm,s,img,output,sttt_in,gttt_in)
type(kmmod2d),    intent(in)    :: kmm
type(shot2d),     intent(in)    :: s
character(len=*), intent(in)    :: output,sttt_in,gttt_in
type(image2d),    intent(inout) :: img
integer                         :: ig
img%img=0.0
call initial_variable(kmm,s)

if (isReadData) open(13,file=output,access="direct",recl=I4*nt,action="read")
if (img%isIllum) img%illum=0.0
if (sort_in.eq."CSG") then
   if (isSaveTTT) then
      if(isSaveBothTTT) then
         call read_ttt(sttt_in,s%sid(1),ttt_s)
      else 
         if (issvs) then                                        ! LSMF 
            call calculate_ttt(kmm%m%vs,ttt_s,s%zs(1),s%xs(1))  ! LSMF
         else                                                   ! LSMF
            call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%xs(1))  
         endif                                                  ! LSMF
      endif
   else
      if (issvs) then                                           ! LSMF
         call calculate_ttt(kmm%m%vs,ttt_s,s%zs(1),s%xs(1))     ! LSMF
      else                                                      ! LSMF
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(1),s%xs(1))
      endif                                                     ! LSMF
   endif
   call dist_cal(dist_s,nz,nx,dz,dx,s%zs(1),s%xs(1))   
   call oblq_cal(dist_s,oblq_s,nz,nx,dx)   
   do ig=1,ng
      if (isReadData) then
         read(13,rec=ig)trace
      else
         trace=s%seis(:,ig)
      endif
      if (isSaveTTT) then
         call read_ttt(gttt_in,s%gid(ig),ttt_g)
      else 
         if (isgvs) then                                        ! LSMF
            call calculate_ttt(kmm%m%vs,ttt_g,s%zg(ig),s%xg(ig))! LSMF
         else
            call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%xg(ig))
         endif                                                  ! LSMF
      endif
      call dist_cal(dist_g,nz,nx,dz,dx,s%zg(ig),s%xg(ig))   
      call oblq_cal(dist_g,oblq_g,nz,nx,dx)   
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,&
            kmm%k%wfft,s%zs(ig),s%xs(ig),s%zg(ig),s%xg(ig))
      else
         call mig_kernal(img%img,img%isIllum,img%illum,kmm%k%wfft)
      endif
   enddo
elseif (sort_in.eq."CRG") then      
   if (isSaveTTT) then
      if(isSaveBothTTT) then
         call read_ttt(gttt_in,s%gid(1),ttt_g)
      else 
         call calculate_ttt(kmm%m%v,ttt_g,s%zg(1),s%xg(1))
      endif
   else
      call calculate_ttt(kmm%m%v,ttt_g,s%zg(1),s%xg(1))
   endif
   call dist_cal(dist_g,nz,nx,dz,dx,s%zg(1),s%xg(1))   
   call oblq_cal(dist_g,oblq_g,nz,nx,dx)   
   do ig=1,ng
      if (isReadData) then
         read(13,rec=ig)trace
      else
         trace=s%seis(:,ig)
      endif
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(ig),ttt_s)
      else 
         call calculate_ttt(kmm%m%v,ttt_s,s%zs(ig),s%xs(ig))
      endif
      call dist_cal(dist_s,nz,nx,dz,dx,s%zs(ig),s%xs(ig))   
      call oblq_cal(dist_s,oblq_s,nz,nx,dz)   
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,&
            kmm%k%wfft,s%zs(ig),s%xs(ig),s%zg(ig),s%xg(ig))
      else
         call mig_kernal(img%img,img%isIllum,img%illum,kmm%k%wfft)
      endif
   enddo
else
   do ig=1,ng
      !write(*,*) ig
      if (isReadData) then
         read(13,rec=ig)trace
      else
         trace=s%seis(:,ig)
      endif
      if (isSaveTTT) then
         call read_ttt(sttt_in,s%sid(ig),ttt_s)
         call read_ttt(gttt_in,s%gid(ig),ttt_g)
      else 
         if (issvs) then
            call calculate_ttt(kmm%m%vs,ttt_s,s%zs(ig),s%xs(ig))
         else
            call calculate_ttt(kmm%m%v,ttt_s,s%zs(ig),s%xs(ig))
         endif
         if (isgvs) then
            call calculate_ttt(kmm%m%vs,ttt_g,s%zg(ig),s%xg(ig))
         else
            call calculate_ttt(kmm%m%v,ttt_g,s%zg(ig),s%xg(ig))
         endif
      endif
      call dist_cal(dist_s,nz,nx,dz,dx,s%zs(ig),s%xs(ig))
      call oblq_cal(dist_s,oblq_s,nz,nx,dz)
      call dist_cal(dist_g,nz,nx,dz,dx,s%zg(ig),s%xg(ig))
      call oblq_cal(dist_g,oblq_g,nz,nx,dx)
      if (isTheta) then
         call mig_kernal_control(img%img,img%isIllum,img%illum,&
            kmm%k%wfft,s%zs(ig),s%xs(ig),s%zg(ig),s%xg(ig))
      else
         !call mig_kernal(img%img,img%isIllum,img%illum,kmm%k%wfft)
         call mig_kernal_zo(img%img,img%isIllum,img%illum,kmm%k%wfft,s%xs(ig),s%xg(ig))
      endif
   enddo
endif
if (isReadData) close(13)
if (img%isIllum) img%illum=merge(1.0,img%illum,img%illum<0.000000001)
call finalize_variable
end subroutine km2d_mig

!------------------------------------------------------------------------

subroutine initial_variable(kmm,s)
type(kmmod2d),  intent(in) :: kmm
type(shot2d),   intent(in) :: s
call copy_variable(kmm%k%nt,nt,kmm%k%nt_conv,nt_conv) 
call copy_variable(kmm%k%isTheta,isTheta,kmm%k%isDip,isDip)
call copy_variable(kmm%k%isWriteData,isWriteData,kmm%k%isReadData,isReadData)
call copy_variable(kmm%k%isSaveTTT,isSaveTTT,kmm%k%isSaveBothTTT,isSaveBothTTT)
call copy_variable(kmm%k%eik_type,eik_type,kmm%k%sort_in,sort_in)
call copy_variable(kmm%k%dt,dt,kmm%m%dx,dx,kmm%m%dz,dz)
call copy_variable(kmm%m%nx,nx,kmm%m%nz,nz,s%ng,ng)
call copy_variable(kmm%k%dfor,dfor)
if (allocated(kmm%m%vs)) then               ! LSMF 
   call copy_variable(kmm%k%issvs,issvs)  ! LSMF
   call copy_variable(kmm%k%isgvs,isgvs)  ! LSMF
else                                      ! LSMF
   issvs=.false.                          ! LSMF
   isgvs=.false.                          ! LSMF
endif                                     ! LSMF



if (dfor.eq."SHORT") call copy_variable(kmm%k%dfor_coe,dfor_coe)
call allocate_and_initial(ttt_s,ttt_g,dist_s,dist_g,oblq_s,oblq_g,nz,nx)
call allocate_and_initial(trace,nt)
call allocate_and_initial(iz0,nx)
iz0=nint(kmm%k%z0_mig/dx)+1
if (isTheta) theta=kmm%k%theta; theta2=(sin(theta/180.0*PI))**2.0
if (isDip) dip=kmm%k%dip; dip2=(sin(dip/180.0*PI))**2.0
end subroutine initial_variable

subroutine finalize_variable
call deallocate_and_free(ttt_s,ttt_g,dist_s,dist_g,oblq_s,oblq_g)
call deallocate_and_free(trace)
call deallocate_and_free(iz0)
end subroutine finalize_variable

subroutine read_ttt(ttt_in,id,ttt)
integer, intent(in) :: id
character(len=*),intent(in) :: ttt_in
real, intent(out)   :: ttt(:,:)
character(len=slen) :: ttt_name
call filename(ttt_name, ttt_in, id, ".bin")
if (dfor.eq."SHORT") then
   call sf90_readself(ttt_name,ttt,nz,nx,dfor_coe)
else
   call sf90_readself(ttt_name,ttt,nz,nx)
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

subroutine mod_kernal_zo(refl,wfft,xs,xg)
real,intent(in)     :: xs,xg
real, intent(in)    :: refl(:,:)
complex, intent(in) :: wfft(:)
real :: trace_temp(nt)
integer             :: ix, iz, it
integer             :: ix_temp
trace_temp=0.0
do ix=nint((xs+xg)/(2*dx))-35,nint((xs+xg)/(2*dx))+35
   ix_temp=ix
   if (ix<1) then
      ix_temp=1
   elseif (ix>nx) then
      ix_temp=nx
   endif
   do iz=iz0(ix_temp),nz
      it=nint((ttt_s(iz,ix_temp)+ttt_g(iz,ix_temp))/dt)*2+1-375
      if (it<1) then
         it=1
      endif
      if (it<=nt) then
         trace_temp(it)=refl(iz,ix_temp)*oblq_s(iz,ix_temp)*oblq_g(iz,ix_temp)&
            /sqrt(dist_s(iz,ix_temp)*dist_g(iz,ix_temp))+trace_temp(it)
      endif
   enddo
enddo
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine mod_kernal_zo

subroutine mod_kernal(refl,wfft)
real, intent(in)    :: refl(:,:)
complex, intent(in) :: wfft(:)
real :: trace_temp(nt)
integer             :: ix, iz, it
trace_temp=0.0
do ix=1,nx
   do iz=iz0(ix),nz
      it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
      if (it<=nt) then
         trace_temp(it)=refl(iz,ix)*oblq_s(iz,ix)*oblq_g(iz,ix)&
            /sqrt(dist_s(iz,ix)*dist_g(iz,ix))+trace_temp(it)
      endif
   enddo
enddo
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine mod_kernal

!=========================================================================
! Add by Bowen Guo
subroutine mig_kernal_zo(img,isIllum,illum,wfft,xs,xg)

real,intent(in)        :: xs,xg
real,    intent(inout) :: img(:,:),illum(:,:)
logical, intent(in)    :: isIllum
complex, intent(in)    :: wfft(:)
integer                :: ix, iz, it!0, it1
real                   :: trace_temp(nt)
integer                :: ix_temp
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
!write(*,*) nint((xs+xg)/(2*dx))
do ix=nint((xs+xg)/(2*dx))-5,nint((xs+xg)/(2*dx))+5
   do iz=iz0(ix),nz
      ix_temp=ix
      if (ix<1) then
         ix_temp=1
      elseif (ix>nx) then
         ix_temp=nx
      endif
      it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)*2+1-375
      if (it<1) then
         it=1
      endif
     ! write(*,*) "it==========",it
     ! write(*,*) "dt==========",dt
      if (it<=nt) then
         img(iz,ix)=trace_temp(it)*oblq_s(iz,ix)*oblq_g(iz,ix)&
            /sqrt(dist_s(iz,ix)*dist_g(iz,ix))+img(iz,ix)
         if (isIllum) then
            illum(iz,ix)=illum(iz,ix)+1.0/sqrt(dist_s(iz,ix)*dist_g(iz,ix))
         endif
      endif 
   enddo
enddo
end subroutine mig_kernal_zo
!==========================================================================



subroutine mig_kernal(img,isIllum,illum,wfft)
real,    intent(inout) :: img(:,:),illum(:,:)
logical, intent(in)    :: isIllum
complex, intent(in)    :: wfft(:)
integer                :: ix, iz, it!0, it1
real :: trace_temp(nt)
!write(*,*) nt_conv
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
do ix=1,nx
   do iz=iz0(ix),nz
     ! it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1-375
      it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
     ! write(*,*) "it==========",it
     ! write(*,*) "dt==========",dt
      if (it<=nt) then
         img(iz,ix)=trace_temp(it)*oblq_s(iz,ix)*oblq_g(iz,ix)&
            /sqrt(dist_s(iz,ix)*dist_g(iz,ix))+img(iz,ix)
         if (isIllum) then
            illum(iz,ix)=illum(iz,ix)+1.0/sqrt(dist_s(iz,ix)*dist_g(iz,ix))
         endif
      endif 
   enddo
enddo
end subroutine mig_kernal

subroutine mod_kernal_control(refl,wfft,zs,xs,zg,xg)
real, intent(in)    :: refl(:,:),zs,xs,zg,xg
complex, intent(in) :: wfft(:)
real :: trace_temp(nt)
integer             :: ix, iz, it
trace_temp=0.0
xm=0.5*(xs+xg);   zm=0.5*(zs+zg)  
do ix=1,nx
   x=(ix-1)*dx; sxx2=(xs-x)**2.0; gxx2=(xg-x)**2.0; xmx2=(xm-x)**2.0
   do iz=iz0(ix),nz
      z=(iz-1)*dz; szz2=(zs-z)**2.0; gzz2=(zg-z)**2.0
      theta_s=sxx2/(sxx2+szz2);   theta_g=gxx2/(gxx2+gzz2)
      if (theta_s<=theta2 .and. theta_g<=theta2) then
         zmz2=(zm-z)**2.0;   dip_sg=xmx2/(xmx2+zmz2)
         if (dip_sg<dip2) then
            it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
            if (it<=nt) then
               trace_temp(it)=refl(iz,ix)*oblq_s(iz,ix)*oblq_g(iz,ix)&
                  /sqrt(dist_s(iz,ix)*dist_g(iz,ix))+trace_temp(it)
            endif
         endif
      endif
   enddo
enddo
call conv_wavelet(trace_temp,trace,wfft,nt,nt_conv)
end subroutine mod_kernal_control

subroutine mig_kernal_control(img,isIllum,illum,wfft,zs,xs,zg,xg)
real, intent(in)    :: zs,xs,zg,xg
real, intent(inout) :: img(:,:),illum(:,:)
logical, intent(in) :: isIllum
complex, intent(in) :: wfft(:)
integer             :: ix,iz,it
real :: trace_temp(nt)
call xcorr_wavelet(trace,trace_temp,wfft,nt,nt_conv)
xm=0.5*(xs+xg);  zm=0.5*(zs+zg)  
do ix=1,nx
   x=(ix-1)*dx; sxx2=(xs-x)**2.0; gxx2=(xg-x)**2.0; xmx2=(xm-x)**2.0
   do iz=iz0(ix),nz
      z=(iz-1)*dz; szz2=(zs-z)**2.0; gzz2=(zg-z)**2.0
      theta_s=sxx2/(sxx2+szz2);   theta_g=gxx2/(gxx2+gzz2)
      if (theta_s<=theta2 .and. theta_g<=theta2) then
         zmz2=(zm-z)**2.0;   dip_sg=xmx2/(xmx2+zmz2)
         if (dip_sg<dip2) then
            it=nint((ttt_s(iz,ix)+ttt_g(iz,ix))/dt)+1
            if (it<=nt) then
               img(iz,ix)=trace_temp(it)*oblq_s(iz,ix)*oblq_g(iz,ix)&
                  /sqrt(dist_s(iz,ix)*dist_g(iz,ix))+img(iz,ix)
               if (isIllum) then
                  illum(iz,ix)=illum(iz,ix)+1.0/sqrt(dist_s(iz,ix)*dist_g(iz,ix))
               endif
            endif
         endif
      endif
   enddo
enddo
end subroutine mig_kernal_control

end module module_km2d
