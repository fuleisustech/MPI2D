! Do marchenko-related project
! @version 1 2014-08-14
! @author Bowen Guo



module module_marchenko

use module_global
use module_utility
use module_fftw
use module_array
use module_parser
use module_io

implicit none


integer,private                ::  is,ig,mute_time,iz,ix,fid_s,fid_g,ipp
integer,private,allocatable    ::  ttt(:,:)
real,private                   ::  img_ip_temp
real,allocatable,private       ::  temp_conv(:),temp(:),temp_stack(:),green_s(:),green_g(:),hann(:)
complex,allocatable,private    ::  G_down_fft(:,:),f0m_fft(:,:),green_s_fft(:)

character(len=slen),private    ::  fn_temp



private initial_variable, finalize_variable


contains

!===============================================================================================
subroutine forward_modeling_green_whole(img,fn_green_s,fn_ttt,data_out,iss,ng,nt,nz,nx,nx_s,nx_e,nz_s,&
                                        nz_e,isgttt)
real,intent(in)      :: img(:,:)
integer,intent(in)   :: iss,ng,nt,nz,nx,nx_s,nx_e,nz_s,nz_e
logical,intent(in)   :: isgttt
character(len=*)     :: fn_green_s,fn_ttt

real,intent(out)     :: data_out(:,:)




call allocate_and_initial(temp_conv,2*nt+1)
call allocate_and_initial(green_g,nt)

fid_g=2000+rank
! Read source side Green's function
if (.not.isgttt) then
   call allocate_and_initial(green_s,nt)
   call allocate_and_initial(green_s_fft,2*nt+1)
   call filename(fn_temp,fn_green_s,iss,".bin")
   fid_s=1000+rank
   open(fid_s,file=fn_temp,access="direct",recl=nt*I4)
else
   call allocate_and_initial(ttt,nz,nx)
   call filename(fn_temp,fn_ttt,iss,".bin")
   call read_binfile(fn_temp,ttt,nz,nx)
endif

do ig=1,ng
   call filename(fn_temp,fn_green_s,ig,".bin")
   open(fid_g,file=fn_temp,access="direct",recl=nt*I4)
   !$OMP PARALLEL DO PRIVATE(ix,iz,ipp)
   do ix=nx_s,nx_e
      do iz=nz_s,nz_e
         ipp=(ix-1)*nz+iz
         read(fid_g,rec=ipp) green_g
         if (.not.isgttt) then
            if (ig.eq.1) then
               read(fid_s,rec=ipp) green_s 
               call fft_wavelet(green_s,green_s_fft,nt,2*nt+1)
            endif
            call conv_wavelet(green_g,temp_conv,green_s_fft,nt,2*nt+1)
         else 
            call shift_backward(green_g,temp_conv,ttt(iz,ix))
         endif
         data_out(:,ig)=temp_conv*img(iz,ix)
      enddo
   enddo
   !$OMP END PARALLEL DO
enddo

if (.not.isgttt) close(fid_s)
close(fid_g)


call deallocate_and_free(temp_conv)
call deallocate_and_free(green_g)

if (.not.isgttt) then
   call deallocate_and_free(green_s)
   call deallocate_and_free(green_s_fft)
else
   call deallocate_and_free(ttt)
endif

end subroutine forward_modeling_green_whole
!======================================================================================================
  
subroutine migrate_green_whole(data_in,fn_green_s,fn_ttt,img_out,iss,ng,nt,nz,nx,nx_s,nx_e,nz_s,&
                               nz_e,isgttt)

real,intent(in)     :: data_in(:,:)
integer,intent(in)  :: iss,ng,nt,nz,nx,nx_s,nx_e,nz_s,nz_e
logical,intent(in)  :: isgttt
character(len=*) :: fn_green_s,fn_ttt

real,intent(out)    :: img_out(:,:)

call allocate_and_initial(temp_conv,2*nt+1)
call allocate_and_initial(green_g,nt)

img_out=0.0
fid_g=rank+4000

! Read source side Green's function
if (.not.isgttt) then
   fid_s=rank+3000
   call allocate_and_initial(green_s,nt)
   call allocate_and_initial(green_s_fft,2*nt+1)
   call filename(fn_temp,fn_green_s,iss,".bin")
   open(fid_s,file=fn_temp,access="direct",recl=nt*I4)
else
   call allocate_and_initial(ttt,nz,nx)
   call filename(fn_temp,fn_ttt,iss,".bin")
   call read_binfile(fn_temp,ttt,nz,nx)
endif

do ig=1,ng
   call filename(fn_temp,fn_green_s,ig,".bin")
   open(fid_g,file=fn_temp,access="direct",recl=nt*I4)
   !$OMP PARALLEL DO PRIVATE(ix,iz,ipp)
   do ix=nx_s,nx_e
      do iz=nz_s,nz_e
         ipp=(ix-1)*nz+iz
         read(fid_g,rec=ipp) green_g
         if (.not.isgttt) then
            if (ig.eq.1) then
               read(fid_s,rec=ipp) green_s
               call fft_wavelet(green_s,green_s_fft,nt,2*nt+1)
            endif
            call conv_wavelet(green_g,temp_conv,green_s_fft,nt,2*nt+1)
         else
            call shift_backward(green_g,temp_conv,ttt(iz,ix))
         endif
         img_out(ix,iz)=img_out(ix,iz)+dot_product(data_in(:,ig),temp_conv)
      enddo
   enddo
   !$OMP END PARALLEL DO
enddo

if (.not.isgttt) close(fid_s)
close(fid_g)
   
call deallocate_and_free(temp_conv)
call deallocate_and_free(green_g)

if (.not.isgttt) then
   call deallocate_and_free(green_s)
   call deallocate_and_free(green_s_fft)
else
   call deallocate_and_free(ttt)
endif

end subroutine migrate_green_whole

!=================================================================================
subroutine forward_green_marchenko(img_ip,G_down,ttt_int,data_out,ns,ng,nt,ntt,map_x,map_z,ip,isgttt)

real,intent(in)     :: G_down(:,:),img_ip
integer,intent(in)  :: ns,ng,ip,nt,ntt,map_x(:),map_z(:),ttt_int(:,:,:)
logical,intent(in)  :: isgttt
real,intent(out)    :: data_out(:,:,:)

call allocate_and_initial(temp_conv,2*nt+1)
if (.not.isgttt) then 
   call allocate_and_initial(G_down_fft,2*nt+1,ns)
   ! fft G_down
   do is=1,ns
      call fft_wavelet(G_down(:,is),G_down_fft(:,is),nt,2*nt+1)
   enddo
endif


do is=1,ns
   do ig=1,ng
      if (isgttt) then
         call shift_backward(G_down(:,is),temp_conv,ttt_int(map_z(ip),map_x(ip),ig))
      else
         call conv_wavelet(G_down(:,is),temp_conv,G_down_fft(:,ig),nt,2*nt+1)
      endif
      data_out(:,ig,is)=temp_conv*img_ip
      data_out(1:nt-1,ig,is)=0.0
   enddo
enddo

call deallocate_and_free(temp_conv)
if (allocated(G_down_fft)) call deallocate_and_free(G_down_fft)

end subroutine forward_green_marchenko
!=================================================================================
subroutine GDM_green_marchenko(img_ip,G_down,ttt_int,data_in,ns,ng,nt,ntt,map_x,map_z,ip,isgttt)
real,intent(in)     :: G_down(:,:),data_in(:,:,:)
integer,intent(in)  :: ns,ng,ip,nt,ntt,map_x(:),map_z(:),ttt_int(:,:,:)
logical,intent(in)  :: isgttt
real,intent(out)    :: img_ip

call allocate_and_initial(temp_conv,2*nt-1)
!if (.not.isgttt) then
   ! call allocate_and_initial(G_down_fft,nt,ns)
   ! fft G_down
   !do is=1,ns
   !   call fft_wavelet(G_down(:,is),G_down_fft(:,is),nt,2*nt+1)
   !enddo
!endif


img_ip=0.0
do is=1,ns
   img_ip_temp=0.0
   do ig=1,ng
      if (isgttt) then
         call shift_backward(G_down(:,is),temp_conv,ttt_int(map_z(ip),map_x(ip),ig))
      else
         call conv_gf(G_down(:,is),G_down(:,ig),temp_conv,nt,2*nt-1)
      endif
      img_ip_temp=img_ip_temp+dot_product(temp_conv,data_in(nt:ntt,ig,is))
    enddo
   img_ip=img_ip+img_ip_temp
enddo
call deallocate_and_free(temp_conv)
!if (allocated(G_down_fft)) call deallocate_and_free(G_down_fft)
end subroutine GDM_green_marchenko

!=====================================================================================
subroutine cal_G_down(G_down,data_in,f0m,ns,ng,nt,ntt,ttt_int,map_x,map_z,ip,period,dt)
real,intent(in)     :: data_in(:,:,:),f0m(:,:),dt
integer,intent(in)  :: ns,ng,nt,ntt,ip,period,map_x(:),map_z(:),ttt_int(:,:,:)
real,intent(out)    :: G_down(:,:)

call initial_variable(hann,period)
call allocate_and_initial(temp,ntt)
call allocate_and_initial(temp_stack,ntt)
!call allocate_and_initial(f0m_fft,ntt,ns)
!do is=1,ns
!   call fft_wavelet(f0m(:,is),f0m_fft(:,is),ntt,2*ntt+1)
!enddo

do is=1,ns
   temp_stack=0.0
   do ig=1,ng
      call xcorr_gf(data_in(:,ig,is),f0m(:,ig),temp,ntt,2*ntt-1)
      temp_stack=temp_stack+temp
   enddo
   call unpad_time_trace(temp_stack,G_down(:,is),nt,ntt)
   mute_time=ttt_int(map_z(ip),map_x(ip),is)-period
  ! mute_time=ttt_int(map_z(ip),map_x(ip),is)-int(period/2*3)
   call mute_before(G_down(:,is),mute_time,hann)
enddo
call deallocate_and_free(temp)
call deallocate_and_free(temp_stack)
call deallocate_and_free(f0m_fft)
call finalize_variable()
end subroutine cal_G_down

!=======================================================================================
subroutine cal_f0m(data_in,f0m,ns,ng,nt,ntt,ttt_int,map_x,map_z,ip,period,dt)
real,intent(in)    :: data_in(:,:,:),dt
real,intent(out)   :: f0m(:,:)
integer,intent(in) :: ns,ng,ip,nt,ntt,period,map_x(:),map_z(:),ttt_int(:,:,:)

call initial_variable(hann,period)
call allocate_and_initial(temp,ntt)
f0m=0.0

do is=1,ns
   do ig=1,ng
      call shift_forward(data_in(:,ig,is),temp,ttt_int(map_z(ip),map_x(ip),ig))
      f0m(:,is)=f0m(:,is)+temp
   enddo
   !mute_time=ttt_int(map_z(ip),map_x(ip),is)+nt-1+int(period/2*3)
   mute_time=ttt_int(map_z(ip),map_x(ip),is)+nt-1+period
   call mute_after(f0m(:,is),mute_time,hann)
enddo
call deallocate_and_free(temp)
call finalize_variable()
end subroutine cal_f0m
!=======================================================================================
subroutine initial_variable(hann,period)

integer,intent(in)        :: period
real,allocatable          :: hann(:)

call allocate_and_initial(hann,int(period/2))
call hanning_window(period,hann)

!call write_binfile("hann.bin",hann,int(period/2))


end subroutine initial_variable

!=======================================================================================
subroutine finalize_variable()

call deallocate_and_free(hann)

end subroutine finalize_variable

!=======================================================================================









end module module_marchenko

















































