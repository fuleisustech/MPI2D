module module_utility
use module_global, only : slen, date, time_now
use module_array
use module_string
use module_io
!use module_source

!use module_fftw       ! Bowen Guo
implicit none
!include "fftw3.f"
include 'fftw3.f'

integer, private :: i1, i2, i3, n11, n12, n21, n22, n31, n32
character(len=slen),private :: output



interface dir_taper   ! Lei Fu
   module procedure dir_taper_remove
end interface dir_taper


interface Matrix_hil   ! Lei Fu
   module procedure Matrix_hil_1d
   module procedure Matrix_hil_2d
end interface Matrix_hil

interface d_integerN   ! Lei Fu
   module procedure d_integer1d_1
   module procedure d_integer1d_N
   module procedure d_integer2d_N
end interface d_integerN

interface bp_filter   ! Lei Fu
   module procedure bp_filter1
!   module procedure wiener_filter1
end interface bp_filter

interface nextpow2   ! Lei Fu
   module procedure nextpow2_1d
end interface nextpow2

interface centroid   ! Lei Fu
   module procedure centroid_1d
end interface centroid

interface mx_peak   ! Lei Fu
   module procedure mx_peak_1d
end interface mx_peak

interface xcorr   ! Lei Fu
   module procedure xcorr_1d
   module procedure xcorr_1dlag
   module procedure xcorr_1d_cor
end interface xcorr

interface dotp   ! Lei Fu
   module procedure dotp_1d
end interface dotp

interface norm_n   ! Lei Fu
   module procedure norm_2
end interface norm_n


interface tau_cal   ! Lei Fu
   module procedure tau_cal_1d
end interface tau_cal

interface res_ewi   ! Lei Fu  data_res for envelope waveform inversion
   module procedure res_ewi_trand
   module procedure res_ewi_trand_ms
end interface res_ewi

interface res_fwi   ! Lei Fu  data_res for fwi
   module procedure res_fwi_wcc
   module procedure res_fwi_trand
   module procedure res_fwi_trand_ms
end interface res_fwi


interface data_res   ! Lei Fu : data_res for wt
   module procedure data_res_wt
   module procedure data_res_empi
   module procedure data_res_phainv
!   module procedure data_res_phainv_wiener ! using wiener filter for bandpass
   module procedure data_res_wt_eik
   module procedure data_res_wt_hilb
   module procedure data_res_wt_multioffset
   module procedure data_res_fti
   module procedure data_res_fti_new
end interface data_res

interface data_der   ! Lei Fu
   module procedure data_der_1d
end interface data_der

interface data_shift   ! Lei Fu
   module procedure data_shift_1d
end interface data_shift

interface trace_dis   ! Lei Fu
   module procedure trace_dis_1d
end interface trace_dis

interface window_shot   ! Lei Fu
   module procedure window_shot_2d
end interface window_shot

interface norm_trace   ! Lei Fu
   module procedure norm_trace_2d
end interface norm_trace

interface norm_shot   ! Lei Fu
   module procedure norm_shot_2d
end interface norm_shot

interface csg_window  ! Lei Fu
   module procedure csg_window_2d
   module procedure csg_win_early_2d
end interface csg_window

interface csg_cut  ! Lei Fu
   module procedure csg_cut_OBC
   module procedure csg_cut_2d
   module procedure csg_cut_2d0
end interface csg_cut

interface res_data_cut  ! Lei Fu
   module procedure res_data_cut_near
end interface res_data_cut

interface csg_cut_streamer  ! Lei Fu
   module procedure csg_cut_2d_redsea
end interface csg_cut_streamer


interface csg_mute   ! Lei Fu
   module procedure csg_mute_2d
   module procedure csg_mute_tt
end interface csg_mute

interface tapering   ! Lei Fu
   module procedure tapering_1d
end interface tapering

interface dist_cal
   module procedure dist_cal_2d
   module procedure dist_cal_3d
end interface dist_cal

interface oblq_cal
   module procedure oblq_cal_2d
   module procedure oblq_cal_3d
end interface oblq_cal

interface cut
   module procedure cut_int_1d
   module procedure cut_int_2d
   module procedure cut_int_3d
   module procedure cut_real_1d
   module procedure cut_real_2d
   module procedure cut_real_3d
end interface cut

interface padarray
  module procedure padarray_real_1d
  module procedure padarray_real_2d
  module procedure padarray_real_3d
  module procedure padarray_int_1d
  module procedure padarray_int_2d
  module procedure padarray_int_3d
end interface padarray

interface smooth
  module procedure smooth_real_1d_v1
  module procedure smooth_real_1d_v2
  module procedure smooth_real_2d_v1
  module procedure smooth_real_2d_v2
  module procedure smooth_real_3d_v1
  module procedure smooth_real_3d_v2
  module procedure smooth_gk_2d_v1
  module procedure smooth_gk_2d_v2
end interface smooth

interface cut_and_pad
   module procedure cut_and_pad_int1d
   module procedure cut_and_pad_int2d
   module procedure cut_and_pad_real2d
   module procedure cut_and_pad_real3d
end interface cut_and_pad

interface copy_variable
   module procedure copy_variable_1i
   module procedure copy_variable_2i
   module procedure copy_variable_3i
   module procedure copy_variable_4i
   module procedure copy_variable_5i
   module procedure copy_variable_1r
   module procedure copy_variable_2r
   module procedure copy_variable_3r
   module procedure copy_variable_4r
   module procedure copy_variable_5r
   module procedure copy_variable_6r
   module procedure copy_variable_1l
   module procedure copy_variable_2l
   module procedure copy_variable_3l
   module procedure copy_variable_1c
   module procedure copy_variable_2c
   module procedure copy_variable_3c
   module procedure copy_variable_1ch
   module procedure copy_variable_2ch
   module procedure copy_variable_3ch
end interface copy_variable

interface wf_assign
   module procedure wf2d_assign_1
   module procedure wf2d_assign_2
   module procedure wf2d_assign_3
   module procedure wf2d_assign_4
   module procedure wf2d_assign_5
   module procedure wf3d_assign_1
   module procedure wf3d_assign_2
   module procedure wf3d_assign_3
   module procedure wf3d_assign_4
   module procedure wf3d_assign_5
end interface wf_assign

interface wf_refresh
   module procedure wf2d_refresh_1
   module procedure wf2d_refresh_2
   module procedure wf2d_refresh_3
   module procedure wf2d_refresh_4
   module procedure wf2d_refresh_5
   module procedure wf3d_refresh_1
   module procedure wf3d_refresh_2
   module procedure wf3d_refresh_3
   module procedure wf3d_refresh_4
   module procedure wf3d_refresh_5
end interface wf_refresh

interface snapshot
   module procedure snapshot2d_wi_bc
   module procedure snapshot2d_wo_bc
   module procedure snapshot3d_wi_bc
   module procedure snapshot3d_wo_bc
end interface snapshot

interface average
   module procedure average_real2
   module procedure average_real3
   module procedure average_real4
   module procedure average_real5
   module procedure average_real8
   module procedure average_1d_real
   module procedure average_2d_real
   module procedure average_3d_real
   module procedure average_4d_real
end interface average

interface accurate_value
   module procedure accurate_value_real_1d
   module procedure accurate_value_real_2d
   module procedure accurate_value_real_3d
end interface accurate_value

interface output_illum
   module procedure output_illum_2d
   module procedure output_illum_3d
end interface output_illum

interface output_prestack_image
   module procedure output_prestack_image_2d
   module procedure output_prestack_image_3d
end interface output_prestack_image

interface shift_forward
   module procedure shift_forward_other
   module procedure shift_forward_self
end interface shift_forward
   
interface shift_backward
   module procedure shift_backward_other
   module procedure shift_backward_self
end interface shift_backward
private cut_core1

contains 


subroutine dir_taper_remove(dpre,dobs,nt,ng,xs,xg,dt,n,T_p)
real                :: dpre(:,:),dobs(:,:)
integer, intent(in) :: nt,ng,n,T_p
real, intent(in)    :: dt
integer             :: ig,it,m
real                :: xs(:,:),xg(:,:),pi
real, allocatable   :: tp0(:),tp1(:),taper(:,:)
integer, allocatable :: tt(:)

pi = 3.1415926
m = 2*n+1
allocate(tt(ng))
allocate(tp0(m))
allocate(tp1(T_p))
allocate(taper(nt,ng))
taper = 1.0
tp1 = 1.0

do it = 1,m
   tp0(it) = 0.5*(1-cos(2*pi*it/m))
enddo

if (m > T_p) then
  if (rank.eq.0) then
    write(*,*) "tapering window parameter error"
  endif
endif

tp1(1:n) = tp0(1:n)
tp1(T_p-n:T_p) = tp0(n+1:m)
do ig=1,ng
  tt(ig)=floor(1.0*abs(xg(ig,1)-xs(1,1))/1500.0/dt)
enddo


if ((T_p+maxval(tt)) < nt ) then
do ig=1,ng
  do it=1,T_p
  taper(it+tt(ig),ig) = 1 - tp1(it) 
  enddo
enddo
endif

dpre = dpre*taper
dobs = dobs*taper

end subroutine dir_taper_remove

subroutine Matrix_hil_1d(xr,nt)
real,intent(inout)           ::  xr(:) 
integer,intent(in)           ::  nt
integer                      ::  kp,nf,k,it
integer(kind=8)              ::  plan, plan1
real, allocatable            ::  h(:)
complex, allocatable  ::  xr_pad(:),Cx(:),Ctmp1(:),trc_work(:)

call nextpow2(nt,kp)
nf = 4*(2**kp);
!nf = 2*nt;

allocate(xr_pad(nf))
xr_pad=0.0
allocate(h(nf))
h=0.0
allocate(Cx(nf))
Cx=0.0
allocate(Ctmp1(nf))
Ctmp1=0.0
allocate(trc_work(nf))
trc_work=0.0

do k=1,nf
  if (k<=nt) then
      xr_pad(k) = cmplx(xr(k),0.0)
  endif
  if (k>nt) then
      xr_pad(k)= cmplx(0.0,0.0);
  endif
enddo

call sfftw_plan_dft_1d(plan,nf,xr_pad,Cx,FFTW_FORWARD,FFTW_ESTIMATE)
call sfftw_execute(plan,xr_pad,Cx)
call sfftw_destroy_plan(plan)

if (mod(nf,2)==0) then
  h(1:nf/2+1) =1
  h(2:nf/2) = 2
else
  h(1)=1
  h(2:(nf+1)/2)=2
endif

Ctmp1 = h*Cx;
call sfftw_plan_dft_1d(plan1,nf,Ctmp1,trc_work,FFTW_BACKWARD,FFTW_ESTIMATE)
call sfftw_execute(plan1,Ctmp1,trc_work)
call sfftw_destroy_plan(plan1)
xr = aimag(trc_work(1:nt))/nf

end subroutine Matrix_hil_1d

subroutine Matrix_hil_2d(d,nz,nx)
real,intent(inout) ::  d(:,:)
integer,intent(in) ::  nz,nx
integer            ::  ix
real               ::  tmp(nx)

do ix=1,nx
   tmp = d(:,ix);
   call Matrix_hil(tmp,nz) 
   d(:,ix) = tmp
enddo

end subroutine Matrix_hil_2d

subroutine d_integer1d_1(d,nt,dout)
integer,intent(in) ::  nt
!real   :: d(:)
real,allocatable   :: d(:), dout(:),tmp(:)
integer :: it
real    :: mean_val

call allocate_and_initial(tmp,nt)

mean_val=sum(abs(d))/nt

if (mean_val>0 )then
 do it=1,nt
  if (it.eq.1) then
     dout(it)=d(it)
  end if
  if (it>1) then
     dout(it) = dout(it-1)+d(it)
  end if
 enddo
tmp = dout*(1.0/(maxval(abs(dout))+1e-36))
dout = tmp
endif

if (mean_val.eq.0) then
  dout=d
endif

end subroutine d_integer1d_1

subroutine d_integer1d_N(d,nt,dout,N)
integer,intent(in) ::  nt
!real   :: d(:)
real,allocatable   :: d(:), dout(:),tmp(:),tmp1(:)
integer :: N,it

call allocate_and_initial(tmp,nt)
call allocate_and_initial(tmp1,nt)

if (N<=0) then
  dout = d;
endif   
if (N>0) then
   do it=1,N
      call d_integer1d_1(d,nt,tmp);
      d=tmp;
   enddo
    dout = tmp;
endif

tmp1 = dout*(1.0/(maxval(abs(dout))+1e-36))
dout = tmp1

end subroutine d_integer1d_N

subroutine d_integer2d_N(d,nt,ng,dout,N)
integer,intent(in) ::  nt,ng
!real   :: d(:,:)
real,allocatable   ::  d(:,:),dout(:,:),trc(:),tmp(:),tmp2(:,:)
integer :: N,it,ig

call allocate_and_initial(trc,nt)
call allocate_and_initial(tmp,nt)
call allocate_and_initial(tmp2,nt,ng)

if (N<=0) then
  dout = d;
endif
if (N>0) then
   do ig=1,ng
      trc = d(:,ig)
      call d_integer1d_N(trc,nt,tmp,N)
      dout(:,ig) = tmp;
   enddo
endif

tmp2 = dout*(1.0/(maxval(abs(dout))+1e-36))
dout = tmp2

end subroutine d_integer2d_N


!subroutine wiener_filter1(F,nt,dt,f_raw,f_target)
!real,allocatable   :: s1(:),s2(:)! s1, raw source wavelet, s2, target  source wavelet
!complex,allocatable :: D(:),otmp(:),trace_work(:)
!integer :: k,it, kp,nf
!real :: pi,dt,f_raw,f_target
!integer :: tmp,nt, i1,i2,i3,i4
!integer(kind=8)   :: plan, plan1
!complex, allocatable :: tracep(:),s1p(:),s2p(:),F(:),Fs1(:),Fs2(:)
!
! pi= 3.1415926535897932384626
!
! call nextpow2(nt,kp)
! nf = 4*(2**kp);
! 
! allocate(tracep(nf))
! allocate(F(nf))
! allocate(Fs1(nf))
! allocate(Fs2(nf))
! 
!
! call ricker_fc(s1,nt,dt,f_raw)
! call ricker_fc(s2,nt,dt,f_target)
!
!! call allocate_and_initial(tracep,nf)
! call allocate_and_initial(s1p,nf)
! call allocate_and_initial(s2p,nf)
! 
!  do k=1,nf
!    if (k<=nt) then
!        !tracep(k) = cmplx(trace(k),0.0)
!        s1p(k) = cmplx(s1(k),0.0)
!        s2p(k) = cmplx(s2(k),0.0)
!    endif
!    if (k>nt) then
!        !tracep(k)= cmplx(0.0,0.0);
!        s1p(k)= cmplx(0.0,0.0);
!        s2p(k)= cmplx(0.0,0.0);
!     endif
!  enddo
!
!! allocate(D(nf))
!! allocate(trace_work(nf))
! call sfftw_plan_dft_1d(plan,nf,s1p,Fs1,FFTW_FORWARD,FFTW_ESTIMATE)
! call sfftw_execute(plan,s1p,Fs1)
! call sfftw_destroy_plan(plan) 
! 
! call sfftw_plan_dft_1d(plan,nf,s2p,Fs2,FFTW_FORWARD,FFTW_ESTIMATE)
! call sfftw_execute(plan,s2p,Fs2)
! call sfftw_destroy_plan(plan) 
!
! F=(Fs1*conjg(Fs2))/(abs(Fs1)**2+0.001*maxval(abs(Fs1))) ! F is the wiener filter
!
!end subroutine wiener_filter1 


subroutine bp_filter1(trace,o,F,nt,dt,f1,f2,f3,f4)
integer,intent(in) ::  nt
real   :: trace(:)
real,intent(out)   :: o(:)
complex,allocatable :: D(:),otmp(:),trace_work(:)
real,allocatable :: x_up(:),x_down(:),up(:),down(:),aux(:),tmp_aux(:),aux2(:),F(:)
integer :: k,it, kp,nf
real :: pi,dt,f1,f2,f3,f4
integer :: tmp, i1,i2,i3,i4
integer(kind=8)   :: plan, plan1
complex, allocatable :: tracep(:)

pi= 3.1415926535897932384626

 call nextpow2(nt,kp)
 nf = 4*(2**kp);
 
 allocate(tracep(nf))
 !call allocate_and_initial(tracep,nf)
  do k=1,nf
    if (k<=nt) then
        tracep(k) = cmplx(trace(k),0.0)
    endif
    if (k>nt) then
        tracep(k)= cmplx(0.0,0.0);
     endif
  enddo

 do k=1,nt
   o(k)=0;
 enddo

 allocate(D(nf))
 allocate(trace_work(nf))
 D=0.0
 trace_work=0.0 

 i1 = floor(nf*f1*dt)+1;
 i2 = floor(nf*f2*dt)+1;
 i3 = floor(nf*f3*dt)+1;
 i4 = floor(nf*f4*dt)+1;

  tmp = i2-i1
  allocate(x_up(tmp))
  do k = 1,(i2-i1),1
   x_up(k) = k
  enddo
  
  tmp=size(x_up,1)
  allocate(up(tmp))
  do k = 1,tmp     
    up(k) = sin(pi*(k-1)/(2*x_up(tmp))) ;     
  enddo

 ! write(*,*) up
  
  tmp = i4-i3
  allocate(x_down(tmp))
  do k = (i4-i3),1,-1
   x_down(k) = i4-i3-k+1
  enddo
!  write(*,*) x_down

  tmp = size(x_down,1)
  allocate(down(tmp))
  do k = 1,tmp;     
    down(k) = cos(pi*(k)/(2*x_down(1))) ;     
  enddo

 ! write(*,*) down
  
  tmp = nf/2+1
  allocate(aux(tmp))
  do k=1,tmp
    if (k<=i1) then
     aux(k)=0
    end if

    if (k>i1.and.k<=i2) then
     aux(k)=up(k-i1)
    end if

    if (k>i2.and.k<=i3) then
     aux(k)=1
    end if

    if (k>i3.and.k<=i4) then
     aux(k)=down(k-i3)
    end if
  
    if (k>i4) then
     aux(k)=0
    end if
  enddo 

  tmp = nf/2-1
  allocate(tmp_aux(tmp))
  do k = 1,tmp
   tmp_aux(k) = aux(k+1)
  enddo 

  allocate(aux2(tmp))
  do k = 1,tmp
   aux2(k) = tmp_aux(tmp-k+1)
  enddo

  !allocate(F(nf))
  do k=1,nf
     if (k <= (nf/2+1)) then 
       F(k) = aux(k)
     endif
     if (k > (nf/2+1)) then
       F(k) = aux2(k-(nf/2+1))
     endif
  enddo

 call sfftw_plan_dft_1d(plan,nf,tracep,D,FFTW_FORWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan,tracep,D)
 call sfftw_destroy_plan(plan) 

 allocate(otmp(nf))
 do k=1,nf
   otmp(k) = D(k)*F(k)
 enddo

 call sfftw_plan_dft_1d(plan1,nf,otmp,trace_work,FFTW_BACKWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan1,otmp,trace_work)
 call sfftw_destroy_plan(plan1) 
 o = real(trace_work(1:nt))/nf

end subroutine bp_filter1 

subroutine nextpow2_1d(ndata,ip)
integer :: ndata
integer,intent(out) :: ip

  if  (ndata > 8192) then
     ip=14
  end if
  if (ndata<8192.and.ndata>4095) then
     ip=13
  end if
  if (ndata<4096.and.ndata>2047) then
     ip=12
  end if
  if (ndata<2048.and.ndata>1023) then
     ip=11
  end if
  if (ndata<1024.and.ndata>511) then
     ip=10
  end if
  if (ndata<512.and.ndata>255) then
     ip=9
  end if
  if (ndata<256.and.ndata>127) then
    ip=8
  end if
  if (ndata<128.and.ndata>63) then
    ip=7
  end if
  if (ndata<64.and.ndata>31) then
    ip=6
  end if
  if  (ndata<32) then
    print *,'Samples Too Small, Check Time Samples...'
    stop
  end if
end subroutine  nextpow2_1d

subroutine dotp_1d(trace1,trace2,nt,norm)
real,intent(in) :: trace1(:),trace2(:)
real,intent(out) :: norm
integer :: it, nt
real :: tmp

tmp = 0.0
do it=1,nt
  tmp = tmp + trace1(it)*trace2(it)
enddo
norm = tmp
end subroutine dotp_1d

subroutine norm_2(trace,nt,norm)
implicit none

real,intent(in) :: trace(:)
real,intent(out) :: norm
integer :: it, nt
real :: tmp

tmp = 0.0
do it=1,nt
  tmp = tmp + trace(it)*trace(it)
enddo

if (tmp.eq.0) then
  norm=0.0
else
  norm = sqrt(tmp)
endif

end subroutine norm_2

subroutine centroid_1d(trace,nt,loc)
real,intent(in) :: trace(:)
integer,intent(out) :: loc
integer :: it, nt
real :: area
real,allocatable :: area_tmp(:)

call allocate_and_initial(area_tmp,nt)

area=0.0
area_tmp=0.0

 do it=1,nt
    area=area+trace(it)
 enddo

 do it=2,nt
    area_tmp(it)=area_tmp(it-1)+trace(it-1)
    if ((area_tmp(it)-0.5*area).LT.0.01) then
     loc=it
    ! exit
    endif
 enddo 

end subroutine centroid_1d


subroutine mx_peak_1d(trace,nt,loc)
real,intent(in) :: trace(:)
integer,intent(out) :: loc
integer :: it, nt
 

  do it=1,nt-1
      if ((trace(it+1).LT.trace(it)).and.(trace(it)>0.8)) then
        loc=it
        exit
      endif
   enddo
end subroutine mx_peak_1d

subroutine xcorr_1d_cor(trace1,trace2,nt,cor)
real,intent(in) :: trace1(:),trace2(:)
real,allocatable :: cor(:)
real :: num,maxnum
integer :: id,it,nt1,nt2, nt,lag,dlag,loc(1)

call allocate_and_initial(cor,2*nt-1)

do id=1,2*nt-1
   lag=id-nt
   num=0.0
   if (lag.eq.0) then
   do it=1,nt
    num=num+trace1(it)*trace2(it)
   enddo
   endif

   if (lag.GT.0) then
   nt1=nt-lag
   do it=1,nt1
   num=num+trace1(lag+it)*trace2(it)
   enddo
   endif

   if (lag.LT.0) then
   nt2=nt+lag
   do it=1,nt2
   num=num+trace1(it)*trace2(it-lag)
   enddo
   endif
   cor(id)=num
enddo
   maxnum=maxval(cor)
   loc= maxloc(cor)
   dlag = loc(1)-nt
   cor = cor/(maxnum+1e-12);

end subroutine xcorr_1d_cor


subroutine xcorr_1dlag(trace1,trace2,nt,dlag)
real,intent(in) :: trace1(:),trace2(:)
real,allocatable :: cor(:)
real :: num,maxnum
integer :: id,it,nt1,nt2, nt,lag,dlag,loc(1)

call allocate_and_initial(cor,2*nt-1)

do id=1,2*nt-1
   lag=id-nt
   num=0.0
   if (lag.eq.0) then
   do it=1,nt
    num=num+trace1(it)*trace2(it)
   enddo
   endif

   if (lag.GT.0) then
   nt1=nt-lag
   do it=1,nt1
   num=num+trace1(lag+it)*trace2(it)
   enddo
   endif

   if (lag.LT.0) then
   nt2=nt+lag
   do it=1,nt2
   num=num+trace1(it)*trace2(it-lag)
   enddo
   endif
   cor(id)=num
enddo
   maxnum=maxval(cor)
   loc= maxloc(cor)
   dlag = loc(1)-nt

end subroutine xcorr_1dlag

subroutine xcorr_1d(trace1,trace2,nt,lag,num)
real,intent(in) :: trace1(:),trace2(:)
real :: num
integer ::it,nt1,nt2, nt,lag
num=0.0
if (lag.eq.0) then
   do it=1,nt
    num=num+trace1(it)*trace2(it)
   enddo
endif

if (lag.GT.0) then
   nt1=nt-lag
   do it=1,nt1
   num=num+trace1(lag+it)*trace2(it)
   enddo
endif

if (lag.LT.0) then
   nt2=nt+lag
   do it=1,nt2
   num=num+trace1(it)*trace2(it-lag)
   enddo
endif
end subroutine xcorr_1d

subroutine tau_cal_1d(lagmax,nlag,dt,tau)
integer,intent(in) :: lagmax,nlag
real, intent(in)   :: dt
real, intent(out)  :: tau(:)
integer :: ilag,lgx

lgx=2*lagmax/(nlag-1)
do ilag=1,nlag
   tau(ilag)=(-lagmax+(ilag-1)*lgx)*dt
enddo
end subroutine tau_cal_1d

subroutine res_ewi_trand_ms(dpre,dobs,nt,ng,res_data,res,ewi_type,dt,f11,f22,f33,f44)
real  :: dpre(:,:),dobs(:,:),dt,f11,f22,f33,f44
integer,intent(in) :: nt,ng,ewi_type ! 1, 2,3
real, intent(out) :: res_data(:,:),res(:)
real,allocatable  :: u(:),u0(:),u_h(:),u0_h(:),u_E(:),u0_E(:),tmp(:),tmp_h(:),u_tmp(:),F(:)
integer :: it,ig,nn,kp
real :: norm_pre,norm_obs,A_max

call allocate_and_initial(u,nt)
call allocate_and_initial(u0,nt)
call allocate_and_initial(u_h,nt)
call allocate_and_initial(u0_h,nt)
call allocate_and_initial(u_E,nt)
call allocate_and_initial(u0_E,nt)
call allocate_and_initial(tmp,nt)
call allocate_and_initial(tmp_h,nt)
call allocate_and_initial(u_tmp,nt)
call norm_shot(dpre,nt,ng)
call norm_shot(dobs,nt,ng)


A_max=maxval(abs(dobs))
call nextpow2(nt,kp)
nn = 4*(2**kp);
allocate(F(nn))
F=0.0

do ig=1,ng
  u_tmp=dpre(:,ig)
  call bp_filter(u_tmp,u,F,nt,dt,f11,f22,f33,f44)
  call norm_n(u,nt,norm_pre)
  u=u/(norm_pre+1e-12)   ! normalized
  u_h =  u
  call Matrix_hil(u_h,nt)
  u_E =  sqrt(u**2+u_h**2)
 
  u_tmp=dobs(:,ig)
  call bp_filter(u_tmp,u0,F,nt,dt,f11,f22,f33,f44)
  call norm_n(u0,nt,norm_obs)
  u0=u0/(norm_obs+1e-12) ! normalized
  u0_h =  u0
  call Matrix_hil(u0_h,nt)
  u0_E =  sqrt(u0**2+u0_h**2)
  if (ewi_type.eq.0) then
   res_data(:,ig)= u - u0;
   res(ig) = sum(res_data(:,ig)**2)
  endif  
  if (ewi_type.eq.1)   then
   tmp = (u_E-u0_E)/(u_E+0.002*maxval(abs(u_E)))
   tmp_h = tmp*u_h
   call Matrix_hil(tmp_h,nt)
   res_data(:,ig)=u*tmp - tmp_h
   res(ig) = sum( (u_E - u0_E)**2)
  endif 
  if (ewi_type.eq.2) then
   tmp = (u_E**2-u0_E**2)
   tmp_h = tmp*u_h
   call Matrix_hil(tmp_h,nt)
   res_data(:,ig)=u*tmp - tmp_h
   res(ig) = sum( (u_E**2 - u0_E**2)**2)
   endif
  if (ewi_type.eq.3) then
   tmp =log(0.001*A_max+u0_E/(u_E+0.002*(maxval(abs(u_E)))))/(u_E**2+0.002*maxval(u_E**2))
   tmp_h = tmp*u_h
   call Matrix_hil(tmp_h,nt)
   res_data(:,ig)= -u*tmp + tmp_h
   res(ig) = sum( (log(0.001*A_max+abs(u_E/(u0_E+0.002*maxval(abs(u0_E))))))**2)
  endif
enddo
end subroutine res_ewi_trand_ms

subroutine res_ewi_trand(dpre,dobs,nt,ng,res_data,res,ewi_type)
real  :: dpre(:,:),dobs(:,:)
integer,intent(in) :: nt,ng,ewi_type ! 1, 2,3
real, intent(out) :: res_data(:,:),res(:)
real,allocatable  :: u(:),u0(:),u_h(:),u0_h(:),u_E(:),u0_E(:),tmp(:),tmp_h(:)
integer :: it,ig
real :: norm_pre,norm_obs,A_max

call allocate_and_initial(u,nt)
call allocate_and_initial(u0,nt)
call allocate_and_initial(u_h,nt)
call allocate_and_initial(u0_h,nt)
call allocate_and_initial(u_E,nt)
call allocate_and_initial(u0_E,nt)
call allocate_and_initial(tmp,nt)
call allocate_and_initial(tmp_h,nt)
call norm_shot(dpre,nt,ng)
call norm_shot(dobs,nt,ng)

A_max=maxval(abs(dobs))

do ig=1,ng
  u=dpre(:,ig)
  call norm_n(u,nt,norm_pre)
  u=u/(norm_pre+1e-12)   ! normalized
  u_h =  u
  call Matrix_hil(u_h,nt)
  u_E =  sqrt(u**2+u_h**2)
 
  u0=dobs(:,ig)           
  call norm_n(u0,nt,norm_obs)
  u0=u0/(norm_obs+1e-12) ! normalized
  u0_h =  u0
  call Matrix_hil(u0_h,nt)
  u0_E =  sqrt(u0**2+u0_h**2)
  if (ewi_type.eq.0) then
   res_data(:,ig)= u - u0;
   res(ig) = sum(res_data(:,ig)**2)
  endif  
  if (ewi_type.eq.1)   then
   tmp = (u_E-u0_E)/(u_E+0.002*maxval(abs(u_E)))
   tmp_h = tmp*u_h
   call Matrix_hil(tmp_h,nt)
   res_data(:,ig)=u*tmp - tmp_h
   res(ig) = sum( (u_E - u0_E)**2)
  endif 
  if (ewi_type.eq.2) then
   tmp = (u_E**2-u0_E**2)
   tmp_h = tmp*u_h
   call Matrix_hil(tmp_h,nt)
   res_data(:,ig)=u*tmp - tmp_h
   res(ig) = sum( (u_E**2 - u0_E**2)**2)
   endif
  if (ewi_type.eq.3) then
   tmp =log(0.01*A_max+u0_E/(u_E+0.002*(maxval(abs(u_E)))))/(u_E**2+0.002*maxval(u_E**2))
   tmp_h = tmp*u_h
   call Matrix_hil(tmp_h,nt)
   res_data(:,ig)= -u*tmp + tmp_h
   res(ig) = sum( (log(0.01*A_max+abs(u_E/(u0_E+0.002*maxval(abs(u0_E))))))**2)
  endif
enddo
end subroutine res_ewi_trand

subroutine res_fwi_trand(dpre,dobs,nt,ng,res_data,res)
real  :: dpre(:,:),dobs(:,:)
integer,intent(in) :: nt,ng
real, intent(out) :: res_data(:,:),res(:)
real,allocatable  :: trace_pre(:),trace_obs(:)
integer :: it,ig
real :: norm_pre,norm_obs

call allocate_and_initial(trace_pre,nt)
call allocate_and_initial(trace_obs,nt)
call norm_shot(dpre,nt,ng)
call norm_shot(dobs,nt,ng)

do ig=1,ng
  trace_pre=dpre(:,ig)
 ! call norm_n(trace_pre,nt,norm_pre)
 ! trace_pre=trace_pre/(norm_pre+1e-12)   ! normalized
  
  trace_obs=dobs(:,ig)           
 ! call norm_n(trace_obs,nt,norm_obs)
 ! trace_obs=trace_obs/(norm_obs+1e-12) ! normalized
  
  res_data(:,ig)= trace_obs - trace_pre;
  res(ig) = sum(res_data(:,ig)**2)
enddo
end subroutine res_fwi_trand

subroutine res_fwi_trand_ms(dpre,dobs,nt,ng,res_data,res,dt,f11,f22,f33,f44)
real  :: dpre(:,:),dobs(:,:),dt,f11,f22,f33,f44
integer,intent(in) :: nt,ng
real, intent(out) :: res_data(:,:),res(:)
real,allocatable  :: trace_pre(:),trace_obs(:),trace_pre_bp(:),trace_obs_bp(:),F(:)
integer :: it,ig,nn,kp
real :: norm_pre,norm_obs

call allocate_and_initial(trace_pre,nt)
call allocate_and_initial(trace_pre_bp,nt)
call allocate_and_initial(trace_obs,nt)
call allocate_and_initial(trace_obs_bp,nt)
call norm_shot(dpre,nt,ng)
call norm_shot(dobs,nt,ng)

call nextpow2(nt,kp)
nn = 4*(2**kp);
allocate(F(nn))
F=0.0

do ig=1,ng
   trace_pre=dpre(:,ig)
 ! call norm_n(trace_pre,nt,norm_pre)
 ! trace_pre=trace_pre/(norm_pre+1e-32)   ! normalized
  call bp_filter(trace_pre,trace_pre_bp,F,nt,dt,f11,f22,f33,f44)
  
  trace_obs=dobs(:,ig)           
 ! call norm_n(trace_obs,nt,norm_obs)
 ! trace_obs=trace_obs/(norm_obs+1e-32) ! normalized
  call bp_filter(trace_obs,trace_obs_bp,F,nt,dt,f11,f22,f33,f44)
  
  res_data(:,ig)= trace_obs_bp - trace_pre_bp;
  res(ig) = sum(res_data(:,ig)**2)
enddo
end subroutine res_fwi_trand_ms


subroutine res_fwi_wcc(dpre,dobs,nt,ng,dt,res_data,res,win,Tmax,lambda)
real  :: dpre(:,:),dobs(:,:)
integer,intent(in) :: nt,ng,win
real, intent(in) :: dt,Tmax,lambda
!real, intent(out) :: res_data(:,:),res(:)
real, allocatable :: res_data(:,:),res(:)
real,allocatable  :: trace_pre(:),trace_pre_A(:),trace_obs(:),trace_obs_A(:),trace(:),cor(:),w1(:)
integer :: t_pre,t_obs,n_shift,ilag,it,ig,shift_max,nlag
real :: dotnum,norm_pre,norm_obs
real, allocatable :: res_data_lag(:,:), res_lag(:)

call allocate_and_initial(trace,nt)
call allocate_and_initial(trace_pre,nt)
call allocate_and_initial(trace_obs,nt)
call allocate_and_initial(trace_pre_A,nt)
call allocate_and_initial(trace_obs_A,nt)

call allocate_and_initial(res_data,nt,ng)
call allocate_and_initial(res,ng)

call csg_mute(dpre,nt,ng,win)
call csg_mute(dobs,nt,ng,win)

shift_max=floor(Tmax/dt)
nlag=2*shift_max-1

call allocate_and_initial(cor,nlag)
call allocate_and_initial(w1,nlag)
call allocate_and_initial(res_data_lag,nt,nlag)
call allocate_and_initial(res_lag,nlag)

do ilag=1,nlag
  w1(ilag)=exp(-(ilag-shift_max-1)**2/(2*(lambda*nlag)**2))
enddo

!if (rank.eq.0) write(*,*) w1

do ig=1,ng
  trace_pre=dpre(:,ig)
  call norm_n(trace_pre,nt,norm_pre)
  trace_pre=trace_pre/(norm_pre+1e-12)   ! normalized
  trace_obs=dobs(:,ig)
  call norm_n(trace_obs,nt,norm_obs)
  trace_obs=trace_obs/(norm_obs+1e-12) ! normalized
  trace_pre_A = trace_pre;
  trace_obs_A = trace_obs;
  call hilbfft(trace_pre_A,nt,1)
  call hilbfft(trace_obs_A,nt,1)
  call xcorr(trace_pre_A,trace_obs_A,nt,cor)
 !if (rank.eq.0) write(*,*), cor

 do ilag=1,nlag
    n_shift=ilag-shift_max-1
    call data_shift(trace_pre,nt,n_shift)
    call dotp(trace_pre,trace_obs,nt,dotnum)
    call norm_n(dpre(:,ig),nt,norm_pre)
    trace = (trace_pre*dotnum-trace_obs)/(norm_pre+1e-12)
    res_data_lag(:,ilag)=trace(:)*w1(ilag)*cor(ilag)
    res_lag(ilag)=w1(ilag)*cor(ilag)
  enddo

  res_data(:,ig)=sum(res_data_lag,2)
  res(ig)=sum(res_lag)
enddo
end subroutine res_fwi_wcc
                                                       
subroutine data_res_fti_new(dpre,dobs,nt,ng,win,dt,Tmax,res_data,res_dt)
! do summation for the data residual and then back propagation
! Tmax : max samples for tau, Tmax = tau_max/dt;
real  :: dpre(:,:),dobs(:,:)
integer,intent(in) :: nt,ng,win,Tmax
real, intent(in) :: dt
real, intent(out) :: res_data(:,:),res_dt(:)
real,allocatable  :: trace(:)
real :: corr,num
real,allocatable :: res_data_tmp(:,:),res_dt_tmp(:)
integer :: n_shift,ig,ilag

call allocate_and_initial(trace,nt)
call allocate_and_initial(res_data_tmp,nt,2*Tmax+1)
call allocate_and_initial(res_dt_tmp,2*Tmax+1)

call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)
call csg_mute(dpre,nt,ng,win) 
call csg_mute(dobs,nt,ng,win) 
call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)

!subroutine xcorr_1d(trace1,trace2,nt,lag,num)
do ig=1,ng
  do ilag=1,2*Tmax+1
    n_shift=ilag-Tmax-1
      call xcorr(dpre(:,ig),dobs(:,ig),nt,n_shift,corr)
      trace=dobs(:,ig)
      call data_shift(trace,nt,n_shift)
      call data_der(trace,nt,dt)
      call trace_dis(trace,nt,num)
      !res_data(:,ig)=-delay*corr**2*trace(:)/(num+1e-9)
      res_data_tmp(:,ilag)=n_shift*dt*corr**2*trace(:)/(num+1e-9)
      !res_data(:,ig)=delay*corr**2*trace(:)/(num)
      res_dt_tmp(ilag)=n_shift*dt*corr
  enddo  
  res_data(:,ig) = sum(res_data_tmp,2);     
  res_dt(ig)=sum(res_dt_tmp(:))

enddo
end subroutine data_res_fti_new

subroutine data_res_fti(dpre,dobs,nt,ng,win,dt,delay,res_data,res_dt)
real  :: dpre(:,:),dobs(:,:)
integer,intent(in) :: nt,ng,win
real, intent(in) :: dt,delay
real, intent(out) :: res_data(:,:),res_dt(:)
real,allocatable  :: trace(:)
real :: corr,num
integer :: n_shift,ig
call allocate_and_initial(trace,nt)
call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)
call csg_mute(dpre,nt,ng,win) 
call csg_mute(dobs,nt,ng,win) 


n_shift=floor(delay/dt)
!if(rank.eq.0) write(*,*) "n_shift", n_shift

do ig=1,ng
   call xcorr(dpre(:,ig),dobs(:,ig),nt,n_shift,corr)
  ! if((rank.eq.0).and.(ig.eq.1)) write(*,*) "corr", corr
  ! if((rank.eq.0).and.(ig.eq.1)) write(*,*) "n_shift", n_shift
  ! if((rank.eq.0).and.(ig.eq.1).and.(n_shift.GT.0)) write(*,*) "dpre", dpre(:,ig)
  ! if((rank.eq.0).and.(ig.eq.1).and.(n_shift.GT.0)) write(*,*) "dobs", dobs(:,ig)
   trace=dobs(:,ig)
   call data_shift(trace,nt,n_shift)
   call data_der(trace,nt,dt)
   call trace_dis(trace,nt,num)
   res_data(:,ig)=-delay*corr**2*trace(:)/(num+1e-9)
   !res_data(:,ig)=delay*corr**2*trace(:)/(num)
   res_dt(ig)=n_shift*dt*corr  
enddo
end subroutine data_res_fti

subroutine data_res_wt_eik(t_pre,t_obs,dobs,nt,ng,win,dt,res_data,res_dt)
real  :: dobs(:,:)
integer,intent(in) :: nt,ng,win
real, intent(in) :: dt,t_pre(:),t_obs(:)
real, intent(out) :: res_data(:,:),res_dt(:)
real,allocatable  :: trace_obs(:),trace(:)
integer :: n_shift,it,ig

call allocate_and_initial(trace,nt)
call allocate_and_initial(trace_obs,nt)

!call norm_trace(dobs,nt,ng)
call csg_mute(dobs,nt,ng,win,t_obs,dt) 
call norm_trace(dobs,nt,ng)
!call norm_trace(dobs,nt,ng)
! write(*,*)  "size ttt= ", size(t_pre)
! write(*,*)  "size tobs = ", size(t_obs)
! write(*,*)  "ng  = ", ng

do ig=1,ng
   n_shift=floor((t_pre(ig)-t_obs(ig))/dt)
!   if(rank.eq.0) write(*,*) n_shift
   trace=dobs(:,ig)
   call data_shift(trace,nt,n_shift)
   call data_der(trace,nt,dt)
   res_data(:,ig)=trace*n_shift*dt
   !res_data(:,ig)=trace
   res_dt(ig)=n_shift*dt  
enddo

end subroutine data_res_wt_eik


!subroutine data_res_phainv_wiener(dpre,dobs,nt,ng,dt,res_adj,res,N,f1,f2,f3,f4,f_raw)
!integer,intent(in) :: nt,ng
!real, intent(in)   :: dt,f1,f2,f3,f4,f_raw
!real               :: dpre(:,:),dobs(:,:),f_target
!real,allocatable   :: res_adj(:,:)
!real,allocatable   :: trace_pre(:),trace_pre_bar(:),trace_pre_barN(:),trace_obs(:),trace_obs_bar(:),trace_obs_barN(:),dp(:),dpN(:),trace(:),res(:),trace_pre_bp(:)
!integer            :: it,k,ig,N,kp,nf
!complex,allocatable :: C_trc_obs(:),C_trc_pre(:),trcp(:),trace_work(:)
!complex,allocatable :: Ctmp1(:),Ctmp2(:),Ctmp3(:)
!integer(kind=8)    :: plan, plan1
!real,allocatable   ::  Mag_trc_obs(:),Mag_trc_pre(:),Pha_trc_obs(:),Pha_trc_pre(:)
!complex :: i
!complex, allocatable :: F(:)
!i = (0.0,1.0);
!
!!allocate(res_adj(nt,ng))
!call allocate_and_initial(trace,nt)
!call allocate_and_initial(trace_pre,nt)
!call allocate_and_initial(trace_pre_bar,nt)
!call allocate_and_initial(trace_pre_barN,nt)
!call allocate_and_initial(trace_pre_bp,nt)
!call allocate_and_initial(trace_obs,nt)
!call allocate_and_initial(trace_obs_bar,nt)
!call allocate_and_initial(trace_obs_barN,nt)
!call allocate_and_initial(dp,nt)
!call allocate_and_initial(dpN,nt)
!call allocate_and_initial(res,ng)
!
!call norm_trace(dpre,nt,ng)
!call norm_trace(dobs,nt,ng)
!
!call nextpow2(nt,kp)
!nf = 4*(2**kp);
!
!allocate(trace_work(nf))
!trace_work=0.0
!allocate(Ctmp1(nf))
!Ctmp1=0.0
!allocate(Ctmp2(nf))
!Ctmp2=0.0
!allocate(Ctmp3(nf))
!Ctmp3=0.0
!allocate(C_trc_obs(nf))
!C_trc_obs=0.0
!allocate(C_trc_pre(nf))
!C_trc_pre=0.0
!allocate(Mag_trc_obs(nf))
!Mag_trc_obs=0.0
!allocate(Mag_trc_pre(nf))
!Mag_trc_pre=0.0
!allocate(Pha_trc_obs(nf))
!Pha_trc_obs=0.0
!allocate(Pha_trc_pre(nf))
!Pha_trc_pre=0.0
!allocate(F(nf))
!F=0.0
!
!!allocate(trcp(nf))
!!do k=1,nf
!!  if (k<=nt) then
!!      trcp(k) = cmplx(trace_obs(k),0.0)
!!  endif
!!  if (k>nt) then
!!      trcp(k)= cmplx(0.0,0.0);
!!  endif
!!enddo
!
! allocate(trcp(nf))
!do ig=1,ng
!! get the magnitude and phase specturm of the observed data
! trace_obs = dobs(:,ig)
! !allocate(trcp(nf))
! do k=1,nf
!  if (k<=nt) then
!      trcp(k) = cmplx(trace_obs(k),0.0)
!  endif
!  if (k>nt) then
!      trcp(k)= cmplx(0.0,0.0);
!  endif
! enddo
!
! call sfftw_plan_dft_1d(plan,nf,trcp,C_trc_obs,FFTW_FORWARD,FFTW_ESTIMATE)
! call sfftw_execute(plan,trcp,C_trc_obs)
! call sfftw_destroy_plan(plan)
! Mag_trc_obs = abs(C_trc_obs)
! Pha_trc_obs = atan2(aimag(C_trc_obs),real(C_trc_obs))
!
!! get the magnitude and phase specturm of the predicted data
! trace_pre = dpre(:,ig)
! do k=1,nf
!  if (k<=nt) then
!      trcp(k) = cmplx(trace_pre(k),0.0)
!  endif
!  if (k>nt) then
!      trcp(k)= cmplx(0.0,0.0);
!  endif
! enddo
!
! call sfftw_plan_dft_1d(plan1,nf,trcp,C_trc_pre,FFTW_FORWARD,FFTW_ESTIMATE)
! call sfftw_execute(plan1,trcp,C_trc_pre)
! call sfftw_destroy_plan(plan1)
! Mag_trc_pre = abs(C_trc_pre)
! Pha_trc_pre = atan2(aimag(C_trc_pre),real(C_trc_pre))
!
!! get the F, lowpass filter, ifft to get trace_obs_bar
!! using wiener filter to get F
! f_target = 0.5*(f2+f3)
! call  bp_filter(F,nt,dt,f_raw,f_target)
!
!! call bp_filter(trace_pre,trace_pre_bp,F,nt,dt,f1,f2,f3,f4)
! Ctmp1=F*Mag_trc_obs*exp(i*Pha_trc_obs)
!
! call sfftw_plan_dft_1d(plan1,nf,Ctmp1,trace_work,FFTW_BACKWARD,FFTW_ESTIMATE)
! call sfftw_execute(plan1,Ctmp1,trace_work)
! call sfftw_destroy_plan(plan1)
! trace_obs_bar = real(trace_work(1:nt))/nf
!
!! lowpass filter, ifft to get trace_pre_bar
! Ctmp2=F*Mag_trc_obs*exp(i*Pha_trc_pre)
!
! call sfftw_plan_dft_1d(plan1,nf,Ctmp2,trace_work,FFTW_BACKWARD,FFTW_ESTIMATE)
! call sfftw_execute(plan1,Ctmp2,trace_work)
! call sfftw_destroy_plan(plan1)
! trace_pre_bar = real(trace_work(1:nt))/nf
!
!
! dp = trace_obs_bar - trace_pre_bar;
!
! call d_integerN(dp,nt,dpN,N-1)
! res_adj(:,ig) = dpN 
! 
! call d_integerN(trace_obs_bar,nt,trace_obs_barN,N)
! call d_integerN(trace_pre_bar,nt,trace_pre_barN,N)
! 
! res(ig)= sum((trace_obs_barN - trace_pre_barN)**2) 
! 
!enddo
! !res = sum(res0)
!
!end subroutine data_res_phainv_wiener

subroutine data_res_empi(dpre,dobs,nt,ng,dt,res_adj,res,N,f1,f2,f3,f4,A_err,fc,N_period)
integer,intent(in) :: nt,ng
real, intent(in)   :: dt,f1,f2,f3,f4
real               :: dpre(:,:),dobs(:,:)
real,allocatable   :: res_adj(:,:)
real,allocatable   :: trace_pre(:),trace_pre_bar(:),trace_pre_barN(:),trace_obs(:),trace_obs_bar(:),trace_obs_barN(:),dp(:),dpN(:),trace(:),res(:),trace_pre_bp(:)
integer            :: it,k,ig,N,kp,nf
complex,allocatable :: C_trc_obs(:),C_trc_pre(:),trcp(:),trace_work(:)
complex,allocatable :: Ctmp1(:),Ctmp2(:),Ctmp3(:)
integer(kind=8)    :: plan, plan1
real,allocatable   ::  Mag_trc_obs(:),Mag_trc_pre(:),Pha_trc_obs(:),Pha_trc_pre(:),F(:)
complex :: i
real               :: A_err,fc,N_period   ! used to extract the first peak in the calculated data

i = (0.0,1.0);

!allocate(res_adj(nt,ng))
call allocate_and_initial(trace,nt)
call allocate_and_initial(trace_pre,nt)
call allocate_and_initial(trace_pre_bar,nt)
call allocate_and_initial(trace_pre_barN,nt)
call allocate_and_initial(trace_pre_bp,nt)
call allocate_and_initial(trace_obs,nt)
call allocate_and_initial(trace_obs_bar,nt)
call allocate_and_initial(trace_obs_barN,nt)
call allocate_and_initial(dp,nt)
call allocate_and_initial(dpN,nt)
call allocate_and_initial(res,ng)


call nextpow2(nt,kp)
nf = 4*(2**kp);

allocate(trace_work(nf))
trace_work=0.0
allocate(Ctmp1(nf))
Ctmp1=0.0
allocate(Ctmp2(nf))
Ctmp2=0.0
allocate(Ctmp3(nf))
Ctmp3=0.0
allocate(C_trc_obs(nf))
C_trc_obs=0.0
allocate(C_trc_pre(nf))
C_trc_pre=0.0
allocate(Mag_trc_obs(nf))
Mag_trc_obs=0.0
allocate(Mag_trc_pre(nf))
Mag_trc_pre=0.0
allocate(Pha_trc_obs(nf))
Pha_trc_obs=0.0
allocate(Pha_trc_pre(nf))
Pha_trc_pre=0.0
allocate(F(nf))
F=0.0


call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)

call csg_window(dobs,nt,ng,dt,A_err,fc,N_period)
call csg_window(dpre,nt,ng,dt,A_err,fc,N_period)


!allocate(trcp(nf))
!do k=1,nf
!  if (k<=nt) then
!      trcp(k) = cmplx(trace_obs(k),0.0)
!  endif
!  if (k>nt) then
!      trcp(k)= cmplx(0.0,0.0);
!  endif
!enddo

 allocate(trcp(nf))
do ig=1,ng
! get the magnitude and phase specturm of the observed data
 trace_obs = dobs(:,ig)
 !allocate(trcp(nf))
 do k=1,nf
  if (k<=nt) then
      trcp(k) = cmplx(trace_obs(k),0.0)
  endif
  if (k>nt) then
      trcp(k)= cmplx(0.0,0.0);
  endif
 enddo

 call sfftw_plan_dft_1d(plan,nf,trcp,C_trc_obs,FFTW_FORWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan,trcp,C_trc_obs)
 call sfftw_destroy_plan(plan)
 Mag_trc_obs = abs(C_trc_obs)
 Pha_trc_obs = atan2(aimag(C_trc_obs),real(C_trc_obs))

! get the magnitude and phase specturm of the predicted data
 trace_pre = dpre(:,ig)
 do k=1,nf
  if (k<=nt) then
      trcp(k) = cmplx(trace_pre(k),0.0)
  endif
  if (k>nt) then
      trcp(k)= cmplx(0.0,0.0);
  endif
 enddo

 call sfftw_plan_dft_1d(plan1,nf,trcp,C_trc_pre,FFTW_FORWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan1,trcp,C_trc_pre)
 call sfftw_destroy_plan(plan1)
 Mag_trc_pre = abs(C_trc_pre)
 Pha_trc_pre = atan2(aimag(C_trc_pre),real(C_trc_pre))

! get the F, lowpass filter, ifft to get trace_obs_bar
 call bp_filter(trace_pre,trace_pre_bp,F,nt,dt,f1,f2,f3,f4)
 Ctmp1=F*Mag_trc_obs*exp(i*Pha_trc_obs)

 call sfftw_plan_dft_1d(plan1,nf,Ctmp1,trace_work,FFTW_BACKWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan1,Ctmp1,trace_work)
 call sfftw_destroy_plan(plan1)
 trace_obs_bar = real(trace_work(1:nt))/nf

! lowpass filter, ifft to get trace_pre_bar
 Ctmp2=F*Mag_trc_obs*exp(i*Pha_trc_pre)

 call sfftw_plan_dft_1d(plan1,nf,Ctmp2,trace_work,FFTW_BACKWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan1,Ctmp2,trace_work)
 call sfftw_destroy_plan(plan1)
 trace_pre_bar = real(trace_work(1:nt))/nf


 dp = trace_obs_bar - trace_pre_bar;

 call d_integerN(dp,nt,dpN,N-1)
 res_adj(:,ig) = dpN 
 
 call d_integerN(trace_obs_bar,nt,trace_obs_barN,N)
 call d_integerN(trace_pre_bar,nt,trace_pre_barN,N)
 
 res(ig)= sum((trace_obs_barN - trace_pre_barN)**2) 
 
enddo
 !res = sum(res0)
 !call deallocate_and_free(dpre)
 !call deallocate_and_free(dobs)

end subroutine data_res_empi


subroutine data_res_phainv(dpre,dobs,nt,ng,dt,res_adj,res,N,f1,f2,f3,f4)
integer,intent(in) :: nt,ng
real, intent(in)   :: dt,f1,f2,f3,f4
real               :: dpre(:,:),dobs(:,:)
real,allocatable   :: res_adj(:,:)
real,allocatable   :: trace_pre(:),trace_pre_bar(:),trace_pre_barN(:),trace_obs(:),trace_obs_bar(:),trace_obs_barN(:),dp(:),dpN(:),trace(:),res(:),trace_pre_bp(:)
integer            :: it,k,ig,N,kp,nf
complex,allocatable :: C_trc_obs(:),C_trc_pre(:),trcp(:),trace_work(:)
complex,allocatable :: Ctmp1(:),Ctmp2(:),Ctmp3(:)
integer(kind=8)    :: plan, plan1
real,allocatable   ::  Mag_trc_obs(:),Mag_trc_pre(:),Pha_trc_obs(:),Pha_trc_pre(:),F(:)
complex :: i

i = (0.0,1.0);

!allocate(res_adj(nt,ng))
call allocate_and_initial(trace,nt)
call allocate_and_initial(trace_pre,nt)
call allocate_and_initial(trace_pre_bar,nt)
call allocate_and_initial(trace_pre_barN,nt)
call allocate_and_initial(trace_pre_bp,nt)
call allocate_and_initial(trace_obs,nt)
call allocate_and_initial(trace_obs_bar,nt)
call allocate_and_initial(trace_obs_barN,nt)
call allocate_and_initial(dp,nt)
call allocate_and_initial(dpN,nt)
call allocate_and_initial(res,ng)

call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)

call nextpow2(nt,kp)
nf = 4*(2**kp);

allocate(trace_work(nf))
trace_work=0.0
allocate(Ctmp1(nf))
Ctmp1=0.0
allocate(Ctmp2(nf))
Ctmp2=0.0
allocate(Ctmp3(nf))
Ctmp3=0.0
allocate(C_trc_obs(nf))
C_trc_obs=0.0
allocate(C_trc_pre(nf))
C_trc_pre=0.0
allocate(Mag_trc_obs(nf))
Mag_trc_obs=0.0
allocate(Mag_trc_pre(nf))
Mag_trc_pre=0.0
allocate(Pha_trc_obs(nf))
Pha_trc_obs=0.0
allocate(Pha_trc_pre(nf))
Pha_trc_pre=0.0
allocate(F(nf))
F=0.0

!allocate(trcp(nf))
!do k=1,nf
!  if (k<=nt) then
!      trcp(k) = cmplx(trace_obs(k),0.0)
!  endif
!  if (k>nt) then
!      trcp(k)= cmplx(0.0,0.0);
!  endif
!enddo

 allocate(trcp(nf))
do ig=1,ng
! get the magnitude and phase specturm of the observed data
 trace_obs = dobs(:,ig)
 !allocate(trcp(nf))
 do k=1,nf
  if (k<=nt) then
      trcp(k) = cmplx(trace_obs(k),0.0)
  endif
  if (k>nt) then
      trcp(k)= cmplx(0.0,0.0);
  endif
 enddo

 call sfftw_plan_dft_1d(plan,nf,trcp,C_trc_obs,FFTW_FORWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan,trcp,C_trc_obs)
 call sfftw_destroy_plan(plan)
 Mag_trc_obs = abs(C_trc_obs)
 Pha_trc_obs = atan2(aimag(C_trc_obs),real(C_trc_obs))

! get the magnitude and phase specturm of the predicted data
 trace_pre = dpre(:,ig)
 do k=1,nf
  if (k<=nt) then
      trcp(k) = cmplx(trace_pre(k),0.0)
  endif
  if (k>nt) then
      trcp(k)= cmplx(0.0,0.0);
  endif
 enddo

 call sfftw_plan_dft_1d(plan1,nf,trcp,C_trc_pre,FFTW_FORWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan1,trcp,C_trc_pre)
 call sfftw_destroy_plan(plan1)
 Mag_trc_pre = abs(C_trc_pre)
 Pha_trc_pre = atan2(aimag(C_trc_pre),real(C_trc_pre))

! get the F, lowpass filter, ifft to get trace_obs_bar
 call bp_filter(trace_pre,trace_pre_bp,F,nt,dt,f1,f2,f3,f4)
 Ctmp1=F*Mag_trc_obs*exp(i*Pha_trc_obs)

 call sfftw_plan_dft_1d(plan1,nf,Ctmp1,trace_work,FFTW_BACKWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan1,Ctmp1,trace_work)
 call sfftw_destroy_plan(plan1)
 trace_obs_bar = real(trace_work(1:nt))/nf

! lowpass filter, ifft to get trace_pre_bar
 Ctmp2=F*Mag_trc_obs*exp(i*Pha_trc_pre)

 call sfftw_plan_dft_1d(plan1,nf,Ctmp2,trace_work,FFTW_BACKWARD,FFTW_ESTIMATE)
 call sfftw_execute(plan1,Ctmp2,trace_work)
 call sfftw_destroy_plan(plan1)
 trace_pre_bar = real(trace_work(1:nt))/nf


 dp = trace_obs_bar - trace_pre_bar;

 call d_integerN(dp,nt,dpN,N-1)
 res_adj(:,ig) = dpN 
 
 call d_integerN(trace_obs_bar,nt,trace_obs_barN,N)
 call d_integerN(trace_pre_bar,nt,trace_pre_barN,N)
 
 res(ig)= sum((trace_obs_barN - trace_pre_barN)**2) 
 
enddo
 !res = sum(res0)
 !call deallocate_and_free(dpre)
 !call deallocate_and_free(dobs)

end subroutine data_res_phainv


subroutine data_res_wt(dpre,dobs,nt,ng,win,dt,res_data,res_dt,flag)
real  :: dpre(:,:),dobs(:,:)
integer,intent(in) :: nt,ng,win
real, intent(in) :: dt
real, intent(out) :: res_data(:,:),res_dt(:)
real,allocatable  :: trace_pre(:),trace_obs(:),trace(:)
integer :: t_pre,t_obs,n_shift,it,ig,flag
integer :: pt1(1),pt2(1)

call allocate_and_initial(trace,nt)
call allocate_and_initial(trace_pre,nt)
call allocate_and_initial(trace_obs,nt)

call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)
call csg_mute(dpre,nt,ng,win) 
call csg_mute(dobs,nt,ng,win) 
call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)

if (flag.eq.1) then
do ig=1,ng
   do it=1,nt
    trace_pre(it)=dpre(it,ig)
    trace_obs(it)=dobs(it,ig)
   enddo
   !call hilbfft(trace_pre,nt,1)
   !call hilbfft(trace_obs,nt,1)
   dpre(:,ig)=trace_pre
   dobs(:,ig)=trace_obs
enddo

do ig=1,ng
   pt1=maxloc(dpre(:,ig))
   pt2=maxloc(dobs(:,ig))
  ! npre(ig)=t1(1)
  ! nobs(ig)=t2(1)
   t_pre=pt1(1)
   t_obs=pt2(1)
!   call centroid(dpre(:,ig),nt,t_pre)
!   call centroid(dobs(:,ig),nt,t_obs)
   n_shift=t_pre-t_obs
!   if(rank.eq.0) write(*,*) ,npre
   trace=dobs(:,ig)
   call data_shift(trace,nt,n_shift)
   call data_der(trace,nt,dt)
   res_data(:,ig)=trace(:)*n_shift*dt
   res_dt(ig)=n_shift*dt  
enddo
endif

if (flag.eq.2) then
do ig=1,ng
   t_pre=0
   t_obs=0
   do it=1,nt-1    
      !if (dpre(it+1,ig)<dpre(it,ig).and.dpre(it,ig)>0.005) then
      if (abs(dpre(it+1,ig))<abs(dpre(it,ig)).and.abs(dpre(it,ig))>0.05) then
        t_pre=it
        exit
      endif
   enddo
   do it=1,nt-1    
     !if (dobs(it+1,ig)<dobs(it,ig).and.dobs(it,ig)>0.005) then
     if (abs(dobs(it+1,ig))<abs(dobs(it,ig)).and.abs(dobs(it,ig))>0.05) then
        t_obs=it
        exit
      endif
   enddo
   n_shift=t_pre-t_obs
   trace=dobs(:,ig)
   call data_shift(trace,nt,n_shift)
   call data_der(trace,nt,dt)
   res_data(:,ig)=trace(:)*n_shift*dt
   res_dt(ig)=n_shift*dt  
enddo
endif
end subroutine data_res_wt

subroutine data_res_wt_multioffset(dpre,dobs,nt,ng,win,dt,res_data,res_dt,flag,is,xs,xg,offset)
real  :: dpre(:,:),dobs(:,:)
integer,intent(in) :: nt,ng,win,is
real, intent(in) :: dt,xs(:),xg(:),offset
real, intent(out) :: res_data(:,:),res_dt(:)
real,allocatable  :: trace_pre(:),trace_obs(:),trace(:)
integer :: t_pre,t_obs,n_shift,it,ig,flag
integer :: pt1(1),pt2(1)

call allocate_and_initial(trace,nt)

call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)
call csg_mute(dpre,nt,ng,win) 
call csg_mute(dobs,nt,ng,win) 

call window_shot(dpre,nt,ng,is,xs,xg,offset)
call window_shot(dobs,nt,ng,is,xs,xg,offset)

if (flag.eq.1) then
do ig=1,ng
   do it=1,nt
    trace_pre(it)=dpre(it,ig)
    trace_obs(it)=dobs(it,ig)
   enddo
   call hilbfft(trace_pre,nt,1)
   call hilbfft(trace_obs,nt,1)
   dpre(:,ig)=trace_pre
   dobs(:,ig)=trace_obs
enddo

do ig=1,ng
   pt1=maxloc(dpre(:,ig))
   pt2=maxloc(dobs(:,ig))
  ! npre(ig)=t1(1)
  ! nobs(ig)=t2(1)
   t_pre=pt1(1)
   t_obs=pt2(1)
!   call centroid(dpre(:,ig),nt,t_pre)
!   call centroid(dobs(:,ig),nt,t_obs)
   n_shift=t_pre-t_obs
!   if(rank.eq.0) write(*,*) ,npre
   trace=dobs(:,ig)
   call data_shift(trace,nt,n_shift)
   call data_der(trace,nt,dt)
   res_data(:,ig)=-trace(:)*n_shift*dt
   res_dt(ig)=n_shift*dt  
enddo
endif

if (flag.eq.2) then
do ig=1,ng
   t_pre=0
   t_obs=0
   do it=1,nt-1    
      if (dpre(it+1,ig)<dpre(it,ig).and.dpre(it,ig)>0.005) then
        t_pre=it
        exit
      endif
   enddo
   do it=1,nt-1    
     if (dobs(it+1,ig)<dobs(it,ig).and.dobs(it,ig)>0.005) then
        t_obs=it
        exit
      endif
   enddo
   n_shift=t_pre-t_obs
   trace=dobs(:,ig)
   call data_shift(trace,nt,n_shift)
   call data_der(trace,nt,dt)
   res_data(:,ig)=-trace(:)*n_shift*dt
   res_dt(ig)=n_shift*dt  
enddo
endif

end subroutine data_res_wt_multioffset

subroutine data_res_wt_hilb(dpre,dobs,nt,ng,win,dt,res_data,res_dt,temp_out)
real  :: dpre(:,:),dobs(:,:)
integer,intent(in) :: nt,ng,win
real, intent(in) :: dt
real, intent(out) :: res_data(:,:),res_dt(:)
real,allocatable  :: trace_pre(:),trace_obs(:),trace(:),temp_out(:,:)
integer :: t_pre,t_obs,n_shift,it,ig
call allocate_and_initial(trace,nt)
call allocate_and_initial(temp_out,nt,ng)

call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)
call csg_mute(dpre,nt,ng,win) 
call csg_mute(dobs,nt,ng,win) 
call norm_trace(dpre,nt,ng)
call norm_trace(dobs,nt,ng)

allocate(trace_pre(1:nt))
allocate(trace_obs(1:nt))

do ig=1,ng
  do it=1,nt
   trace_pre(it)=dpre(it,ig)
   trace_obs(it)=dobs(it,ig)
  enddo
   call hilbfft(trace_pre,nt,1)
   call hilbfft(trace_obs,nt,1)
   !temp_out(:,ig)=trace_obs
   temp_out(:,ig)=trace_pre
  
   call xcorr(trace_pre,trace_obs,nt,n_shift)
    
  ! call mx_peak(trace_pre,nt,t_pre)
  ! call mx_peak(trace_obs,nt,t_obs)
  ! if(rank.eq.0.) write(*,*) "maxxx" ,  n_shift 
  ! t_pre=maxloc(trace_pre)
  ! t_obs=maxloc(trace_obs)

   ! n_shift=t_pre-t_obs
   !temp_out(:,ig)=n_shift
   
   trace=dobs(:,ig)

   call data_shift(trace,nt,n_shift)
  ! call data_der(trace,nt,dt)
   res_data(:,ig)=trace(:)*n_shift*dt
   res_dt(ig)=n_shift*dt  
enddo
end subroutine data_res_wt_hilb

subroutine data_der_1d(trace,nt,dt)
real ,intent(inout) :: trace(:)
real , intent(in)  :: dt
integer , intent(in) :: nt
integer :: it
do it=1,nt-1
   trace(it)=(trace(it+1)-trace(it))/dt
enddo 
end subroutine data_der_1d

subroutine data_shift_1d(trace,nt,n_shift)
real ,intent(inout) :: trace(:)
integer,intent(in)  :: nt,n_shift
real ,allocatable   ::  temp(:)
integer :: lag

allocate(temp(nt))

if (n_shift.eq.0) then 
   temp=trace
elseif (n_shift.LT.0) then 
   lag=-n_shift
   temp(1:nt-lag+1)=trace(lag:nt)
elseif (n_shift.GT.0) then
   lag=n_shift
   temp(1:lag)=0
   temp(lag+1:nt)=trace(1:nt-lag)
endif
   trace=temp
end subroutine data_shift_1d

subroutine trace_dis_1d(trace,nt,num)
integer,intent(in) :: nt
real,intent(in) :: trace(:)
integer            :: it
real :: num
num=0.0
do it=1,nt
   num=num+trace(it)**2
enddo
end subroutine trace_dis_1d

subroutine window_shot_2d(csg,nt,ng,is,xs,xg,offset)
integer,intent(in) :: nt,ng,is
real ,intent(in)   :: xs(:),xg(:),offset
real,intent(inout) :: csg(:,:)
integer :: ig
do ig=1,ng
   if (abs(xs(is)-xg(ig)).GT.offset) then
    csg(:,ig)=0
   endif
enddo
end subroutine window_shot_2d

subroutine norm_trace_2d(csg,nt,ng)
integer,intent(in) :: nt,ng
real,intent(inout) :: csg(:,:)
integer            :: ig
real               :: eps,mean_val

do ig=1,ng
   mean_val=sum(csg(:,ig))/nt
   if (mean_val>0) then
      csg(:,ig)=csg(:,ig)/(maxval(csg(:,ig))+0.0001*maxval(csg(:,ig)))
   endif
enddo
end subroutine norm_trace_2d

subroutine norm_shot_2d(csg,nt,ng)
integer,intent(in) :: nt,ng
real,intent(inout) :: csg(:,:)
integer            :: ig

do ig=1,ng
   csg(:,ig)=csg(:,ig)/(maxval(csg(:,:))+1e-36)
enddo
end subroutine norm_shot_2d

subroutine csg_window_2d(dpre,dobs,nt,ng,dt,Tmax,n_taper)
implicit none
real  :: dpre(:,:),dobs(:,:)
integer, intent(in) :: nt,ng
real   , intent(in) :: dt,Tmax
integer :: it,ig,n0,n_taper
real, allocatable   :: taper(:) , taper0(:)

n0=int(Tmax/dt)

allocate(taper0(n_taper))
taper0=0.0
allocate(taper(nt))
taper=0.0

do it=1,n_taper
   taper0(it)=0.5*(1 + cos(6.283*(it-1)/(2.0*n_taper)))
enddo

do it=1,nt
   if (it<n0) then
    taper(it)=1.0
  endif

  if( (it>n0).and.(it<=n0+n_taper) )then
    taper(it) = taper0(it-n0)  
  endif
 
  if (it>n0+n_taper) then
   taper(it)=0.0
  endif
enddo


do ig=1,ng
   dpre(:,ig)=dpre(:,ig)*taper 
   dobs(:,ig)=dobs(:,ig)*taper
enddo

deallocate(taper,taper0)

end subroutine csg_window_2d

subroutine csg_win_early_2d(dpre,nt,ng,dt,A_err,fc,N_period)
implicit none
real  :: dpre(:,:)
integer, intent(in) :: nt,ng
real   , intent(in) :: dt
integer             :: it,ig,n0,wd,gap,b1,b2,b3,b4
real, allocatable   :: tau(:),taper(:),mask(:,:)
real                :: A_err,fc,N_period,T0

T0=1.0/(1.0*fc)
wd= int((T0*N_period)/dt)
gap=int(0.2*wd)


allocate(tau(ng))
tau=0.0
allocate(mask(nt,ng))
mask=0.0
allocate(taper(nt))
taper=1.0

! extract the first peak for early arrival , for predicted data
do ig=1,ng
 do it=1,nt-1
    if (abs(dpre(it,ig))>A_err) then
     tau(ig)=it
     exit
    endif
 enddo
enddo

! mask data
do ig =1,ng
  if (tau(ig)-gap>1) then
   b1=tau(ig)-gap
  else
   b1=1
  endif
  b2 = tau(ig)
  b3 = tau(ig)+wd-gap
  if (tau(ig)+wd>nt) then  
   b4 = nt
  else
   b4 = tau(ig)+wd
  endif

 do it=1,nt
  if (it<=b1) then
   taper(it)=0.0
  elseif ((it>b1).and.(it<b2)) then
   taper(it)=0.5*(1-cos(6.283*(it-b1)/(2*gap)))
  elseif ((it>=b2).and.(it<=b3)) then
   taper(it)=1.0  
  elseif ((it>b3).and.(it<b4)) then
   taper(it)=0.5*(1+cos(6.283*(it-b1)/(2*gap)))
  else
   taper(it)=0.0
  endif

 enddo
 mask(:,ig)=taper;
enddo
dpre=dpre*mask

deallocate(tau,taper,mask)

end subroutine csg_win_early_2d

subroutine csg_cut_2d_redsea(dpre,dobs,nt,ng,offset_cut,offset_cut2)
real  :: dpre(:,:),dobs(:,:)
integer, intent(in) :: nt,ng,offset_cut,offset_cut2
integer :: ig
do ig=1,ng
   if (abs(1-ig).LT.offset_cut) then
   dpre(:,ig)=0.0  
   dobs(:,ig)=0.0  
   endif
enddo
do ig=1,ng
   if (abs(ng-ig).LT.offset_cut2) then
   dpre(:,ig)=0.0  
   dobs(:,ig)=0.0  
   endif
enddo
end subroutine csg_cut_2d_redsea


subroutine csg_cut_2d0(dpre,dobs,nt,ng,offset_cut)
real  :: dpre(:,:),dobs(:,:)
integer, intent(in) :: nt,ng,offset_cut
integer :: ig
do ig=1,ng
   if (abs(1-ig).LT.offset_cut) then
   dpre(:,ig)=0.0  
   dobs(:,ig)=0.0  
   endif
enddo
end subroutine csg_cut_2d0

subroutine csg_cut_2d(dpre,dobs,nt,ng,is,offset_cut)
real  :: dpre(:,:),dobs(:,:)
integer, intent(in) :: nt,ng,is,offset_cut
integer :: ig
do ig=1,ng
   if (abs(is-ig).LT.offset_cut) then
   dpre(:,ig)=0.0  
   dobs(:,ig)=0.0  
   endif
enddo
end subroutine csg_cut_2d

subroutine csg_cut_OBC(dpre,dobs,nt,ng,xs,xg,near_cut,far_cut)
real                :: dpre(:,:),dobs(:,:)
integer, intent(in) :: nt,ng
real, intent(in)    :: near_cut,far_cut
integer             :: ig
real                :: xs(:,:),xg(:,:)    
do ig=1,ng
   !if (abs(xs(1,1)-xg(ig,1)).GT.offset_cut) then
   if (abs(xs(1,1)-xg(ig,1)).LT.near_cut) then
   dpre(:,ig)=0.0  
   dobs(:,ig)=0.0  
   endif
   if (abs(xs(1,1)-xg(ig,1)).GT.far_cut) then
   dpre(:,ig)=0.0  
   dobs(:,ig)=0.0  
   endif
enddo
end subroutine csg_cut_OBC

subroutine res_data_cut_near(res_data,nt,ng,xs,xg,it,vmin,freq)
real                :: res_data(:,:)
integer, intent(in) :: nt,ng,it
integer             :: ig
real                :: xs(:,:),xg(:,:),vmin,freq,lamd

lamd = vmin/freq
!write(*,*) "lambda = ",lamd

if (mod(it,2).eq.0) then 
do ig=1,ng
   if (abs(xs(1,1)-xg(ig,1)).GT.lamd) then
   res_data(:,ig)=0.0  
   endif
enddo
endif

end subroutine res_data_cut_near


subroutine csg_mute_tt(csg,nt,ng,win,tobs,dt)
integer, intent(in) :: nt,ng,win
real,intent(in)     :: tobs(:),dt
real, intent(inout) :: csg(:,:)
integer  :: pwin,pmax,it,ig
real  :: b1,b2
real, allocatable   :: csg_temp(:,:),tapering(:)

allocate(csg_temp(nt,ng))
call allocate_and_initial(tapering,nt)

do ig=1,ng
   !csg_temp(:,ig)=csg(:,ig)/maxval(csg(:,ig)+1e-9);
   pmax=tobs(ig)/dt
   !do it=1,nt-1    
   ! ! if (abs(csg_temp(it+1,ig))<abs(csg_temp(it,ig)).and.abs(csg_temp(it,ig))>0.005) then
   !  if (csg_temp(it+1,ig)<csg_temp(it,ig).and.csg_temp(it,ig)>0.005) then
   !     pmax=it
   !     exit
   !  endif
   !enddo
   b1=pmax-0.2*win
   b2=pmax+0.8*win
   !pwin=pmax+win
   call tapering_1d(tapering,b1,b2,nt)
   csg(:,ig)=csg(:,ig)*tapering
!   csg(pwin:nt,ig)=0
enddo
end subroutine csg_mute_tt

subroutine csg_mute_2d(csg,nt,ng,win)
integer, intent(in) :: nt,ng,win
real, intent(inout) :: csg(:,:)
integer  :: pwin,pmax,it,ig
real  :: b1,b2
real, allocatable   :: csg_temp(:,:),tapering(:)

allocate(csg_temp(nt,ng))
call allocate_and_initial(tapering,nt)

do ig=1,ng
   csg_temp(:,ig)=csg(:,ig)/maxval(csg(:,ig)+1e-9);
   pmax=0.0
   do it=1,nt-1    
    ! if (abs(csg_temp(it+1,ig))<abs(csg_temp(it,ig)).and.abs(csg_temp(it,ig))>0.005) then
     if (csg_temp(it+1,ig)<csg_temp(it,ig).and.csg_temp(it,ig)>0.005) then
        pmax=it
        exit
     endif
   enddo
   b1=pmax-0.2*win
   b2=pmax+0.8*win
   pwin=pmax+win
   call tapering_1d(tapering,b1,b2,nt)
   csg(:,ig)=csg(:,ig)*tapering
!   csg(pwin:nt,ig)=0
enddo

end subroutine csg_mute_2d


subroutine tapering_1d(tapering,b1,b2,nt)
integer,intent(in) :: nt
real,intent(in) :: b1,b2
real,allocatable   :: tapering(:)
real :: sigma
integer :: it
sigma =abs( 0.3*(b2-b1));
call allocate_and_initial(tapering,nt)   
do it=1,nt
  if (it.LT.b1) then
     tapering(it)=exp(-((it-b1)**2)/(2*sigma**2))
  elseif ((it.GE.b1).and.(it.LE.b2)) then
     tapering(it)=1
  else
     tapering(it)=exp(-((it-b2)**2)/(2*sigma**2))
  endif
enddo

end subroutine tapering_1d

subroutine dist_cal_2d(dist,nz,nx,dz,dx,zs,xs)
real, intent(in) :: dz,dx,zs,xs
integer, intent(in) :: nz,nx
real, intent(out) :: dist(:,:)
real              :: x,z,dist_x2,dist_z2
integer           :: ix, iz
do ix=1,nx
   x=(ix-1)*dx;   dist_x2=(xs-x)**2.0
   do iz=1,nz
      z=(iz-1)*dz;   dist_z2=(zs-z)**2.0
      dist(iz,ix)=sqrt(dist_x2+dist_z2)
      if (dist(iz,ix) < dx/2.0) then
         dist(iz,ix)=dx
      endif
   enddo
enddo
end subroutine dist_cal_2d

subroutine dist_cal_3d(dist,nz,ny,nx,dz,dy,dx,zs,ys,xs)
real, intent(in) :: dz,dy,dx,zs,ys,xs
integer, intent(in) :: nz,ny,nx
real, intent(out) :: dist(:,:,:)
real              :: x,y,z,dist_x2,dist_y2,dist_z2,dist_xy2
integer           :: ix, iy, iz
do ix=1,nx
   x=(ix-1)*dx;   dist_x2=(xs-x)**2.0
   do iy=1,ny
      y=(iy-1)*dy;   dist_y2=(ys-y)**2.0
      dist_xy2=dist_x2+dist_y2
      do iz=1,nz
         z=(iz-1)*dz;   dist_z2=(zs-z)**2.0
         dist(iz,iy,ix)=sqrt(dist_xy2+dist_z2)
         if (dist(iz,iy,ix) < dx/2.0) then
            dist(iz,iy,ix)=dx
         endif
      enddo
   enddo
enddo
end subroutine dist_cal_3d

subroutine oblq_cal_2d(dist,oblq,nz,nx,dz)
real, intent(in)  ::  dist(:,:),dz
real, intent(out) ::  oblq(:,:)
integer, intent(in) :: nz,nx
integer           ::  iz,ix
do ix=1,nx
   do iz=1,nz
      if (dist(iz,ix) < dz+0.0005) then
         !oblq(iz,ix)=1.0
         oblq(iz,ix)=0.0
      else
         oblq(iz,ix)=abs(iz-1)*dz/dist(iz,ix)
      endif
   enddo
enddo
end subroutine oblq_cal_2d

subroutine oblq_cal_3d(dist,oblq,nz,ny,nx,dz)
real, intent(in)  ::  dist(:,:,:),dz
real, intent(out) ::  oblq(:,:,:)
integer, intent(in) :: nz,ny,nx
integer           ::  iz,iy,ix
do ix=1,nx
   do iy=1,ny
      do iz=1,nz
         if (dist(iz,iy,ix) < dz+0.0005) then
            !oblq(iz,ix)=1.0
            oblq(iz,iy,ix)=0.0
         else
            oblq(iz,iy,ix)=abs(iz-1)*dz/dist(iz,iy,ix)
         endif
      enddo
   enddo
enddo
end subroutine oblq_cal_3d

!=====================================================================

!subroutine hanning(taper,n)
!integer, intent(in)  :: n
!real,    intent(out) :: taper(:)
!integer              :: i
!real                 :: m
!m=2.0*n+1.0
!do i=1,n
!  taper(i)=0.5*(1-cos(2*PI*i/m))
!enddo
!end subroutine hanning

!=====================================================================

subroutine cut_core1(n1_in,i1,i2,n11,n12,n1_out)
integer, intent(in)  :: n1_in, i1, i2
integer, intent(out) :: n11, n12, n1_out
if (i1<=0) then
   n11=-i1+1
else
   n11=0
endif
if (i2>n1_in) then
   n12=i2-n1_in
else
   n12=0
endif
n1_out=n1_in+n11+n12
end subroutine cut_core1

!=====================================================================

subroutine cut_int_1d(data_in,n1_in,i1,i2,n1_out,data_out)
integer, intent(in)  :: n1_in,i1,i2,data_in(:)
integer, intent(out) :: n1_out,data_out(:)
integer, allocatable :: data_tmp(:)
integer              :: i, n11, n12
call cut_core1(n1_in,i1,i2,n11,n12,n1_out)
allocate(data_tmp(n1_out))
call padarray_int_1d(data_in,n1_in,n11,n12,data_tmp,"replicate")
if (i1<=0) then
   do i=1,i2-i1+1
      data_out(i)=data_tmp(i)
   enddo
else
   do i=i1,i2
      data_out(i-i1+1)=data_tmp(i)
   enddo
endif
deallocate(data_tmp)
end subroutine cut_int_1d

!=====================================================================

subroutine cut_real_1d(data_in,n1_in,i1,i2,n1_out,data_out)
integer, intent(in)  :: n1_in,i1,i2
real,    intent(in)  :: data_in(:)
integer, intent(out) :: n1_out
real,    intent(out) :: data_out(:)
real,    allocatable :: data_tmp(:)
integer :: i, n11, n12
call cut_core1(n1_in,i1,i2,n11,n12,n1_out)
allocate(data_tmp(n1_out))
call padarray_real_1d(data_in,n1_in,n11,n12,data_tmp,"replicate")
if (i1<=0) then
   do i=i1-i1+1,i2-i1+1
      data_out(i)=data_tmp(i)
   enddo
else
   do i=i1,i2
      data_out(i-i1+1)=data_tmp(i)
   enddo
endif
deallocate(data_tmp)
end subroutine cut_real_1d

!=====================================================================

subroutine cut_int_2d(data_in,n1_in,i1,i2,n1_out,n2_in,j1,j2,n2_out,data_out)
integer, intent(in)  :: n1_in,i1,i2,n2_in,j1,j2,data_in(:,:)
integer, intent(out) :: n1_out,n2_out,data_out(:,:)
integer, allocatable :: data_tmp(:,:)
integer :: i, n11, n12, j, n21, n22
call cut_core1(n1_in,i1,i2,n11,n12,n1_out)
call cut_core1(n2_in,j1,j2,n21,n22,n2_out)
allocate(data_tmp(n1_out,n2_out))
call padarray_int_2d(data_in,n1_in,n11,n12,n2_in,n21,n22,data_tmp,"replicate")
if (j1<=0) then
   if (i1<=0) then
      do j=1,j2-j1+1
         do i=1,i2-i1+1
            data_out(i,j)=data_tmp(i,j)
         enddo
      enddo
   else
      do j=1,j2-j1+1
         do i=i1,i2
            data_out(i-i1+1,j)=data_tmp(i,j)
         enddo
      enddo
   endif
else
   if (i1<=0) then
      do j=j1,j2
         do i=1,i2-i1+1
            data_out(i,j-j1+1)=data_tmp(i,j)
         enddo
      enddo
   else
      do j=j1,j2
         do i=i1,i2
            data_out(i-i1+1,j-j1+1)=data_tmp(i,j)
         enddo
      enddo
   endif
endif
deallocate(data_tmp)
end subroutine cut_int_2d

!====================================================================

subroutine cut_real_2d(data_in,n1_in,i1,i2,n1_out,n2_in,j1,j2,n2_out,data_out)
integer, intent(in)  :: n1_in,i1,i2,n2_in,j1,j2
real,    intent(in)  :: data_in(:,:)
integer, intent(out) :: n1_out,n2_out
real,    intent(out) :: data_out(:,:)
real,    allocatable :: data_tmp(:,:)
integer :: i, n11, n12, j, n21, n22
call cut_core1(n1_in,i1,i2,n11,n12,n1_out)
call cut_core1(n2_in,j1,j2,n21,n22,n2_out)
allocate(data_tmp(n1_out,n2_out))
call padarray_real_2d(data_in,n1_in,n11,n12,n2_in,n21,n22,data_tmp,"replicate")
if (j1<=0) then
   if (i1<=0) then
      do j=1,j2-j1+1
         do i=1,i2-i1+1
            data_out(i,j)=data_tmp(i,j)
         enddo
      enddo
   else
      do j=1,j2-j1+1
         do i=i1,i2
            data_out(i-i1+1,j)=data_tmp(i,j)
         enddo
      enddo
   endif
else
   if (i1<=0) then
      do j=j1,j2
         do i=1,i2-i1+1
            data_out(i,j-j1+1)=data_tmp(i,j)
         enddo
      enddo
   else
      do j=j1,j2
         do i=i1,i2
            data_out(i-i1+1,j-j1+1)=data_tmp(i,j)
         enddo
      enddo
   endif
endif
deallocate(data_tmp)
end subroutine cut_real_2d

!====================================================================

subroutine cut_int_3d(data_in,n1_in,i1,i2,n1_out,n2_in,j1,j2,n2_out,&
   n3_in,k1,k2,n3_out,data_out)
integer, intent(in)  :: n1_in,i1,i2,n2_in,j1,j2,n3_in,k1,k2,data_in(:,:,:)
integer, intent(out) :: n1_out,n2_out,n3_out,data_out(:,:,:)
integer, allocatable :: data_tmp(:,:,:)
integer :: i, n11, n12, j, n21, n22, k, n31, n32
call cut_core1(n1_in,i1,i2,n11,n12,n1_out)
call cut_core1(n2_in,j1,j2,n21,n22,n2_out)
call cut_core1(n3_in,k1,k2,n31,n32,n3_out)
allocate(data_tmp(n1_out,n2_out,n3_out))
call padarray_int_3d(data_in,n1_in,n11,n12,n2_in,n21,n22,n3_in,n31,&
                     n32,data_tmp,"replicate")
if (k1<=0) then
   if (j1<=0) then
      if (i1<=0) then
         do k=1,k2-k1+1
            do j=1,j2-j1+1
               do i=1,i2-i1+1
                  data_out(i,j,k)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      else
         do k=1,k2-k1+1
            do j=1,j2-j1+1
               do i=i1,i2
                  data_out(i-i1+1,j,k)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      endif
   else
      if (i1<=0) then
         do k=1,k2-k1+1
            do j=j1,j2
               do i=1,i2-i1+1
                  data_out(i,j-j1+1,k)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      else
         do k=1,k2-k1+1
            do j=j1,j2
               do i=i1,i2
                  data_out(i-i1+1,j-j1+1,k)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      endif
   endif
else
   if (j1<=0) then
      if (i1<=0) then
         do k=k1,k2
            do j=1,j2-j1+1
               do i=1,i2-i1+1
                  data_out(i,j,k-k1+1)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      else
         do k=k1,k2
            do j=1,j2-j1+1
               do i=i1,i2
                  data_out(i-i1+1,j,k-k1+1)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      endif
   else
      if (i1<=0) then
         do k=k1,k2
            do j=j1,j2
               do i=1,i2-i1+1
                  data_out(i,j-j1+1,k-k1+1)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      else
         do k=k1,k2
            do j=j1,j2
               do i=i1,i2
                  data_out(i-i1+1,j-j1+1,k-k1+1)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      endif
   endif
endif
deallocate(data_tmp)
end subroutine cut_int_3d

!====================================================================

subroutine cut_real_3d(data_in,n1_in,i1,i2,n1_out,n2_in,j1,j2,n2_out,&
   n3_in,k1,k2,n3_out,data_out)
integer, intent(in)  :: n1_in,i1,i2,n2_in,j1,j2,n3_in,k1,k2
real,    intent(in)  :: data_in(:,:,:)
integer, intent(out) :: n1_out,n2_out,n3_out
real,    intent(out) :: data_out(:,:,:)
real,    allocatable :: data_tmp(:,:,:)
integer :: i, n11, n12, j, n21, n22, k, n31, n32
call cut_core1(n1_in,i1,i2,n11,n12,n1_out)
call cut_core1(n2_in,j1,j2,n21,n22,n2_out)
call cut_core1(n3_in,k1,k2,n31,n32,n3_out)
allocate(data_tmp(n1_out,n2_out,n3_out))
call padarray_real_3d(data_in,n1_in,n11,n12,n2_in,n21,n22,n3_in,n31,n32,&
                      data_tmp,"replicate")
if (k1<=0) then
   if (j1<=0) then
      if (i1<=0) then
         do k=1,k2-k1+1
            do j=1,j2-j1+1
               do i=1,i2-i1+1
                  data_out(i,j,k)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      else
         do k=1,k2-k1+1
            do j=1,j2-j1+1
               do i=i1,i2
                  data_out(i-i1+1,j,k)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      endif
   else
      if (i1<=0) then
         do k=1,k2-k1+1
            do j=j1,j2
               do i=1,i2-i1+1
                  data_out(i,j-j1+1,k)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      else
         do k=1,k2-k1+1
            do j=j1,j2
               do i=i1,i2
                  data_out(i-i1+1,j-j1+1,k)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      endif
   endif
else
   if (j1<=0) then
      if (i1<=0) then
         do k=k1,k2
            do j=1,j2-j1+1
               do i=1,i2-i1+1
                  data_out(i,j,k-k1+1)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      else
         do k=k1,k2
            do j=1,j2-j1+1
               do i=i1,i2
                  data_out(i-i1+1,j,k-k1+1)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      endif
   else
      if (i1<=0) then
         do k=k1,k2
            do j=j1,j2
               do i=1,i2-i1+1
                  data_out(i,j-j1+1,k-k1+1)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      else
         do k=k1,k2
            do j=j1,j2
               do i=i1,i2
                  data_out(i-i1+1,j-j1+1,k-k1+1)=data_tmp(i,j,k)
               enddo
            enddo
         enddo
      endif
   endif
endif
deallocate(data_tmp)
end subroutine cut_real_3d

!=====================================================================

subroutine padarray_real_1d(d_in,n1,n1_l,n1_r,d_out,pad_type,dc)
integer, intent(in)  :: n1, n1_l, n1_r
real,    intent(in)  :: d_in(:)
real,intent(in),optional ::dc
real,    intent(out) :: d_out(:)
character(len=*),intent(in) :: pad_type
integer              :: n11, n12, i1
n11=n1_l+1;   n12=n1_l+n1
d_out(n11:n12)=d_in
if (pad_type.eq."replicate") then
   do i1=1,n1_l
      d_out(i1) = d_out(n11)
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = d_out(n12)
   enddo
elseif (pad_type.eq."mirror") then
   do i1=1,n1_l
      d_out(i1)=d_out(n11+n1_l+1-i1)
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = d_out(n12-i1)
   enddo
elseif (pad_type.eq."symmetric") then
   do i1=1,n1_l
      d_out(i1)=d_out(n11+n1_l-i1)
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = d_out(n12-i1+1)
   enddo
elseif (pad_type.eq."constant") then
   do i1=1,n1_l
      d_out(i1)=dc
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = dc
   enddo  
elseif (pad_type.eq."coefficient") then
   do i1=1,n1_l
      d_out(i1)=d_out(n11)*dc
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = d_out(n12)*dc
   enddo  
endif   
end subroutine padarray_real_1d

!=====================================================================

subroutine padarray_real_2d(d_in,n1,n1_u,n1_d,n2,n2_l,n2_r,d_out,pad_type,dc)
integer, intent(in)  :: n1, n1_u, n1_d, n2, n2_l, n2_r
real,    intent(in)  :: d_in(:,:)
real,intent(in),optional :: dc
real,    intent(out) :: d_out(:,:)
character(len=*),intent(in) :: pad_type
integer              :: n11, n12, i1, n21, n22, i2
n11=n1_u+1;   n12=n1_u+n1;   n21=n2_l+1;   n22=n2_l+n2
d_out(n11:n12,n21:n22)=d_in
!write(*,*)"pad 2d",n11,n12,n21,n22
do i1=n11,n12
   call padarray_real_1d(d_in(i1-n11+1,:),n2,n2_l,n2_r,d_out(i1,:),pad_type,dc)
enddo
do i2=1,n2+n2_l+n2_r
   call padarray_real_1d(d_out(n11:n12,i2),n1,n1_u,n1_d,d_out(:,i2),pad_type,dc)
enddo
end subroutine padarray_real_2d

!=====================================================================

subroutine padarray_real_3d(d_in,n1,n1_u,n1_d,n2,n2_b,n2_f,n3,n3_l,n3_r,d_out,pad_type,dc)
integer, intent(in)  :: n1,n1_u,n1_d,n2,n2_f,n2_b,n3,n3_l,n3_r
real,    intent(in)  :: d_in(:,:,:)
real,intent(in),optional :: dc
real,    intent(out) :: d_out(:,:,:)
character(len=*),intent(in) :: pad_type
integer              :: n11, n12, i1, n21, n22, i2, n31, n32, i3
n11=n1_u+1;   n12=n1_u+n1;   n21=n2_b+1
n22=n2_b+n2;  n31=n3_l+1;    n32=n3_l+n3
d_out(n11:n12,n21:n22,n31:n32)=d_in
do i1=n11,n12
   call padarray_real_2d(d_in(i1-n11+1,:,:),n2,n2_b,n2_f,n3,n3_l,n3_r,&
                         d_out(i1,:,:),pad_type,dc)
enddo
do i3=1,n3+n3_l+n3_r
   do i2=1,n2+n2_f+n2_b
      call padarray_real_1d(d_out(n11:n12,i2,i3),n1,n1_u,n1_d,&
                           d_out(:,i2,i3),pad_type,dc)
   enddo
enddo
end subroutine padarray_real_3d

!=====================================================================

subroutine padarray_int_1d(d_in,n1,n1_l,n1_r,d_out,pad_type,dc)
integer, intent(in)  :: n1, n1_l, n1_r
integer, intent(in)  :: d_in(:)
integer,intent(in),optional :: dc
integer, intent(out) :: d_out(:)
character(len=*),intent(in) :: pad_type
integer              :: n11, n12, i1
n11=n1_l+1;   n12=n1_l+n1
d_out(n11:n12)=d_in
if (pad_type.eq."replicate") then
   do i1=1,n1_l
      d_out(i1) = d_out(n11)
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = d_out(n12)
   enddo
elseif (pad_type.eq."mirror") then
   do i1=1,n1_l
      d_out(i1)=d_out(n11+n1_l+1-i1)
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = d_out(n12-i1)
   enddo
elseif (pad_type.eq."symmetric") then
   do i1=1,n1_l
      d_out(i1)=d_out(n11+n1_l-i1)
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = d_out(n12-i1+1)
   enddo
elseif (pad_type.eq."constant") then
   do i1=1,n1_l
      d_out(i1)=dc
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = dc
   enddo
elseif (pad_type.eq."coefficient") then
   do i1=1,n1_l
      d_out(i1)=d_out(n11)*dc
   enddo
   do i1=1,n1_r
      d_out(n12+i1) = d_out(n12)*dc
   enddo
endif   
end subroutine padarray_int_1d

!=====================================================================

subroutine padarray_int_2d(d_in,n1,n1_u,n1_d,n2,n2_l,n2_r,d_out,pad_type,dc)
integer, intent(in)  :: n1, n1_u, n1_d, n2, n2_l, n2_r
integer, intent(in)  :: d_in(:,:)
integer,intent(in),optional :: dc
integer, intent(out) :: d_out(:,:)
character(len=*),intent(in) :: pad_type
integer              :: n11, n12, i1, n21, n22, i2
n11=n1_u+1;   n12=n1_u+n1;   n21=n2_l+1;   n22=n2_l+n2
d_out(n11:n12,n21:n22)=d_in
do i1=n11,n12
   call padarray_int_1d(d_in(i1-n11+1,:),n2,n2_l,n2_r,d_out(i1,:),pad_type,dc)
enddo
do i2=1,n2+n2_l+n2_r
   call padarray_int_1d(d_out(n11:n12,i2),n1,n1_u,n1_d,d_out(:,i2),pad_type,dc)
enddo
end subroutine padarray_int_2d

!=====================================================================

subroutine padarray_int_3d(d_in,n1,n1_u,n1_d,n2,n2_b,n2_f,n3,n3_l,n3_r,d_out,pad_type,dc)
integer, intent(in)  :: n1,n1_u,n1_d,n2,n2_f,n2_b,n3,n3_l,n3_r
integer, intent(in)  :: d_in(:,:,:)
integer,intent(in),optional :: dc
integer, intent(out) :: d_out(:,:,:)
character(len=*),intent(in) :: pad_type
integer              :: n11, n12, i1, n21, n22, i2, n31, n32, i3
n11=n1_u+1;n12=n1_u+n1;n21=n2_b+1;n22=n2_b+n2;n31=n3_l+1;n32=n3_l+n3
d_out(n11:n12,n21:n22,n31:n32)=d_in
do i1=n11,n12
   call padarray_int_2d(d_in(i1-n11+1,:,:),n2,n2_b,n2_f,n3,n3_l,n3_r,&
                        d_out(i1,:,:),pad_type,dc)
enddo
do i3=1,n3+n3_l+n3_r
   do i2=1,n2+n2_f+n2_b
      call padarray_int_1d(d_out(n11:n12,i2,i3),n1,n1_u,n1_d,&
                           d_out(:,i2,i3),pad_type,dc)
   enddo
enddo
end subroutine padarray_int_3d

!=====================================================================

subroutine smooth_real_1d_v1(d_in,n1,n1p,d_out,smooth_type)
integer, intent(in)  :: n1, n1p 
real,    intent(in)  :: d_in(:)
real,    intent(out) :: d_out(:)
character(len=*),intent(in) :: smooth_type
real,allocatable     :: d_tmp(:)
integer              :: i1,hn1p
hn1p=(n1p-1)/2
allocate(d_tmp(n1+2*hn1p))
call padarray(d_in,n1,hn1p,hn1p,d_tmp,"replicate")
if (smooth_type.eq."moving") then
   do i1=hn1p+1,n1+hn1p
      d_out(i1-hn1p)=average(d_tmp(i1-hn1p:i1+hn1p))
   enddo
endif
deallocate(d_tmp)
end subroutine smooth_real_1d_v1

!=====================================================================

subroutine smooth_real_1d_v2(d,n1,n1p,smooth_type)
integer, intent(in)    :: n1, n1p 
real,    intent(inout) :: d(:)
character(len=*),intent(in) :: smooth_type
real,allocatable       :: d_tmp(:)
integer                :: i1,hn1p
hn1p=(n1p-1)/2
allocate(d_tmp(n1+2*hn1p))
call padarray(d,n1,hn1p,hn1p,d_tmp,"replicate")
if (smooth_type.eq."moving") then
   do i1=hn1p+1,n1+hn1p
      d(i1-hn1p)=average(d_tmp(i1-hn1p:i1+hn1p))
   enddo
endif
deallocate(d_tmp)
end subroutine smooth_real_1d_v2

!=====================================================================

subroutine smooth_real_2d_v1(d_in,n1,n2,n1p,n2p,d_out,smooth_type)
integer, intent(in)  :: n1, n2, n1p, n2p
real,    intent(in)  :: d_in(:,:)
real,    intent(out) :: d_out(:,:)
character(len=*),intent(in) :: smooth_type
real,allocatable     :: d_tmp(:,:)
integer              :: i1, i2, hn1p, hn2p
hn1p=(n1p-1)/2;   hn2p=(n2p-1)/2
allocate(d_tmp(n1+2*hn1p,n2+2*hn2p))
call padarray(d_in,n1,hn1p,hn1p,n2,hn2p,hn2p,d_tmp,"replicate")
if (smooth_type.eq."moving") then
   do i2=hn2p+1,n2+hn2p
      do i1=hn1p+1,n1+hn1p
         d_out(i1-hn1p,i2-hn2p)=average(d_tmp(i1-hn1p:i1+hn1p,i2-hn2p:i2+hn2p))
      enddo
   enddo
endif
deallocate(d_tmp)
end subroutine smooth_real_2d_v1

!=====================================================================

subroutine smooth_real_2d_v2(d,n1,n2,n1p,n2p,smooth_type)
integer, intent(in)  :: n1, n2, n1p, n2p
real,    intent(inout)  :: d(:,:)
character(len=*),intent(in) :: smooth_type
real,allocatable     :: d_tmp(:,:)
integer              :: i1, i2, hn1p, hn2p
hn1p=(n1p-1)/2;   hn2p=(n2p-1)/2
allocate(d_tmp(n1+2*hn1p,n2+2*hn2p))
call padarray(d,n1,hn1p,hn1p,n2,hn2p,hn2p,d_tmp,"replicate")
if (smooth_type.eq."moving") then
   do i2=hn2p+1,n2+hn2p
      do i1=hn1p+1,n1+hn1p
         d(i1-hn1p,i2-hn2p)=average(d_tmp(i1-hn1p:i1+hn1p,i2-hn2p:i2+hn2p))
      enddo
   enddo
endif
deallocate(d_tmp)
end subroutine smooth_real_2d_v2


!------------------------------------------------------------------------------
! Square smoothing filter
subroutine smooth_gk_2d_v1(x, n1, n2, x1, hw)

integer, intent(in) :: n1, n2, hw
real, intent(in) :: x(:,:)
real, intent(out) :: x1(:,:)
integer             :: i, j, k, i1, j1
real                :: sumx
!real, allocatable   :: temp(:,:)

!allocate(temp(n1,n2))
!temp = x
do i=1,n1
  do j=1,n2
    k = 0
    sumx = 0.0
    do i1=max(1,i-hw),min(n1,i+hw)
      do j1=max(1,j-hw),min(n2,j+hw)
        k = k+1
        sumx = sumx + x(i1,j1)
      enddo
    enddo
!    if (k > 0) then
      x1(i,j) = sumx/real(k)
!    endif
  enddo
enddo
!deallocate(temp)

end subroutine smooth_gk_2d_v1

!================================================================
! Square smoothing filter
subroutine smooth_gk_2d_v2(x, n1, n2, hw)

integer, intent(in) :: n1, n2, hw
real, intent(inout) :: x(:,:)
integer             :: i, j, k, i1, j1
real                :: sumx
real, allocatable   :: temp(:,:)

allocate(temp(n1,n2))
temp = x
do i=1,n1
  do j=1,n2
    k = 0
    sumx = 0.0
    do i1=max(1,i-hw),min(n1,i+hw)
      do j1=max(1,j-hw),min(n2,j+hw)
        k = k+1
        sumx = sumx + temp(i1,j1)
      enddo
    enddo
!    if (k > 0) then
      x(i,j) = sumx/real(k)
!    endif
  enddo
enddo
deallocate(temp)

end subroutine smooth_gk_2d_v2

!=====================================================================

subroutine smooth_real_3d_v1(d_in,n1,n2,n3,n1p,n2p,n3p,d_out,smooth_type)
integer, intent(in)  :: n1, n2, n3, n1p, n2p, n3p
real,    intent(in)  :: d_in(:,:,:)
real,    intent(out) :: d_out(:,:,:)
character(len=*),intent(in) :: smooth_type
real,allocatable     :: d_tmp(:,:,:)
integer              :: i1, i2, i3, hn1p, hn2p, hn3p
hn1p=(n1p-1)/2;   hn2p=(n2p-1)/2;   hn3p=(n3p-1)/2
allocate(d_tmp(n1+2*hn1p,n2+2*hn2p,n3+2*hn3p))
call padarray(d_in,n1,hn1p,hn1p,n2,hn2p,hn2p,n3,hn3p,hn3p,d_tmp,"replicate")
if (smooth_type.eq."moving") then
   do i3=hn3p+1,n3+hn3p
      do i2=hn2p+1,n2+hn2p
         do i1=hn1p+1,n1+hn1p
            d_out(i1-hn1p,i2-hn2p,i3-hn3p)&
           =average(d_tmp(i1-hn1p:i1+hn1p,i2-hn2p:i2+hn2p,i3-hn3p:i3+hn3p))
         enddo
      enddo
   enddo
endif
deallocate(d_tmp)
end subroutine smooth_real_3d_v1

!=====================================================================

subroutine smooth_real_3d_v2(d,n1,n2,n3,n1p,n2p,n3p,smooth_type)
integer, intent(in)  :: n1, n2, n3, n1p, n2p, n3p
real,    intent(inout)  :: d(:,:,:)
character(len=*),intent(in) :: smooth_type
real,allocatable     :: d_tmp(:,:,:)
integer              :: i1, i2, i3, hn1p, hn2p, hn3p
hn1p=(n1p-1)/2;   hn2p=(n2p-1)/2;   hn3p=(n3p-1)/2
allocate(d_tmp(n1+2*hn1p,n2+2*hn2p,n3+2*hn3p))
call padarray(d,n1,hn1p,hn1p,n2,hn2p,hn2p,n3,hn3p,hn3p,d_tmp,"replicate")
if (smooth_type.eq."moving") then
   do i3=hn3p+1,n3+hn3p
      do i2=hn2p+1,n2+hn2p
         do i1=hn1p+1,n1+hn1p
            d(i1-hn1p,i2-hn2p,i3-hn3p)&
           =average(d_tmp(i1-hn1p:i1+hn1p,i2-hn2p:i2+hn2p,i3-hn3p:i3+hn3p))
         enddo
      enddo
   enddo
endif
deallocate(d_tmp)
end subroutine smooth_real_3d_v2

!=====================================================================

subroutine cut_and_pad_int1d(d_in,n1_in,i1,i2,npml,n1,d_out)
integer,   intent(in) :: d_in(:),n1_in,i1,i2,npml,n1
integer,  allocatable :: d_temp(:)
integer,allocatable,intent(out) ::d_out(:)
integer               :: n1_out, n1_tmp
call allocate_and_initial(d_temp,n1)
n1_out=n1+2*npml
call allocate_and_initial(d_out,n1_out)
!call cut(d_in,n1_in,i1,i2,n1,d_temp)
call cut(d_in,n1_in,i1,i2,n1_tmp,d_temp)
call padarray(d_temp,n1,npml,npml,d_out,"replicate")
call deallocate_and_free(d_temp)
end subroutine cut_and_pad_int1d

!====================================================================

subroutine cut_and_pad_int2d(d_in,n1_in,n2_in,i1,i2,j1,j2,npml,n1,n2,d_out)
integer,   intent(in) :: d_in(:,:),n1_in,n2_in,i1,i2,j1,j2,npml,n1,n2
integer,  allocatable :: d_temp(:,:)
integer,allocatable,intent(out) :: d_out(:,:)
integer               :: n1_out,n2_out,n1_tmp,n2_tmp
call allocate_and_initial(d_temp,n1,n2)
n1_out=n1+2*npml;   n2_out=n2+2*npml
call allocate_and_initial(d_out,n1_out,n2_out)
call cut(d_in,n1_in,i1,i2,n1_tmp,n2_in,j1,j2,n2_tmp,d_temp)
call padarray(d_temp,n1,npml,npml,n2,npml,npml,d_out,"replicate")
call deallocate_and_free(d_temp)
end subroutine cut_and_pad_int2d

!====================================================================

subroutine cut_and_pad_real2d(d_in,n1_in,n2_in,i1,i2,j1,j2,npml,n1,n2,d_out)
integer,   intent(in) :: n1_in,n2_in,i1,i2,j1,j2,npml,n1,n2
real,      intent(in) :: d_in(:,:)
real,     allocatable :: d_temp(:,:)
real,allocatable,intent(out) :: d_out(:,:)
integer               :: n1_out,n2_out,n1_tmp,n2_tmp
call allocate_and_initial(d_temp,n1,n2)
n1_out=n1+2*npml;   n2_out=n2+2*npml
call allocate_and_initial(d_out,n1_out,n2_out)
call cut(d_in,n1_in,i1,i2,n1_tmp,n2_in,j1,j2,n2_tmp,d_temp)
call padarray(d_temp,n1,npml,npml,n2,npml,npml,d_out,"replicate")
call deallocate_and_free(d_temp)
end subroutine cut_and_pad_real2d

!====================================================================

subroutine cut_and_pad_real3d(d_in,n1_in,n2_in,n3_in,i1,i2,j1,j2,&
               k1,k2,npml,n1,n2,n3,d_out)
integer,   intent(in) :: n1_in,n2_in,n3_in,i1,i2,j1,j2,k1,k2,npml,n1,n2,n3
real,      intent(in) :: d_in(:,:,:)
real,     allocatable :: d_temp(:,:,:)
real,allocatable,intent(out) :: d_out(:,:,:)
integer               :: n1_out,n2_out,n3_out,n1_tmp,n2_tmp,n3_tmp
call allocate_and_initial(d_temp,n1,n2,n3)
n1_out=n1+2*npml;   n2_out=n2+2*npml;   n3_out=n3+2*npml
call allocate_and_initial(d_out,n1_out,n2_out,n3_out)
call cut(d_in,n1_in,i1,i2,n1_tmp,n2_in,j1,j2,n2_tmp,n3_in,k1,k2,n3_tmp,d_temp)
call padarray(d_temp,n1,npml,npml,n2,npml,npml,n3,npml,npml,d_out,"replicate")
call deallocate_and_free(d_temp)
end subroutine cut_and_pad_real3d

!====================================================================

subroutine copy_variable_1i(i1,oi1)
integer, intent(in)  :: i1
integer, intent(out) :: oi1
oi1=i1
end subroutine copy_variable_1i

!====================================================================

subroutine copy_variable_2i(i1,oi1,i2,oi2)
integer, intent(in)  :: i1,i2
integer, intent(out) :: oi1,oi2
oi1=i1;   oi2=i2
end subroutine copy_variable_2i

!====================================================================

subroutine copy_variable_3i(i1,oi1,i2,oi2,i3,oi3)
integer, intent(in)  :: i1,i2,i3
integer, intent(out) :: oi1,oi2,oi3
oi1=i1;   oi2=i2;   oi3=i3
end subroutine copy_variable_3i

!====================================================================

subroutine copy_variable_4i(i1,oi1,i2,oi2,i3,oi3,i4,oi4)
integer, intent(in)  :: i1,i2,i3,i4
integer, intent(out) :: oi1,oi2,oi3,oi4
oi1=i1;   oi2=i2;   oi3=i3;   oi4=i4
end subroutine copy_variable_4i

!====================================================================

subroutine copy_variable_5i(i1,oi1,i2,oi2,i3,oi3,i4,oi4,i5,oi5)
integer, intent(in)  :: i1,i2,i3,i4,i5
integer, intent(out) :: oi1,oi2,oi3,oi4,oi5
oi1=i1;   oi2=i2;   oi3=i3;   oi4=i4;   oi5=i5
end subroutine copy_variable_5i

!====================================================================

subroutine copy_variable_1r(r1,or1)
real, intent(in)  :: r1
real, intent(out) :: or1
or1=r1
end subroutine copy_variable_1r

!====================================================================

subroutine copy_variable_2r(r1,or1,r2,or2)
real, intent(in)  :: r1,r2
real, intent(out) :: or1,or2
or1=r1;   or2=r2
end subroutine copy_variable_2r

!====================================================================

subroutine copy_variable_3r(r1,or1,r2,or2,r3,or3)
real, intent(in)  :: r1,r2,r3
real, intent(out) :: or1,or2,or3
or1=r1;   or2=r2;   or3=r3
end subroutine copy_variable_3r

!====================================================================

subroutine copy_variable_4r(r1,or1,r2,or2,r3,or3,r4,or4)
real, intent(in)  :: r1,r2,r3,r4
real, intent(out) :: or1,or2,or3,or4
or1=r1;   or2=r2;   or3=r3;   or4=r4
end subroutine copy_variable_4r

!====================================================================

subroutine copy_variable_5r(r1,or1,r2,or2,r3,or3,r4,or4,r5,or5)
real, intent(in)  :: r1,r2,r3,r4,r5
real, intent(out) :: or1,or2,or3,or4,or5
or1=r1;   or2=r2;   or3=r3;   or4=r4;   or5=r5
end subroutine copy_variable_5r

!====================================================================

subroutine copy_variable_6r(r1,or1,r2,or2,r3,or3,r4,or4,r5,or5,r6,or6)
real, intent(in)  :: r1,r2,r3,r4,r5,r6
real, intent(out) :: or1,or2,or3,or4,or5,or6
or1=r1;   or2=r2;   or3=r3;   or4=r4;   or5=r5;   or6=r6
end subroutine copy_variable_6r

!====================================================================

subroutine copy_variable_1l(l1,ol1)
logical, intent(in)  :: l1
logical, intent(out) :: ol1
ol1=l1
end subroutine copy_variable_1l

!====================================================================

subroutine copy_variable_2l(l1,ol1,l2,ol2)
logical, intent(in)  :: l1,l2
logical, intent(out) :: ol1,ol2
ol1=l1;   ol2=l2
end subroutine copy_variable_2l

!====================================================================

subroutine copy_variable_3l(l1,ol1,l2,ol2,l3,ol3)
logical, intent(in)  :: l1,l2,l3
logical, intent(out) :: ol1,ol2,ol3
ol1=l1;   ol2=l2;   ol3=l3
end subroutine copy_variable_3l
!====================================================================

subroutine copy_variable_1c(c1,oc1)
complex, intent(in)  :: c1
complex, intent(out) :: oc1
oc1=c1
end subroutine copy_variable_1c

!====================================================================

subroutine copy_variable_2c(c1,oc1,c2,oc2)
complex, intent(in)  :: c1,c2
complex, intent(out) :: oc1,oc2
oc1=c1;   oc2=c2
end subroutine copy_variable_2c

!====================================================================

subroutine copy_variable_3c(c1,oc1,c2,oc2,c3,oc3)
complex, intent(in)  :: c1,c2,c3
complex, intent(out) :: oc1,oc2,oc3
oc1=c1;   oc2=c2;   oc3=c3
end subroutine copy_variable_3c

!====================================================================

subroutine copy_variable_1ch(ch1,och1)
character(len=*),intent(in)  :: ch1
character(len=*),intent(out) :: och1
och1=ch1
end subroutine copy_variable_1ch

!====================================================================

subroutine copy_variable_2ch(ch1,och1,ch2,och2)
character(len=*),intent(in)  :: ch1, ch2
character(len=*),intent(out) :: och1,och2
och1=ch1;   och2=ch2
end subroutine copy_variable_2ch

!====================================================================

subroutine copy_variable_3ch(ch1,och1,ch2,och2,ch3,och3)
character(len=*),intent(in)  :: ch1,  ch2,  ch3
character(len=*),intent(out) :: och1, och2, och3
och1=ch1;   och2=ch2;   och3=ch3
end subroutine copy_variable_3ch

!====================================================================

!====================================================================

subroutine wf2d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
real,target, intent(in) :: d1_0(:,:),d1_1(:,:),d1(:,:)
real,pointer    :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
d1p_0=>d1_0;   d1p_1=>d1_1;   d1p=>d1
end subroutine wf2d_assign_1

!====================================================================

subroutine wf2d_assign_2(d1_0,d1p_0,d1_1,d1p_1,d1,d1p,&
                         d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
real, intent(in) :: d1_0(:,:),d1_1(:,:),d1(:,:)
real, intent(in) :: d2_0(:,:),d2_1(:,:),d2(:,:)
real, pointer    :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
real, pointer    :: d2p_0(:,:),d2p_1(:,:),d2p(:,:)
call wf2d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
call wf2d_assign_1(d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
end subroutine wf2d_assign_2

!====================================================================

subroutine wf2d_assign_3(d1_0,d1p_0,d1_1,d1p_1,d1,d1p,&
                         d2_0,d2p_0,d2_1,d2p_1,d2,d2p,&
                         d3_0,d3p_0,d3_1,d3p_1,d3,d3p)
real, intent(in) :: d1_0(:,:),d1_1(:,:),d1(:,:)
real, intent(in) :: d2_0(:,:),d2_1(:,:),d2(:,:)
real, intent(in) :: d3_0(:,:),d3_1(:,:),d3(:,:)
real, pointer    :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
real, pointer    :: d2p_0(:,:),d2p_1(:,:),d2p(:,:)
real, pointer    :: d3p_0(:,:),d3p_1(:,:),d3p(:,:)
call wf2d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
call wf2d_assign_1(d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
call wf2d_assign_1(d3_0,d3p_0,d3_1,d3p_1,d3,d3p)
end subroutine wf2d_assign_3

!====================================================================

subroutine wf2d_assign_4(d1_0,d1p_0,d1_1,d1p_1,d1,d1p,&
                         d2_0,d2p_0,d2_1,d2p_1,d2,d2p,&
                         d3_0,d3p_0,d3_1,d3p_1,d3,d3p,&
                         d4_0,d4p_0,d4_1,d4p_1,d4,d4p)
real, intent(in) :: d1_0(:,:),d1_1(:,:),d1(:,:)
real, intent(in) :: d2_0(:,:),d2_1(:,:),d2(:,:)
real, intent(in) :: d3_0(:,:),d3_1(:,:),d3(:,:)
real, intent(in) :: d4_0(:,:),d4_1(:,:),d4(:,:)
real, pointer    :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
real, pointer    :: d2p_0(:,:),d2p_1(:,:),d2p(:,:)
real, pointer    :: d3p_0(:,:),d3p_1(:,:),d3p(:,:)
real, pointer    :: d4p_0(:,:),d4p_1(:,:),d4p(:,:)
call wf2d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
call wf2d_assign_1(d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
call wf2d_assign_1(d3_0,d3p_0,d3_1,d3p_1,d3,d3p)
call wf2d_assign_1(d4_0,d4p_0,d4_1,d4p_1,d4,d4p)
end subroutine wf2d_assign_4

!====================================================================

subroutine wf2d_assign_5(d1_0,d1p_0,d1_1,d1p_1,d1,d1p,&
                         d2_0,d2p_0,d2_1,d2p_1,d2,d2p,&
                         d3_0,d3p_0,d3_1,d3p_1,d3,d3p,&
                         d4_0,d4p_0,d4_1,d4p_1,d4,d4p,&
                         d5_0,d5p_0,d5_1,d5p_1,d5,d5p)
real, intent(in) :: d1_0(:,:),d1_1(:,:),d1(:,:)
real, intent(in) :: d2_0(:,:),d2_1(:,:),d2(:,:)
real, intent(in) :: d3_0(:,:),d3_1(:,:),d3(:,:)
real, intent(in) :: d4_0(:,:),d4_1(:,:),d4(:,:)
real, intent(in) :: d5_0(:,:),d5_1(:,:),d5(:,:)
real, pointer    :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
real, pointer    :: d2p_0(:,:),d2p_1(:,:),d2p(:,:)
real, pointer    :: d3p_0(:,:),d3p_1(:,:),d3p(:,:)
real, pointer    :: d4p_0(:,:),d4p_1(:,:),d4p(:,:)
real, pointer    :: d5p_0(:,:),d5p_1(:,:),d5p(:,:)
call wf2d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
call wf2d_assign_1(d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
call wf2d_assign_1(d3_0,d3p_0,d3_1,d3p_1,d3,d3p)
call wf2d_assign_1(d4_0,d4p_0,d4_1,d4p_1,d4,d4p)
call wf2d_assign_1(d5_0,d5p_0,d5_1,d5p_1,d5,d5p)
end subroutine wf2d_assign_5

!====================================================================

subroutine wf3d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
real,target,intent(in) :: d1_0(:,:,:),d1_1(:,:,:),d1(:,:,:)
real, pointer    :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
d1p_0=>d1_0;   d1p_1=>d1_1;   d1p=>d1
end subroutine wf3d_assign_1

!====================================================================

subroutine wf3d_assign_2(d1_0,d1p_0,d1_1,d1p_1,d1,d1p,&
                         d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
real, intent(in) :: d1_0(:,:,:),d1_1(:,:,:),d1(:,:,:)
real, intent(in) :: d2_0(:,:,:),d2_1(:,:,:),d2(:,:,:)
real, pointer    :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
real, pointer    :: d2p_0(:,:,:),d2p_1(:,:,:),d2p(:,:,:)
call wf3d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
call wf3d_assign_1(d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
end subroutine wf3d_assign_2

!====================================================================

subroutine wf3d_assign_3(d1_0,d1p_0,d1_1,d1p_1,d1,d1p,&
                         d2_0,d2p_0,d2_1,d2p_1,d2,d2p,&
                         d3_0,d3p_0,d3_1,d3p_1,d3,d3p)
real, intent(in) :: d1_0(:,:,:),d1_1(:,:,:),d1(:,:,:)
real, intent(in) :: d2_0(:,:,:),d2_1(:,:,:),d2(:,:,:)
real, intent(in) :: d3_0(:,:,:),d3_1(:,:,:),d3(:,:,:)
real, pointer    :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
real, pointer    :: d2p_0(:,:,:),d2p_1(:,:,:),d2p(:,:,:)
real, pointer    :: d3p_0(:,:,:),d3p_1(:,:,:),d3p(:,:,:)
call wf3d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
call wf3d_assign_1(d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
call wf3d_assign_1(d3_0,d3p_0,d3_1,d3p_1,d3,d3p)
end subroutine wf3d_assign_3

!====================================================================

subroutine wf3d_assign_4(d1_0,d1p_0,d1_1,d1p_1,d1,d1p,&
                         d2_0,d2p_0,d2_1,d2p_1,d2,d2p,&
                         d3_0,d3p_0,d3_1,d3p_1,d3,d3p,&
                         d4_0,d4p_0,d4_1,d4p_1,d4,d4p)
real, intent(in) :: d1_0(:,:,:),d1_1(:,:,:),d1(:,:,:)
real, intent(in) :: d2_0(:,:,:),d2_1(:,:,:),d2(:,:,:)
real, intent(in) :: d3_0(:,:,:),d3_1(:,:,:),d3(:,:,:)
real, intent(in) :: d4_0(:,:,:),d4_1(:,:,:),d4(:,:,:)
real, pointer    :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
real, pointer    :: d2p_0(:,:,:),d2p_1(:,:,:),d2p(:,:,:)
real, pointer    :: d3p_0(:,:,:),d3p_1(:,:,:),d3p(:,:,:)
real, pointer    :: d4p_0(:,:,:),d4p_1(:,:,:),d4p(:,:,:)
call wf3d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
call wf3d_assign_1(d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
call wf3d_assign_1(d3_0,d3p_0,d3_1,d3p_1,d3,d3p)
call wf3d_assign_1(d4_0,d4p_0,d4_1,d4p_1,d4,d4p)
end subroutine wf3d_assign_4

!====================================================================

subroutine wf3d_assign_5(d1_0,d1p_0,d1_1,d1p_1,d1,d1p,&
                         d2_0,d2p_0,d2_1,d2p_1,d2,d2p,&
                         d3_0,d3p_0,d3_1,d3p_1,d3,d3p,&
                         d4_0,d4p_0,d4_1,d4p_1,d4,d4p,&
                         d5_0,d5p_0,d5_1,d5p_1,d5,d5p)
real, intent(in) :: d1_0(:,:,:),d1_1(:,:,:),d1(:,:,:)
real, intent(in) :: d2_0(:,:,:),d2_1(:,:,:),d2(:,:,:)
real, intent(in) :: d3_0(:,:,:),d3_1(:,:,:),d3(:,:,:)
real, intent(in) :: d4_0(:,:,:),d4_1(:,:,:),d4(:,:,:)
real, intent(in) :: d5_0(:,:,:),d5_1(:,:,:),d5(:,:,:)
real, pointer    :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
real, pointer    :: d2p_0(:,:,:),d2p_1(:,:,:),d2p(:,:,:)
real, pointer    :: d3p_0(:,:,:),d3p_1(:,:,:),d3p(:,:,:)
real, pointer    :: d4p_0(:,:,:),d4p_1(:,:,:),d4p(:,:,:)
real, pointer    :: d5p_0(:,:,:),d5p_1(:,:,:),d5p(:,:,:)
call wf3d_assign_1(d1_0,d1p_0,d1_1,d1p_1,d1,d1p)
call wf3d_assign_1(d2_0,d2p_0,d2_1,d2p_1,d2,d2p)
call wf3d_assign_1(d3_0,d3p_0,d3_1,d3p_1,d3,d3p)
call wf3d_assign_1(d4_0,d4p_0,d4_1,d4p_1,d4,d4p)
call wf3d_assign_1(d5_0,d5p_0,d5_1,d5p_1,d5,d5p)
end subroutine wf3d_assign_5

!====================================================================

subroutine wf2d_refresh_1(d1p_0,d1p_1,d1p)
real,pointer :: d1p_0(:,:),d1p_1(:,:),d1p(:,:),d1p_t(:,:)
d1p_t=>d1p_0;   d1p_0=>d1p_1;   d1p_1=>d1p;   d1p=>d1p_t
end subroutine wf2d_refresh_1

!====================================================================

subroutine wf2d_refresh_2(d1p_0,d1p_1,d1p,d2p_0,d2p_1,d2p)
real,pointer :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
real,pointer :: d2p_0(:,:),d2p_1(:,:),d2p(:,:)
call wf2d_refresh_1(d1p_0,d1p_1,d1p)
call wf2d_refresh_1(d2p_0,d2p_1,d2p)
end subroutine wf2d_refresh_2

!====================================================================

subroutine wf2d_refresh_3(d1p_0,d1p_1,d1p,d2p_0,d2p_1,d2p,d3p_0,d3p_1,d3p)
real,pointer :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
real,pointer :: d2p_0(:,:),d2p_1(:,:),d2p(:,:)
real,pointer :: d3p_0(:,:),d3p_1(:,:),d3p(:,:)
call wf2d_refresh_1(d1p_0,d1p_1,d1p)
call wf2d_refresh_1(d2p_0,d2p_1,d2p)
call wf2d_refresh_1(d3p_0,d3p_1,d3p)
end subroutine wf2d_refresh_3

!====================================================================

subroutine wf2d_refresh_4(d1p_0,d1p_1,d1p,d2p_0,d2p_1,d2p,&
   d3p_0,d3p_1,d3p,d4p_0,d4p_1,d4p)
real,pointer :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
real,pointer :: d2p_0(:,:),d2p_1(:,:),d2p(:,:)
real,pointer :: d3p_0(:,:),d3p_1(:,:),d3p(:,:)
real,pointer :: d4p_0(:,:),d4p_1(:,:),d4p(:,:)
call wf2d_refresh_1(d1p_0,d1p_1,d1p)
call wf2d_refresh_1(d2p_0,d2p_1,d2p)
call wf2d_refresh_1(d3p_0,d3p_1,d3p)
call wf2d_refresh_1(d4p_0,d4p_1,d4p)
end subroutine wf2d_refresh_4

!====================================================================

subroutine wf2d_refresh_5(d1p_0,d1p_1,d1p,d2p_0,d2p_1,d2p,&
   d3p_0,d3p_1,d3p,d4p_0,d4p_1,d4p,d5p_0,d5p_1,d5p)
real,pointer :: d1p_0(:,:),d1p_1(:,:),d1p(:,:)
real,pointer :: d2p_0(:,:),d2p_1(:,:),d2p(:,:)
real,pointer :: d3p_0(:,:),d3p_1(:,:),d3p(:,:)
real,pointer :: d4p_0(:,:),d4p_1(:,:),d4p(:,:)
real,pointer :: d5p_0(:,:),d5p_1(:,:),d5p(:,:)
call wf2d_refresh_1(d1p_0,d1p_1,d1p)
call wf2d_refresh_1(d2p_0,d2p_1,d2p)
call wf2d_refresh_1(d3p_0,d3p_1,d3p)
call wf2d_refresh_1(d4p_0,d4p_1,d4p)
call wf2d_refresh_1(d5p_0,d5p_1,d5p)
end subroutine wf2d_refresh_5

!====================================================================

subroutine wf3d_refresh_1(d1p_0,d1p_1,d1p)
real,pointer :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:),d1p_t(:,:,:)
d1p_t=>d1p_0;   d1p_0=>d1p_1;   d1p_1=>d1p;   d1p=>d1p_t
end subroutine wf3d_refresh_1

!====================================================================

subroutine wf3d_refresh_2(d1p_0,d1p_1,d1p,d2p_0,d2p_1,d2p)
real, pointer :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
real, pointer :: d2p_0(:,:,:),d2p_1(:,:,:),d2p(:,:,:)
call wf3d_refresh_1(d1p_0,d1p_1,d1p)
call wf3d_refresh_1(d2p_0,d2p_1,d2p)
end subroutine wf3d_refresh_2

!====================================================================

subroutine wf3d_refresh_3(d1p_0,d1p_1,d1p,d2p_0,d2p_1,d2p,d3p_0,d3p_1,d3p)
real, pointer :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
real, pointer :: d2p_0(:,:,:),d2p_1(:,:,:),d2p(:,:,:)
real, pointer :: d3p_0(:,:,:),d3p_1(:,:,:),d3p(:,:,:)
call wf3d_refresh_1(d1p_0,d1p_1,d1p)
call wf3d_refresh_1(d2p_0,d2p_1,d2p)
call wf3d_refresh_1(d3p_0,d3p_1,d3p)
end subroutine wf3d_refresh_3

!====================================================================

subroutine wf3d_refresh_4(d1p_0,d1p_1,d1p,d2p_0,d2p_1,d2p,&
   d3p_0,d3p_1,d3p,d4p_0,d4p_1,d4p)
real, pointer :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
real, pointer :: d2p_0(:,:,:),d2p_1(:,:,:),d2p(:,:,:)
real, pointer :: d3p_0(:,:,:),d3p_1(:,:,:),d3p(:,:,:)
real, pointer :: d4p_0(:,:,:),d4p_1(:,:,:),d4p(:,:,:)
call wf3d_refresh_1(d1p_0,d1p_1,d1p)
call wf3d_refresh_1(d2p_0,d2p_1,d2p)
call wf3d_refresh_1(d3p_0,d3p_1,d3p)
call wf3d_refresh_1(d4p_0,d4p_1,d4p)
end subroutine wf3d_refresh_4

!====================================================================

subroutine wf3d_refresh_5(d1p_0,d1p_1,d1p,d2p_0,d2p_1,d2p,&
   d3p_0,d3p_1,d3p,d4p_0,d4p_1,d4p,d5p_0,d5p_1,d5p)
real, pointer :: d1p_0(:,:,:),d1p_1(:,:,:),d1p(:,:,:)
real, pointer :: d2p_0(:,:,:),d2p_1(:,:,:),d2p(:,:,:)
real, pointer :: d3p_0(:,:,:),d3p_1(:,:,:),d3p(:,:,:)
real, pointer :: d4p_0(:,:,:),d4p_1(:,:,:),d4p(:,:,:)
real, pointer :: d5p_0(:,:,:),d5p_1(:,:,:),d5p(:,:,:)
call wf3d_refresh_1(d1p_0,d1p_1,d1p)
call wf3d_refresh_1(d2p_0,d2p_1,d2p)
call wf3d_refresh_1(d3p_0,d3p_1,d3p)
call wf3d_refresh_1(d4p_0,d4p_1,d4p)
call wf3d_refresh_1(d5p_0,d5p_1,d5p)
end subroutine wf3d_refresh_5

!====================================================================

subroutine snapshot2d_wi_bc(field_name,wf,it,nz,nx)
use module_string
use module_io
character(len=*), intent(in) :: field_name
integer,   intent(in) :: it, nz, nx
real,      intent(in) :: wf(:,:)
call filename(output,field_name,it,".bin")
call write_binfile(output,wf,nz,nx)
end subroutine snapshot2d_wi_bc

!====================================================================

subroutine snapshot2d_wo_bc(field_name,wf,it,nz,nx,nbc)
use module_string
use module_io
character(len=*), intent(in) :: field_name
integer,   intent(in) :: it, nz, nx, nbc
real,      intent(in) :: wf(:,:)
call filename(output,field_name,it,".bin")
call write_binfile(output,wf(nbc+1:nz-nbc,nbc+1:nx-nbc))
end subroutine snapshot2d_wo_bc

!====================================================================

subroutine snapshot3d_wi_bc(field_name,wf,it,nz,ny,nx)
use module_string
use module_io
character(len=*), intent(in) :: field_name
integer,   intent(in) :: it, nz, ny, nx
real,      intent(in) :: wf(:,:,:)
call filename(output,field_name,it,".bin")
call write_binfile(output,wf,nz,ny,nx)
end subroutine snapshot3d_wi_bc

!====================================================================

subroutine snapshot3d_wo_bc(field_name,wf,it,nz,ny,nx,nbc)
use module_string
use module_io
character(len=*), intent(in) :: field_name
integer,   intent(in) :: it, nz, ny, nx, nbc
real,      intent(in) :: wf(:,:,:)
call filename(output,field_name,it,".bin")
call write_binfile(output,wf(nbc+1:nz-nbc,nbc+1:ny-nbc,nbc+1:nx-nbc))
end subroutine snapshot3d_wo_bc

!====================================================================

real function average_real2(d1,d2)
real, intent(in) :: d1, d2
average_real2=(d1+d2)/2.0
end function average_real2

!====================================================================

real function average_real3(d1,d2,d3)
real, intent(in) :: d1, d2, d3
average_real3=(d1+d2+d3)/3.0
end function average_real3

!====================================================================

real function average_real4(d1,d2,d3,d4)
real, intent(in) :: d1, d2, d3, d4
average_real4=(d1+d2+d3+d4)/4.0
end function average_real4

!====================================================================

real function average_real5(d1,d2,d3,d4,d5)
real, intent(in) :: d1, d2, d3, d4, d5
average_real5=(d1+d2+d3+d4+d5)/5.0
end function average_real5

!====================================================================

real function average_real8(d1,d2,d3,d4,d5,d6,d7,d8)
real, intent(in) :: d1, d2, d3, d4, d5, d6, d7, d8
average_real8=(d1+d2+d3+d4+d5+d6+d7+d8)/8.0
end function average_real8

!====================================================================

real function average_1d_real(d)
real, intent(in) :: d(:)
average_1d_real=sum(d)/real(size(d,1))
end function average_1d_real

!====================================================================

real function average_2d_real(d)
real, intent(in) :: d(:,:)
average_2d_real=sum(d)/real(size(d,1)*size(d,2))
end function average_2d_real

!====================================================================

real function average_3d_real(d)
real, intent(in) :: d(:,:,:)
average_3d_real=sum(d)/real(size(d,1)*size(d,2)*size(d,3))
end function average_3d_real

!====================================================================

real function average_4d_real(d)
real, intent(in) :: d(:,:,:,:)
average_4d_real=sum(d)/real(size(d,1)*size(d,2)*size(d,3)*size(d,4))
end function average_4d_real

!====================================================================

real function accurate_value_real_1d(d_in,x1,dx)
real,    intent(in) :: d_in(:),x1,dx
integer :: i0,i2
real    :: x0, x2
i0=floor(x1/dx)+1;   i2=i0+1
x0=1.0-(x1-(i0-1)*dx)/dx;   x2=1.0-x0
accurate_value_real_1d=d_in(i0)*x0+d_in(i2)*x2
end function accurate_value_real_1d

real function accurate_value_real_2d(d_in,z1,x1,dz,dx)
real,    intent(in) :: d_in(:,:),z1,x1,dz,dx
integer :: i0, i2, j0, j2
real    :: z0, z2, x0, x2
i0=floor(z1/dz)+1;   i2=i0+1
j0=floor(x1/dx)+1;   j2=j0+1
x0=1.0-(x1-(j0-1)*dx)/dx;   x2=1.0-x0
z0=1.0-(z1-(i0-1)*dz)/dz;   z2=1.0-z0
accurate_value_real_2d=d_in(i0,j0)*z0*x0+d_in(i0,j2)*z0*x2 &
                      +d_in(i2,j0)*z2*x0+d_in(i2,j2)*z2*x2
end function accurate_value_real_2d

real function accurate_value_real_3d(d_in,z1,y1,x1,dz,dy,dx)
real,    intent(in) :: d_in(:,:,:),z1,y1,x1,dz,dy,dx
integer :: i0, i2, j0, j2, k0, k2
real    :: z0, z2, y0, y2, x0, x2
i0=floor(z1/dz)+1;   i2=i0+1
j0=floor(y1/dy)+1;   j2=j0+1
k0=floor(x1/dx)+1;   k2=k0+1
z0=1.0-(z1-(i0-1)*dz)/dz;   z2=1.0-z0
y0=1.0-(y1-(j0-1)*dy)/dy;   y2=1.0-y0
x0=1.0-(x1-(k0-1)*dx)/dx;   x2=1.0-x0
accurate_value_real_3d=d_in(i0,j0,k0)*z0*y0*x0+d_in(i0,j0,k2)*z0*y0*x2 &
                      +d_in(i2,j0,k0)*z2*y0*x0+d_in(i2,j0,k2)*z2*y0*x2 &
                      +d_in(i0,j2,k0)*z0*y2*x0+d_in(i0,j2,k2)*z0*y2*x2 &
                      +d_in(i2,j2,k0)*z2*y2*x0+d_in(i2,j2,k2)*z2*y2*x2
end function accurate_value_real_3d

subroutine output_illum_2d(isIllum,isSaveIllum,pre_illum_out,proid,&
   illum,ix1,ix2,nz,nx)
logical, intent(in) :: isIllum,isSaveIllum
character(len=*), intent(in) :: pre_illum_out
integer, intent(in) :: proid,ix1,ix2,nz,nx
real,    intent(in) :: illum(:,:)
character(len=slen) :: output_temp1,output_temp2,output
if (isIllum) then
   if (isSaveIllum) then
      call filename(output_temp1,pre_illum_out,proid,"_ix_")
      call filename(output_temp2,output_temp1,ix1,"_")
      call filename(output, output_temp2,ix2,".bin")
      call write_binfile(output,illum,nz,nx)
   endif
endif
end subroutine output_illum_2d

subroutine output_illum_3d(isIllum,isSaveIllum,pre_illum_out,proid,&
   illum,ix1,ix2,iy1,iy2,nz,ny,nx)
logical, intent(in) :: isIllum,isSaveIllum
character(len=*), intent(in) :: pre_illum_out
integer, intent(in) :: proid,ix1,ix2,iy1,iy2,nz,ny,nx
real,    intent(in) :: illum(:,:,:)
character(len=slen) :: output_temp1,output_temp2,output
if (isIllum) then
   if (isSaveIllum) then
      call filename(output_temp1,pre_illum_out,proid,"_ix_")
      call filename(output_temp2,output_temp1,ix1,"_")
      call filename(output_temp1,output_temp2,ix2,"_iy_")
      call filename(output_temp2,output_temp1,iy1,"_iy_")
      call filename(output,output_temp2,iy2,".bin")
      call write_binfile(output,illum,nz,ny,nx)
   endif
endif
end subroutine output_illum_3d

subroutine output_prestack_image_2d(isPreStackImg,pre_illum_img_out,&
   proid,img,ix1,ix2,nz,nx)
logical, intent(in) :: isPreStackImg
character(len=*), intent(in) :: pre_illum_img_out
real, intent(in) :: img(:,:)
integer, intent(in) :: proid,ix1,ix2,nz,nx
character(len=slen) :: output_temp1,output_temp2,output
if (isPreStackImg) then
   call filename(output,pre_illum_img_out,proid,".bin")
!   call filename(output_temp2,output_temp1,ix1,"_")
!   call filename(output,output_temp2,ix2,".bin")
   call write_binfile(output,img,nz,nx)
endif
end subroutine output_prestack_image_2d

subroutine output_prestack_image_3d(isPreStackImg,pre_illum_img_out,&
   proid,img,ix1,ix2,iy1,iy2,nz,ny,nx)
logical, intent(in) :: isPreStackImg
character(len=*), intent(in) :: pre_illum_img_out
real, intent(in) :: img(:,:,:)
integer, intent(in) :: proid,ix1,ix2,iy1,iy2,nz,ny,nx
character(len=slen) :: output_temp1,output_temp2,output
if (isPreStackImg) then
   call filename(output_temp1,pre_illum_img_out,proid,"_ix_")
   call filename(output_temp2,output_temp1,ix1,"_")
   call filename(output_temp1,output_temp2,ix2,"_iy_")
   call filename(output_temp2,output_temp1,iy1,"_iy_")
   call filename(output,output_temp2,iy2,".bin")
   call write_binfile(output,img,nz,ny,nx)
endif
end subroutine output_prestack_image_3d

subroutine log_end(logid)
integer,intent(in) :: logid
call date_and_time(date,time_now)
write(logid,*)"Date ",date
write(logid,*)"Time ",time_now
write(logid,*)"All finished!"
close(logid)
end subroutine log_end

integer function umod(n,nmod)
integer, intent(in) :: n, nmod
if (mod(n,nmod)==0) then
   umod=int(n/nmod)
else
   umod=int(n/nmod)+1
endif
end function umod

integer function lmod(n,nmod)
integer, intent(in) :: n, nmod
lmod=int(n/nmod)
end function lmod




!===========================================================================
! The following subroutines are added by Bowen Guo
subroutine shift_forward_other(trace_in,trace_out,tt)
real,intent(in)    :: trace_in(:)
integer,intent(in) :: tt
real,intent(out)   :: trace_out(:)
integer            :: it,nt
nt=size(trace_in,1)
do it=1+tt,nt
   trace_out(it-tt)=trace_in(it)
enddo
trace_out(nt-tt+1:nt)=0.0
end subroutine shift_forward_other
!===========================================================================

subroutine shift_forward_self(trace,tt)
real,intent(inout)   :: trace(:)
integer,intent(in)   :: tt

integer              :: it,nt
real,allocatable     :: temp(:)

nt=size(trace,1)
allocate(temp(nt))
do it=1+tt,nt
   temp(it-tt)=trace(it)
enddo
temp(nt-tt+1:nt)=0.0

trace=temp
end subroutine shift_forward_self
!===========================================================================


subroutine shift_backward_other(trace_in,trace_out,tt)
real,intent(in)     :: trace_in(:)
integer,intent(in)  :: tt
real,intent(out)    :: trace_out(:)
integer             :: it,nt

nt=size(trace_in,1)
do it=1,nt-tt
   trace_out(it+tt)=trace_in(it)
enddo
trace_out(1:tt)=0.0

end subroutine shift_backward_other
!=============================================================================
subroutine shift_backward_self(trace,tt)
real,intent(inout)   :: trace(:)
integer,intent(in)   :: tt

integer              :: it,nt
real,allocatable     :: temp(:)

nt=size(trace,1)
allocate(temp(nt))
do it=1,nt-tt
   temp(it+tt)=trace(it)
enddo
temp(1:tt)=0.0

trace=temp

end subroutine shift_backward_self
!===========================================================================
subroutine mute_after(trace,mute_time,hann)
! Not cut suddenly, but to use a hanning window
real,intent(in)      :: hann(:)
real,intent(inout)   :: trace(:)
integer,intent(in)   :: mute_time

integer              :: hann_length,it

hann_length=size(hann)

do it=1,hann_length
   trace(mute_time+it)=trace(mute_time+it)*hann(it)
enddo

trace(mute_time+hann_length+1:size(trace))=0.0


end subroutine mute_after
!===============================================================================
subroutine mute_before(trace,mute_time,hann)
! Not cut suddenly, but to use a hanning window
real,intent(in)      :: hann(:)
real,intent(inout)   :: trace(:)
integer,intent(in)   :: mute_time

integer              :: hann_length,it

hann_length=size(hann)
do it=1,hann_length
   trace(mute_time-it)=trace(mute_time-it)*hann(it)
enddo
trace(1:mute_time-hann_length-1)=0.0


end subroutine mute_before

!==============================================================================
subroutine pad_time_gather(seis_in,seis_out,nt,ntt)
real,intent(in)    :: seis_in(:,:)
integer,intent(in) :: nt,ntt
real,intent(out)   :: seis_out(:,:)
integer            :: ig,ng

ng=size(seis_in,2)
do ig=1,ng
   seis_out(nt:ntt,ig)=seis_in(:,ig)
enddo

end subroutine pad_time_gather
!==================================================================================

subroutine unpad_time_trace(trace_in,trace_out,nt,ntt)
real,intent(in)    :: trace_in(:)
integer,intent(in) :: nt,ntt
real,intent(out)   :: trace_out(:)

trace_out=trace_in(1:nt)
end subroutine unpad_time_trace



!==================================================================================

subroutine adjust_ttt(ttt_in,ttt_out,dt,nz,nx)
real,intent(in)     :: ttt_in(:,:),dt
integer,intent(in)  :: nz,nx
integer,intent(out) :: ttt_out(:,:)
integer             :: is,iz,ix
do ix=1,nx
   do iz=1,nz
      ttt_out(iz,ix)=floor(ttt_in(iz,ix)/dt)+1
   enddo
enddo


end subroutine adjust_ttt


!========================================================================
integer function find_index(v_in,aim)
real,intent(in)      :: v_in(:)
real,intent(in)      :: aim

integer              :: i
do i=1,size(v_in,1)
   if (v_in(i).eq.aim) then
      find_index=i
      return
   endif
enddo
end function find_index
!=================================================================================
! Add by Bowen Guo @2015-02-23  TO control the matrix between max and min
subroutine control_value(m,vmax,vmin)
real,intent(inout)    :: m(:,:)
real,intent(in)       :: vmax,vmin
integer               :: nz,nx,ix,iz

nz=size(m,1);nx=size(m,2)
do ix=1,nx
   do iz=1,nz
      if (m(iz,ix)>vmax) then
         m(iz,ix)=vmax
      endif
      if (m(iz,ix)<vmin) then
         m(iz,ix)=vmin
      endif
   enddo
enddo

end subroutine control_value
!==================================================================================
!subroutine FourierShift(trace_in,trace_out,shift)
!real,intent(in)        :: shift
!real,intent(in)        :: trace_in(:)
!real,intent(out)       :: trace_out(:)
!integer                :: nt,it
!complex,allocatable    :: trace_in_fft(:),trace_out_fft(:),temp(:),temp_fft(:)
!
!nt=size(trace_in)
!!do it=1:int(nt/2)
!!   temp(it)=it-1
!!enddo
!!do it=int(nt/2)+1:nt
!!   temp(it)=
!temp=[0:int(nt/2)-1,int((-nt)/2):-1]/nt
!
!call allocate_and_initial(trace_in_fft,nt)
!call allocate_and_initial(temp,nt)
!call allocate_and_initial(temp_fft,nt)
!call fft_wavelet(trace_in,trace_in_fft,nt,nt)
!
!temp_fft=exp(-1*J*2*PI*shift*temp)
!if (mod(nt,2).eq.0) then
!   temp_fft(nt/2+1)=real(temp_fft(nt/2+1))
!endif
!
!trace_out_fft=trace_in_fft*temp_fft
!call ifft_wavelet(trace_out_fft,trace_out,nt)
!
!if (shift<0) then
!   trace_out(1:shift)=0.0
!else
!   trace_out(nt-shift+1:nt)=0.0
!endif
!
!
!end subroutine FourierShift



end module module_utility
