module module_fftw
use module_sf90_mpi
use module_array
implicit none
include "fftw3.f"
integer(kind=8)   :: plan, plan1, plan2, plan3

contains

subroutine fft_wavelet(w,wf,nw,nt_conv)
integer, intent(in) :: nw,nt_conv
real,    intent(in) :: w(nw)
complex, intent(out):: wf(nt_conv)
real,   allocatable :: wp(:) ! w_padding
call allocate_and_initial(wp,nt_conv)
wp=0.0
wp(1:nw)=w
call sfftw_plan_dft_r2c_1d(plan,nt_conv,wp,wf,FFTW_ESTIMATE)
call sfftw_execute(plan,wp,wf)
call sfftw_destroy_plan(plan)
end subroutine fft_wavelet

! ifft_wavelet is written by Bowen Guo
!subroutine ifft_wavelet(wf,w,nt)
!integer, intent(in) :: nt
!complex,intent(in)  :: wf(:)
!real,intent(out)    :: w(:)

!call sfftw_plan_dft_c2r_1d(plan1,nt_conv,w,wf,FFTW_ESTIMATE)
!call sfftw_execute_dft_c2r(plan1,wf,w)
!call sfftw_destroy_plan(plan1)
!w=w(1:nt)/nt
!end subroutine ifft_wavelet


subroutine conv_wavelet(trace_in,trace_out,wf,nt,nt_conv)
integer, intent(in)        :: nt,nt_conv
complex, intent(in)        :: wf(nt_conv)
real,    intent(in)        :: trace_in(nt)
real,    intent(out)       :: trace_out(nt)
real                       :: trace_work(nt_conv)
complex                    :: trace_in_fft(nt_conv),trace_out_fft(nt_conv)
trace_work(1:nt)=trace_in
trace_work(nt+1:nt_conv)=0.0
call sfftw_plan_dft_r2c_1d(plan1,nt_conv,trace_work,trace_in_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_r2c(plan1,trace_work,trace_in_fft)
call sfftw_destroy_plan(plan1)
trace_out_fft=trace_in_fft*wf
call sfftw_plan_dft_c2r_1d(plan2,nt_conv,trace_work,trace_out_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_c2r(plan2,trace_out_fft,trace_work)
call sfftw_destroy_plan(plan2)
!trace_out=trace_work(1:nt)/nt_conv
trace_out=trace_work(1:nt)/nt_conv
end subroutine conv_wavelet

subroutine xcorr_wavelet(trace_in,trace_out,wf,nt,nt_conv)
integer, intent(in)    :: nt,nt_conv
complex, intent(in)    :: wf(nt_conv)
real,    intent(in)    :: trace_in(nt)
real,    intent(out)   :: trace_out(nt)
real                   :: trace_work(nt_conv)
complex                :: trace_in_fft(nt_conv),trace_out_fft(nt_conv)
trace_work(1:nt)=trace_in
trace_work(nt+1:nt_conv)=0.0
call sfftw_plan_dft_r2c_1d(plan1,nt_conv,trace_work,trace_in_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_r2c(plan1,trace_work,trace_in_fft)
call sfftw_destroy_plan(plan1)
trace_out_fft=trace_in_fft*conjg(wf)
call sfftw_plan_dft_c2r_1d(plan2,nt_conv,trace_out_fft,trace_work,FFTW_ESTIMATE)
call sfftw_execute_dft_c2r(plan2,trace_out_fft,trace_work)
call sfftw_destroy_plan(plan2)
trace_out=trace_work(1:nt)/nt_conv
end subroutine xcorr_wavelet

subroutine conv_gf(gf_s,gf_g,gf_sg,nt,nt_conv)
integer, intent(in)        :: nt,nt_conv
real,    intent(in)        :: gf_s(nt),gf_g(nt)
real,    intent(out)       :: gf_sg(nt)
real                       :: gf_s_work(nt_conv),gf_g_work(nt_conv),gf_sg_work(nt_conv)
complex                    :: gf_s_fft(nt_conv),gf_g_fft(nt_conv),gf_sg_fft(nt_conv)
gf_s_work(1:nt)=gf_s
gf_s_work(nt+1:nt_conv)=0.0
gf_g_work(1:nt)=gf_g
gf_g_work(nt+1:nt_conv)=0.0
call sfftw_plan_dft_r2c_1d(plan1,nt_conv,gf_s_work,gf_s_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_r2c(plan1,gf_s_work,gf_s_fft)
call sfftw_destroy_plan(plan1)
call sfftw_plan_dft_r2c_1d(plan2,nt_conv,gf_g_work,gf_g_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_r2c(plan2,gf_g_work,gf_g_fft)
call sfftw_destroy_plan(plan2)
gf_sg_fft=gf_s_fft*gf_g_fft
call sfftw_plan_dft_c2r_1d(plan3,nt_conv,gf_sg_work,gf_sg_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_c2r(plan3,gf_sg_fft,gf_sg_work)
call sfftw_destroy_plan(plan3)
gf_sg=gf_sg_work(1:nt)/nt_conv
end subroutine conv_gf
!========================================================================
! Added by Bowen Guo
subroutine xcorr_gf(gf_s,gf_g,gf_sg,nt,nt_conv)
integer, intent(in)        :: nt,nt_conv
real,    intent(in)        :: gf_s(nt),gf_g(nt)
real,    intent(out)       :: gf_sg(nt)
real                       :: gf_s_work(nt_conv),gf_g_work(nt_conv),gf_sg_work(nt_conv)
complex                    :: gf_s_fft(nt_conv),gf_g_fft(nt_conv),gf_sg_fft(nt_conv)
gf_s_work(1:nt)=gf_s
gf_s_work(nt+1:nt_conv)=0.0
gf_g_work(1:nt)=gf_g
gf_g_work(nt+1:nt_conv)=0.0
call sfftw_plan_dft_r2c_1d(plan1,nt_conv,gf_s_work,gf_s_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_r2c(plan1,gf_s_work,gf_s_fft)
call sfftw_destroy_plan(plan1)
call sfftw_plan_dft_r2c_1d(plan2,nt_conv,gf_g_work,gf_g_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_r2c(plan2,gf_g_work,gf_g_fft)
call sfftw_destroy_plan(plan2)
gf_sg_fft=gf_s_fft*conjg(gf_g_fft)
call sfftw_plan_dft_c2r_1d(plan3,nt_conv,gf_sg_work,gf_sg_fft,FFTW_ESTIMATE)
call sfftw_execute_dft_c2r(plan3,gf_sg_fft,gf_sg_work)
call sfftw_destroy_plan(plan3)
gf_sg=gf_sg_work(1:nt)/nt_conv
end subroutine xcorr_gf







end module module_fftw
