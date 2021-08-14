module module_source
use module_global, only : PI, I4
use module_datatype
use module_sf90_mpi
implicit none

contains

!=====================================================================
subroutine ricker_bak(s, nt, f, dt)
! Creating a Ricker source wavelet with the same length as the traces

real,    parameter   :: timeshift = 1.41421356237309504880
integer, intent(in)  :: nt
real,    intent(in)  :: f, dt
real,    intent(out) :: s(:)
real,    allocatable :: s1(:)
integer              :: it, np
real                 :: pi2, const, b, tshift, tim1, tim2, u, amp, smax

np = floor(timeshift/(f*dt)+0.5)*2+1
allocate(s1(np))
tshift = timeshift/f
!tshift = timeshift/f + 29.0*dt ! For BP data, fpeak = 3.5 Hz
pi2 = sqrt(PI)/2.0
b = sqrt(6.0)/(PI*f)
const = 2.0*sqrt(6.0)/b

smax = 0.0
do it=1,np
  tim1 = real(it-1)*dt
  tim2 = tim1 - tshift
  u = const*tim2
  amp = ((u*u)/4.0-0.5)*pi2*(exp(-u*u/4.0))
  s1(it) = -amp
  if (smax .lt. abs(amp)) then
    smax = abs(amp)
  endif
enddo

smax = smax*2
do it=1,np
  s1(it) = s1(it)/smax
enddo

do it=1,nt
   if(it.le.np) then
      s(it)=s1(it)
   else
      s(it)=0.0
   endif
enddo

s=s*1.0e-3
end subroutine ricker_bak


subroutine ricker(s)
! Creaing a Ricker source wavelet 
type(source), intent(inout) :: s
real,    parameter   :: timeshift = 1.41421356237309504880

integer :: it
real    :: pi2, const, b, tshift, tim1, tim2, u, amp, smax

tshift = timeshift/s%freq
pi2 = sqrt(PI)/2.0
b = sqrt(6.0)/(PI*s%freq)
const = 2.0*sqrt(6.0)/b

smax = 0.0
do it=1,s%nw
  tim1 = real(it-1)*s%dt
  tim2 = tim1 - tshift
  u = const*tim2
  amp = ((u*u)/4.0-0.5)*pi2*(exp(-u*u/4.0))
  s%sou(it) = -amp
  if (smax .lt. abs(amp)) then
    smax = abs(amp)
  endif
enddo

smax = smax*2
do it=1,s%nw
  s%sou(it) = s%sou(it)/smax
enddo

end subroutine ricker


subroutine ricker_fc(s,nt,dt,fc)
! Creaing a Ricker source wavelet 
real,    parameter   :: timeshift = 1.41421356237309504880
integer :: it,nt
real    :: pi2, const, b, tshift, tim1, tim2, u, amp, smax,fc,dt
real,allocatable :: s(:)

allocate(s(nt))

tshift = timeshift/fc
pi2 = sqrt(PI)/2.0
b = sqrt(6.0)/(PI*fc)
const = 2.0*sqrt(6.0)/b

smax = 0.0
do it=1,nt
  tim1 = real(it-1)*dt
  tim2 = tim1 - tshift
  u = const*tim2
  amp = ((u*u)/4.0-0.5)*pi2*(exp(-u*u/4.0))
  s(it) = -amp
  if (smax .lt. abs(amp)) then
    smax = abs(amp)
  endif
enddo

smax = smax*2
do it=1,nt
  s(it) = s(it)/smax
enddo

end subroutine ricker_fc

!====================================================================

subroutine wavelet_expand(s_in,s_out,nt)
real, intent(in) :: s_in(:)
real, intent(out):: s_out(:)
integer, intent(in):: nt
integer :: nsource
! setup source
s_out=0.0
nsource=size(s_in,1)
if (nsource .le. nt) then
   s_out(1:nsource)=s_in
else
   s_out=s_in(1:nsource)
endif
end subroutine wavelet_expand
!-----------------------------------------------------------------------------------------

subroutine read_sourcefile_mpi(source_file,s,nt_source)
character(len=*),intent(in)  :: source_file
real,            intent(out) :: s(:)
integer,         intent(in)  :: nt_source
integer                      :: nt
s=0.0
nt=size(s,1)
if (rank.eq.0) then
   if (nt_source.le.nt) then
      open(10,file=source_file,access="direct",recl=I4*nt_source,status='old')
      read(10,rec=1) s(1:nt_source)
      close(10)
   else
      open(10,file=source_file,access="direct",recl=I4*nt,status='old')
      read(10,rec=1) s(1:nt)
      close(10)
  endif
endif
call MPI_BCAST(s,nt,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine read_sourcefile_mpi
!===========================================================================================
! Added by Bowen Guo

subroutine deter_source_shift(s,source_shift)
real,intent(in)         :: s(:)
integer,intent(out)     :: source_shift

integer                 :: i,nt
real                    :: max_value

nt=size(s,1)
max_value=s(1)
source_shift=1

do i=2,nt
   if (max_value<s(i)) then 
      max_value=s(i)
      source_shift=i
   endif
enddo
end subroutine deter_source_shift
!===========================================================================================


















end module module_source
