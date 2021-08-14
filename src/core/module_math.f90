module module_math
use module_global, only : PI

implicit none

interface signn
   module procedure signn_i1
   module procedure signn_i2
   module procedure signn_i4
   module procedure signn_i8
   module procedure signn_r4
   module procedure signn_r8
end interface signn

interface local_min
   module procedure local_min_r4_1d
   !module procedure local_min_r4_1d
end interface local_min

interface local_max
   module procedure local_max_r4_1d
   !module procedure local_max_r4_1d
end interface local_max

interface integrate                        ! modified by Bowen Guo
   module procedure inde_integrate_r4_1d       ! modified by Bowen Guo
   module procedure inde_integrate_r8_1d       ! modified by Bowen Guo
end interface integrate               ! modified by Bowen Guo

interface data_type_transfer          ! modified by Bowen Guo
   module procedure real4to8_2d          ! modified by Bowen Guo
end interface data_type_transfer      ! modified by Bowen Guo

contains

subroutine real4to8_2d(d_4,d_8)         ! modified by Bowen Guo
real(kind=4),intent(in)  :: d_4(:,:)    ! modified by Bowen Guo
real(kind=8),intent(out) :: d_8(:,:)    ! modified by Bowen Guo
d_8=d_4                                 ! modified by Bowen Guo
end subroutine real4to8_2d                 ! modified by Bowen Guo




subroutine  inde_integrate_r4_1d(d,n)               ! modified by Bowen Guo
real,intent(inout)  :: d(:)
integer,intent(in)  :: n                ! modified by Bowen Guo
integer             :: i,j
real                :: d_int(n)             ! modified by Bowen Guo
d_int=0.0                                   ! modified by Bowen Guo
do i=1,n                                    ! modified by Bowen Guo
   do j=1,i                                 ! modified by Bowen Guo
       d_int(i)=d_int(i)+d(j)               ! modified by Bowen Guo
   enddo                                    ! modified by Bowen Guo
enddo                                       ! modified by Bowen GUo
d=d_int                                     ! modified by Bowen Guo
end subroutine  inde_integrate_r4_1d              ! modified by Bowen Guo

subroutine  inde_integrate_r8_1d(d,n)               ! modified by Bowen Guo
real(kind=8),intent(inout)  :: d(:)
integer,intent(in)          :: n                ! modified by Bowen Guo
integer                     :: i,j
real(kind=8)                :: d_int(n)             ! modified by Bowen Guo
d_int=0.0                                   ! modified by Bowen Guo
do i=1,n                                    ! modified by Bowen Guo
   do j=1,i                                 ! modified by Bowen Guo
       d_int(i)=d_int(i)+d(j)               ! modified by Bowen Guo
   enddo                                    ! modified by Bowen Guo
enddo                                       ! modified by Bowen GUo
d=d_int                                     ! modified by Bowen Guo
end subroutine  inde_integrate_r8_1d              ! modified by Bowen Guo

real function signn_i1(d)
integer(kind=1), intent(in) :: d
if (d<0) then;signn_i1=-1.0;else;signn_i1=1.0;endif
end function signn_i1

real function signn_i2(d)
integer(kind=2), intent(in) :: d
if (d<0) then;signn_i2=-1.0;else;signn_i2=1.0;endif
end function signn_i2

real function signn_i4(d)
integer(kind=4), intent(in) :: d
if (d<0) then;signn_i4=-1.0;else;signn_i4=1.0;endif
end function signn_i4

real function signn_i8(d)
integer(kind=8), intent(in) :: d
if (d<0) then;signn_i8=-1.0;else;signn_i8=1.0;endif
end function signn_i8

real function signn_r4(d)
real(kind=4), intent(in) :: d
if (d<0) then;signn_r4=-1.0;else;signn_r4=1.0;endif
end function signn_r4

real function signn_r8(d)
real(kind=8), intent(in) :: d
if (d<0) then;signn_r8=-1.0;else;signn_r8=1.0;endif
end function signn_r8

function hanning(n)
real :: hanning(n)
integer,intent(in) :: n
integer            :: i
real               :: m
m=2.0*n+1.0
do i=1,n
  hanning(i)=0.5*(1-cos(2*PI*i/m))
enddo
end function hanning

function reverse_hanning(n)
real :: reverse_hanning(n)
integer,intent(in) :: n
integer            :: i
real               :: m
real,allocatable   :: h_temp(:)
allocate(h_temp(n))
m=2.0*n+1.0
do i=1,n
   h_temp(i)=0.5*(1-cos(2*PI*i/m))
enddo
do i=1,n
   reverse_hanning(i)=h_temp(n-i+1)
enddo
deallocate(h_temp)
end function reverse_hanning

FUNCTION ran1(idum)
INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
REAL ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
if (idum.le.0.or.iy.eq.0) then
   idum=max(-idum,1)
   do j=NTAB+8,1,-1
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      if (j.le.NTAB) iv(j)=idum
   enddo
   iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
END FUNCTION ran1

FUNCTION gasdev(idum)
INTEGER idum
REAL gasdev
!    USES ran1
INTEGER iset
REAL fac,gset,rsq,v1,v2!,ran1
SAVE iset,gset
DATA iset/0/
if (iset.eq.0) then
1  v1=2.*ran1(idum)-1.
   v2=2.*ran1(idum)-1.
   rsq=v1**2+v2**2
   if (rsq.ge.1..or.rsq.eq.0.) goto 1
   fac=sqrt(-2.*log(rsq)/rsq)
   gset=v1*fac
   gasdev=v2*fac
   iset=1
else
   gasdev=gset
   iset=0
endif
return
END FUNCTION gasdev

subroutine local_min_r4_1d(d,n,min_number,min_loc)
real,    intent(in) :: d(:)
integer, intent(in) :: n
integer,intent(out) :: min_number, min_loc(:)
integer             :: i
min_number=0;   min_loc=0
! deal with the left end
if ( d(1)<=d(2) ) then
   min_number=min_number+1;   min_loc(min_number)=1
endif
do i=2,n-1
   if ( d(i)<= d(i-1) .and. d(i)<= d(i+1)) then
      min_number=min_number+1;   min_loc(min_number)=i
   endif
enddo

! deal with the right end
if (d(n)<=d(n-1) ) then
   min_number=min_number+1;   min_loc(min_number)=n
endif
end subroutine local_min_r4_1d

subroutine local_max_r4_1d(d,n,max_number,max_loc)
real,    intent(in) :: d(:)
integer, intent(in) :: n
integer,intent(out) :: max_number, max_loc(:)
integer            :: i
max_number=0;   max_loc=0
! deal with the left end
if ( d(1)>=d(2) ) then
   max_number=max_number+1;   max_loc(max_number)=1
endif
do i=2,n-1
   if ( d(i) >= d(i-1) .and. d(i) >= d(i+1)) then
      max_number=max_number+1;   max_loc(max_number)=i
   endif
enddo
! deal with the right end
if (d(n) >= d(n-1) ) then
   max_number=max_number+1;   max_loc(max_number)=n
endif
end subroutine local_max_r4_1d

!=================================================================================================
! Create Hanning window Written by Bowen Guo
subroutine hanning_window(period,hann)
integer, intent(in)      :: period
real, intent(out)        :: hann(:)

integer                  :: it

do it=1,int(period/2)
   hann(it)=0.5*(cos(2*PI/real(period)*real(it))+1)
enddo

end subroutine hanning_window

!================================================================================================












end module module_math
