module module_io
use module_global, only : I4
use module_sf90_mpi
use module_math
!use module_utility
implicit none
integer, private :: n1_get,n2_get,n3_get,n4_get,i1, i2, i3, ii4, it, n
real,    private :: df_sign

interface read_asciifile
   module procedure read_asciifile_i4_1d
   module procedure read_asciifile_i4_2d
   module procedure read_asciifile_r4_1d
   module procedure read_asciifile_r4_2d
   module procedure read_asciifile_r8_1d
   module procedure read_asciifile_r8_2d
end interface read_asciifile

interface read_asciifile_mpi
   module procedure read_asciifile_mpi_i4_1d
   module procedure read_asciifile_mpi_i4_2d
   module procedure read_asciifile_mpi_r4_1d
   module procedure read_asciifile_mpi_r4_2d
   module procedure read_asciifile_mpi_r8_1d
   module procedure read_asciifile_mpi_r8_2d
end interface read_asciifile_mpi

interface write_asciifile
   module procedure write_asciifile_i4_1d
   module procedure write_asciifile_i4_2d
   module procedure write_asciifile_r4_1d
   module procedure write_asciifile_r4_2d
   module procedure write_asciifile_r8_1d
   module procedure write_asciifile_r8_2d
   module procedure write_asciifile_i4_1d_wo_dim
   module procedure write_asciifile_i4_2d_wo_dim
   module procedure write_asciifile_r4_1d_wo_dim
   module procedure write_asciifile_r4_2d_wo_dim
   module procedure write_asciifile_r8_1d_wo_dim
   module procedure write_asciifile_r8_2d_wo_dim
end interface write_asciifile

interface read_binfile
   module procedure read_binfile_i2_1d
   module procedure read_binfile_i2_2d
   module procedure read_binfile_i2_3d
   module procedure read_binfile_i2_4d
   module procedure read_binfile_i4_1d
   module procedure read_binfile_i4_2d
   module procedure read_binfile_i4_3d
   module procedure read_binfile_i4_4d
   module procedure read_binfile_r4_1d
   module procedure read_binfile_r4_2d
   module procedure read_binfile_r4_3d
   module procedure read_binfile_r4_4d
   module procedure read_binfile_r8_1d
   module procedure read_binfile_r8_2d
   module procedure read_binfile_r8_3d
   module procedure read_binfile_r8_4d
   module procedure read_binfile_dfor_i2tor4_1d  
   module procedure read_binfile_dfor_i2tor4_2d  
   module procedure read_binfile_dfor_i2tor4_3d  
   module procedure read_binfile_dfor_i2tor4_4d  
end interface read_binfile

interface read_binfile_mpi
   module procedure read_binfile_mpi_i2_1d
   module procedure read_binfile_mpi_i2_2d
   module procedure read_binfile_mpi_i2_3d
   module procedure read_binfile_mpi_i2_4d
   module procedure read_binfile_mpi_i4_1d
   module procedure read_binfile_mpi_i4_2d
   module procedure read_binfile_mpi_i4_3d
   module procedure read_binfile_mpi_i4_4d
   module procedure read_binfile_mpi_r4_1d
   module procedure read_binfile_mpi_r4_2d
   module procedure read_binfile_mpi_r4_3d
   module procedure read_binfile_mpi_r4_4d
   module procedure read_binfile_mpi_r8_1d
   module procedure read_binfile_mpi_r8_2d
   module procedure read_binfile_mpi_r8_3d
   module procedure read_binfile_mpi_r8_4d
end interface read_binfile_mpi

interface write_binfile
   module procedure write_binfile_i2_1d
   module procedure write_binfile_i2_1d_wo_dim
   module procedure write_binfile_i2_2d
   module procedure write_binfile_i2_2d_wo_dim
   module procedure write_binfile_i2_3d
   module procedure write_binfile_i2_3d_wo_dim
   module procedure write_binfile_i2_4d
   module procedure write_binfile_i2_4d_wo_dim
   module procedure write_binfile_i4_1d
   module procedure write_binfile_i4_1d_wo_dim
   module procedure write_binfile_i4_2d
   module procedure write_binfile_i4_2d_wo_dim
   module procedure write_binfile_i4_3d
   module procedure write_binfile_i4_3d_wo_dim
   module procedure write_binfile_i4_4d
   module procedure write_binfile_i4_4d_wo_dim
   module procedure write_binfile_r4_1d
   module procedure write_binfile_r4_1d_wo_dim
   module procedure write_binfile_r4_2d
   module procedure write_binfile_r4_2d_wo_dim
   module procedure write_binfile_r4_3d
   module procedure write_binfile_r4_3d_wo_dim
   module procedure write_binfile_r4_4d
   module procedure write_binfile_r4_4d_wo_dim
   module procedure write_binfile_r8_1d
   module procedure write_binfile_r8_1d_wo_dim
   module procedure write_binfile_r8_2d
   module procedure write_binfile_r8_2d_wo_dim
   module procedure write_binfile_r8_3d
   module procedure write_binfile_r8_3d_wo_dim
   module procedure write_binfile_r8_4d
   module procedure write_binfile_r8_4d_wo_dim
   module procedure write_binfile_dfor_r4toi2_1d  
   module procedure write_binfile_dfor_r4toi2_2d  
   module procedure write_binfile_dfor_r4toi2_3d  
   module procedure write_binfile_dfor_r4toi2_4d  
end interface write_binfile

interface read_binfile_delete
  module procedure read_binfile_delete_r4_2d
end interface read_binfile_delete

contains

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_i4_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
integer,          intent(out) :: d(:)
integer,          intent(in)  :: n1
open (11,file=fn,status="old",action="read")
do i1=1,n1;   read(11,*)d(i1);   enddo;close(11)
end subroutine read_asciifile_i4_1d

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_mpi_i4_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
integer,          intent(out) :: d(:)
integer,          intent(in)  :: n1
if (rank.eq.0) call read_asciifile_i4_1d(fn,d,n1)
call sf90_bcast(d,n1)
end subroutine read_asciifile_mpi_i4_1d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_i4_1d(fn,d,n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: d(:)
integer,          intent(in) :: n1
open (11,file=fn,action="write",form="formatted")
do i1=1,n1;   write(11,*)d(i1);   enddo;close(11)
end subroutine write_asciifile_i4_1d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_i4_1d_wo_dim(fn,d)
character(len=*), intent(in) :: fn
integer,          intent(in) :: d(:)
n1_get=size(d);open (11,file=fn,action="write",form="formatted")
do i1=1,n1_get;   write(11,*)d(i1);   enddo;   close(11)
end subroutine write_asciifile_i4_1d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_r4_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
real,             intent(out) :: d(:)
integer,          intent(in)  :: n1
open (11,file=fn,status="old",action="read")
do i1=1,n1;   read(11,*)d(i1);   enddo;close(11)
end subroutine read_asciifile_r4_1d

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_mpi_r4_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
real,             intent(out) :: d(:)
integer,          intent(in)  :: n1
if (rank.eq.0) call read_asciifile_r4_1d(fn,d,n1)
call sf90_bcast(d,n1)
end subroutine read_asciifile_mpi_r4_1d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_r4_1d(fn,d,n1)
character(len=*), intent(in) :: fn
real,             intent(in) :: d(:)
integer,          intent(in) :: n1
open (11,file=fn,action="write",form="formatted")
do i1=1,n1;   write(11,*)d(i1);   enddo;close(11)
end subroutine write_asciifile_r4_1d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_r4_1d_wo_dim(fn,d)
character(len=*), intent(in) :: fn
real, intent(in) :: d(:)
n1_get=size(d);open (11,file=fn,action="write",form="formatted")
do i1=1,n1_get;   write(11,*)d(i1);   enddo;close(11)
end subroutine write_asciifile_r4_1d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_r8_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
double precision, intent(out) :: d(:)
integer,          intent(in)  :: n1
open (11,file=fn,status="old",action="read")
do i1=1,n1;   read(11,*)d(i1);   enddo;close(11)
end subroutine read_asciifile_r8_1d

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_mpi_r8_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
double precision, intent(out) :: d(:)
integer,          intent(in)  :: n1
if (rank.eq.0) call read_asciifile_r8_1d(fn,d,n1)
call sf90_bcast(d,n1)
end subroutine read_asciifile_mpi_r8_1d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_r8_1d(fn,d,n1)
character(len=*), intent(in) :: fn
double precision, intent(in) :: d(:)
integer,          intent(in) :: n1
open (11,file=fn,action="write",form="formatted")
do i1=1,n1;   write(11,*)d(i1);   enddo;close(11)
end subroutine write_asciifile_r8_1d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_r8_1d_wo_dim(fn,d)
character(len=*), intent(in) :: fn
double precision, intent(in) :: d(:)
open (11,file=fn,action="write",form="formatted")
n1_get=size(d);do i1=1,n1_get;write(11,*)d(i1);enddo;close(11)
end subroutine write_asciifile_r8_1d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_i4_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
integer,          intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
open (11,file=fn,status="old",action="read")
do i2=1,n2;   read(11,*)(d(i1,i2),i1=1,n1);   enddo;close(11)
end subroutine read_asciifile_i4_2d

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_mpi_i4_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
integer,          intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
if (rank.eq.0) call read_asciifile_i4_2d(fn,d,n1,n2)
call sf90_bcast(d,n1,n2)
end subroutine read_asciifile_mpi_i4_2d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_i4_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: d(:,:)
integer,          intent(in) :: n1,n2
open (11,file=fn,action="write",form="formatted")
do i2=1,n2;   write(11,*)(d(i1,i2),i1=1,n1);   enddo;   close(11)
end subroutine write_asciifile_i4_2d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_i4_2d_wo_dim(fn,d)
character(len=*), intent(in) :: fn
integer,          intent(in) :: d(:,:)
n1_get=size(d,1);   n2_get=size(d,2);
open (11,file=fn,action="write",form="formatted")
do i2=1,n2_get;write(11,*)(d(i1,i2),i1=1,n1_get);enddo;close(11)
end subroutine write_asciifile_i4_2d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_r4_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
real,             intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
open (11,file=fn,status="old",action="read")
do i2=1,n2;read(11,*)(d(i1,i2),i1=1,n1);enddo;close(11)
end subroutine read_asciifile_r4_2d

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_mpi_r4_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
real,             intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
if (rank.eq.0) call read_asciifile_r4_2d(fn,d,n1,n2)
call sf90_bcast(d,n1,n2)
end subroutine read_asciifile_mpi_r4_2d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_r4_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
real,             intent(in) :: d(:,:)
integer,          intent(in) :: n1,n2
open (11,file=fn,action="write",form="formatted")
do i2=1,n2;   write(11,*)(d(i1,i2),i1=1,n1);   enddo;   close(11)
end subroutine write_asciifile_r4_2d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_r4_2d_wo_dim(fn,d)
character(len=*), intent(in) :: fn
real,             intent(in) :: d(:,:)
n1_get=size(d,1);   n2_get=size(d,2)
open (11,file=fn,action="write",form="formatted")
do i2=1,n2_get;write(11,*)(d(i1,i2),i1=1,n1_get);enddo;close(11)
end subroutine write_asciifile_r4_2d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_r8_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
double precision, intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
open (11,file=fn,status="old",action="read")
do i2=1,n2;read(11,*)(d(i1,i2),i1=1,n1);enddo;close(11)
end subroutine read_asciifile_r8_2d

!-----------------------------------------------------------------------------------------
subroutine read_asciifile_mpi_r8_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
double precision, intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
if (rank.eq.0) call read_asciifile_r8_2d(fn,d,n1,n2)
call sf90_bcast(d,n1,n2)
end subroutine read_asciifile_mpi_r8_2d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_r8_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
double precision, intent(in) :: d(:,:)
integer,          intent(in) :: n1,n2
open (11,file=fn,action="write",form="formatted")
do i2=1,n2;write(11,*)(d(i1,i2),i1=1,n1);enddo;close(11)
end subroutine write_asciifile_r8_2d

!-----------------------------------------------------------------------------------------
subroutine write_asciifile_r8_2d_wo_dim(fn,d)
character(len=*), intent(in) :: fn
double precision, intent(in) :: d(:,:)
n1_get=size(d,1);   n2_get=size(d,2)
open (11,file=fn,action="write",form="formatted")
do i2=1,n2_get;write(11,*)(d(i1,i2),i1=1,n1_get);enddo;close(11)
end subroutine write_asciifile_r8_2d_wo_dim

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_binfile_i2_1d(fn,d,n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1
integer(kind=2),  intent(out):: d(:)
open(10,file=fn,access="direct",recl=n1*I4/2,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_i2_1d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_i2_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
integer(kind=2),  intent(out) :: d(:)
integer,          intent(in)  :: n1
if (rank.eq.0) call read_binfile_i2_1d(fn,d,n1)
call sf90_bcast(d,n1)
end subroutine read_binfile_mpi_i2_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i2_1d(fn, d, n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1
integer(kind=2),  intent(in) :: d(:)
open(10,file=fn,access="direct",recl=n1*I4/2,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i2_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i2_1d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
integer(kind=2),  intent(in) :: d(:)
n1_get=size(d)
open(10,file=fn,access="direct",recl=n1_get*I4/2,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i2_1d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_i2_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2
integer(kind=2),  intent(out):: d(:,:)
open(10,file=fn,access="direct",recl=n1*n2*I4/2,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_i2_2d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_i2_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
integer(kind=2),  intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
if (rank.eq.0) call read_binfile_i2_2d(fn,d,n1,n2)
call sf90_bcast(d,n1,n2)
end subroutine read_binfile_mpi_i2_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i2_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2
integer(kind=2),  intent(in) :: d(:,:)
open(10,file=fn,access="direct",recl=n1*n2*I4/2,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i2_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i2_2d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
integer(kind=2),  intent(in) :: d(:,:)
n1_get=size(d,1); n2_get=size(d,2)
open(10,file=fn,access="direct",recl=n1_get*n2_get*I4/2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i2_2d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_i2_3d(fn,d,n1,n2,n3)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3
integer(kind=2),  intent(out):: d(:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*I4/2,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_i2_3d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_i2_3d(fn,d,n1,n2,n3)
character(len=*), intent(in)  :: fn
integer(kind=2),  intent(out) :: d(:,:,:)
integer,          intent(in)  :: n1,n2,n3
if (rank.eq.0) call read_binfile_i2_3d(fn,d,n1,n2,n3)
call sf90_bcast(d,n1,n2,n3)
end subroutine read_binfile_mpi_i2_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i2_3d(fn,d,n1,n2,n3)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3
integer(kind=2),  intent(in) :: d(:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*I4/2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i2_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i2_3d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
integer(kind=2),  intent(in) :: d(:,:,:)
n1_get=size(d,1); n2_get=size(d,2); n3_get=size(d,3)
open(10,file=fn,access="direct",recl=n1_get*n2_get*n3_get*I4/2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i2_3d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_i2_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3,n4
integer(kind=2),  intent(out):: d(:,:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*n4*I4/2,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_i2_4d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_i2_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in)  :: fn
integer(kind=2),  intent(out) :: d(:,:,:,:)
integer,          intent(in)  :: n1,n2,n3,n4
if (rank.eq.0) call read_binfile_i2_4d(fn,d,n1,n2,n3,n4)
call sf90_bcast(d,n1,n2,n3,n4)
end subroutine read_binfile_mpi_i2_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i2_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3,n4
integer(kind=2),  intent(in) :: d(:,:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*n4*I4/2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i2_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i2_4d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
integer(kind=2),  intent(in) :: d(:,:,:,:)
n1_get=size(d,1); n2_get=size(d,2); n3_get=size(d,3); n4_get=size(d,4)
open(10,file=fn,access="direct",recl=n1_get*n2_get*n3_get*n4_get*I4/2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i2_4d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_i4_1d(fn,d,n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1
integer(kind=4),  intent(out):: d(:)
open(10,file=fn,access="direct",recl=n1*I4,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_i4_1d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_i4_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
integer(kind=4),  intent(out) :: d(:)
integer,          intent(in)  :: n1
if (rank.eq.0) call read_binfile_i4_1d(fn,d,n1)
call sf90_bcast(d,n1)
end subroutine read_binfile_mpi_i4_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i4_1d(fn, d, n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1
integer(kind=4),  intent(in) :: d(:)
open(10,file=fn,access="direct",recl=n1*I4,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i4_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i4_1d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
integer(kind=4),  intent(in) :: d(:)
n1_get=size(d)
open(10,file=fn,access="direct",recl=n1_get*I4,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i4_1d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_i4_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2
integer(kind=4),  intent(out):: d(:,:)
open(10,file=fn,access="direct",recl=n1*n2*I4,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_i4_2d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_i4_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
integer(kind=4),  intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
if (rank.eq.0) call read_binfile_i4_2d(fn,d,n1,n2)
call sf90_bcast(d,n1,n2)
end subroutine read_binfile_mpi_i4_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i4_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2
integer(kind=4),  intent(in) :: d(:,:)
open(10,file=fn,access="direct",recl=n1*n2*I4,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i4_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i4_2d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
integer(kind=4),  intent(in) :: d(:,:)
n1_get=size(d,1); n2_get=size(d,2)
open(10,file=fn,access="direct",recl=n1_get*n2_get*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i4_2d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_i4_3d(fn,d,n1,n2,n3)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3
integer(kind=4),  intent(out):: d(:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*I4,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_i4_3d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_i4_3d(fn,d,n1,n2,n3)
character(len=*), intent(in)  :: fn
integer(kind=4),  intent(out) :: d(:,:,:)
integer,          intent(in)  :: n1,n2,n3
if (rank.eq.0) call read_binfile_i4_3d(fn,d,n1,n2,n3)
call sf90_bcast(d,n1,n2,n3)
end subroutine read_binfile_mpi_i4_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i4_3d(fn,d,n1,n2,n3)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3
integer(kind=4),  intent(in) :: d(:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i4_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i4_3d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
integer(kind=4),  intent(in) :: d(:,:,:)
n1_get=size(d,1); n2_get=size(d,2); n3_get=size(d,3)
open(10,file=fn,access="direct",recl=n1_get*n2_get*n3_get*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i4_3d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_i4_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3,n4
integer(kind=4),  intent(out):: d(:,:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*n4*I4,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_i4_4d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_i4_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in)  :: fn
integer(kind=4),  intent(out) :: d(:,:,:,:)
integer,          intent(in)  :: n1,n2,n3,n4
if (rank.eq.0) call read_binfile_i4_4d(fn,d,n1,n2,n3,n4)
call sf90_bcast(d,n1,n2,n3,n4)
end subroutine read_binfile_mpi_i4_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i4_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3,n4
integer(kind=4),  intent(in) :: d(:,:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*n4*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i4_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_i4_4d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
integer(kind=4),  intent(in) :: d(:,:,:,:)
n1_get=size(d,1); n2_get=size(d,2); n3_get=size(d,3); n4_get=size(d,4)
open(10,file=fn,access="direct",recl=n1_get*n2_get*n3_get*n4_get*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_i4_4d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_r4_1d(fn,d,n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1
real(kind=4),     intent(out):: d(:)
open(10,file=fn,access="direct",recl=n1*I4,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_r4_1d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_r4_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
real(kind=4),     intent(out) :: d(:)
integer,          intent(in)  :: n1
if (rank.eq.0) call read_binfile_r4_1d(fn,d,n1)
call sf90_bcast(d,n1)
end subroutine read_binfile_mpi_r4_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r4_1d(fn, d, n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1
real(kind=4),  intent(in) :: d(:)
open(10,file=fn,access="direct",recl=n1*I4,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r4_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r4_1d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
real(kind=4),  intent(in) :: d(:)
n1_get=size(d)
open(10,file=fn,access="direct",recl=n1_get*I4,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r4_1d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_r4_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2
real(kind=4),  intent(out):: d(:,:)
open(10,file=fn,access="direct",recl=n1*n2*I4,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_r4_2d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_r4_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
real(kind=4),  intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
if (rank.eq.0) call read_binfile_r4_2d(fn,d,n1,n2)
call sf90_bcast(d,n1,n2)
end subroutine read_binfile_mpi_r4_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r4_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2
real(kind=4),  intent(in) :: d(:,:)
open(10,file=fn,access="direct",recl=n1*n2*I4,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r4_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r4_2d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
real(kind=4),  intent(in) :: d(:,:)
n1_get=size(d,1); n2_get=size(d,2)
open(10,file=fn,access="direct",recl=n1_get*n2_get*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r4_2d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_r4_3d(fn,d,n1,n2,n3)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3
real(kind=4),  intent(out):: d(:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*I4,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_r4_3d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_r4_3d(fn,d,n1,n2,n3)
character(len=*), intent(in)  :: fn
real(kind=4),  intent(out) :: d(:,:,:)
integer,          intent(in)  :: n1,n2,n3
if (rank.eq.0) call read_binfile_r4_3d(fn,d,n1,n2,n3)
call sf90_bcast(d,n1,n2,n3)
end subroutine read_binfile_mpi_r4_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r4_3d(fn,d,n1,n2,n3)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3
real(kind=4),  intent(in) :: d(:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r4_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r4_3d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
real(kind=4),  intent(in) :: d(:,:,:)
n1_get=size(d,1); n2_get=size(d,2); n3_get=size(d,3)
open(10,file=fn,access="direct",recl=n1_get*n2_get*n3_get*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r4_3d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_r4_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3,n4
real(kind=4),  intent(out):: d(:,:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*n4*I4,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_r4_4d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_r4_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in)  :: fn
real(kind=4),  intent(out) :: d(:,:,:,:)
integer,          intent(in)  :: n1,n2,n3,n4
if (rank.eq.0) call read_binfile_r4_4d(fn,d,n1,n2,n3,n4)
call sf90_bcast(d,n1,n2,n3,n4)
end subroutine read_binfile_mpi_r4_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r4_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3,n4
real(kind=4),  intent(in) :: d(:,:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*n4*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r4_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r4_4d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
real(kind=4),  intent(in) :: d(:,:,:,:)
n1_get=size(d,1); n2_get=size(d,2); n3_get=size(d,3); n4_get=size(d,4)
open(10,file=fn,access="direct",recl=n1_get*n2_get*n3_get*n4_get*I4,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r4_4d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_r8_1d(fn,d,n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1
real(kind=8),     intent(out):: d(:)
open(10,file=fn,access="direct",recl=n1*I4*2,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_r8_1d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_r8_1d(fn,d,n1)
character(len=*), intent(in)  :: fn
real(kind=8),     intent(out) :: d(:)
integer,          intent(in)  :: n1
if (rank.eq.0) call read_binfile_r8_1d(fn,d,n1)
call sf90_bcast(d,n1)
end subroutine read_binfile_mpi_r8_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r8_1d(fn, d, n1)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1
real(kind=8),  intent(in) :: d(:)
open(10,file=fn,access="direct",recl=n1*I4*2,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r8_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r8_1d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
real(kind=8),  intent(in) :: d(:)
n1_get=size(d)
open(10,file=fn,access="direct",recl=n1_get*I4*2,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r8_1d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_r8_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2
real(kind=8),  intent(out):: d(:,:)
open(10,file=fn,access="direct",recl=n1*n2*I4*2,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_r8_2d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_r8_2d(fn,d,n1,n2)
character(len=*), intent(in)  :: fn
real(kind=8),  intent(out) :: d(:,:)
integer,          intent(in)  :: n1,n2
if (rank.eq.0) call read_binfile_r8_2d(fn,d,n1,n2)
call sf90_bcast(d,n1,n2)
end subroutine read_binfile_mpi_r8_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r8_2d(fn,d,n1,n2)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2
real(kind=8),  intent(in) :: d(:,:)
open(10,file=fn,access="direct",recl=n1*n2*I4*2,status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r8_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r8_2d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
real(kind=8),  intent(in) :: d(:,:)
n1_get=size(d,1); n2_get=size(d,2)
open(10,file=fn,access="direct",recl=n1_get*n2_get*I4*2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r8_2d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_r8_3d(fn,d,n1,n2,n3)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3
real(kind=8),  intent(out):: d(:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*I4*2,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_r8_3d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_r8_3d(fn,d,n1,n2,n3)
character(len=*), intent(in)  :: fn
real(kind=8),  intent(out) :: d(:,:,:)
integer,          intent(in)  :: n1,n2,n3
if (rank.eq.0) call read_binfile_r8_3d(fn,d,n1,n2,n3)
call sf90_bcast(d,n1,n2,n3)
end subroutine read_binfile_mpi_r8_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r8_3d(fn,d,n1,n2,n3)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3
real(kind=8),  intent(in) :: d(:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*I4*2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r8_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r8_3d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
real(kind=8),  intent(in) :: d(:,:,:)
n1_get=size(d,1); n2_get=size(d,2); n3_get=size(d,3)
open(10,file=fn,access="direct",recl=n1_get*n2_get*n3_get*I4*2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r8_3d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_r8_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3,n4
real(kind=8),  intent(out):: d(:,:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*n4*I4*2,status="old",action="read")
read(10, rec=1) d;   close(10)
end subroutine read_binfile_r8_4d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_mpi_r8_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in)  :: fn
real(kind=8),  intent(out) :: d(:,:,:,:)
integer,          intent(in)  :: n1,n2,n3,n4
if (rank.eq.0) call read_binfile_r8_4d(fn,d,n1,n2,n3,n4)
call sf90_bcast(d,n1,n2,n3,n4)
end subroutine read_binfile_mpi_r8_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r8_4d(fn,d,n1,n2,n3,n4)
character(len=*), intent(in) :: fn
integer,          intent(in) :: n1,n2,n3,n4
real(kind=8),  intent(in) :: d(:,:,:,:)
open(10,file=fn,access="direct",recl=n1*n2*n3*n4*I4*2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r8_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_r8_4d_wo_dim(fn, d)
character(len=*), intent(in) :: fn
real(kind=8),  intent(in) :: d(:,:,:,:)
n1_get=size(d,1); n2_get=size(d,2); n3_get=size(d,3); n4_get=size(d,4)
open(10,file=fn,access="direct",recl=n1_get*n2_get*n3_get*n4_get*I4*2,&
   status="replace",action="write")
write(10,rec=1) d;   close(10)
end subroutine write_binfile_r8_4d_wo_dim

!-----------------------------------------------------------------------------------------
subroutine read_binfile_delete_r4_2d(fn, d, n1, n2)
character(len=*), intent(in) :: fn
integer, intent(in) :: n1, n2
real                :: d(:,:)
open(10,file=fn,access="direct",recl=n1*I4,status="old",action="read")
do i2 = 1, n2
   read(10, rec=i2)(d(i1,i2),i1=1,n1)
enddo
close(10,status="delete")
end subroutine read_binfile_delete_r4_2d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_dfor_i2tor4_1d(fn,d,n1,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),    intent(out) :: d(:)
real,             intent(in) :: dfor_coe
integer,          intent(in) :: n1
integer(kind=2) :: d1(n1)
! if df_sign>0 -- read, then divide   10e(abs(dfor_coe))
! if df_sign<0 -- read, then multiple 10e(abs(dfor_coe))
call read_binfile_i2_1d(fn,d1,n1)
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d=d1/(10**abs(dfor_coe))
else
   d=d1*(10**abs(dfor_coe))
endif
end subroutine read_binfile_dfor_i2tor4_1d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_dfor_r4toi2_1d(fn,d,n1,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),     intent(in) :: d(:)
real,             intent(in) :: dfor_coe
integer,          intent(in) :: n1
integer(kind=2) :: d1(n1)
! if df_sign>0 -- multiple 10e(abs(dfor_coe)), then write
! if df_sign<0 -- divide   10e(abs(dfor_coe)), then write
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d1=d*(10**abs(dfor_coe))
else
   d1=d/(10**abs(dfor_coe))
endif
call write_binfile_i2_1d(fn,d1,n1)
end subroutine write_binfile_dfor_r4toi2_1d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_dfor_i2tor4_2d(fn,d,n1,n2,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),    intent(out) :: d(:,:)
real,             intent(in) :: dfor_coe
integer,          intent(in) :: n1,n2
integer(kind=2) :: d1(n1,n2)
! if df_sign>0 -- read, then divide   10e(abs(dfor_coe))
! if df_sign<0 -- read, then multiple 10e(abs(dfor_coe))
call read_binfile_i2_2d(fn,d1,n1,n2)
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d=d1/(10**abs(dfor_coe))
else
   d=d1*(10**abs(dfor_coe))
endif
end subroutine read_binfile_dfor_i2tor4_2d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_dfor_r4toi2_2d(fn,d,n1,n2,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),     intent(in) :: d(:,:)
real,             intent(in) :: dfor_coe
integer,          intent(in) :: n1,n2
integer(kind=2) :: d1(n1,n2)
! if df_sign>0 -- multiple 10e(abs(dfor_coe)), then write
! if df_sign<0 -- divide   10e(abs(dfor_coe)), then write
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d1=d*(10**abs(dfor_coe))
else
   d1=d/(10**abs(dfor_coe))
endif
call write_binfile_i2_2d(fn,d1,n1,n2)
end subroutine write_binfile_dfor_r4toi2_2d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_dfor_i2tor4_3d(fn,d,n1,n2,n3,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),    intent(out) :: d(:,:,:)
real,             intent(in) :: dfor_coe
integer,          intent(in) :: n1,n2,n3
integer(kind=2) :: d1(n1,n2,n3)
! if df_sign>0 -- read, then divide   10e(abs(dfor_coe))
! if df_sign<0 -- read, then multiple 10e(abs(dfor_coe))
call read_binfile_i2_3d(fn,d1,n1,n2,n3)
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d=d1/(10**abs(dfor_coe))
else
   d=d1*(10**abs(dfor_coe))
endif
end subroutine read_binfile_dfor_i2tor4_3d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_dfor_r4toi2_3d(fn,d,n1,n2,n3,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),     intent(in) :: d(:,:,:)
real,             intent(in) :: dfor_coe
integer,          intent(in) :: n1,n2,n3
integer(kind=2) :: d1(n1,n2,n3)
! if df_sign>0 -- multiple 10e(abs(dfor_coe)), then write
! if df_sign<0 -- divide   10e(abs(dfor_coe)), then write
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d1=d*(10**abs(dfor_coe))
else
   d1=d/(10**abs(dfor_coe))
endif
call write_binfile_i2_3d(fn,d1,n1,n2,n3)
end subroutine write_binfile_dfor_r4toi2_3d

!-----------------------------------------------------------------------------------------
subroutine read_binfile_dfor_i2tor4_4d(fn,d,n1,n2,n3,n4,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),    intent(out) :: d(:,:,:,:)
real,             intent(in) :: dfor_coe
integer,          intent(in) :: n1,n2,n3,n4
integer(kind=2) :: d1(n1,n2,n3,n4)
! if df_sign>0 -- read, then divide   10e(abs(dfor_coe))
! if df_sign<0 -- read, then multiple 10e(abs(dfor_coe))
call read_binfile_i2_4d(fn,d1,n1,n2,n3,n4)
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d=d1/(10**abs(dfor_coe))
else
   d=d1*(10**abs(dfor_coe))
endif
end subroutine read_binfile_dfor_i2tor4_4d

!-----------------------------------------------------------------------------------------
subroutine write_binfile_dfor_r4toi2_4d(fn,d,n1,n2,n3,n4,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),     intent(in) :: d(:,:,:,:)
real,             intent(in) :: dfor_coe
integer,          intent(in) :: n1,n2,n3,n4
integer(kind=2) :: d1(n1,n2,n3,n4)
! if df_sign>0 -- multiple 10e(abs(dfor_coe)), then write
! if df_sign<0 -- divide   10e(abs(dfor_coe)), then write
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d1=d*(10**abs(dfor_coe))
else
   d1=d/(10**abs(dfor_coe))
endif
call write_binfile_i2_4d(fn,d1,n1,n2,n3,n4)
end subroutine write_binfile_dfor_r4toi2_4d
!====================================================================================
!! Following is added by Bowen Guo
!! To read the whole shot gather data
!! Memory expensive
!subroutine read_seis_all(seis_all,fn_head,c,a,m,o)
!character(len=*),intent(in)               :: fn_head
!type(coord2d),intent(in)                  :: c
!type(aperture),intent(in)                 :: a
!type(model2d),intent(in)                  :: m
!type(other),intent(in)                    :: o

!integer    :: nt,ng,ns

!nt=o%nt_out
!ng=c%ng(is)

















end module module_io
