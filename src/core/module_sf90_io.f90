module module_sf90_io
use module_sf90_mpi
use module_io
use module_math
implicit none

integer, private :: fid
real,    private :: df_sign

interface sf90_readself
   module procedure mpi_readself_i2_1d_1i
   module procedure mpi_readself_i4_1d_1i
   module procedure mpi_readself_r4_1d_1r
   module procedure mpi_readself_r8_1d_1r
   module procedure mpi_readself_i2_2d_1i
   module procedure mpi_readself_i4_2d_1i
   module procedure mpi_readself_r4_2d_1r
   module procedure mpi_readself_r8_2d_1r
   module procedure mpi_readself_i2_3d_1i
   module procedure mpi_readself_i4_3d_1i
   module procedure mpi_readself_r4_3d_1r
   module procedure mpi_readself_r8_3d_1r
   module procedure mpi_readself_i2_4d_1i
   module procedure mpi_readself_i4_4d_1i
   module procedure mpi_readself_r4_4d_1r
   module procedure mpi_readself_r8_4d_1r
   module procedure mpi_readself_dfor_i2tor4_1d_1r
   module procedure mpi_readself_dfor_i2tor4_2d_1r
   module procedure mpi_readself_dfor_i2tor4_3d_1r
   module procedure mpi_readself_dfor_i2tor4_4d_1r
end interface sf90_readself

contains

!---------------------------------------------------------------------------------------
subroutine mpi_readself_i2_1d_1i(fn,d1,n1)
character(len=*), intent(in)  :: fn
integer(kind=2),  intent(out) :: d1(:)
integer,          intent(in)  :: n1
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_i2_1d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_readself_i4_1d_1i(fn,d1,n1)
character(len=*), intent(in)  :: fn
integer,          intent(out) :: d1(:)
integer,          intent(in)  :: n1
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_i4_1d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_readself_r4_1d_1r(fn,d1,n1)
character(len=*), intent(in)  :: fn
real,             intent(out) :: d1(:)
integer,          intent(in)  :: n1
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1,MPI_REAL,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_r4_1d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_r8_1d_1r(fn,d1,n1)
character(len=*), intent(in)  :: fn
real(kind=8),     intent(out) :: d1(:)
integer,          intent(in)  :: n1
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_r8_1d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_i2_2d_1i(fn,d1,n1,n2)
character(len=*), intent(in)  :: fn
integer(kind=2),  intent(out) :: d1(:,:)
integer,          intent(in)  :: n1,n2
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_i2_2d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_readself_i4_2d_1i(fn,d1,n1,n2)
character(len=*), intent(in)  :: fn
integer,          intent(out) :: d1(:,:)
integer,          intent(in)  :: n1,n2
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_i4_2d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_readself_r4_2d_1r(fn,d1,n1,n2)
character(len=*), intent(in)  :: fn
real,             intent(out) :: d1(:,:)
integer,          intent(in)  :: n1,n2
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2,MPI_REAL,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_r4_2d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_r8_2d_1r(fn,d1,n1,n2)
character(len=*), intent(in)  :: fn
real(kind=8),     intent(out) :: d1(:,:)
integer,          intent(in)  :: n1,n2
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_r8_2d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_i2_3d_1i(fn,d1,n1,n2,n3)
character(len=*), intent(in)  :: fn
integer(kind=2),  intent(out) :: d1(:,:,:)
integer,          intent(in)  :: n1,n2,n3
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2*n3,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_i2_3d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_readself_i4_3d_1i(fn,d1,n1,n2,n3)
character(len=*), intent(in)  :: fn
integer,          intent(out) :: d1(:,:,:)
integer,          intent(in)  :: n1,n2,n3
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2*n3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_i4_3d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_readself_r4_3d_1r(fn,d1,n1,n2,n3)
character(len=*), intent(in)  :: fn
real,             intent(out) :: d1(:,:,:)
integer,          intent(in)  :: n1,n2,n3
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2*n3,MPI_REAL,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_r4_3d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_r8_3d_1r(fn,d1,n1,n2,n3)
character(len=*), intent(in)  :: fn
real(kind=8),     intent(out) :: d1(:,:,:)
integer,          intent(in)  :: n1,n2,n3
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2*n3,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_r8_3d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_i2_4d_1i(fn,d1,n1,n2,n3,n4)
character(len=*), intent(in)  :: fn
integer(kind=2),  intent(out) :: d1(:,:,:,:)
integer,          intent(in)  :: n1,n2,n3,n4
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2*n3*n4,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_i2_4d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_readself_i4_4d_1i(fn,d1,n1,n2,n3,n4)
character(len=*), intent(in)  :: fn
integer,          intent(out) :: d1(:,:,:,:)
integer,          intent(in)  :: n1,n2,n3,n4
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2*n3*n4,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_i4_4d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_readself_r4_4d_1r(fn,d1,n1,n2,n3,n4)
character(len=*), intent(in)  :: fn
real,             intent(out) :: d1(:,:,:,:)
integer,          intent(in)  :: n1,n2,n3,n4
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2*n3*n4,MPI_REAL,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_r4_4d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_r8_4d_1r(fn,d1,n1,n2,n3,n4)
character(len=*), intent(in)  :: fn
real(kind=8),     intent(out) :: d1(:,:,:,:)
integer,          intent(in)  :: n1,n2,n3,n4
call MPI_FILE_OPEN(MPI_COMM_SELF, fn, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fid, ierr)
call MPI_FILE_READ(fid,d1,n1*n2*n3*n4,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(fid,ierr)
end subroutine mpi_readself_r8_4d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_dfor_i2tor4_1d_1r(fn,d1,n1,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),    intent(out) :: d1(:)
integer,          intent(in) :: n1
real,             intent(in) :: dfor_coe
integer(kind=2) :: d(n1)
call mpi_readself_i2_1d_1i(fn,d,n1)
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d1=d/(10**abs(dfor_coe))
else
   d1=d*(10**abs(dfor_coe))
endif
end subroutine mpi_readself_dfor_i2tor4_1d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_dfor_i2tor4_2d_1r(fn,d1,n1,n2,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),    intent(out) :: d1(n1,n2)
integer,          intent(in) :: n1,n2
real,             intent(in) :: dfor_coe
integer(kind=2) :: d(n1,n2)
call mpi_readself_i2_2d_1i(fn,d,n1,n2)
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d1=d/(10**abs(dfor_coe))
else
   d1=d*(10**abs(dfor_coe))
endif
end subroutine mpi_readself_dfor_i2tor4_2d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_dfor_i2tor4_3d_1r(fn,d1,n1,n2,n3,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),    intent(out) :: d1(n1,n2,n3)
integer,          intent(in) :: n1,n2,n3
real,             intent(in) :: dfor_coe
integer(kind=2) :: d(n1,n2,n3)
call mpi_readself_i2_3d_1i(fn,d,n1,n2,n3)
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d1=d/(10**abs(dfor_coe))
else
   d1=d*(10**abs(dfor_coe))
endif
end subroutine mpi_readself_dfor_i2tor4_3d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_readself_dfor_i2tor4_4d_1r(fn,d1,n1,n2,n3,n4,dfor_coe)
character(len=*), intent(in) :: fn
real(kind=4),    intent(out) :: d1(n1,n2,n3,n4)
integer,          intent(in) :: n1,n2,n3,n4
real,             intent(in) :: dfor_coe
integer(kind=2) :: d(n1,n2,n3,n4)
call mpi_readself_i2_4d_1i(fn,d,n1,n2,n3,n4)
df_sign=signn(dfor_coe)
if (df_sign>0.0) then
   d1=d/(10**abs(dfor_coe))
else
   d1=d*(10**abs(dfor_coe))
endif
end subroutine mpi_readself_dfor_i2tor4_4d_1r

end module module_sf90_io
