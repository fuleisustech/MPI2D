!
!  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
!  King Abdullah University of Science and Technology, All rights reserved.
!
!  Sponsors of CSIM are granted a non-exclusive, irrevocable royalty free
!  world-wide license to use this software and associated documentation files
!  (the "Software"), in their business, including the rights to modify and
!  distribute the Software to their affiliates, partners, clients or consultants
!  as necessary in the conduct of the sponsors business, all without accounting
!  to the King Abdullah University of Science and Technology, subject to the
!  following conditions:
!
!  The above copyright notice and this permission notice shall be included
!  in all copies or substantial portions of the Software.
!
!  Warranty Disclaimer:
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!  DEALINGS IN THE SOFTWARE.
!

!
!  file:     sf90_mpi.f90
!  @version: 1.0
!  @author:  Xin Wang
!  email:    xin.wang@kaust.edu.sa
!  date:     December 2013
!  purpose:  module for MPI
!

module module_sf90_mpi
implicit none
include 'mpif.h'
integer  :: root, rank, nsize, tag, status(MPI_STATUS_SIZE), &
            ierr, sendtag, recvtag, currentShot, &
            source_tag, data_tag, terminate_tag, result_tag, &
            slave, count, req(12), st(MPI_STATUS_SIZE,12)

interface sf90_sum
  module procedure mpi_sum_i4
  module procedure mpi_sum_i8
  module procedure mpi_sum_i4_1d
  module procedure mpi_sum_i4_2d
  module procedure mpi_sum_i4_3d
  module procedure mpi_sum_i4_4d
  module procedure mpi_sum_r4
  module procedure mpi_sum_r8
  module procedure mpi_sum_r4_1d
  module procedure mpi_sum_r4_2d
  module procedure mpi_sum_r4_3d
  module procedure mpi_sum_r4_4d
  module procedure mpi_sum_r8_1d
  module procedure mpi_sum_r8_2d
end interface sf90_sum

interface sf90_allsum
  module procedure mpi_allsum_i4
  module procedure mpi_allsum_i8
  module procedure mpi_allsum_i4_1d
  module procedure mpi_allsum_i4_2d
  module procedure mpi_allsum_i4_3d
  module procedure mpi_allsum_i4_4d
  module procedure mpi_allsum_r4
  module procedure mpi_allsum_r8
  module procedure mpi_allsum_r4_1d
  module procedure mpi_allsum_r4_2d
  module procedure mpi_allsum_r4_3d
  module procedure mpi_allsum_r4_4d
  module procedure mpi_allsum_r8_1d
  module procedure mpi_allsum_r8_2d
end interface sf90_allsum

interface sf90_bcast
   module procedure mpi_bcast_i2_1i
   module procedure mpi_bcast_i2_2i
   module procedure mpi_bcast_i2_3i
   module procedure mpi_bcast_i2_4i
   module procedure mpi_bcast_i4_1i
   module procedure mpi_bcast_i4_2i
   module procedure mpi_bcast_i4_3i
   module procedure mpi_bcast_i4_4i
   module procedure mpi_bcast_i8_1i
   module procedure mpi_bcast_i8_2i
   module procedure mpi_bcast_i8_3i
   module procedure mpi_bcast_i8_4i
   module procedure mpi_bcast_r4_1r
   module procedure mpi_bcast_r4_2r
   module procedure mpi_bcast_r4_3r
   module procedure mpi_bcast_r4_4r
   module procedure mpi_bcast_r8_1r
   module procedure mpi_bcast_r8_2r
   module procedure mpi_bcast_r8_3r
   module procedure mpi_bcast_r8_4r
   module procedure mpi_bcast_i2_1d_1i
   module procedure mpi_bcast_i2_1d_2i
   module procedure mpi_bcast_i2_1d_3i
   module procedure mpi_bcast_i2_1d_4i
   module procedure mpi_bcast_i4_1d_1i
   module procedure mpi_bcast_i4_1d_2i
   module procedure mpi_bcast_i4_1d_3i
   module procedure mpi_bcast_i4_1d_4i
   module procedure mpi_bcast_i8_1d_1i
   module procedure mpi_bcast_i8_1d_2i
   module procedure mpi_bcast_i8_1d_3i
   module procedure mpi_bcast_i8_1d_4i
   module procedure mpi_bcast_r4_1d_1r
   module procedure mpi_bcast_r4_1d_2r
   module procedure mpi_bcast_r4_1d_3r
   module procedure mpi_bcast_r4_1d_4r
   module procedure mpi_bcast_r4_1d_5r
   module procedure mpi_bcast_r4_1d_7r
   module procedure mpi_bcast_r8_1d_1r
   module procedure mpi_bcast_r8_1d_2r
   module procedure mpi_bcast_r8_1d_3r
   module procedure mpi_bcast_r8_1d_4r
   module procedure mpi_bcast_i2_2d_1i
   module procedure mpi_bcast_i2_2d_2i
   module procedure mpi_bcast_i2_2d_3i
   module procedure mpi_bcast_i2_2d_4i
   module procedure mpi_bcast_i4_2d_1i
   module procedure mpi_bcast_i4_2d_2i
   module procedure mpi_bcast_i4_2d_3i
   module procedure mpi_bcast_i4_2d_4i
   module procedure mpi_bcast_i8_2d_1i
   module procedure mpi_bcast_i8_2d_2i
   module procedure mpi_bcast_i8_2d_3i
   module procedure mpi_bcast_i8_2d_4i
   module procedure mpi_bcast_r4_2d_1r
   module procedure mpi_bcast_r4_2d_2r
   module procedure mpi_bcast_r4_2d_3r
   module procedure mpi_bcast_r4_2d_4r
   module procedure mpi_bcast_r4_2d_5r
   module procedure mpi_bcast_r4_2d_6r
   module procedure mpi_bcast_r4_2d_7r
   module procedure mpi_bcast_r8_2d_1r
   module procedure mpi_bcast_r8_2d_2r
   module procedure mpi_bcast_r8_2d_3r
   module procedure mpi_bcast_r8_2d_4r
   module procedure mpi_bcast_i2_3d_1i
   module procedure mpi_bcast_i2_3d_2i
   module procedure mpi_bcast_i2_3d_3i
   module procedure mpi_bcast_i2_3d_4i
   module procedure mpi_bcast_i4_3d_1i
   module procedure mpi_bcast_i4_3d_2i
   module procedure mpi_bcast_i4_3d_3i
   module procedure mpi_bcast_i4_3d_4i
   module procedure mpi_bcast_i8_3d_1i
   module procedure mpi_bcast_i8_3d_2i
   module procedure mpi_bcast_i8_3d_3i
   module procedure mpi_bcast_i8_3d_4i
   module procedure mpi_bcast_r4_3d_1r
   module procedure mpi_bcast_r4_3d_2r
   module procedure mpi_bcast_r4_3d_3r
   module procedure mpi_bcast_r4_3d_4r
   module procedure mpi_bcast_r8_3d_1r
   module procedure mpi_bcast_r8_3d_2r
   module procedure mpi_bcast_r8_3d_3r
   module procedure mpi_bcast_r8_3d_4r
   module procedure mpi_bcast_i2_4d_1i
   module procedure mpi_bcast_i2_4d_2i
   module procedure mpi_bcast_i2_4d_3i
   module procedure mpi_bcast_i2_4d_4i
   module procedure mpi_bcast_i4_4d_1i
   module procedure mpi_bcast_i4_4d_2i
   module procedure mpi_bcast_i4_4d_3i
   module procedure mpi_bcast_i4_4d_4i
   module procedure mpi_bcast_i8_4d_1i
   module procedure mpi_bcast_i8_4d_2i
   module procedure mpi_bcast_i8_4d_3i
   module procedure mpi_bcast_i8_4d_4i
   module procedure mpi_bcast_r4_4d_1r
   module procedure mpi_bcast_r4_4d_2r
   module procedure mpi_bcast_r4_4d_3r
   module procedure mpi_bcast_r4_4d_4r
   module procedure mpi_bcast_r8_4d_1r
   module procedure mpi_bcast_r8_4d_2r
   module procedure mpi_bcast_r8_4d_3r
   module procedure mpi_bcast_r8_4d_4r
end interface sf90_bcast

interface message
   module procedure message_c
   module procedure message_c2
   module procedure message_c1i1
   module procedure message_c1i1c1
   module procedure message_cicdp
   module procedure message_cicr
   module procedure message_c1long1
   module procedure message_c1long2
   module procedure message_c1r1
   module procedure message_c1r1c1dp1
   module procedure message_c1d1
   module procedure message_c1l1
   module procedure message_c2i2
   module procedure message_c3i3
   module procedure message_c2r2
   module procedure message_c2d2
   module procedure message_c2l2
   module procedure message_c3r3
   module procedure message_c3d3
   module procedure message_c3l3
   module procedure message_i3
   module procedure message_r3
end interface message

interface message_all
   module procedure message_all_c
   module procedure message_all_c2
   module procedure message_all_c1i1
   module procedure message_all_c1long1
   module procedure message_all_c1long2
   module procedure message_all_c1r1
   module procedure message_all_c1d1
   module procedure message_all_c1l1
   module procedure message_all_c2i2
   module procedure message_all_c2r2
   module procedure message_all_c2d2
   module procedure message_all_c2l2
   module procedure message_all_c3i3
   module procedure message_all_c3r3
   module procedure message_all_c3d3
   module procedure message_all_c3l3
   module procedure message_all_c_i_c_c_i_c
end interface message_all

contains

!------------------------------------------------------------------------------
subroutine init_mpi

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nsize,ierr)

end subroutine init_mpi

!------------------------------------------------------------------------------
subroutine stop_mpi

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_Finalize(ierr)

end subroutine stop_mpi

!-----------------------------------------------------------------------------------------
! Static load balancing
!
! Purpose: Assign starting and ending shot numbers to an MPI process
!
! Modified from the original version written by Jianming Sheng
!
! is1 = starting shot number of the original list of shot numbers
! is2 = ending shot number of the original list of shot numbers
! is_begin = starting shot number for an MPI process
! is2_end = ending shot number for an MPI process
!
subroutine get_assigned_continue(ns, ns_me, is_me)

integer, intent(in)  :: ns
integer, intent(out) :: ns_me,is_me(:)
integer              :: is,ns_node, ns_left, ns_assigned1, nsize_assigned1,&
                        ns_assigned2, nsize_assigned2, ns_assigned

! ns # less than nsize
if (ns<nsize) then
   if (rank.eq.0) write(*,*)"WARNING: # of shots ",ns," is less than # of cores",nsize
   !call flush(6)
   if (rank<=ns-1) then
      ns_me=1;   is_me=rank+1
   else
      ns_me=0
   endif
   return
endif

! one cpu case
if (nsize == 1) then
  ns_me=ns
  do is=1,ns_me
     is_me(is)=is
  enddo
  return
endif

ns_node = int(ns/nsize)

!if (ns_node == 0) then
!  ns_me=ns_node
!  sid_me=sid(rank*ns_me+1:(rank+1)*ns_me)
!  return
!endif

ns_left = ns - ns_node*nsize
if ( rank < (nsize-ns_left) ) then
  ns_me=ns_node
  nsize_assigned1=rank
  ns_assigned1=ns_me*nsize_assigned1
  ns_assigned=ns_assigned1
else 
  ns_me=ns_node+1
  nsize_assigned1=nsize-ns_left
  ns_assigned1=ns_node*nsize_assigned1
  nsize_assigned2=rank-(nsize-ns_left)
  ns_assigned2=ns_me*nsize_assigned2
  ns_assigned=ns_assigned1+ns_assigned2
endif
do is=1,ns_me
   is_me(is)=ns_assigned+is
enddo

end subroutine get_assigned_continue

!-----------------------------------------------------------------------------------------
subroutine get_assigned_evenly(ns,ns_me,is_me)

integer, intent(in)  :: ns
integer, intent(out) :: ns_me,is_me(:)
integer              :: is,ns_node, ns_left

! ns # less than nsize
if (ns<nsize) then
   if (rank.eq.0) write(*,*)"WARNING: # of shots ",ns," is less than # of cores",nsize
   !call flush(6)
   if (rank<=ns-1) then
      ns_me=1;   is_me=rank+1
   else
      ns_me=0
   endif
   return
endif

! one cpu case
if (nsize == 1) then
  ns_me=ns
  do is=1,ns_me
     is_me(is)=is
  enddo
  return
endif

ns_node = int(ns/nsize)

ns_left = ns - ns_node*nsize
if ( rank < (nsize-ns_left) ) then
  ns_me=ns_node
else
  ns_me=ns_node+1
endif
if ( rank < (nsize-ns_left) ) then
   do is=1,ns_me
     is_me(is)=(is-1)*nsize+rank+1
   enddo
else
   do is=1,ns_me-1
      is_me(is)=(is-1)*nsize+rank+1
   enddo
   is_me(ns_me)=(ns_me-1)*nsize+rank+1-(nsize-ns_left)
endif   
end subroutine get_assigned_evenly

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i4(d1,d2)
integer, intent(in) :: d1
integer,intent(out) :: d2
call MPI_REDUCE(d1,d2,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i4

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i8(d1,d2)
integer(kind=8), intent(in) :: d1
integer(kind=8),intent(out) :: d2
call MPI_REDUCE(d1,d2,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i8

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r4(d1,d2)
real, intent(in)  :: d1
real, intent(out) :: d2
call MPI_REDUCE(d1,d2,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r4

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r8(d1,d2)
real(kind=8), intent(in)  :: d1
real(kind=8), intent(out) :: d2
call MPI_REDUCE(d1,d2,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r8

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i4_1d(d1,d2,n1)
integer, intent(in) :: d1(:),n1
integer,intent(out) :: d2(:)
call MPI_REDUCE(d1,d2,n1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i4_1d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i8_1d(d1,d2,n1)
integer(kind=8), intent(in) :: d1(:)
integer,         intent(in) :: n1
integer(kind=8),intent(out) :: d2(:)
call MPI_REDUCE(d1,d2,n1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i8_1d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r4_1d(d1,d2,n1)
real,  intent(in)  :: d1(:)
integer,intent(in) :: n1
real,  intent(out) :: d2(:)
call MPI_REDUCE(d1,d2,n1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r4_1d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r8_1d(d1,d2,n1)
real(kind=8), intent(in)  :: d1(:)
integer,      intent(in)  :: n1
real(kind=8), intent(out) :: d2(:)
call MPI_REDUCE(d1,d2,n1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r8_1d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i4_2d(d1,d2,n1,n2)
integer, intent(in) :: d1(:,:),n1,n2
integer,intent(out) :: d2(:,:)
call MPI_REDUCE(d1,d2,n1*n2,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i4_2d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i8_2d(d1,d2,n1,n2)
integer(kind=8), intent(in) :: d1(:,:)
integer,         intent(in) :: n1,n2
integer(kind=8),intent(out) :: d2(:,:)
call MPI_REDUCE(d1,d2,n1*n2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i8_2d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r4_2d(d1,d2,n1,n2)
real,  intent(in)  :: d1(:,:)
integer,intent(in) :: n1,n2
real,  intent(out) :: d2(:,:)
call MPI_REDUCE(d1,d2,n1*n2,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r4_2d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r8_2d(d1,d2,n1,n2)
real(kind=8), intent(in)  :: d1(:,:)
integer,      intent(in)  :: n1,n2
real(kind=8), intent(out) :: d2(:,:)
call MPI_REDUCE(d1,d2,n1*n2,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r8_2d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i4_3d(d1,d2,n1,n2,n3)
integer, intent(in) :: d1(:,:,:),n1,n2,n3
integer,intent(out) :: d2(:,:,:)
call MPI_REDUCE(d1,d2,n1*n2*n3,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i4_3d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i8_3d(d1,d2,n1,n2,n3)
integer(kind=8), intent(in) :: d1(:,:,:)
integer,         intent(in) :: n1,n2,n3
integer(kind=8),intent(out) :: d2(:,:,:)
call MPI_REDUCE(d1,d2,n1*n2*n3,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i8_3d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r4_3d(d1,d2,n1,n2,n3)
real,  intent(in)  :: d1(:,:,:)
integer,intent(in) :: n1,n2,n3
real,  intent(out) :: d2(:,:,:)
call MPI_REDUCE(d1,d2,n1*n2*n3,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r4_3d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r8_3d(d1,d2,n1,n2,n3)
real(kind=8), intent(in)  :: d1(:,:,:)
integer,      intent(in)  :: n1,n2,n3
real(kind=8), intent(out) :: d2(:,:,:)
call MPI_REDUCE(d1,d2,n1*n2*n3,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r8_3d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i4_4d(d1,d2,n1,n2,n3,n4)
integer, intent(in) :: d1(:,:,:,:),n1,n2,n3,n4
integer,intent(out) :: d2(:,:,:,:)
call MPI_REDUCE(d1,d2,n1*n2*n3*n4,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i4_4d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_i8_4d(d1,d2,n1,n2,n3,n4)
integer(kind=8), intent(in) :: d1(:,:,:,:)
integer,         intent(in) :: n1,n2,n3,n4
integer(kind=8),intent(out) :: d2(:,:,:,:)
call MPI_REDUCE(d1,d2,n1*n2*n3*n4,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_i8_4d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r4_4d(d1,d2,n1,n2,n3,n4)
real,  intent(in)  :: d1(:,:,:,:)
integer,intent(in) :: n1,n2,n3,n4
real,  intent(out) :: d2(:,:,:,:)
call MPI_REDUCE(d1,d2,n1*n2*n3*n4,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r4_4d

!---------------------------------------------------------------------------------------
subroutine mpi_sum_r8_4d(d1,d2,n1,n2,n3,n4)
real(kind=8), intent(in)  :: d1(:,:,:,:)
integer,      intent(in)  :: n1,n2,n3,n4
real(kind=8), intent(out) :: d2(:,:,:,:)
call MPI_REDUCE(d1,d2,n1*n2*n3*n4,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_sum_r8_4d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i4(d1,d2)
integer, intent(in) :: d1
integer,intent(out) :: d2
call MPI_ALLREDUCE(d1,d2,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i4

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i8(d1,d2)
integer(kind=8), intent(in) :: d1
integer(kind=8),intent(out) :: d2
call MPI_ALLREDUCE(d1,d2,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i8

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r4(d1,d2)
real, intent(in)  :: d1
real, intent(out) :: d2
call MPI_ALLREDUCE(d1,d2,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r4

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r8(d1,d2)
real(kind=8), intent(in)  :: d1
real(kind=8), intent(out) :: d2
call MPI_ALLREDUCE(d1,d2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r8

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i4_1d(d1,d2,n1)
integer, intent(in) :: d1(:),n1
integer,intent(out) :: d2(:)
call MPI_ALLREDUCE(d1,d2,n1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i4_1d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i8_1d(d1,d2,n1)
integer(kind=8), intent(in) :: d1(:)
integer,         intent(in) :: n1
integer(kind=8),intent(out) :: d2(:)
call MPI_ALLREDUCE(d1,d2,n1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i8_1d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r4_1d(d1,d2,n1)
real,  intent(in)  :: d1(:)
integer,intent(in) :: n1
real,  intent(out) :: d2(:)
call MPI_ALLREDUCE(d1,d2,n1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r4_1d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r8_1d(d1,d2,n1)
real(kind=8), intent(in)  :: d1(:)
integer,      intent(in)  :: n1
real(kind=8), intent(out) :: d2(:)
call MPI_ALLREDUCE(d1,d2,n1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r8_1d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i4_2d(d1,d2,n1,n2)
integer, intent(in) :: d1(:,:),n1,n2
integer,intent(out) :: d2(:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i4_2d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i8_2d(d1,d2,n1,n2)
integer(kind=8), intent(in) :: d1(:,:)
integer,         intent(in) :: n1,n2
integer(kind=8),intent(out) :: d2(:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i8_2d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r4_2d(d1,d2,n1,n2)
real,  intent(in)  :: d1(:,:)
integer,intent(in) :: n1,n2
real,  intent(out) :: d2(:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r4_2d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r8_2d(d1,d2,n1,n2)
real(kind=8), intent(in)  :: d1(:,:)
integer,      intent(in)  :: n1,n2
real(kind=8), intent(out) :: d2(:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r8_2d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i4_3d(d1,d2,n1,n2,n3)
integer, intent(in) :: d1(:,:,:),n1,n2,n3
integer,intent(out) :: d2(:,:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2*n3,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i4_3d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i8_3d(d1,d2,n1,n2,n3)
integer(kind=8), intent(in) :: d1(:,:,:)
integer,         intent(in) :: n1,n2,n3
integer(kind=8),intent(out) :: d2(:,:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2*n3,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i8_3d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r4_3d(d1,d2,n1,n2,n3)
real,  intent(in)  :: d1(:,:,:)
integer,intent(in) :: n1,n2,n3
real,  intent(out) :: d2(:,:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2*n3,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r4_3d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r8_3d(d1,d2,n1,n2,n3)
real(kind=8), intent(in)  :: d1(:,:,:)
integer,      intent(in)  :: n1,n2,n3
real(kind=8), intent(out) :: d2(:,:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2*n3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r8_3d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i4_4d(d1,d2,n1,n2,n3,n4)
integer, intent(in) :: d1(:,:,:,:),n1,n2,n3,n4
integer,intent(out) :: d2(:,:,:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2*n3*n4,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i4_4d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_i8_4d(d1,d2,n1,n2,n3,n4)
integer(kind=8), intent(in) :: d1(:,:,:,:)
integer,         intent(in) :: n1,n2,n3,n4
integer(kind=8),intent(out) :: d2(:,:,:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2*n3*n4,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_i8_4d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r4_4d(d1,d2,n1,n2,n3,n4)
real,  intent(in)  :: d1(:,:,:,:)
integer,intent(in) :: n1,n2,n3,n4
real,  intent(out) :: d2(:,:,:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2*n3*n4,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r4_4d

!---------------------------------------------------------------------------------------
subroutine mpi_allsum_r8_4d(d1,d2,n1,n2,n3,n4)
real(kind=8), intent(in)  :: d1(:,:,:,:)
integer,      intent(in)  :: n1,n2,n3,n4
real(kind=8), intent(out) :: d2(:,:,:,:)
call MPI_ALLREDUCE(d1,d2,n1*n2*n3*n4,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine mpi_allsum_r8_4d

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_1i(i1)
integer(kind=2), intent(inout) :: i1
call MPI_BCAST(i1,1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_2i(i1,i2)
integer(kind=2), intent(inout) :: i1,i2
call MPI_BCAST(i1,1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i2,1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_3i(i1,i2,i3)
integer(kind=2), intent(inout) :: i1,i2,i3
call mpi_bcast_i2_2i(i1,i2)
call mpi_bcast_i2_1i(i3)
end subroutine mpi_bcast_i2_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_4i(i1,i2,i3,i4)
integer(kind=2), intent(inout) :: i1,i2,i3,i4
call mpi_bcast_i2_2i(i1,i2)
call mpi_bcast_i2_2i(i3,i4)
end subroutine mpi_bcast_i2_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_1i(i1)
integer, intent(inout) :: i1
call MPI_BCAST(i1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_2i(i1,i2)
integer, intent(inout) :: i1,i2
call MPI_BCAST(i1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_3i(i1,i2,i3)
integer, intent(inout) :: i1,i2,i3
call mpi_bcast_i4_2i(i1,i2)
call mpi_bcast_i4_1i(i3)
end subroutine mpi_bcast_i4_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_4i(i1,i2,i3,i4)
integer, intent(inout) :: i1,i2,i3,i4
call mpi_bcast_i4_2i(i1,i2)
call mpi_bcast_i4_2i(i3,i4)
end subroutine mpi_bcast_i4_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_1i(i1)
integer(kind=8), intent(inout) :: i1
call MPI_BCAST(i1,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_2i(i1,i2)
integer(kind=8), intent(inout) :: i1,i2
call MPI_BCAST(i1,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i2,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_3i(i1,i2,i3)
integer(kind=8), intent(inout) :: i1,i2,i3
call mpi_bcast_i8_2i(i1,i2)
call mpi_bcast_i8_1i(i3)
end subroutine mpi_bcast_i8_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_4i(i1,i2,i3,i4)
integer(kind=8), intent(inout) :: i1,i2,i3,i4
call mpi_bcast_i8_2i(i1,i2)
call mpi_bcast_i8_2i(i3,i4)
end subroutine mpi_bcast_i8_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_1r(r1)
real, intent(inout) :: r1
call MPI_BCAST(r1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_2r(r1,r2)
real, intent(inout) :: r1,r2
call MPI_BCAST(r1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(r2,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_3r(r1,r2,r3)
real, intent(inout) :: r1,r2,r3
call mpi_bcast_r4_2r(r1,r2)
call mpi_bcast_r4_1r(r3)
end subroutine mpi_bcast_r4_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_4r(r1,r2,r3,r4)
real, intent(inout) :: r1,r2,r3,r4
call mpi_bcast_r4_2r(r1,r2)
call mpi_bcast_r4_2r(r3,r4)
end subroutine mpi_bcast_r4_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_1r(r1)
real(kind=8), intent(inout) :: r1
call MPI_BCAST(r1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_2r(r1,r2)
real(kind=8), intent(inout) :: r1,r2
call MPI_BCAST(r1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(r2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_3r(r1,r2,r3)
real(kind=8), intent(inout) :: r1,r2,r3
call mpi_bcast_r8_2r(r1,r2)
call mpi_bcast_r8_1r(r3)
end subroutine mpi_bcast_r8_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_4r(r1,r2,r3,r4)
real(kind=8), intent(inout) :: r1,r2,r3,r4
call mpi_bcast_r8_2r(r1,r2)
call mpi_bcast_r8_2r(r3,r4)
end subroutine mpi_bcast_r8_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_1d_1i(d1,n1)
integer(kind=2),intent(inout) :: d1(:)
integer,        intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_1d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_1d_2i(d1,d2,n1)
integer(kind=2),intent(inout) :: d1(:),d2(:)
integer,        intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_1d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_1d_3i(d1,d2,d3,n1)
integer(kind=2),intent(inout) :: d1(:),d2(:),d3(:)
integer,        intent(in)    :: n1
call mpi_bcast_i2_1d_2i(d1,d2,n1)
call mpi_bcast_i2_1d_1i(d3,n1)
end subroutine mpi_bcast_i2_1d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_1d_4i(d1,d2,d3,d4,n1)
integer(kind=2),intent(inout) :: d1(:),d2(:),d3(:),d4(:)
integer,        intent(in)    :: n1
call mpi_bcast_i2_1d_2i(d1,d2,n1)
call mpi_bcast_i2_1d_2i(d3,d4,n1)
end subroutine mpi_bcast_i2_1d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_2d_1i(d1,n1,n2)
integer(kind=2),intent(inout) :: d1(:,:)
integer,        intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_2d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_2d_2i(d1,d2,n1,n2)
integer(kind=2),intent(inout) :: d1(:,:),d2(:,:)
integer,        intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_2d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_2d_3i(d1,d2,d3,n1,n2)
integer(kind=2),intent(inout) :: d1(:,:),d2(:,:),d3(:,:)
integer,        intent(in)    :: n1,n2
call mpi_bcast_i2_2d_2i(d1,d2,n1,n2)
call mpi_bcast_i2_2d_1i(d3,n1,n2)
end subroutine mpi_bcast_i2_2d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_2d_4i(d1,d2,d3,d4,n1,n2)
integer(kind=2),intent(inout) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:)
integer,        intent(in)    :: n1,n2
call mpi_bcast_i2_2d_2i(d1,d2,n1,n2)
call mpi_bcast_i2_2d_2i(d3,d4,n1,n2)
end subroutine mpi_bcast_i2_2d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_3d_1i(d1,n1,n2,n3)
integer(kind=2),intent(inout) :: d1(:,:,:)
integer,        intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_3d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_3d_2i(d1,d2,n1,n2,n3)
integer(kind=2),intent(inout) :: d1(:,:,:),d2(:,:,:)
integer,        intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_3d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_3d_3i(d1,d2,d3,n1,n2,n3)
integer(kind=2),intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:)
integer,        intent(in)    :: n1,n2,n3
call mpi_bcast_i2_3d_2i(d1,d2,n1,n2,n3)
call mpi_bcast_i2_3d_1i(d3,n1,n2,n3)
end subroutine mpi_bcast_i2_3d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_3d_4i(d1,d2,d3,d4,n1,n2,n3)
integer(kind=2),intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),d4(:,:,:)
integer,        intent(in)    :: n1,n2,n3
call mpi_bcast_i2_3d_2i(d1,d2,n1,n2,n3)
call mpi_bcast_i2_3d_2i(d3,d4,n1,n2,n3)
end subroutine mpi_bcast_i2_3d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_4d_1i(d1,n1,n2,n3,n4)
integer(kind=2),intent(inout) :: d1(:,:,:,:)
integer,        intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_4d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_4d_2i(d1,d2,n1,n2,n3,n4)
integer(kind=2),intent(inout) :: d1(:,:,:,:),d2(:,:,:,:)
integer,        intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3*n4,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i2_4d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_4d_3i(d1,d2,d3,n1,n2,n3,n4)
integer(kind=2),intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:)
integer,        intent(in)    :: n1,n2,n3,n4
call mpi_bcast_i2_4d_2i(d1,d2,n1,n2,n3,n4)
call mpi_bcast_i2_4d_1i(d3,n1,n2,n3,n4)
end subroutine mpi_bcast_i2_4d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i2_4d_4i(d1,d2,d3,d4,n1,n2,n3,n4)
integer(kind=2),intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:),d4(:,:,:,:)
integer,        intent(in)    :: n1,n2,n3,n4
call mpi_bcast_i2_4d_2i(d1,d2,n1,n2,n3,n4)
call mpi_bcast_i2_4d_2i(d3,d4,n1,n2,n3,n4)
end subroutine mpi_bcast_i2_4d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_1d_1i(d1,n1)
integer,intent(inout) :: d1(:)
integer,intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_1d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_1d_2i(d1,d2,n1)
integer,intent(inout) :: d1(:),d2(:)
integer,intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_1d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_1d_3i(d1,d2,d3,n1)
integer,intent(inout) :: d1(:),d2(:),d3(:)
integer,intent(in)    :: n1
call mpi_bcast_i4_1d_2i(d1,d2,n1)
call mpi_bcast_i4_1d_1i(d3,n1)
end subroutine mpi_bcast_i4_1d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_1d_4i(d1,d2,d3,d4,n1)
integer,intent(inout) :: d1(:),d2(:),d3(:),d4(:)
integer,intent(in)    :: n1
call mpi_bcast_i4_1d_2i(d1,d2,n1)
call mpi_bcast_i4_1d_2i(d3,d4,n1)
end subroutine mpi_bcast_i4_1d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_2d_1i(d1,n1,n2)
integer,intent(inout) :: d1(:,:)
integer,intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_2d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_2d_2i(d1,d2,n1,n2)
integer,intent(inout) :: d1(:,:),d2(:,:)
integer,intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_2d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_2d_3i(d1,d2,d3,n1,n2)
integer,intent(inout) :: d1(:,:),d2(:,:),d3(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_i4_2d_2i(d1,d2,n1,n2)
call mpi_bcast_i4_2d_1i(d3,n1,n2)
end subroutine mpi_bcast_i4_2d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_2d_4i(d1,d2,d3,d4,n1,n2)
integer,intent(inout) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_i4_2d_2i(d1,d2,n1,n2)
call mpi_bcast_i4_2d_2i(d3,d4,n1,n2)
end subroutine mpi_bcast_i4_2d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_3d_1i(d1,n1,n2,n3)
integer,intent(inout) :: d1(:,:,:)
integer,intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_3d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_3d_2i(d1,d2,n1,n2,n3)
integer,intent(inout) :: d1(:,:,:),d2(:,:,:)
integer,intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_3d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_3d_3i(d1,d2,d3,n1,n2,n3)
integer,intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:)
integer,intent(in)    :: n1,n2,n3
call mpi_bcast_i4_3d_2i(d1,d2,n1,n2,n3)
call mpi_bcast_i4_3d_1i(d3,n1,n2,n3)
end subroutine mpi_bcast_i4_3d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_3d_4i(d1,d2,d3,d4,n1,n2,n3)
integer,intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),d4(:,:,:)
integer,intent(in)    :: n1,n2,n3
call mpi_bcast_i4_3d_2i(d1,d2,n1,n2,n3)
call mpi_bcast_i4_3d_2i(d3,d4,n1,n2,n3)
end subroutine mpi_bcast_i4_3d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_4d_1i(d1,n1,n2,n3,n4)
integer,intent(inout) :: d1(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_4d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_4d_2i(d1,d2,n1,n2,n3,n4)
integer,intent(inout) :: d1(:,:,:,:),d2(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3*n4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i4_4d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_4d_3i(d1,d2,d3,n1,n2,n3,n4)
integer,intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call mpi_bcast_i4_4d_2i(d1,d2,n1,n2,n3,n4)
call mpi_bcast_i4_4d_1i(d3,n1,n2,n3,n4)
end subroutine mpi_bcast_i4_4d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i4_4d_4i(d1,d2,d3,d4,n1,n2,n3,n4)
integer,intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:),d4(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call mpi_bcast_i4_4d_2i(d1,d2,n1,n2,n3,n4)
call mpi_bcast_i4_4d_2i(d3,d4,n1,n2,n3,n4)
end subroutine mpi_bcast_i4_4d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_1d_1i(d1,n1)
integer(kind=8),intent(inout) :: d1(:)
integer,        intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_1d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_1d_2i(d1,d2,n1)
integer(kind=8),intent(inout) :: d1(:),d2(:)
integer,        intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_1d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_1d_3i(d1,d2,d3,n1)
integer(kind=8),intent(inout) :: d1(:),d2(:),d3(:)
integer,        intent(in)    :: n1
call mpi_bcast_i8_1d_2i(d1,d2,n1)
call mpi_bcast_i8_1d_1i(d3,n1)
end subroutine mpi_bcast_i8_1d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_1d_4i(d1,d2,d3,d4,n1)
integer(kind=8),intent(inout) :: d1(:),d2(:),d3(:),d4(:)
integer,        intent(in)    :: n1
call mpi_bcast_i8_1d_2i(d1,d2,n1)
call mpi_bcast_i8_1d_2i(d3,d4,n1)
end subroutine mpi_bcast_i8_1d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_2d_1i(d1,n1,n2)
integer(kind=8),intent(inout) :: d1(:,:)
integer,        intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_2d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_2d_2i(d1,d2,n1,n2)
integer(kind=8),intent(inout) :: d1(:,:),d2(:,:)
integer,        intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_2d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_2d_3i(d1,d2,d3,n1,n2)
integer(kind=8),intent(inout) :: d1(:,:),d2(:,:),d3(:,:)
integer,        intent(in)    :: n1,n2
call mpi_bcast_i8_2d_2i(d1,d2,n1,n2)
call mpi_bcast_i8_2d_1i(d3,n1,n2)
end subroutine mpi_bcast_i8_2d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_2d_4i(d1,d2,d3,d4,n1,n2)
integer(kind=8),intent(inout) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:)
integer,        intent(in)    :: n1,n2
call mpi_bcast_i8_2d_2i(d1,d2,n1,n2)
call mpi_bcast_i8_2d_2i(d3,d4,n1,n2)
end subroutine mpi_bcast_i8_2d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_3d_1i(d1,n1,n2,n3)
integer(kind=8),intent(inout) :: d1(:,:,:)
integer,        intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_3d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_3d_2i(d1,d2,n1,n2,n3)
integer(kind=8),intent(inout) :: d1(:,:,:),d2(:,:,:)
integer,        intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_3d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_3d_3i(d1,d2,d3,n1,n2,n3)
integer(kind=8),intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:)
integer,        intent(in)    :: n1,n2,n3
call mpi_bcast_i8_3d_2i(d1,d2,n1,n2,n3)
call mpi_bcast_i8_3d_1i(d3,n1,n2,n3)
end subroutine mpi_bcast_i8_3d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_3d_4i(d1,d2,d3,d4,n1,n2,n3)
integer(kind=8),intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),d4(:,:,:)
integer,        intent(in)    :: n1,n2,n3
call mpi_bcast_i8_3d_2i(d1,d2,n1,n2,n3)
call mpi_bcast_i8_3d_2i(d3,d4,n1,n2,n3)
end subroutine mpi_bcast_i8_3d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_4d_1i(d1,n1,n2,n3,n4)
integer(kind=8),intent(inout) :: d1(:,:,:,:)
integer,        intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_4d_1i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_4d_2i(d1,d2,n1,n2,n3,n4)
integer(kind=8),intent(inout) :: d1(:,:,:,:),d2(:,:,:,:)
integer,        intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3*n4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_i8_4d_2i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_4d_3i(d1,d2,d3,n1,n2,n3,n4)
integer(kind=8),intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:)
integer,        intent(in)    :: n1,n2,n3,n4
call mpi_bcast_i8_4d_2i(d1,d2,n1,n2,n3,n4)
call mpi_bcast_i8_4d_1i(d3,n1,n2,n3,n4)
end subroutine mpi_bcast_i8_4d_3i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_i8_4d_4i(d1,d2,d3,d4,n1,n2,n3,n4)
integer(kind=8),intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:),d4(:,:,:,:)
integer,        intent(in)    :: n1,n2,n3,n4
call mpi_bcast_i8_4d_2i(d1,d2,n1,n2,n3,n4)
call mpi_bcast_i8_4d_2i(d3,d4,n1,n2,n3,n4)
end subroutine mpi_bcast_i8_4d_4i

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_1d_1r(d1,n1)
real,   intent(inout) :: d1(:)
integer,intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_1d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_1d_2r(d1,d2,n1)
real,   intent(inout) :: d1(:),d2(:)
integer,intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_1d_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_1d_3r(d1,d2,d3,n1)
real,   intent(inout) :: d1(:),d2(:),d3(:)
integer,intent(in)    :: n1
call mpi_bcast_r4_1d_2r(d1,d2,n1)
call mpi_bcast_r4_1d_1r(d3,n1)
end subroutine mpi_bcast_r4_1d_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_1d_4r(d1,d2,d3,d4,n1)
real,   intent(inout) :: d1(:),d2(:),d3(:),d4(:)
integer,intent(in)    :: n1
call mpi_bcast_r4_1d_2r(d1,d2,n1)
call mpi_bcast_r4_1d_2r(d3,d4,n1)
end subroutine mpi_bcast_r4_1d_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_1d_5r(d1,d2,d3,d4,d5,n1)
real,   intent(inout) :: d1(:),d2(:),d3(:),d4(:),d5(:)
integer,intent(in)    :: n1
call mpi_bcast_r4_1d_3r(d1,d2,d3,n1)
call mpi_bcast_r4_1d_2r(d4,d5,n1)
end subroutine mpi_bcast_r4_1d_5r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_1d_7r(d1,d2,d3,d4,d5,d6,d7,n1)
real,   intent(inout) :: d1(:),d2(:),d3(:),d4(:),d5(:),d6(:),d7(:)
integer,intent(in)    :: n1
call mpi_bcast_r4_1d_4r(d1,d2,d3,d4,n1)
call mpi_bcast_r4_1d_3r(d5,d6,d7,n1)
end subroutine mpi_bcast_r4_1d_7r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_2d_1r(d1,n1,n2)
real,   intent(inout) :: d1(:,:)
integer,intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_2d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_2d_2r(d1,d2,n1,n2)
real,   intent(inout) :: d1(:,:),d2(:,:)
integer,intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_2d_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_2d_3r(d1,d2,d3,n1,n2)
real,   intent(inout) :: d1(:,:),d2(:,:),d3(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_r4_2d_2r(d1,d2,n1,n2)
call mpi_bcast_r4_2d_1r(d3,n1,n2)
end subroutine mpi_bcast_r4_2d_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_2d_4r(d1,d2,d3,d4,n1,n2)
real,   intent(inout) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_r4_2d_2r(d1,d2,n1,n2)
call mpi_bcast_r4_2d_2r(d3,d4,n1,n2)
end subroutine mpi_bcast_r4_2d_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_2d_5r(d1,d2,d3,d4,d5,n1,n2)
real,   intent(inout) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),d5(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_r4_2d_3r(d1,d2,d3,n1,n2)
call mpi_bcast_r4_2d_2r(d4,d5,n1,n2)
end subroutine mpi_bcast_r4_2d_5r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_2d_6r(d1,d2,d3,d4,d5,d6,n1,n2)
real,   intent(inout) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),d5(:,:),d6(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_r4_2d_3r(d1,d2,d3,n1,n2)
call mpi_bcast_r4_2d_3r(d4,d5,d6,n1,n2)
end subroutine mpi_bcast_r4_2d_6r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_2d_7r(d1,d2,d3,d4,d5,d6,d7,n1,n2)
real,   intent(inout) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),d5(:,:),d6(:,:),d7(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_r4_2d_4r(d1,d2,d3,d4,n1,n2)
call mpi_bcast_r4_2d_3r(d5,d6,d7,n1,n2)
end subroutine mpi_bcast_r4_2d_7r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_3d_1r(d1,n1,n2,n3)
real,   intent(inout) :: d1(:,:,:)
integer,intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_3d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_3d_2r(d1,d2,n1,n2,n3)
real,   intent(inout) :: d1(:,:,:),d2(:,:,:)
integer,intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_3d_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_3d_3r(d1,d2,d3,n1,n2,n3)
real,   intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:)
integer,intent(in)    :: n1,n2,n3
call mpi_bcast_r4_3d_2r(d1,d2,n1,n2,n3)
call mpi_bcast_r4_3d_1r(d3,n1,n2,n3)
end subroutine mpi_bcast_r4_3d_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_3d_4r(d1,d2,d3,d4,n1,n2,n3)
real,   intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),d4(:,:,:)
integer,intent(in)    :: n1,n2,n3
call mpi_bcast_r4_3d_2r(d1,d2,n1,n2,n3)
call mpi_bcast_r4_3d_2r(d3,d4,n1,n2,n3)
end subroutine mpi_bcast_r4_3d_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_4d_1r(d1,n1,n2,n3,n4)
real,   intent(inout) :: d1(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_4d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_4d_2r(d1,d2,n1,n2,n3,n4)
real,   intent(inout) :: d1(:,:,:,:),d2(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3*n4,MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r4_4d_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_4d_3r(d1,d2,d3,n1,n2,n3,n4)
real,   intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call mpi_bcast_r4_4d_2r(d1,d2,n1,n2,n3,n4)
call mpi_bcast_r4_4d_1r(d3,n1,n2,n3,n4)
end subroutine mpi_bcast_r4_4d_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r4_4d_4r(d1,d2,d3,d4,n1,n2,n3,n4)
real,   intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:),d4(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call mpi_bcast_r4_4d_2r(d1,d2,n1,n2,n3,n4)
call mpi_bcast_r4_4d_2r(d3,d4,n1,n2,n3,n4)
end subroutine mpi_bcast_r4_4d_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_1d_1r(d1,n1)
real(kind=8),intent(inout) :: d1(:)
integer,intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_1d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_1d_2r(d1,d2,n1)
real(kind=8),   intent(inout) :: d1(:),d2(:)
integer,intent(in)    :: n1
call MPI_BCAST(d1,n1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_1d_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_1d_3r(d1,d2,d3,n1)
real(kind=8),   intent(inout) :: d1(:),d2(:),d3(:)
integer,intent(in)    :: n1
call mpi_bcast_r8_1d_2r(d1,d2,n1)
call mpi_bcast_r8_1d_1r(d3,n1)
end subroutine mpi_bcast_r8_1d_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_1d_4r(d1,d2,d3,d4,n1)
real(kind=8),   intent(inout) :: d1(:),d2(:),d3(:),d4(:)
integer,intent(in)    :: n1
call mpi_bcast_r8_1d_2r(d1,d2,n1)
call mpi_bcast_r8_1d_2r(d3,d4,n1)
end subroutine mpi_bcast_r8_1d_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_2d_1r(d1,n1,n2)
real(kind=8),   intent(inout) :: d1(:,:)
integer,intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_2d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_2d_2r(d1,d2,n1,n2)
real(kind=8),   intent(inout) :: d1(:,:),d2(:,:)
integer,intent(in)    :: n1,n2
call MPI_BCAST(d1,n1*n2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_2d_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_2d_3r(d1,d2,d3,n1,n2)
real(kind=8),   intent(inout) :: d1(:,:),d2(:,:),d3(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_r8_2d_2r(d1,d2,n1,n2)
call mpi_bcast_r8_2d_1r(d3,n1,n2)
end subroutine mpi_bcast_r8_2d_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_2d_4r(d1,d2,d3,d4,n1,n2)
real(kind=8),   intent(inout) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:)
integer,intent(in)    :: n1,n2
call mpi_bcast_r8_2d_2r(d1,d2,n1,n2)
call mpi_bcast_r8_2d_2r(d3,d4,n1,n2)
end subroutine mpi_bcast_r8_2d_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_3d_1r(d1,n1,n2,n3)
real(kind=8),   intent(inout) :: d1(:,:,:)
integer,intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_3d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_3d_2r(d1,d2,n1,n2,n3)
real(kind=8),   intent(inout) :: d1(:,:,:),d2(:,:,:)
integer,intent(in)    :: n1,n2,n3
call MPI_BCAST(d1,n1*n2*n3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_3d_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_3d_3r(d1,d2,d3,n1,n2,n3)
real(kind=8),   intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:)
integer,intent(in)    :: n1,n2,n3
call mpi_bcast_r8_3d_2r(d1,d2,n1,n2,n3)
call mpi_bcast_r8_3d_1r(d3,n1,n2,n3)
end subroutine mpi_bcast_r8_3d_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_3d_4r(d1,d2,d3,d4,n1,n2,n3)
real(kind=8),   intent(inout) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),d4(:,:,:)
integer,intent(in)    :: n1,n2,n3
call mpi_bcast_r8_3d_2r(d1,d2,n1,n2,n3)
call mpi_bcast_r8_3d_2r(d3,d4,n1,n2,n3)
end subroutine mpi_bcast_r8_3d_4r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_4d_1r(d1,n1,n2,n3,n4)
real(kind=8),   intent(inout) :: d1(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_4d_1r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_4d_2r(d1,d2,n1,n2,n3,n4)
real(kind=8),   intent(inout) :: d1(:,:,:,:),d2(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call MPI_BCAST(d1,n1*n2*n3*n4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(d2,n1*n2*n3*n4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine mpi_bcast_r8_4d_2r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_4d_3r(d1,d2,d3,n1,n2,n3,n4)
real(kind=8),   intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call mpi_bcast_r8_4d_2r(d1,d2,n1,n2,n3,n4)
call mpi_bcast_r8_4d_1r(d3,n1,n2,n3,n4)
end subroutine mpi_bcast_r8_4d_3r

!---------------------------------------------------------------------------------------
subroutine mpi_bcast_r8_4d_4r(d1,d2,d3,d4,n1,n2,n3,n4)
real(kind=8),   intent(inout) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:),d4(:,:,:,:)
integer,intent(in)    :: n1,n2,n3,n4
call mpi_bcast_r8_4d_2r(d1,d2,n1,n2,n3,n4)
call mpi_bcast_r8_4d_2r(d3,d4,n1,n2,n3,n4)
end subroutine mpi_bcast_r8_4d_4r

!---------------------------------------------------------------------------------------
subroutine message_c(rank0,msg)
integer,      intent(in) :: rank0
character(*), intent(in) :: msg
if (rank == rank0) then 
   write(*,*) msg
   !call flush(6)
endif 
end subroutine message_c

!---------------------------------------------------------------------------------------

subroutine message_c2(rank0,msg1,msg2)
integer,      intent(in) :: rank0
character(*), intent(in) :: msg1,msg2
if (rank == rank0) then 
   write(*,*) msg1," ",msg2
   !call flush(6)
endif 
end subroutine message_c2

!---------------------------------------------------------------------------------------

subroutine message_c1i1(rank0,char1,int1)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1
integer,      intent(in) :: int1
if (rank == rank0) then 
   write(*,*) char1,int1
   !call flush(6)
endif 
end subroutine message_c1i1

!---------------------------------------------------------------------------------------

subroutine message_c1i1c1(rank0,char1,int1,char2)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1,char2
integer,      intent(in) :: int1
if (rank == rank0) then 
   write(*,*) char1,int1,char2
   !call flush(6)
endif 
end subroutine message_c1i1c1

!---------------------------------------------------------------------------------------

subroutine message_cicdp(rank0,char1,int1,char2,dp1)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1,char2
integer,      intent(in) :: int1
double precision, intent(in) :: dp1
if (rank == rank0) then 
   write(*,*) char1,int1,char2,dp1
   !call flush(6)
endif 
end subroutine message_cicdp

!---------------------------------------------------------------------------------------

subroutine message_cicr(rank0,char1,int1,char2,r1)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1,char2
integer,      intent(in) :: int1
real,         intent(in) :: r1
if (rank == rank0) then 
   write(*,*) char1,int1,char2,r1
   !call flush(6)
endif 
end subroutine message_cicr

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine message_c1long1(rank0,char1,int1)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1
integer(8),   intent(in) :: int1
if (rank == rank0) then 
   write(*,*) char1,int1
   !call flush(6)
endif 
end subroutine message_c1long1

!---------------------------------------------------------------------------------------

subroutine message_c1long2(rank0,char1,int1,int2)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1
integer(8),   intent(in) :: int1,int2
if (rank == rank0) then 
   write(*,*) char1,int1," ",int2
   !call flush(6)
endif 
end subroutine message_c1long2

!---------------------------------------------------------------------------------------

subroutine message_c2i2(rank0,char1,int1,char2,int2)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1, char2
integer,      intent(in) :: int1, int2
if (rank == rank0) then 
   write(*,*) char1,int1,char2,int2
   !call flush(6)
endif 
end subroutine message_c2i2

!---------------------------------------------------------------------------------------

subroutine message_c3i3(rank0,char1,int1,char2,int2,char3,int3)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1, char2, char3
integer,      intent(in) :: int1, int2, int3
if (rank == rank0) then 
   write(*,*) char1,int1,char2,int2,char3,int3
   !call flush(6)
endif 
end subroutine message_c3i3

!---------------------------------------------------------------------------------------

subroutine message_c1r1(rank0,char1,real1)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1
real,         intent(in) :: real1
if (rank == rank0) then 
   write(*,*) char1,real1
   !call flush(6)
endif 
end subroutine message_c1r1

!---------------------------------------------------------------------------------------

subroutine message_c1r1c1dp1(rank0,char1,real1,char2,dp1)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1,char2
real,         intent(in) :: real1
double precision,intent(in) :: dp1
if (rank == rank0) then 
   write(*,*) char1,real1,char2,dp1
   !call flush(6)
endif 
end subroutine message_c1r1c1dp1

!---------------------------------------------------------------------------------------

subroutine message_c2r2(rank0,char1,real1,char2,real2)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1, char2
real,         intent(in) :: real1, real2
if (rank == rank0) then 
   write(*,*) char1,real1,char2,real2
   !call flush(6)
endif 
end subroutine message_c2r2

!---------------------------------------------------------------------------------------

subroutine message_c3r3(rank0,char1,real1,char2,real2,char3,real3)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1, char2, char3
real,         intent(in) :: real1, real2, real3
if (rank == rank0) then 
   write(*,*) char1,real1,char2,real2,char3,real3
   !call flush(6)
endif 
end subroutine message_c3r3

!---------------------------------------------------------------------------------------

subroutine message_c1d1(rank0,char1,double1)
integer,           intent(in) :: rank0
character(*),      intent(in) :: char1
double precision,  intent(in) :: double1
if (rank == rank0) then 
   write(*,*) char1,double1
   !call flush(6)
endif 
end subroutine message_c1d1

!---------------------------------------------------------------------------------------

subroutine message_c2d2(rank0,char1,double1,char2,double2)
integer,          intent(in) :: rank0
character(*),     intent(in) :: char1, char2
double precision, intent(in) :: double1, double2
if (rank == rank0) then 
   write(*,*) char1,double1,char2,double2
   !call flush(6)
endif 
end subroutine message_c2d2

!---------------------------------------------------------------------------------------

subroutine message_c3d3(rank0,char1,double1,char2,double2,char3,double3)
integer,          intent(in) :: rank0
character(*),     intent(in) :: char1, char2, char3
double precision, intent(in) :: double1, double2, double3
if (rank == rank0) then 
   write(*,*) char1,double1,char2,double2,char3,double3
   !call flush(6)
endif 
end subroutine message_c3d3

!---------------------------------------------------------------------------------------

subroutine message_c1l1(rank0,char1,logic1)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1
logical,      intent(in) :: logic1
if (rank == rank0) then 
   write(*,*) char1,logic1
   !call flush(6)
endif 
end subroutine message_c1l1

!---------------------------------------------------------------------------------------

subroutine message_c2l2(rank0,char1,logic1,char2,logic2)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1, char2
logical,      intent(in) :: logic1, logic2
if (rank == rank0) then 
   write(*,*) char1,logic1,char2,logic2
   !call flush(6)
endif 
end subroutine message_c2l2

!---------------------------------------------------------------------------------------

subroutine message_c3l3(rank0,char1,logic1,char2,logic2,char3,logic3)
integer,      intent(in) :: rank0
character(*), intent(in) :: char1, char2, char3
logical,      intent(in) :: logic1, logic2, logic3
if (rank == rank0) then 
   write(*,*) char1,logic1,char2,logic2,char3,logic3
   !call flush(6)
endif 
end subroutine message_c3l3

!---------------------------------------------------------------------------------------

subroutine message_i3(rank0,int1,int2,int3)
integer,      intent(in) :: rank0, int1, int2, int3
if (rank == rank0) then 
   write(*,*) int1, int2, int3
   !call flush(6)
endif 
end subroutine message_i3

!---------------------------------------------------------------------------------------

subroutine message_r3(rank0,real1,real2,real3)
integer,      intent(in) :: rank0
real,         intent(in) :: real1, real2, real3
if (rank == rank0) then 
   write(*,*) real1, real2, real3
   !call flush(6)
endif 
end subroutine message_r3

!---------------------------------------------------------------------------------------

subroutine message_all_c(msg)

character(*), intent(in) :: msg

write(*,*) msg
!call flush(6)

end subroutine message_all_c

!---------------------------------------------------------------------------------------

subroutine message_all_c2(msg1,msg2)

character(*), intent(in) :: msg1,msg2

write(*,*) msg1,msg2
!call flush(6)

end subroutine message_all_c2

!---------------------------------------------------------------------------------------

subroutine message_all_c1i1(char1,int1)
character(*), intent(in) :: char1
integer,      intent(in) :: int1
write(*,*) char1,int1
!call flush(6)
end subroutine message_all_c1i1

!---------------------------------------------------------------------------------------

subroutine message_all_c1long1(char1,int1)
character(*), intent(in) :: char1
integer(8),   intent(in) :: int1
write(*,*) char1,int1
!call flush(6)
end subroutine message_all_c1long1

!---------------------------------------------------------------------------------------

subroutine message_all_c1long2(char1,int1,int2)
character(*), intent(in) :: char1
integer(8),   intent(in) :: int1,int2
write(*,*) char1,int1," ",int2
!call flush(6)
end subroutine message_all_c1long2

!---------------------------------------------------------------------------------------



subroutine message_all_c2i2(char1,int1,char2,int2)
character(*), intent(in) :: char1, char2
integer,      intent(in) :: int1, int2
write(*,*) char1,int1,char2,int2
!call flush(6)
end subroutine message_all_c2i2

!---------------------------------------------------------------------------------------

subroutine message_all_c3i3(char1,int1,char2,int2,char3,int3)
character(*), intent(in) :: char1, char2, char3
integer,      intent(in) :: int1, int2, int3
write(*,*) char1,int1,char2,int2,char3,int3
!call flush(6)
end subroutine message_all_c3i3

!---------------------------------------------------------------------------------------

subroutine message_all_c1r1(char1,real1)
character(*), intent(in) :: char1
real,         intent(in) :: real1
write(*,*) char1,real1
!call flush(6)
end subroutine message_all_c1r1

!---------------------------------------------------------------------------------------

subroutine message_all_c2r2(char1,real1,char2,real2)
character(*), intent(in) :: char1, char2
real,         intent(in) :: real1, real2
write(*,*) char1,real1,char2,real2
!call flush(6)
end subroutine message_all_c2r2

!---------------------------------------------------------------------------------------

subroutine message_all_c3r3(char1,real1,char2,real2,char3,real3)
character(*), intent(in) :: char1, char2, char3
real,         intent(in) :: real1, real2, real3
write(*,*) char1,real1,char2,real2,char3,real3
!call flush(6)
end subroutine message_all_c3r3

!---------------------------------------------------------------------------------------

subroutine message_all_c1d1(char1,double1)
character(*),      intent(in) :: char1
double precision,  intent(in) :: double1
write(*,*) char1,double1
!call flush(6)
end subroutine message_all_c1d1

!---------------------------------------------------------------------------------------

subroutine message_all_c2d2(char1,double1,char2,double2)
character(*),     intent(in) :: char1, char2
double precision, intent(in) :: double1, double2
write(*,*) char1,double1,char2,double2
!call flush(6)
end subroutine message_all_c2d2

!---------------------------------------------------------------------------------------

subroutine message_all_c3d3(char1,double1,char2,double2,char3,double3)
character(*),     intent(in) :: char1, char2, char3
double precision, intent(in) :: double1, double2, double3
write(*,*) char1,double1,char2,double2,char3,double3
!call flush(6)
end subroutine message_all_c3d3

!---------------------------------------------------------------------------------------

subroutine message_all_c1l1(char1,logic1)
character(*),      intent(in) :: char1
logical,           intent(in) :: logic1
write(*,*) char1,logic1
!call flush(6)
end subroutine message_all_c1l1

!---------------------------------------------------------------------------------------

subroutine message_all_c2l2(char1,logic1,char2,logic2)
character(*),     intent(in) :: char1, char2
logical,          intent(in) :: logic1, logic2
write(*,*) char1,logic1,char2,logic2
!call flush(6)
end subroutine message_all_c2l2

!---------------------------------------------------------------------------------------

subroutine message_all_c3l3(char1,logic1,char2,logic2,char3,logic3)
character(*),     intent(in) :: char1, char2, char3
logical,          intent(in) :: logic1, logic2, logic3
write(*,*) char1,logic1,char2,logic2,char3,logic3
!call flush(6)
end subroutine message_all_c3l3

!---------------------------------------------------------------------------------------

subroutine message_all_c_i_c_c_i_c(c1,i1,c2,c3,i2,c4)
character(len=*),intent(in) :: c1, c2, c3, c4
integer,         intent(in) :: i1, i2
write(*,*)trim(c1)," ",i1," ",trim(c2)," ",trim(c3)," ",i2," ",trim(c4)
end subroutine message_all_c_i_c_c_i_c
!---------------------------------------------------------------------------------------

! Following subroutines are written by Bowen Guo
subroutine map_model(nx_s,nx_e,nz_s,nz_e,np,map_x,map_z)
integer, intent(in)    :: nx_s,nx_e,nz_s,nz_e
integer, intent(out)   :: np,map_x(:),map_z(:)
integer                :: ix,iz,i,npx,npz
! The target area has to be a square
npx=nx_e-nx_s+1
npz=nz_e-nz_s+1
np=npx*npz
if (rank.eq.0) write(*,*) 'The total number of points ==',np
!call allocate_and_initial(map_x,np)
!call allocate_and_initial(map_z,np)
do ix=1,npx
   do iz=1,npz
      i=(ix-1)*npz+iz
      map_x(i)=nx_s+ix-1
      map_z(i)=nz_s+iz-1
   enddo
enddo
end subroutine map_model
!====================================================================================
subroutine assign_model(np,np_me,ip_base)
integer, intent(in)   :: np
integer, intent(out)  :: np_me,ip_base
integer               :: np_left

! number of points in the model less than nsize
if (np<nsize) then
   if (rank.eq.0) then
      write(*,*) "WARNING: # of points ",np,"is less than # of cores",nsize
   endif
   if (rank<=np-1) then
      np_me=1
      ip_base=rank
   else
      np_me=0
   endif
   return
endif

! Only one CPU case
if (nsize == 1) then
   np_me=np
   ip_base=0
   return
endif

! Normal case
np_me = int(np/nsize)
np_left = np-np_me*nsize


if (rank.le.np_left-1) then
   np_me=np_me+1
   ip_base=rank*np_me
else
   ip_base=rank*np_me+np_left
endif

end subroutine assign_model
!==================================================================================
























end module module_sf90_mpi
