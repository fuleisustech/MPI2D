module module_array
use module_global, only : MEM_UNIT
use module_array_compose
implicit none

interface allocate_and_initial
   module procedure allocate_and_initial_1d_1i
   module procedure allocate_and_initial_1d_2i
   module procedure allocate_and_initial_1d_3i
   module procedure allocate_and_initial_1d_4i
   module procedure allocate_and_initial_1d_6i
   module procedure allocate_and_initial_1d_1r
   module procedure allocate_and_initial_1d_2r
   module procedure allocate_and_initial_1d_3r
   module procedure allocate_and_initial_1d_4r
   module procedure allocate_and_initial_1d_5r
   module procedure allocate_and_initial_1d_6r
   module procedure allocate_and_initial_1d_7r
   module procedure allocate_and_initial_1d_1c
   module procedure allocate_and_initial_1d_2c
   module procedure allocate_and_initial_1d_3c
   module procedure allocate_and_initial_2d_1c ! Added by Bowen Guo
   module procedure allocate_and_initial_2d_1i
   module procedure allocate_and_initial_2d_2i
   module procedure allocate_and_initial_2d_1r
   module procedure allocate_and_initial_2d_2r
   module procedure allocate_and_initial_2d_3r
   module procedure allocate_and_initial_2d_4r
   module procedure allocate_and_initial_2d_5r
   module procedure allocate_and_initial_2d_6r
   module procedure allocate_and_initial_2d_7r
   module procedure allocate_and_initial_2d_9r
   module procedure allocate_and_initial_2d_12r
   module procedure allocate_and_initial_2d_15r
   module procedure allocate_and_initial_2d_18r
   module procedure allocate_and_initial_1d_1dp
   module procedure allocate_and_initial_2d_1dp
   module procedure allocate_and_initial_2d_2dp
   module procedure allocate_and_initial_3d_1i  ! modified by Bowen Guo
   module procedure allocate_and_initial_3d_1r
   module procedure allocate_and_initial_3d_2r
   module procedure allocate_and_initial_3d_3r
   module procedure allocate_and_initial_3d_4r
   module procedure allocate_and_initial_3d_5r
   module procedure allocate_and_initial_3d_6r
   module procedure allocate_and_initial_4d_1r
   module procedure allocate_and_initial_4d_2r
   module procedure allocate_and_initial_4d_3r
   module procedure allocate_and_initial_4d_6r
end interface allocate_and_initial

interface deallocate_and_free
   module procedure deallocate_and_free_1d_1i
   module procedure deallocate_and_free_1d_2i
   module procedure deallocate_and_free_1d_3i
   module procedure deallocate_and_free_1d_4i
   module procedure deallocate_and_free_1d_6i
   module procedure deallocate_and_free_1d_1r
   module procedure deallocate_and_free_1d_2r
   module procedure deallocate_and_free_1d_3r
   module procedure deallocate_and_free_1d_4r
   module procedure deallocate_and_free_1d_5r
   module procedure deallocate_and_free_1d_6r
   module procedure deallocate_and_free_1d_7r
   module procedure deallocate_and_free_1d_1c
   module procedure deallocate_and_free_1d_2c
   module procedure deallocate_and_free_2d_1c  ! added by Bowen Guo
   module procedure deallocate_and_free_2d_1i
   module procedure deallocate_and_free_2d_2i
   module procedure deallocate_and_free_2d_1r
   module procedure deallocate_and_free_2d_2r
   module procedure deallocate_and_free_2d_3r
   module procedure deallocate_and_free_2d_4r
   module procedure deallocate_and_free_2d_5r
   module procedure deallocate_and_free_2d_6r
   module procedure deallocate_and_free_2d_7r
   module procedure deallocate_and_free_2d_9r
   module procedure deallocate_and_free_2d_12r
   module procedure deallocate_and_free_2d_15r
   module procedure deallocate_and_free_2d_18r
   module procedure deallocate_and_free_1d_1dp
   module procedure deallocate_and_free_2d_1dp
   module procedure deallocate_and_free_2d_2dp
   module procedure deallocate_and_free_3d_1i  ! Added by Bowen Guo
   module procedure deallocate_and_free_3d_1r
   module procedure deallocate_and_free_3d_2r
   module procedure deallocate_and_free_3d_3r
   module procedure deallocate_and_free_3d_4r
   module procedure deallocate_and_free_3d_5r
   module procedure deallocate_and_free_3d_6r
   module procedure deallocate_and_free_4d_1r
   module procedure deallocate_and_free_4d_2r
   module procedure deallocate_and_free_4d_3r
   module procedure deallocate_and_free_4d_6r
end interface deallocate_and_free

interface initial
   module procedure initial_1d_1i
   module procedure initial_1d_2i
   module procedure initial_1d_3i
   module procedure initial_1d_1r
   module procedure initial_1d_2r
   module procedure initial_1d_3r
   module procedure initial_1d_6r
   module procedure initial_2d_1r
   module procedure initial_2d_2r
   module procedure initial_2d_3r
   module procedure initial_2d_4r
   module procedure initial_2d_5r
   module procedure initial_2d_6r
   module procedure initial_3d_1r
   module procedure initial_3d_2r
   module procedure initial_3d_3r
   module procedure initial_3d_4r
   module procedure initial_3d_6r
   module procedure initial_4d_1r
   module procedure initial_4d_2r
   module procedure initial_4d_3r
end interface initial

interface copy_array
   module procedure copy_1d_int_array1
   module procedure copy_1d_int_array2
   module procedure copy_1d_int_array3
   module procedure copy_2d_int_array1
   module procedure copy_2d_int_array2
   module procedure copy_2d_int_array3
   module procedure copy_3d_int_array1
   module procedure copy_3d_int_array2
   module procedure copy_3d_int_array3
   module procedure copy_1d_real_array1
   module procedure copy_1d_real_array2
   module procedure copy_1d_real_array3
   module procedure copy_2d_real_array1
   module procedure copy_2d_real_array2
   module procedure copy_2d_real_array3
   module procedure copy_3d_real_array1
   module procedure copy_3d_real_array2
   module procedure copy_3d_real_array3
   module procedure copy_1d_double_array1
   module procedure copy_1d_double_array2
   module procedure copy_1d_double_array3
   module procedure copy_2d_double_array1
   module procedure copy_2d_double_array2
   module procedure copy_2d_double_array3
   module procedure copy_3d_double_array1
   module procedure copy_3d_double_array2
   module procedure copy_3d_double_array3
   module procedure copy_1d_complex_array1
   module procedure copy_1d_complex_array2
   module procedure copy_1d_complex_array3
   module procedure copy_2d_complex_array1
   module procedure copy_2d_complex_array2
   module procedure copy_2d_complex_array3
   module procedure copy_3d_complex_array1
   module procedure copy_3d_complex_array2
   module procedure copy_3d_complex_array3
end interface copy_array

interface sub2ind
   module procedure sub2ind_2d
   module procedure sub2ind_3d
   module procedure sub2ind_4d
end interface sub2ind

interface ind2sub
   module procedure ind2sub_2d
   module procedure ind2sub_3d
   module procedure ind2sub_4d
end interface ind2sub

interface permute
   module procedure permute2d_real
   module procedure permute3d_real
   module procedure permute2d_int
   module procedure permute3d_int
end interface permute

interface arr2vec
   module procedure arr2vec_real_2d
   module procedure arr2vec_real_3d
   module procedure arr2vec_real_4d
   module procedure arr2vec_int_2d
   module procedure arr2vec_int_3d
   module procedure arr2vec_int_4d
   module procedure arr2vec_double_2d
   module procedure arr2vec_double_3d
   module procedure arr2vec_double_4d
!   module procedure arr2vec_complex_2d
!   module procedure arr2vec_complex_3d
!   module procedure arr2vec_complex_4d
end interface arr2vec

interface vec2arr
   module procedure vec2arr_real_2d
   module procedure vec2arr_real_3d
   module procedure vec2arr_real_4d
   module procedure vec2arr_int_2d
   module procedure vec2arr_int_3d
   module procedure vec2arr_int_4d
   module procedure vec2arr_double_2d
   module procedure vec2arr_double_3d
   module procedure vec2arr_double_4d
!   module procedure vec2arr_complex_2d
!   module procedure vec2arr_complex_3d
!   module procedure vec2arr_complex_4d
end interface vec2arr

interface fliplr
   module procedure flipud_int_1d
   module procedure flipud_real_1d
   module procedure flipud_double_1d
   module procedure fliplr_int_2d
   module procedure fliplr_real_2d
   module procedure fliplr_double_2d
end interface fliplr

interface flipud
   module procedure flipud_int_1d
   module procedure flipud_real_1d
   module procedure flipud_double_1d
   module procedure flipud_int_2d
   module procedure flipud_real_2d
   module procedure flipud_double_2d
end interface flipud

interface rotate
   module procedure rotate_int_2d
   module procedure rotate_real_2d
   module procedure rotate_double_2d
end interface rotate

!interface locate
!   module procedure locate_int_1d
!   module procedure locate_int_2d
!end interface

contains

!=====================================================================

subroutine allocate_and_initial_1d_1i(d,n1,mem)
integer, allocatable, intent(out) :: d(:)
integer, intent(in)     :: n1
integer(8),optional,  intent(inout) :: mem
allocate(d(n1))
d=0
if (present(mem)) mem=mem+n1/MEM_UNIT
end subroutine allocate_and_initial_1d_1i

!=====================================================================

subroutine deallocate_and_free_1d_1i(d,mem)
integer,allocatable,  intent(inout) :: d(:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1
n1=size(d,1)
if (present(mem)) mem=mem-n1/MEM_UNIT
deallocate(d)
end subroutine deallocate_and_free_1d_1i

!=====================================================================

subroutine allocate_and_initial_1d_2i(d1,d2,n1,mem)
integer, allocatable, intent(out) :: d1(:),d2(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_1i(d1,n1,mem)
call allocate_and_initial_1d_1i(d2,n1,mem)
end subroutine allocate_and_initial_1d_2i

!=====================================================================

subroutine deallocate_and_free_1d_2i(d1,d2,mem)
integer, allocatable, intent(inout) :: d1(:), d2(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_1i(d1,mem)
call deallocate_and_free_1d_1i(d2,mem)
end subroutine deallocate_and_free_1d_2i

!=====================================================================

subroutine allocate_and_initial_1d_3i(d1,d2,d3,n1,mem)
integer, allocatable, intent(out) :: d1(:),d2(:),d3(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_2i(d1,d2,n1,mem)
call allocate_and_initial_1d_1i(d3,n1,mem)
end subroutine allocate_and_initial_1d_3i

!=====================================================================

subroutine deallocate_and_free_1d_3i(d1,d2,d3,mem)
integer, allocatable, intent(inout) :: d1(:), d2(:),d3(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_2i(d1,d2,mem)
call deallocate_and_free_1d_1i(d3,mem)
end subroutine deallocate_and_free_1d_3i

!=====================================================================

subroutine allocate_and_initial_1d_4i(d1,d2,d3,d4,n1,mem)
integer, allocatable, intent(out) :: d1(:),d2(:),d3(:),d4(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_2i(d1,d2,n1,mem)
call allocate_and_initial_1d_2i(d3,d4,n1,mem)
end subroutine allocate_and_initial_1d_4i

!=====================================================================

subroutine deallocate_and_free_1d_4i(d1,d2,d3,d4,mem)
integer, allocatable, intent(inout) :: d1(:), d2(:),d3(:),d4(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_2i(d1,d2,mem)
call deallocate_and_free_1d_2i(d3,d4,mem)
end subroutine deallocate_and_free_1d_4i

!=====================================================================

subroutine allocate_and_initial_1d_6i(d1,d2,d3,d4,d5,d6,n1,mem)
integer, allocatable, intent(out) :: d1(:),d2(:),d3(:),d4(:),d5(:),d6(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_3i(d1,d2,d3,n1,mem)
call allocate_and_initial_1d_3i(d4,d5,d6,n1,mem)
end subroutine allocate_and_initial_1d_6i

!=====================================================================

subroutine deallocate_and_free_1d_6i(d1,d2,d3,d4,d5,d6,mem)
integer, allocatable, intent(inout) :: d1(:), d2(:),d3(:),d4(:),d5(:),d6(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_3i(d1,d2,d3,mem)
call deallocate_and_free_1d_3i(d4,d5,d6,mem)
end subroutine deallocate_and_free_1d_6i

!=====================================================================

subroutine allocate_and_initial_1d_1r(d,n1,mem)
real, allocatable, intent(out) :: d(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
allocate(d(n1))
d=0.0
if (present(mem)) mem=mem+n1/MEM_UNIT
end subroutine allocate_and_initial_1d_1r

!=====================================================================

subroutine deallocate_and_free_1d_1r(d,mem)
real, allocatable,  intent(inout) :: d(:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1
n1=size(d,1)
if (present(mem)) mem=mem-n1/MEM_UNIT
deallocate(d)
end subroutine deallocate_and_free_1d_1r

!=====================================================================

subroutine allocate_and_initial_1d_2r(d1,d2,n1,mem)
real, allocatable, intent(out) :: d1(:),d2(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_1r(d1,n1,mem)
call allocate_and_initial_1d_1r(d2,n1,mem)
end subroutine allocate_and_initial_1d_2r

!=====================================================================

subroutine deallocate_and_free_1d_2r(d1,d2,mem)
real, allocatable, intent(inout) :: d1(:), d2(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_1r(d1,mem)
call deallocate_and_free_1d_1r(d2,mem)
end subroutine deallocate_and_free_1d_2r

!=====================================================================

subroutine allocate_and_initial_1d_3r(d1,d2,d3,n1,mem)
real, allocatable, intent(out) :: d1(:),d2(:),d3(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_2r(d1,d2,n1,mem)
call allocate_and_initial_1d_1r(d3,n1,mem)
end subroutine allocate_and_initial_1d_3r

!=====================================================================

subroutine deallocate_and_free_1d_3r(d1,d2,d3,mem)
real, allocatable, intent(inout) :: d1(:), d2(:), d3(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_2r(d1,d2,mem)
call deallocate_and_free_1d_1r(d3,mem)
end subroutine deallocate_and_free_1d_3r

!=====================================================================

subroutine allocate_and_initial_1d_4r(d1,d2,d3,d4,n1,mem)
real, allocatable, intent(out) :: d1(:),d2(:),d3(:),d4(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_2r(d1,d2,n1,mem)
call allocate_and_initial_1d_2r(d3,d4,n1,mem)
end subroutine allocate_and_initial_1d_4r

!=====================================================================

subroutine deallocate_and_free_1d_4r(d1,d2,d3,d4,mem)
real, allocatable, intent(inout) :: d1(:), d2(:), d3(:), d4(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_2r(d1,d2,mem)
call deallocate_and_free_1d_2r(d3,d4,mem)
end subroutine deallocate_and_free_1d_4r

!=====================================================================

subroutine allocate_and_initial_1d_5r(d1,d2,d3,d4,d5,n1,mem)
real, allocatable, intent(out) :: d1(:),d2(:),d3(:),d4(:),d5(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_3r(d1,d2,d3,n1,mem)
call allocate_and_initial_1d_2r(d4,d5,n1,mem)
end subroutine allocate_and_initial_1d_5r

!=====================================================================

subroutine deallocate_and_free_1d_5r(d1,d2,d3,d4,d5,mem)
real, allocatable, intent(inout) :: d1(:), d2(:), d3(:), d4(:),d5(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_3r(d1,d2,d3,mem)
call deallocate_and_free_1d_2r(d4,d5,mem)
end subroutine deallocate_and_free_1d_5r

!=====================================================================

subroutine allocate_and_initial_1d_6r(d1,d2,d3,d4,d5,d6,n1,mem)
real, allocatable, intent(out) :: d1(:),d2(:),d3(:),d4(:),d5(:),d6(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_3r(d1,d2,d3,n1,mem)
call allocate_and_initial_1d_3r(d4,d5,d6,n1,mem)
end subroutine allocate_and_initial_1d_6r

!=====================================================================

subroutine deallocate_and_free_1d_6r(d1,d2,d3,d4,d5,d6,mem)
real, allocatable, intent(inout) :: d1(:), d2(:), d3(:), d4(:),d5(:),d6(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_3r(d1,d2,d3,mem)
call deallocate_and_free_1d_3r(d4,d5,d6,mem)
end subroutine deallocate_and_free_1d_6r

!=====================================================================

subroutine allocate_and_initial_1d_7r(d1,d2,d3,d4,d5,d6,d7,n1,mem)
real, allocatable, intent(out) :: d1(:),d2(:),d3(:),d4(:),d5(:),d6(:),d7(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_4r(d1,d2,d3,d4,n1,mem)
call allocate_and_initial_1d_3r(d5,d6,d7,n1,mem)
end subroutine allocate_and_initial_1d_7r

!=====================================================================

subroutine deallocate_and_free_1d_7r(d1,d2,d3,d4,d5,d6,d7,mem)
real, allocatable, intent(inout) :: d1(:), d2(:), d3(:), d4(:),d5(:),d6(:),d7(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_4r(d1,d2,d3,d4,mem)
call deallocate_and_free_1d_3r(d5,d6,d7,mem)
end subroutine deallocate_and_free_1d_7r

!=====================================================================
subroutine allocate_and_initial_1d_1c(d,n1,mem)
complex, allocatable, intent(out) :: d(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
allocate(d(n1))
d=0.0
if (present(mem)) mem=n1*2/MEM_UNIT
end subroutine allocate_and_initial_1d_1c

!=====================================================================

subroutine deallocate_and_free_1d_1c(d,mem)
complex, allocatable, intent(inout) :: d(:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1
n1=size(d,1)
if (present(mem)) mem=mem-n1*2/MEM_UNIT
deallocate(d)
end subroutine deallocate_and_free_1d_1c

!=====================================================================

subroutine allocate_and_initial_1d_2c(d1,d2,n1,mem)
complex, allocatable, intent(out) :: d1(:),d2(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_1d_1c(d1,n1,mem)
call allocate_and_initial_1d_1c(d2,n1,mem)
end subroutine allocate_and_initial_1d_2c

!=====================================================================

subroutine deallocate_and_free_1d_2c(d1,d2,mem)
complex, allocatable, intent(inout) :: d1(:), d2(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_1c(d1,mem)
call deallocate_and_free_1d_1c(d2,mem)
end subroutine deallocate_and_free_1d_2c
!======================================================================

subroutine deallocate_and_free_2d_1c(d1,mem)
complex, allocatable, intent(inout) :: d1(:,:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1,n2
n1=size(d1,1);n2=size(d1,2)
if (present(mem)) mem=mem-n1*n2*2/MEM_UNIT
deallocate(d1)
end subroutine deallocate_and_free_2d_1c
!=====================================================================

subroutine allocate_and_initial_1d_3c(d1,d2,d3,n1,mem)
complex, allocatable, intent(out) :: d1(:),d2(:),d3(:)
integer(8),optional,  intent(inout) :: mem
integer, intent(in) :: n1
call allocate_and_initial_1d_2c(d1,d2,n1,mem)
call allocate_and_initial_1d_1c(d3,n1,mem)
end subroutine allocate_and_initial_1d_3c

!=====================================================================
! Added by Bowen Guo
subroutine allocate_and_initial_2d_1c(d1,n1,n2,mem)
complex, allocatable, intent(out) :: d1(:,:)
integer(8),optional,  intent(inout) :: mem
integer, intent(in) :: n1,n2
allocate(d1(n1,n2))
d1=0.0
end subroutine allocate_and_initial_2d_1c

!=====================================================================


subroutine deallocate_and_free_1d_3c(d1,d2,d3,mem)
complex, allocatable, intent(inout) :: d1(:), d2(:),d3(:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_1d_2c(d1,d2,mem)
call deallocate_and_free_1d_1c(d3,mem)
end subroutine deallocate_and_free_1d_3c

!=====================================================================

subroutine allocate_and_initial_2d_1i(d1,n1,n2,mem)
integer, allocatable, intent(out) :: d1(:,:)
integer, intent(in) :: n1,n2
integer(8),optional,  intent(inout) :: mem
allocate(d1(n1,n2))
d1=0
if (present(mem)) mem=mem+n1*n2/MEM_UNIT
end subroutine allocate_and_initial_2d_1i

!=====================================================================

subroutine deallocate_and_free_2d_1i(d,mem)
integer, allocatable, intent(inout) :: d(:,:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1,n2
n1=size(d,1)
n2=size(d,2)
if (present(mem)) mem=mem-n1*n2/MEM_UNIT
deallocate(d)
end subroutine deallocate_and_free_2d_1i

!=====================================================================

subroutine allocate_and_initial_2d_2i(d1,d2,n1,n2,mem)
integer, allocatable, intent(out) :: d1(:,:),d2(:,:)
integer, intent(in) :: n1,n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_1i(d1,n1,n2,mem)
call allocate_and_initial_2d_1i(d2,n1,n2,mem)
end subroutine allocate_and_initial_2d_2i

!=====================================================================

subroutine deallocate_and_free_2d_2i(d1,d2,mem)
integer, allocatable, intent(inout) :: d1(:,:),d2(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_1i(d1,mem)
call deallocate_and_free_2d_1i(d2,mem)
end subroutine deallocate_and_free_2d_2i

!=====================================================================

subroutine allocate_and_initial_2d_1r(d,n1,n2,mem)
real, allocatable, intent(out) :: d(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
allocate(d(n1,n2))
d=0.0
if (present(mem)) mem=mem+n1*n2/MEM_UNIT
end subroutine allocate_and_initial_2d_1r

!=====================================================================

subroutine deallocate_and_free_2d_1r(d,mem)
real,    allocatable, intent(inout) :: d(:,:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1,n2
n1=size(d,1)
n2=size(d,2)
if (present(mem)) mem=mem-n1*n2/MEM_UNIT
deallocate(d)
end subroutine deallocate_and_free_2d_1r

!=====================================================================

subroutine allocate_and_initial_2d_2r(d1,d2,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_1r(d1,n1,n2,mem)
call allocate_and_initial_2d_1r(d2,n1,n2,mem)
end subroutine allocate_and_initial_2d_2r

!=====================================================================

subroutine deallocate_and_free_2d_2r(d1,d2,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_1r(d1,mem)
call deallocate_and_free_2d_1r(d2,mem)
end subroutine deallocate_and_free_2d_2r

!=====================================================================

subroutine allocate_and_initial_2d_3r(d1,d2,d3,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_2r(d1,d2,n1,n2,mem)
call allocate_and_initial_2d_1r(d3,n1,n2,mem)
end subroutine allocate_and_initial_2d_3r

!=====================================================================

subroutine deallocate_and_free_2d_3r(d1,d2,d3,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_2r(d1,d2,mem)
call deallocate_and_free_2d_1r(d3,mem)
end subroutine deallocate_and_free_2d_3r

!=====================================================================

subroutine allocate_and_initial_2d_4r(d1,d2,d3,d4,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_2r(d1,d2,n1,n2,mem)
call allocate_and_initial_2d_2r(d3,d4,n1,n2,mem)
end subroutine allocate_and_initial_2d_4r

!=====================================================================

subroutine deallocate_and_free_2d_4r(d1,d2,d3,d4,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:),d4(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_2r(d1,d2,mem)
call deallocate_and_free_2d_2r(d3,d4,mem)
end subroutine deallocate_and_free_2d_4r

!=====================================================================

subroutine allocate_and_initial_2d_5r(d1,d2,d3,d4,d5,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),d5(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_3r(d1,d2,d3,n1,n2,mem)
call allocate_and_initial_2d_2r(d4,d5,n1,n2,mem)
end subroutine allocate_and_initial_2d_5r

!=====================================================================

subroutine deallocate_and_free_2d_5r(d1,d2,d3,d4,d5,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:),d4(:,:),d5(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_3r(d1,d2,d3,mem)
call deallocate_and_free_2d_2r(d4,d5,mem)
end subroutine deallocate_and_free_2d_5r

!=====================================================================

subroutine allocate_and_initial_2d_6r(d1,d2,d3,d4,d5,d6,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:),&
                                  d4(:,:),d5(:,:),d6(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_3r(d1,d2,d3,n1,n2,mem)
call allocate_and_initial_2d_3r(d4,d5,d6,n1,n2,mem)
end subroutine allocate_and_initial_2d_6r

!=====================================================================

subroutine deallocate_and_free_2d_6r(d1,d2,d3,d4,d5,d6,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_3r(d1,d2,d3,mem)
call deallocate_and_free_2d_3r(d4,d5,d6,mem)
end subroutine deallocate_and_free_2d_6r

!=====================================================================

subroutine allocate_and_initial_2d_7r(d1,d2,d3,d4,d5,d6,d7,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:),&
                                  d4(:,:),d5(:,:),d6(:,:),d7(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_4r(d1,d2,d3,d4,n1,n2,mem)
call allocate_and_initial_2d_3r(d5,d6,d7,n1,n2,mem)
end subroutine allocate_and_initial_2d_7r

!=====================================================================

subroutine deallocate_and_free_2d_7r(d1,d2,d3,d4,d5,d6,d7,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:),d7(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_4r(d1,d2,d3,d4,mem)
call deallocate_and_free_2d_3r(d5,d6,d7,mem)
end subroutine deallocate_and_free_2d_7r

!=====================================================================

subroutine allocate_and_initial_2d_9r(d1,d2,d3,d4,d5,d6,d7,d8,d9,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:),&
   d4(:,:),d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_6r(d1,d2,d3,d4,d5,d6,n1,n2,mem)
call allocate_and_initial_2d_3r(d7,d8,d9,n1,n2,mem)
end subroutine allocate_and_initial_2d_9r

!=====================================================================

subroutine deallocate_and_free_2d_9r(d1,d2,d3,d4,d5,d6,d7,d8,d9,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_6r(d1,d2,d3,d4,d5,d6,mem)
call deallocate_and_free_2d_3r(d7,d8,d9,mem)
end subroutine deallocate_and_free_2d_9r

!=====================================================================

subroutine allocate_and_initial_2d_12r(d1,d2,d3,d4,d5,d6,d7,d8,d9,&
   d10,d11,d12,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:),d10(:,:),d11(:,:),d12(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_6r(d1,d2,d3,d4,d5,d6,n1,n2,mem)
call allocate_and_initial_2d_6r(d7,d8,d9,d10,d11,d12,n1,n2,mem)
end subroutine allocate_and_initial_2d_12r

!=====================================================================

subroutine deallocate_and_free_2d_12r(d1,d2,d3,d4,d5,d6,d7,d8,d9,&
   d10,d11,d12,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:),d10(:,:),d11(:,:),d12(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_6r(d1,d2,d3,d4,d5,d6,mem)
call deallocate_and_free_2d_6r(d7,d8,d9,d10,d11,d12,mem)
end subroutine deallocate_and_free_2d_12r

!=====================================================================

subroutine allocate_and_initial_2d_15r(d1,d2,d3,d4,d5,d6,d7,d8,d9,&
   d10,d11,d12,d13,d14,d15,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:),d10(:,:),d11(:,:),d12(:,:),&
   d13(:,:),d14(:,:),d15(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_9r(d1,d2,d3,d4,d5,d6,d7,d8,d9,n1,n2,mem)
call allocate_and_initial_2d_6r(d10,d11,d12,d13,d14,d15,n1,n2,mem)
end subroutine allocate_and_initial_2d_15r

!=====================================================================

subroutine deallocate_and_free_2d_15r(d1,d2,d3,d4,d5,d6,d7,d8,d9,&
   d10,d11,d12,d13,d14,d15,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:),d10(:,:),d11(:,:),d12(:,:),&
   d13(:,:),d14(:,:),d15(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_9r(d1,d2,d3,d4,d5,d6,d7,d8,d9,mem)
call deallocate_and_free_2d_6r(d10,d11,d12,d13,d14,d15,mem)
end subroutine deallocate_and_free_2d_15r

!=====================================================================

subroutine allocate_and_initial_2d_18r(d1,d2,d3,d4,d5,d6,d7,d8,d9,&
   d10,d11,d12,d13,d14,d15,d16,d17,d18,n1,n2,mem)
real, allocatable, intent(out) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:),d10(:,:),d11(:,:),d12(:,:),&
   d13(:,:),d14(:,:),d15(:,:),d16(:,:),d17(:,:),d18(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_9r(d1,d2,d3,d4,d5,d6,d7,d8,d9,n1,n2,mem)
call allocate_and_initial_2d_9r(d10,d11,d12,d13,d14,d15,&
   d16,d17,d18,n1,n2,mem)
end subroutine allocate_and_initial_2d_18r

!=====================================================================

subroutine deallocate_and_free_2d_18r(d1,d2,d3,d4,d5,d6,d7,d8,d9,&
   d10,d11,d12,d13,d14,d15,d16,d17,d18,mem)
real, allocatable, intent(inout) :: d1(:,:), d2(:,:),d3(:,:),d4(:,:),&
   d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:),d10(:,:),d11(:,:),d12(:,:),&
   d13(:,:),d14(:,:),d15(:,:),d16(:,:),d17(:,:),d18(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_9r(d1,d2,d3,d4,d5,d6,d7,d8,d9,mem)
call deallocate_and_free_2d_9r(d10,d11,d12,d13,d14,d15,&
   d16,d17,d18,mem)
end subroutine deallocate_and_free_2d_18r

!=====================================================================

subroutine allocate_and_initial_1d_1dp(d,n1,mem)
double precision, allocatable, intent(out) :: d(:)
integer, intent(in) :: n1
integer(8),optional,  intent(inout) :: mem
allocate(d(n1))
d=0.0
if (present(mem)) mem=mem+n1/MEM_UNIT*2
end subroutine allocate_and_initial_1d_1dp

!=====================================================================

subroutine deallocate_and_free_1d_1dp(d,mem)
double precision,allocatable, intent(inout) :: d(:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1
n1=size(d)
if (present(mem)) mem=mem-n1/MEM_UNIT*2
deallocate(d)
end subroutine deallocate_and_free_1d_1dp

!=====================================================================


subroutine allocate_and_initial_2d_1dp(d,n1,n2,mem)
double precision, allocatable, intent(out) :: d(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
allocate(d(n1,n2))
d=0.0
if (present(mem)) mem=mem+n1*n2/MEM_UNIT*2
end subroutine allocate_and_initial_2d_1dp

!=====================================================================

subroutine deallocate_and_free_2d_1dp(d,mem)
double precision,allocatable, intent(inout) :: d(:,:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1,n2
n1=size(d,1)
n2=size(d,2)
if (present(mem)) mem=mem-n1*n2/MEM_UNIT*2
deallocate(d)
end subroutine deallocate_and_free_2d_1dp

!=====================================================================

subroutine allocate_and_initial_2d_2dp(d1,d2,n1,n2,mem)
double precision, allocatable, intent(out) :: d1(:,:),d2(:,:)
integer, intent(in) :: n1, n2
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_2d_1dp(d1,n1,n2,mem)
call allocate_and_initial_2d_1dp(d2,n1,n2,mem)
end subroutine allocate_and_initial_2d_2dp

!=====================================================================

subroutine deallocate_and_free_2d_2dp(d1,d2,mem)
double precision,allocatable, intent(inout) :: d1(:,:),d2(:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_2d_1dp(d1,mem)
call deallocate_and_free_2d_1dp(d2,mem)
end subroutine deallocate_and_free_2d_2dp

!=====================================================================
! Added by Bowen Guo
subroutine allocate_and_initial_3d_1i(d,n1,n2,n3,mem)
integer,allocatable,intent(out)   :: d(:,:,:)
integer,intent(in)                :: n1,n2,n3
integer(8),optional,intent(inout) :: mem
allocate(d(n1,n2,n3))
d=0.0
if (present(mem)) mem=mem+n1*n2*n3/MEM_UNIT
end subroutine allocate_and_initial_3d_1i
!=====================================================================

subroutine allocate_and_initial_3d_1r(d,n1,n2,n3,mem)
real, allocatable, intent(out) :: d(:,:,:)
integer, intent(in) :: n1, n2, n3
integer(8),optional,  intent(inout) :: mem
allocate(d(n1,n2,n3))
d=0.0
if (present(mem)) mem=mem+n1*n2*n3/MEM_UNIT
end subroutine allocate_and_initial_3d_1r

!=====================================================================
! Added by Bowen Guo
subroutine deallocate_and_free_3d_1i(d,mem)
integer,    allocatable, intent(inout) :: d(:,:,:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1,n2,n3
n1=size(d,1)
n2=size(d,2)
n3=size(d,3)
if (present(mem)) mem=mem-n1*n2*n3/MEM_UNIT
deallocate(d)
end subroutine deallocate_and_free_3d_1i

!========================================================================

subroutine deallocate_and_free_3d_1r(d,mem)
real,    allocatable, intent(inout) :: d(:,:,:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1,n2,n3
n1=size(d,1)
n2=size(d,2)
n3=size(d,3)
if (present(mem)) mem=mem-n1*n2*n3/MEM_UNIT
deallocate(d)
end subroutine deallocate_and_free_3d_1r

!=====================================================================

subroutine allocate_and_initial_3d_2r(d1,d2,n1,n2,n3,mem)
real, allocatable, intent(out) :: d1(:,:,:),d2(:,:,:)
integer, intent(in) :: n1, n2, n3
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_3d_1r(d1,n1,n2,n3,mem)
call allocate_and_initial_3d_1r(d2,n1,n2,n3,mem)
end subroutine allocate_and_initial_3d_2r

!=====================================================================

subroutine deallocate_and_free_3d_2r(d1,d2,mem)
real, allocatable, intent(inout) :: d1(:,:,:), d2(:,:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_3d_1r(d1,mem)
call deallocate_and_free_3d_1r(d2,mem)
end subroutine deallocate_and_free_3d_2r

!=====================================================================

subroutine allocate_and_initial_3d_3r(d1,d2,d3,n1,n2,n3,mem)
real, allocatable, intent(out) :: d1(:,:,:),d2(:,:,:),d3(:,:,:)
integer, intent(in) :: n1, n2, n3
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_3d_2r(d1,d2,n1,n2,n3,mem)
call allocate_and_initial_3d_1r(d3,n1,n2,n3,mem)
end subroutine allocate_and_initial_3d_3r

!=====================================================================

subroutine deallocate_and_free_3d_3r(d1,d2,d3,mem)
real, allocatable, intent(inout) :: d1(:,:,:), d2(:,:,:), d3(:,:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_3d_2r(d1,d2,mem)
call deallocate_and_free_3d_1r(d3,mem)
end subroutine deallocate_and_free_3d_3r

!=====================================================================

subroutine allocate_and_initial_3d_4r(d1,d2,d3,d4,n1,n2,n3,mem)
real, allocatable, intent(out) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),&
   d4(:,:,:)
integer, intent(in) :: n1, n2, n3
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_3d_2r(d1,d2,n1,n2,n3,mem)
call allocate_and_initial_3d_2r(d3,d4,n1,n2,n3,mem)
end subroutine allocate_and_initial_3d_4r

!=====================================================================

subroutine deallocate_and_free_3d_4r(d1,d2,d3,d4,mem)
real, allocatable, intent(inout) :: d1(:,:,:), d2(:,:,:), d3(:,:,:),&
   d4(:,:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_3d_2r(d1,d2,mem)
call deallocate_and_free_3d_2r(d3,d4,mem)
end subroutine deallocate_and_free_3d_4r

!=====================================================================

subroutine allocate_and_initial_3d_5r(d1,d2,d3,d4,d5,n1,n2,n3,mem)
real, allocatable, intent(out) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),&
   d4(:,:,:),d5(:,:,:)
integer, intent(in) :: n1, n2, n3
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_3d_3r(d1,d2,d3,n1,n2,n3,mem)
call allocate_and_initial_3d_2r(d4,d5,n1,n2,n3,mem)
end subroutine allocate_and_initial_3d_5r

!=====================================================================

subroutine deallocate_and_free_3d_5r(d1,d2,d3,d4,d5,mem)
real, allocatable, intent(inout) :: d1(:,:,:), d2(:,:,:), d3(:,:,:),&
   d4(:,:,:),d5(:,:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_3d_3r(d1,d2,d3,mem)
call deallocate_and_free_3d_2r(d4,d5,mem)
end subroutine deallocate_and_free_3d_5r

!=====================================================================

subroutine allocate_and_initial_3d_6r(d1,d2,d3,d4,d5,d6,n1,n2,n3,mem)
real, allocatable, intent(out) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),&
   d4(:,:,:),d5(:,:,:),d6(:,:,:)
integer, intent(in) :: n1, n2, n3
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_3d_3r(d1,d2,d3,n1,n2,n3,mem)
call allocate_and_initial_3d_3r(d4,d5,d6,n1,n2,n3,mem)
end subroutine allocate_and_initial_3d_6r

!=====================================================================

subroutine deallocate_and_free_3d_6r(d1,d2,d3,d4,d5,d6,mem)
real, allocatable, intent(inout) :: d1(:,:,:), d2(:,:,:), d3(:,:,:),&
   d4(:,:,:),d5(:,:,:),d6(:,:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_3d_3r(d1,d2,d3,mem)
call deallocate_and_free_3d_3r(d4,d5,d6,mem)
end subroutine deallocate_and_free_3d_6r

!=====================================================================

subroutine allocate_and_initial_4d_1r(d,n1,n2,n3,n4,mem)
real, allocatable, intent(out) :: d(:,:,:,:)
integer, intent(in) :: n1, n2, n3, n4
integer(8),optional,  intent(inout) :: mem
allocate(d(n1,n2,n3,n4))
d=0.0
if (present(mem)) mem=mem+n1*n2*n3*n4/MEM_UNIT
end subroutine allocate_and_initial_4d_1r

!=====================================================================

subroutine deallocate_and_free_4d_1r(d,mem)
real,    allocatable, intent(inout) :: d(:,:,:,:)
integer(8),optional,  intent(inout) :: mem
integer                             :: n1,n2,n3,n4
n1=size(d,1)
n2=size(d,2)
n3=size(d,3)
n4=size(d,4)
if (present(mem)) mem=mem-n1*n2*n3*n4/MEM_UNIT
deallocate(d)
end subroutine deallocate_and_free_4d_1r

!=====================================================================

subroutine allocate_and_initial_4d_2r(d1,d2,n1,n2,n3,n4,mem)
real, allocatable, intent(out) :: d1(:,:,:,:),d2(:,:,:,:)
integer, intent(in) :: n1, n2, n3, n4
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_4d_1r(d1,n1,n2,n3,n4,mem)
call allocate_and_initial_4d_1r(d2,n1,n2,n3,n4,mem)
end subroutine allocate_and_initial_4d_2r

!=====================================================================

subroutine deallocate_and_free_4d_2r(d1,d2,mem)
real, allocatable, intent(inout) :: d1(:,:,:,:), d2(:,:,:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_4d_1r(d1,mem)
call deallocate_and_free_4d_1r(d2,mem)
end subroutine deallocate_and_free_4d_2r

!=====================================================================

subroutine allocate_and_initial_4d_3r(d1,d2,d3,n1,n2,n3,n4,mem)
real, allocatable, intent(out) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:)
integer, intent(in) :: n1, n2, n3, n4
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_4d_2r(d1,d2,n1,n2,n3,n4,mem)
call allocate_and_initial_4d_1r(d3,n1,n2,n3,n4,mem)
end subroutine allocate_and_initial_4d_3r

!=====================================================================

subroutine deallocate_and_free_4d_3r(d1,d2,d3,mem)
real, allocatable, intent(inout) :: d1(:,:,:,:), d2(:,:,:,:), d3(:,:,:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_4d_2r(d1,d2,mem)
call deallocate_and_free_4d_1r(d3,mem)
end subroutine deallocate_and_free_4d_3r

!=====================================================================

subroutine allocate_and_initial_4d_6r(d1,d2,d3,d4,d5,d6,n1,n2,n3,n4,mem)
real, allocatable, intent(out) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:),&
   d4(:,:,:,:),d5(:,:,:,:),d6(:,:,:,:)
integer, intent(in) :: n1, n2, n3, n4
integer(8),optional,  intent(inout) :: mem
call allocate_and_initial_4d_3r(d1,d2,d3,n1,n2,n3,n4,mem)
call allocate_and_initial_4d_3r(d4,d5,d6,n1,n2,n3,n4,mem)
end subroutine allocate_and_initial_4d_6r

!=====================================================================

subroutine deallocate_and_free_4d_6r(d1,d2,d3,d4,d5,d6,mem)
real, allocatable, intent(inout) :: d1(:,:,:,:), d2(:,:,:,:), d3(:,:,:,:),&
   d4(:,:,:,:),d5(:,:,:,:),d6(:,:,:,:)
integer(8),optional,  intent(inout) :: mem
call deallocate_and_free_4d_3r(d1,d2,d3,mem)
call deallocate_and_free_4d_3r(d4,d5,d6,mem)
end subroutine deallocate_and_free_4d_6r

!=====================================================================

subroutine initial_1d_1i(d)
integer, intent(out) :: d(:)
d=0
end subroutine initial_1d_1i

!=====================================================================

subroutine initial_1d_2i(d1,d2)
integer, intent(out) :: d1(:),d2(:)
call initial_1d_1i(d1)
call initial_1d_1i(d2)
end subroutine initial_1d_2i

!=====================================================================

subroutine initial_1d_3i(d1,d2,d3)
integer, intent(out) :: d1(:),d2(:),d3(:)
call initial_1d_2i(d1,d2)
call initial_1d_1i(d3)
end subroutine initial_1d_3i

!=====================================================================

subroutine initial_1d_1r(d)
real,    intent(out) :: d(:)
d=0.0
end subroutine initial_1d_1r

!=====================================================================

subroutine initial_1d_2r(d1,d2)
real,    intent(out) :: d1(:),d2(:)
call initial_1d_1r(d1)
call initial_1d_1r(d2)
end subroutine initial_1d_2r

!=====================================================================

subroutine initial_1d_3r(d1,d2,d3)
real, intent(out) :: d1(:),d2(:),d3(:)
call initial_1d_2r(d1,d2)
call initial_1d_1r(d3)
end subroutine initial_1d_3r

!=====================================================================

subroutine initial_1d_6r(d1,d2,d3,d4,d5,d6)
real,    intent(out) :: d1(:),d2(:),d3(:),d4(:),d5(:),d6(:)
call initial_1d_3r(d1,d2,d3)
call initial_1d_3r(d4,d5,d6)
end subroutine initial_1d_6r

!=====================================================================

subroutine initial_2d_1r(d)
real,    intent(out) :: d(:,:)
d=0.0
end subroutine initial_2d_1r

!=====================================================================

subroutine initial_2d_2r(d1,d2)
real,    intent(out) :: d1(:,:),d2(:,:)
call initial_2d_1r(d1)
call initial_2d_1r(d2)
end subroutine initial_2d_2r

!=====================================================================

subroutine initial_2d_3r(d1,d2,d3)
real,    intent(out) :: d1(:,:),d2(:,:),d3(:,:)
call initial_2d_2r(d1,d2)
call initial_2d_1r(d3)
end subroutine initial_2d_3r

!=====================================================================

subroutine initial_2d_4r(d1,d2,d3,d4)
real,    intent(out) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:)
call initial_2d_2r(d1,d2)
call initial_2d_2r(d3,d4)
end subroutine initial_2d_4r

!=====================================================================

subroutine initial_2d_5r(d1,d2,d3,d4,d5)
real,    intent(out) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),d5(:,:)
call initial_2d_3r(d1,d2,d3)
call initial_2d_2r(d4,d5)
end subroutine initial_2d_5r

!=====================================================================

subroutine initial_2d_6r(d1,d2,d3,d4,d5,d6)
real,intent(out) :: d1(:,:),d2(:,:),d3(:,:),d4(:,:),d5(:,:),d6(:,:)
call initial_2d_3r(d1,d2,d3)
call initial_2d_3r(d4,d5,d6)
end subroutine initial_2d_6r

!=====================================================================

subroutine initial_3d_1r(d)
real,    intent(out) :: d(:,:,:)
d=0.0
end subroutine initial_3d_1r

!=====================================================================

subroutine initial_3d_2r(d1,d2)
real,   intent(out) :: d1(:,:,:),d2(:,:,:)
call initial_3d_1r(d1)
call initial_3d_1r(d2)
end subroutine initial_3d_2r

!=====================================================================

subroutine initial_3d_3r(d1,d2,d3)
real,    intent(out) :: d1(:,:,:),d2(:,:,:),d3(:,:,:)
call initial_3d_2r(d1,d2)
call initial_3d_1r(d3)
end subroutine initial_3d_3r

!=====================================================================

subroutine initial_3d_4r(d1,d2,d3,d4)
real, intent(out) :: d1(:,:,:),d2(:,:,:),d3(:,:,:),d4(:,:,:)
call initial_3d_2r(d1,d2)
call initial_3d_2r(d3,d4)
end subroutine initial_3d_4r

!=====================================================================

subroutine initial_3d_6r(d1,d2,d3,d4,d5,d6)
real,intent(out)::d1(:,:,:),d2(:,:,:),d3(:,:,:),d4(:,:,:),d5(:,:,:),d6(:,:,:)
call initial_3d_3r(d1,d2,d3)
call initial_3d_3r(d4,d5,d6)
end subroutine initial_3d_6r

!=====================================================================

subroutine initial_4d_1r(d)
real,    intent(out) :: d(:,:,:,:)
d=0.0
end subroutine initial_4d_1r

!=====================================================================

subroutine initial_4d_2r(d1,d2)
real,    intent(out) :: d1(:,:,:,:),d2(:,:,:,:)
call initial_4d_1r(d1)
call initial_4d_1r(d2)
end subroutine initial_4d_2r

!=====================================================================

subroutine initial_4d_3r(d1,d2,d3)
real,intent(out) :: d1(:,:,:,:),d2(:,:,:,:),d3(:,:,:,:)
call initial_4d_2r(d1,d2)
call initial_4d_1r(d3)
end subroutine initial_4d_3r

!=====================================================================

subroutine copy_1d_int_array1(i1,oi1)
integer, intent(in)  :: i1(:)
integer, intent(out) :: oi1(:)
oi1=i1
end subroutine copy_1d_int_array1

!====================================================================

subroutine copy_1d_int_array2(i1,oi1,i2,oi2)
integer, intent(in)  :: i1(:),i2(:)
integer, intent(out) :: oi1(:),oi2(:)
call copy_1d_int_array1(i1,oi1)
call copy_1d_int_array1(i2,oi2)
end subroutine copy_1d_int_array2

!====================================================================

subroutine copy_1d_int_array3(i1,oi1,i2,oi2,i3,oi3)
integer, intent(in)  :: i1(:),i2(:),i3(:)
integer, intent(out) :: oi1(:),oi2(:),oi3(:)
call copy_1d_int_array2(i1,oi1,i2,oi2)
call copy_1d_int_array1(i3,oi3)
end subroutine copy_1d_int_array3

!====================================================================

subroutine copy_2d_int_array1(i1,oi1)
integer, intent(in)  :: i1(:,:)
integer, intent(out) :: oi1(:,:)
oi1=i1
end subroutine copy_2d_int_array1

!====================================================================

subroutine copy_2d_int_array2(i1,oi1,i2,oi2)
integer, intent(in)  :: i1(:,:),i2(:,:)
integer, intent(out) :: oi1(:,:),oi2(:,:)
call copy_2d_int_array1(i1,oi1)
call copy_2d_int_array1(i2,oi2)
end subroutine copy_2d_int_array2

!====================================================================

subroutine copy_2d_int_array3(i1,oi1,i2,oi2,i3,oi3)
integer, intent(in)  :: i1(:,:),i2(:,:),i3(:,:)
integer, intent(out) :: oi1(:,:),oi2(:,:),oi3(:,:)
call copy_2d_int_array2(i1,oi1,i2,oi2)
call copy_2d_int_array1(i3,oi3)
end subroutine copy_2d_int_array3

!====================================================================

subroutine copy_3d_int_array1(i1,oi1)
integer, intent(in)  :: i1(:,:,:)
integer, intent(out) :: oi1(:,:,:)
oi1=i1
end subroutine copy_3d_int_array1

!====================================================================

subroutine copy_3d_int_array2(i1,oi1,i2,oi2)
integer, intent(in)  :: i1(:,:,:),i2(:,:,:)
integer, intent(out) :: oi1(:,:,:),oi2(:,:,:)
call copy_3d_int_array1(i1,oi1)
call copy_3d_int_array1(i2,oi2)
end subroutine copy_3d_int_array2

!====================================================================

subroutine copy_3d_int_array3(i1,oi1,i2,oi2,i3,oi3)
integer, intent(in)  :: i1(:,:,:),i2(:,:,:),i3(:,:,:)
integer, intent(out) :: oi1(:,:,:),oi2(:,:,:),oi3(:,:,:)
call copy_3d_int_array2(i1,oi1,i2,oi2)
call copy_3d_int_array1(i3,oi3)
end subroutine copy_3d_int_array3

!====================================================================

subroutine copy_1d_real_array1(r1,or1)
real, intent(in)  :: r1(:)
real, intent(out) :: or1(:)
or1=r1
end subroutine copy_1d_real_array1

!====================================================================

subroutine copy_1d_real_array2(r1,or1,r2,or2)
real, intent(in)  :: r1(:),r2(:)
real, intent(out) :: or1(:),or2(:)
call copy_1d_real_array1(r1,or1)
call copy_1d_real_array1(r2,or2)
end subroutine copy_1d_real_array2

!====================================================================

subroutine copy_1d_real_array3(r1,or1,r2,or2,r3,or3)
real, intent(in)  :: r1(:),r2(:),r3(:)
real, intent(out) :: or1(:),or2(:),or3(:)
call copy_1d_real_array2(r1,or1,r2,or2)
call copy_1d_real_array1(r3,or3)
end subroutine copy_1d_real_array3

!====================================================================

subroutine copy_2d_real_array1(r1,or1)
real, intent(in)  :: r1(:,:)
real, intent(out) :: or1(:,:)
or1=r1
end subroutine copy_2d_real_array1

!====================================================================

subroutine copy_2d_real_array2(r1,or1,r2,or2)
real, intent(in)  :: r1(:,:),r2(:,:)
real, intent(out) :: or1(:,:),or2(:,:)
call copy_2d_real_array1(r1,or1)
call copy_2d_real_array1(r2,or2)
end subroutine copy_2d_real_array2

!====================================================================

subroutine copy_2d_real_array3(r1,or1,r2,or2,r3,or3)
real, intent(in)  :: r1(:,:),r2(:,:),r3(:,:)
real, intent(out) :: or1(:,:),or2(:,:),or3(:,:)
call copy_2d_real_array2(r1,or1,r2,or2)
call copy_2d_real_array1(r3,or3)
end subroutine copy_2d_real_array3

!====================================================================

subroutine copy_3d_real_array1(r1,or1)
real, intent(in)  :: r1(:,:,:)
real, intent(out) :: or1(:,:,:)
or1=r1
end subroutine copy_3d_real_array1

!====================================================================

subroutine copy_3d_real_array2(r1,or1,r2,or2)
real, intent(in)  :: r1(:,:,:),r2(:,:,:)
real, intent(out) :: or1(:,:,:),or2(:,:,:)
call copy_3d_real_array1(r1,or1)
call copy_3d_real_array1(r2,or2)
end subroutine copy_3d_real_array2

!====================================================================

subroutine copy_3d_real_array3(r1,or1,r2,or2,r3,or3)
real, intent(in)  :: r1(:,:,:),r2(:,:,:),r3(:,:,:)
real, intent(out) :: or1(:,:,:),or2(:,:,:),or3(:,:,:)
call copy_3d_real_array2(r1,or1,r2,or2)
call copy_3d_real_array1(r3,or3)
end subroutine copy_3d_real_array3

!====================================================================

subroutine copy_1d_double_array1(d1,od1)
double precision, intent(in)  :: d1(:)
double precision, intent(out) :: od1(:)
od1=d1
end subroutine copy_1d_double_array1

!====================================================================

subroutine copy_1d_double_array2(d1,od1,d2,od2)
double precision, intent(in)  :: d1(:),d2(:)
double precision, intent(out) :: od1(:),od2(:)
call copy_1d_double_array1(d1,od1)
call copy_1d_double_array1(d2,od2)
end subroutine copy_1d_double_array2

!====================================================================

subroutine copy_1d_double_array3(d1,od1,d2,od2,d3,od3)
double precision, intent(in)  :: d1(:),d2(:),d3(:)
double precision, intent(out) :: od1(:),od2(:),od3(:)
call copy_1d_double_array2(d1,od1,d2,od2)
call copy_1d_double_array1(d3,od3)
end subroutine copy_1d_double_array3

!====================================================================

subroutine copy_2d_double_array1(d1,od1)
double precision, intent(in)  :: d1(:,:)
double precision, intent(out) :: od1(:,:)
od1=d1
end subroutine copy_2d_double_array1

!====================================================================

subroutine copy_2d_double_array2(d1,od1,d2,od2)
double precision, intent(in)  :: d1(:,:),d2(:,:)
double precision, intent(out) :: od1(:,:),od2(:,:)
call copy_2d_double_array1(d1,od1)
call copy_2d_double_array1(d2,od2)
end subroutine copy_2d_double_array2

!====================================================================

subroutine copy_2d_double_array3(d1,od1,d2,od2,d3,od3)
double precision, intent(in)  :: d1(:,:),d2(:,:),d3(:,:)
double precision, intent(out) :: od1(:,:),od2(:,:),od3(:,:)
call copy_2d_double_array2(d1,od1,d2,od2)
call copy_2d_double_array1(d3,od3)
end subroutine copy_2d_double_array3

!====================================================================

subroutine copy_3d_double_array1(d1,od1)
double precision, intent(in)  :: d1(:,:,:)
double precision, intent(out) :: od1(:,:,:)
od1=d1
end subroutine copy_3d_double_array1

!====================================================================

subroutine copy_3d_double_array2(d1,od1,d2,od2)
double precision, intent(in)  :: d1(:,:,:),d2(:,:,:)
double precision, intent(out) :: od1(:,:,:),od2(:,:,:)
call copy_3d_double_array1(d1,od1)
call copy_3d_double_array1(d2,od2)
end subroutine copy_3d_double_array2

!====================================================================

subroutine copy_3d_double_array3(d1,od1,d2,od2,d3,od3)
double precision, intent(in)  :: d1(:,:,:),d2(:,:,:),d3(:,:,:)
double precision, intent(out) :: od1(:,:,:),od2(:,:,:),od3(:,:,:)
call copy_3d_double_array2(d1,od1,d2,od2)
call copy_3d_double_array1(d3,od3)
end subroutine copy_3d_double_array3

!====================================================================

subroutine copy_1d_complex_array1(c1,oc1)
complex, intent(in)  :: c1(:)
complex, intent(out) :: oc1(:)
oc1=c1
end subroutine copy_1d_complex_array1

!====================================================================

subroutine copy_1d_complex_array2(c1,oc1,c2,oc2)
complex, intent(in)  :: c1(:),c2(:)
complex, intent(out) :: oc1(:),oc2(:)
call copy_1d_complex_array1(c1,oc1)
call copy_1d_complex_array1(c2,oc2)
end subroutine copy_1d_complex_array2

!====================================================================

subroutine copy_1d_complex_array3(c1,oc1,c2,oc2,c3,oc3)
complex, intent(in)  :: c1(:),c2(:),c3(:)
complex, intent(out) :: oc1(:),oc2(:),oc3(:)
call copy_1d_complex_array2(c1,oc1,c2,oc2)
call copy_1d_complex_array1(c3,oc3)
end subroutine copy_1d_complex_array3

!====================================================================

subroutine copy_2d_complex_array1(c1,oc1)
complex, intent(in)  :: c1(:,:)
complex, intent(out) :: oc1(:,:)
oc1=c1
end subroutine copy_2d_complex_array1

!====================================================================

subroutine copy_2d_complex_array2(c1,oc1,c2,oc2)
complex, intent(in)  :: c1(:,:),c2(:,:)
complex, intent(out) :: oc1(:,:),oc2(:,:)
call copy_2d_complex_array1(c1,oc1)
call copy_2d_complex_array1(c2,oc2)
end subroutine copy_2d_complex_array2

!====================================================================

subroutine copy_2d_complex_array3(c1,oc1,c2,oc2,c3,oc3)
complex, intent(in)  :: c1(:,:),c2(:,:),c3(:,:)
complex, intent(out) :: oc1(:,:),oc2(:,:),oc3(:,:)
call copy_2d_complex_array2(c1,oc1,c2,oc2)
call copy_2d_complex_array1(c3,oc3)
end subroutine copy_2d_complex_array3

!====================================================================

subroutine copy_3d_complex_array1(c1,oc1)
complex, intent(in)  :: c1(:,:,:)
complex, intent(out) :: oc1(:,:,:)
oc1=c1
end subroutine copy_3d_complex_array1

!====================================================================

subroutine copy_3d_complex_array2(c1,oc1,c2,oc2)
complex, intent(in)  :: c1(:,:,:),c2(:,:,:)
complex, intent(out) :: oc1(:,:,:),oc2(:,:,:)
call copy_3d_complex_array1(c1,oc1)
call copy_3d_complex_array1(c2,oc2)
end subroutine copy_3d_complex_array2

!====================================================================

subroutine copy_3d_complex_array3(c1,oc1,c2,oc2,c3,oc3)
complex, intent(in)  :: c1(:,:,:),c2(:,:,:),c3(:,:,:)
complex, intent(out) :: oc1(:,:,:),oc2(:,:,:),oc3(:,:,:)
call copy_3d_complex_array2(c1,oc1,c2,oc2)
call copy_3d_complex_array1(c3,oc3)
end subroutine copy_3d_complex_array3

!====================================================================

subroutine sub2ind_2d(n1,i1,i2,ind)
integer, intent(in) :: n1,i1,i2
integer, intent(out):: ind
ind=(i2-1)*n1+i1
end subroutine sub2ind_2d

!====================================================================

subroutine sub2ind_3d(n1,n2,i1,i2,i3,ind)
integer, intent(in) :: n1,n2,i1,i2,i3
integer, intent(out):: ind
ind=(i3-1)*(n1*n2)+(i2-1)*n1+i1
end subroutine sub2ind_3d

!====================================================================

subroutine sub2ind_4d(n1,n2,n3,i1,i2,i3,i4,ind)
integer, intent(in) :: n1,n2,n3,i1,i2,i3,i4
integer, intent(out):: ind
ind=(i4-1)*(n1*n2*n3)+(i3-1)*(n1*n2)+(i2-1)*n1+i1
end subroutine sub2ind_4d

!====================================================================

subroutine ind2sub_2d(n1,ind,i1,i2)
integer, intent(in) :: n1,ind
integer, intent(out):: i1,i2
if (mod(ind,n1).ne.0) then
   i2=int(ind/n1)+1
   i1=ind-int(ind/n1)*n1
else
   i2=int(ind/n1)
   i1=n1
endif
end subroutine ind2sub_2d

!====================================================================

subroutine ind2sub_3d(n1,n2,ind,i1,i2,i3)
integer, intent(in) :: n1,n2,ind
integer, intent(out):: i1,i2,i3
integer             :: ind2,n12
n12=n1*n2
if (mod(ind,n12).ne.0) then
   i3=int(ind/n12)+1
   ind2=ind-n12*(i3-1)
   if (mod(ind2,n1).ne.0) then
      i2=int(ind2/n1)+1
      i1=ind2-int(ind2/n1)*n1
   else
      i2=int(ind2/n1)
      i1=n1
   endif
else
   i3=int(ind/n12)
   i2=n2
   i1=n1
endif
end subroutine ind2sub_3d

!====================================================================

subroutine ind2sub_4d(n1,n2,n3,ind,i1,i2,i3,i4)
integer, intent(in) :: n1,n2,n3,ind
integer, intent(out):: i1,i2,i3,i4
integer             :: ind2,ind3,n123,n12
n123=n1*n2*n3
n12=n1*n2
if (mod(ind,n123).ne.0) then
   i4=int(ind/n123)+1
   ind3=ind-n123*(i4-1)
   if (mod(ind3,n12).ne.0) then
      i3=int(ind3/n12)+1
      ind2=ind3-n12*(i3-1)
      if (mod(ind2,n1).ne.0) then
         i2=int(ind2/n1)+1
         i1=ind2-int(ind2/n1)*n1
      else
         i2=int(ind2/n1)
         i1=n1
      endif
   else
      i3=int(ind3/n12)
      i2=n2
      i1=n1
   endif
else
   i3=int(ind/n12)
   i2=n2
   i1=n1
endif
end subroutine ind2sub_4d

!====================================================================

subroutine permute2d_real(d_in,d_out)
real, intent(in) :: d_in(:,:)
real, intent(out):: d_out(:,:)
d_out=transpose(d_in)
end subroutine permute2d_real

!====================================================================

subroutine permute2d_int(d_in,d_out)
integer, intent(in)  :: d_in(:,:)
integer, intent(out) :: d_out(:,:)
d_out=transpose(d_in)
end subroutine permute2d_int

!====================================================================

subroutine permute3d_real(d_in,d_out,o1,o2,o3) 
real, intent(in) :: d_in(:,:,:)
integer, intent(in) :: o1,o2,o3
real, intent(out) :: d_out(:,:,:)
integer i1, i2, i3, n1, n2, n3
if (o1*o2*o3/=6 .or. (o1+o2+o3)/=6) then
   write(*,*)"o1, o2 and o3 should be 3 different number of 1,2 and 3!"
   stop
endif
n1=size(d_in,1)
n2=size(d_in,2)
n3=size(d_in,3)
if (o1==1) then
   if (o2==2) then
      d_out=d_in 
   elseif (o2==3) then
      forall (i1=1:n1,i2=1:n3,i3=1:n2)   
         d_out(i1,i2,i3)=d_in(i1,i3,i2)
      end forall
   endif
elseif (o1==2) then
   if (o2==1) then
      forall (i1=1:n2,i2=1:n1,i3=1:n3)   
         d_out(i1,i2,i3)=d_in(i2,i1,i3)
      end forall
   elseif (o2==3) then
      forall (i1=1:n2,i2=1:n3,i3=1:n1)   
         d_out(i1,i2,i3)=d_in(i2,i3,i1)
      end forall
   endif
elseif (o1==3) then
   if (o2==1) then
      forall (i1=1:n3,i2=1:n1,i3=1:n2)   
         d_out(i1,i2,i3)=d_in(i3,i1,i2)
      end forall
   else
      forall (i1=1:n3,i2=1:n2,i3=1:n1)   
         d_out(i1,i2,i3)=d_in(i3,i2,i1)
      end forall
   endif
endif
end subroutine permute3d_real

!====================================================================

subroutine permute3d_int(d_in,d_out,o1,o2,o3) 
integer, intent(in) :: d_in(:,:,:)
integer, intent(in) :: o1,o2,o3
integer, intent(out) :: d_out(:,:,:)
integer i1, i2, i3, n1, n2, n3
if (o1*o2*o3/=6 .or. (o1+o2+o3)/=6) then
   write(*,*)"o1, o2 and o3 should be 3 different number of 1,2 and 3!"
   stop
endif
n1=size(d_in,1)
n2=size(d_in,2)
n3=size(d_in,3)
if (o1==1) then
   if (o2==2) then
      d_out=d_in 
   elseif (o2==3) then
      forall (i1=1:n1,i2=1:n3,i3=1:n2)   
         d_out(i1,i2,i3)=d_in(i1,i3,i2)
      end forall
   endif
elseif (o1==2) then
   if (o2==1) then
      forall (i1=1:n2,i2=1:n1,i3=1:n3)   
         d_out(i1,i2,i3)=d_in(i2,i1,i3)
      end forall
   elseif (o2==3) then
      forall (i1=1:n2,i2=1:n3,i3=1:n1)   
         d_out(i1,i2,i3)=d_in(i2,i3,i1)
      end forall
   endif
elseif (o1==3) then
   if (o2==1) then
      forall (i1=1:n3,i2=1:n1,i3=1:n2)   
         d_out(i1,i2,i3)=d_in(i3,i1,i2)
      end forall
   else
      forall (i1=1:n3,i2=1:n2,i3=1:n1)   
         d_out(i1,i2,i3)=d_in(i3,i2,i1)
      end forall
   endif
endif
end subroutine permute3d_int

!====================================================================

function arr2vec_real_2d(d_in,n1,n2)
real :: arr2vec_real_2d(n1*n2)
integer, intent(in) :: n1, n2
real,    intent(in) :: d_in(:,:)
integer :: i1, i2, ind, ind2
do i2=1,n2
   ind2=(i2-1)*n1
   do i1=1,n1
      ind=ind2+i1
      arr2vec_real_2d(ind)=d_in(i1,i2)
   enddo
enddo
end function arr2vec_real_2d

!====================================================================

function arr2vec_real_3d(d_in,n1,n2,n3)
integer, intent(in) :: n1,n2,n3
real :: arr2vec_real_3d(n1*n2*n3)
real, intent(in) :: d_in(:,:,:)
integer :: i1, i2, i3, ind, ind3, ind2, ind32, n12
n12=n1*n2
do i3=1,n3
   ind3=(i3-1)*n12
   do i2=1,n2
      ind2=(i2-1)*n1
      ind32=ind3+ind2
      do i1=1,n1
         ind=ind32+i1
         arr2vec_real_3d(ind)=d_in(i1,i2,i3)
      enddo
   enddo
enddo
end function arr2vec_real_3d

!====================================================================

function arr2vec_real_4d(d_in,n1,n2,n3,n4)
integer, intent(in) :: n1,n2,n3,n4
real :: arr2vec_real_4d(n1*n2*n3*n4)
real, intent(in) :: d_in(:,:,:,:)
integer :: i1, i2, i3, i4, ind, ind4, ind3, ind2, ind43, ind432, n12, n123
n123=n1*n2*n3
n12=n1*n2
do i4=1,n4
   ind4=(i4-1)*n123
   do i3=1,n3
      ind3=(i3-1)*n12
      ind43=ind4+ind3
      do i2=1,n2
         ind2=(i2-1)*n1
         ind432=ind43+ind2
         do i1=1,n1
            ind=ind432+i1
            arr2vec_real_4d(ind)=d_in(i1,i2,i3,i4)
         enddo
      enddo
   enddo
enddo
end function arr2vec_real_4d

!====================================================================

function arr2vec_int_2d(d_in,n1,n2)
integer :: arr2vec_int_2d(n1*n2)
integer, intent(in) :: n1, n2
integer,    intent(in) :: d_in(:,:)
integer :: i1, i2, ind, ind2
do i2=1,n2
   ind2=(i2-1)*n1
   do i1=1,n1
      ind=ind2+i1
      arr2vec_int_2d(ind)=d_in(i1,i2)
   enddo
enddo
end function arr2vec_int_2d

!====================================================================

function arr2vec_int_3d(d_in,n1,n2,n3)
integer, intent(in) :: n1,n2,n3
integer :: arr2vec_int_3d(n1*n2*n3)
integer, intent(in) :: d_in(:,:,:)
integer :: i1, i2, i3, ind, ind3, ind2, ind32, n12
n12=n1*n2
do i3=1,n3
   ind3=(i3-1)*n12
   do i2=1,n2
      ind2=(i2-1)*n1
      ind32=ind3+ind2
      do i1=1,n1
         ind=ind32+i1
         arr2vec_int_3d(ind)=d_in(i1,i2,i3)
      enddo
   enddo
enddo
end function arr2vec_int_3d

!====================================================================

function arr2vec_int_4d(d_in,n1,n2,n3,n4)
integer, intent(in) :: n1,n2,n3,n4
integer :: arr2vec_int_4d(n1*n2*n3*n4)
integer, intent(in) :: d_in(:,:,:,:)
integer :: i1, i2, i3, i4, ind, ind4, ind3, ind2, ind43, ind432, n12, n123
n123=n1*n2*n3
n12=n1*n2
do i4=1,n4
   ind4=(i4-1)*n123
   do i3=1,n3
      ind3=(i3-1)*n12
      ind43=ind4+ind3
      do i2=1,n2
         ind2=(i2-1)*n1
         ind432=ind43+ind2
         do i1=1,n1
            ind=ind432+i1
            arr2vec_int_4d(ind)=d_in(i1,i2,i3,i4)
         enddo
      enddo
   enddo
enddo
end function arr2vec_int_4d

!====================================================================

function arr2vec_double_2d(d_in,n1,n2)
double precision :: arr2vec_double_2d(n1*n2)
integer, intent(in) :: n1, n2
double precision,intent(in) :: d_in(:,:)
integer :: i1, i2, ind, ind2
do i2=1,n2
   ind2=(i2-1)*n1
   do i1=1,n1
      ind=ind2+i1
      arr2vec_double_2d(ind)=d_in(i1,i2)
   enddo
enddo
end function arr2vec_double_2d

!====================================================================

function arr2vec_double_3d(d_in,n1,n2,n3)
integer, intent(in) :: n1,n2,n3
double precision :: arr2vec_double_3d(n1*n2*n3)
double precision, intent(in) :: d_in(:,:,:)
integer :: i1, i2, i3, ind, ind3, ind2, ind32, n12
n12=n1*n2
do i3=1,n3
   ind3=(i3-1)*n12
   do i2=1,n2
      ind2=(i2-1)*n1
      ind32=ind3+ind2
      do i1=1,n1
         ind=ind32+i1
         arr2vec_double_3d(ind)=d_in(i1,i2,i3)
      enddo
   enddo
enddo
end function arr2vec_double_3d

!====================================================================

function arr2vec_double_4d(d_in,n1,n2,n3,n4)
integer, intent(in) :: n1,n2,n3,n4
double precision :: arr2vec_double_4d(n1*n2*n3*n4)
double precision, intent(in) :: d_in(:,:,:,:)
integer :: i1, i2, i3, i4, ind, ind4, ind3, ind2, ind43, ind432, n12, n123
n123=n1*n2*n3
n12=n1*n2
do i4=1,n4
   ind4=(i4-1)*n123
   do i3=1,n3
      ind3=(i3-1)*n12
      ind43=ind4+ind3
      do i2=1,n2
         ind2=(i2-1)*n1
         ind432=ind43+ind2
         do i1=1,n1
            ind=ind432+i1
            arr2vec_double_4d(ind)=d_in(i1,i2,i3,i4)
         enddo
      enddo
   enddo
enddo
end function arr2vec_double_4d

!====================================================================

function vec2arr_real_2d(d_in,n1,n2)
real :: vec2arr_real_2d(n1,n2)
integer, intent(in) :: n1, n2
real,    intent(in) :: d_in(:)
integer :: i1, i2, ind, ind2
do i2=1,n2
   ind2=(i2-1)*n1
   do i1=1,n1
      ind=ind2+i1
      vec2arr_real_2d(i1,i2)=d_in(ind)
   enddo
enddo
end function vec2arr_real_2d

!====================================================================

function vec2arr_real_3d(d_in,n1,n2,n3)
integer, intent(in) :: n1,n2,n3
real :: vec2arr_real_3d(n1,n2,n3)
real, intent(in) :: d_in(:)
integer :: i1, i2, i3, ind, ind3, ind2, ind32, n12
n12=n1*n2
do i3=1,n3
   ind3=(i3-1)*n12
   do i2=1,n2
      ind2=(i2-1)*n1
      ind32=ind3+ind2
      do i1=1,n1
         ind=ind32+i1
         vec2arr_real_3d(i1,i2,i3)=d_in(ind)
      enddo
   enddo
enddo
end function vec2arr_real_3d

!====================================================================

function vec2arr_real_4d(d_in,n1,n2,n3,n4)
integer, intent(in) :: n1,n2,n3,n4
real :: vec2arr_real_4d(n1,n2,n3,n4)
real, intent(in) :: d_in(:)
integer :: i1, i2, i3, i4, ind, ind4, ind3, ind2, ind43, ind432, n12, n123
n123=n1*n2*n3
n12=n1*n2
do i4=1,n4
   ind4=(i4-1)*n123
   do i3=1,n3
      ind3=(i3-1)*n12
      ind43=ind4+ind3
      do i2=1,n2
         ind2=(i2-1)*n1
         ind432=ind43+ind2
         do i1=1,n1
            ind=ind432+i1
            vec2arr_real_4d(i1,i2,i3,i4)=d_in(ind)
         enddo
      enddo
   enddo
enddo
end function vec2arr_real_4d

!====================================================================

function vec2arr_int_2d(d_in,n1,n2)
integer :: vec2arr_int_2d(n1,n2)
integer, intent(in) :: n1, n2
integer, intent(in) :: d_in(:)
integer :: i1, i2, ind, ind2
do i2=1,n2
   ind2=(i2-1)*n1
   do i1=1,n1
      ind=ind2+i1
      vec2arr_int_2d(i1,i2)=d_in(ind)
   enddo
enddo
end function vec2arr_int_2d

!====================================================================

function vec2arr_int_3d(d_in,n1,n2,n3)
integer, intent(in) :: n1,n2,n3
integer :: vec2arr_int_3d(n1,n2,n3)
integer, intent(in) :: d_in(:)
integer :: i1, i2, i3, ind, ind3, ind2, ind32,n12
n12=n1*n2
do i3=1,n3
   ind3=(i3-1)*n12
   do i2=1,n2
      ind2=(i2-1)*n1
      ind32=ind3+ind2
      do i1=1,n1
         ind=ind32+i1
         vec2arr_int_3d(i1,i2,i3)=d_in(ind)
      enddo
   enddo
enddo
end function vec2arr_int_3d

!====================================================================

function vec2arr_int_4d(d_in,n1,n2,n3,n4)
integer, intent(in) :: n1,n2,n3,n4
integer :: vec2arr_int_4d(n1,n2,n3,n4)
integer, intent(in) :: d_in(:)
integer :: i1, i2, i3, i4, ind, ind4, ind3, ind2, ind43, ind432, n12, n123
n123=n1*n2*n3
n12=n1*n2
do i4=1,n4
   ind4=(i4-1)*n123
   do i3=1,n3
      ind3=(i3-1)*n12
      ind43=ind4+ind3
      do i2=1,n2
         ind2=(i2-1)*n1
         ind432=ind43+ind2
         do i1=1,n1
            ind=ind432+i1
            vec2arr_int_4d(i1,i2,i3,i4)=d_in(ind)
         enddo
      enddo
   enddo
enddo
end function vec2arr_int_4d

!====================================================================

function vec2arr_double_2d(d_in,n1,n2)
double precision :: vec2arr_double_2d(n1,n2)
integer, intent(in) :: n1, n2
double precision, intent(in) :: d_in(:)
integer :: i1, i2, ind, ind2
do i2=1,n2
   ind2=(i2-1)*n1
   do i1=1,n1
      ind=ind2+i1
      vec2arr_double_2d(i1,i2)=d_in(ind)
   enddo
enddo
end function vec2arr_double_2d

!====================================================================

function vec2arr_double_3d(d_in,n1,n2,n3)
integer, intent(in) :: n1,n2,n3
double precision :: vec2arr_double_3d(n1,n2,n3)
double precision, intent(in) :: d_in(:)
integer :: i1, i2, i3, ind, ind3, ind2, ind32,n12
n12=n1*n2
do i3=1,n3
   ind3=(i3-1)*n12
   do i2=1,n2
      ind2=(i2-1)*n1
      ind32=ind3+ind2
      do i1=1,n1
         ind=ind32+i1
         vec2arr_double_3d(i1,i2,i3)=d_in(ind)
      enddo
   enddo
enddo
end function vec2arr_double_3d

!====================================================================

function vec2arr_double_4d(d_in,n1,n2,n3,n4)
integer, intent(in) :: n1,n2,n3,n4
double precision :: vec2arr_double_4d(n1,n2,n3,n4)
double precision, intent(in) :: d_in(:)
integer :: i1, i2, i3, i4, ind, ind4, ind3, ind2, ind43, ind432, n12, n123
n123=n1*n2*n3
n12=n1*n2
do i4=1,n4
   ind4=(i4-1)*n123
   do i3=1,n3
      ind3=(i3-1)*n12
      ind43=ind4+ind3
      do i2=1,n2
         ind2=(i2-1)*n1
         ind432=ind43+ind2
         do i1=1,n1
            ind=ind432+i1
            vec2arr_double_4d(i1,i2,i3,i4)=d_in(ind)
         enddo
      enddo
   enddo
enddo
end function vec2arr_double_4d

!====================================================================

function flipud_int_1d(d_in,n1)
integer :: flipud_int_1d(n1) 
integer,         intent(in) :: d_in(:),n1
integer                     :: i1
do i1=1,n1
   flipud_int_1d(i1)=d_in(n1-i1+1)
enddo
end function flipud_int_1d

!====================================================================

function flipud_real_1d(d_in,n1)
real :: flipud_real_1d(n1) 
real,     intent(in) :: d_in(:)
integer,  intent(in) :: n1
integer              :: i1
do i1=1,n1
   flipud_real_1d(i1)=d_in(n1-i1+1)
enddo
end function flipud_real_1d

!====================================================================

function flipud_double_1d(d_in,n1)
double precision :: flipud_double_1d(n1) 
double precision,  intent(in) :: d_in(:)
integer,  intent(in) :: n1
integer              :: i1
do i1=1,n1
   flipud_double_1d(i1)=d_in(n1-i1+1)
enddo
end function flipud_double_1d

!====================================================================

function flipud_int_2d(d_in,n1,n2)
integer :: flipud_int_2d(n1,n2)
integer,   intent(in) :: d_in(:,:),n1,n2
integer               :: i1
do i1=1,n1
   flipud_int_2d(i1,:)=d_in(n1-i1+1,:)
enddo
end function flipud_int_2d

!====================================================================

function flipud_real_2d(d_in,n1,n2)
real :: flipud_real_2d(n1,n2)
real,   intent(in) :: d_in(:,:)
integer,intent(in) :: n1,n2
integer            :: i1
do i1=1,n1
   flipud_real_2d(i1,:)=d_in(n1-i1+1,:)
enddo
end function flipud_real_2d

!====================================================================

function flipud_double_2d(d_in,n1,n2)
double precision :: flipud_double_2d(n1,n2)
double precision,intent(in) :: d_in(:,:)
integer,intent(in) :: n1,n2
integer            :: i1
do i1=1,n1
   flipud_double_2d(i1,:)=d_in(n1-i1+1,:)
enddo
end function flipud_double_2d

!====================================================================

function fliplr_int_2d(d_in,n1,n2)
integer :: fliplr_int_2d(n1,n2)
integer,intent(in) :: d_in(:,:),n1,n2
integer            :: i2
do i2=1,n2
   fliplr_int_2d(:,i2)=d_in(:,n2-i2+1)
enddo
end function fliplr_int_2d

!====================================================================

function fliplr_real_2d(d_in,n1,n2)
real :: fliplr_real_2d(n1,n2)
real   ,intent(in) :: d_in(:,:)
integer,intent(in) :: n1, n2
integer            :: i2
do i2=1,n2
   fliplr_real_2d(:,i2)=d_in(:,n2-i2+1)
enddo
end function fliplr_real_2d

!====================================================================

function fliplr_double_2d(d_in,n1,n2)
double precision :: fliplr_double_2d(n1,n2)
double precision,intent(in) :: d_in(:,:)
integer,intent(in) :: n1, n2
integer            :: i2
do i2=1,n2
   fliplr_double_2d(:,i2)=d_in(:,n2-i2+1)
enddo
end function fliplr_double_2d

!====================================================================

!subroutine rotate_int_2d(d_in,n1,n2,o1,o2,deg,d_out)
subroutine rotate_int_2d(d_in,n1,n2,o2,deg,d_out)
!integer,     intent(in)  :: d_in(:,:), n1,n2,o1,o2,deg
integer,     intent(in)  :: d_in(:,:), n1,n2,o2,deg
integer,     intent(out) :: d_out(:,:)
if (o2.eq.1) then
   if (deg==1) then
      d_out=fliplr(transpose(d_in),n2,n1)
   elseif (deg==2) then
      d_out=fliplr(flipud(d_in,n1,n2),n1,n2)
   elseif (deg==3) then
      d_out=flipud(transpose(d_in),n2,n1)
   endif
elseif (o2.eq.n2) then
   if (deg==1) then
      d_out=transpose(fliplr(d_in,n1,n2))
   elseif (deg==2) then
      d_out=flipud(fliplr(d_in,n1,n2),n1,n2)
   elseif (deg==3) then
      d_out=flipud(transpose(d_in),n2,n1)
   endif
endif
end subroutine rotate_int_2d

!====================================================================

!subroutine rotate_real_2d(d_in,n1,n2,o1,o2,deg,d_out)
subroutine rotate_real_2d(d_in,n1,n2,o2,deg,d_out)
real,        intent(in)  :: d_in(:,:)
!integer,     intent(in)  :: n1,n2,o1,o2,deg
integer,     intent(in)  :: n1,n2,o2,deg
real,        intent(out) :: d_out(:,:)
if (o2.eq.1) then
   if (deg==1) then
      d_out=fliplr(transpose(d_in),n2,n1)
   elseif (deg==2) then
      d_out=fliplr(flipud(d_in,n1,n2),n1,n2)
   elseif (deg==3) then
      d_out=flipud(transpose(d_in),n2,n1)
   endif
elseif (o2.eq.n2) then
   if (deg==1) then
      d_out=transpose(fliplr(d_in,n1,n2))
   elseif (deg==2) then
      d_out=flipud(fliplr(d_in,n1,n2),n1,n2)
   elseif (deg==3) then
      d_out=flipud(transpose(d_in),n2,n1)
   endif
endif
end subroutine rotate_real_2d

!====================================================================

!subroutine rotate_double_2d(d_in,n1,n2,o1,o2,deg,d_out)
subroutine rotate_double_2d(d_in,n1,n2,o2,deg,d_out)
double precision, intent(in)  :: d_in(:,:)
!integer,          intent(in)  :: n1,n2,o1,o2,deg
integer,          intent(in)  :: n1,n2,o2,deg
double precision, intent(out) :: d_out(:,:)
if (o2.eq.1) then
   if (deg==1) then
      d_out=fliplr(transpose(d_in),n2,n1)
   elseif (deg==2) then
      d_out=fliplr(flipud(d_in,n1,n2),n1,n2)
   elseif (deg==3) then
      d_out=flipud(transpose(d_in),n2,n1)
   endif
elseif (o2.eq.n2) then
   if (deg==1) then
      d_out=transpose(fliplr(d_in,n1,n2))
   elseif (deg==2) then
      d_out=flipud(fliplr(d_in,n1,n2),n1,n2)
   elseif (deg==3) then
      d_out=flipud(transpose(d_in),n2,n1)
   endif
endif
end subroutine rotate_double_2d

!====================================================================

!subroutine locate_int_1d(d_in,di,n1,i1)
!integer, intent(in)  :: d_in(:),di,n1
!integer, intent(out) :: i1
!do i=1,n1
!   if (d_in(i).eq.di) then
!      i1=i
!      goto 100
!enddo
!100 continue
!return
!end subroutine locate_int_1d
      
!====================================================================

SUBROUTINE indexx(n,arr,indx)
INTEGER n,indx(n),M,NSTACK
REAL arr(n)
PARAMETER (M=7,NSTACK=50)
INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
REAL a
do j=1,n
   indx(j)=j
enddo
jstack=0
l=1
ir=n
1 continue
if (ir-l.lt.M) then
   do j=l+1,ir
      indxt=indx(j)
      a=arr(indxt)
      do i=j-1,1,-1
         if (arr(indx(i)).le.a) goto 2
         indx(i+1)=indx(i)
      enddo
      i=0
2  continue
      indx(i+1)=indxt
   enddo
   if (jstack.eq.0) return
   ir=istack(jstack)
   l=istack(jstack-1)
   jstack=jstack-2
else
   k=(l+ir)/2
   itemp=indx(k)
   indx(k)=indx(l+1)
   indx(l+1)=itemp
   if (arr(indx(l+1)).gt.arr(indx(ir))) then
      itemp=indx(l+1)
      indx(l+1)=indx(ir)
      indx(ir)=itemp
   endif
   if (arr(indx(l)).gt.arr(indx(ir))) then
      itemp=indx(l)
      indx(l)=indx(ir)
      indx(ir)=itemp
   endif
   if (arr(indx(l+1)).gt.arr(indx(l))) then
      itemp=indx(l+1)
      indx(l+1)=indx(l)
      indx(l)=itemp
   endif
   i=l+1
   j=ir
   indxt=indx(l)
   a=arr(indxt)
3  continue
   i=i+1
   if (arr(indx(i)).lt.a) goto 3
4  continue
   j=j-1
   if (arr(indx(j)).gt.a) goto 4
   if (j.lt.i) goto 5
   itemp=indx(i)
   indx(i)=indx(j)
   indx(j)=itemp
   goto 3
5  indx(l)=indx(j)
   indx(j)=indxt
   jstack=jstack+2
   if (jstack.gt.NSTACK) then
      write(*,*)"NSTACK too small in indexx"
      stop
   endif
   if (ir-i+1.ge.j-l) then
      istack(jstack)=ir
      istack(jstack-1)=i
      ir=j-1
   else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
   endif
endif
goto 1
END SUBROUTINE indexx

end module module_array
