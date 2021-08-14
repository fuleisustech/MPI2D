module module_array_compose
implicit none

interface compose_funs
   module procedure compose_2v_to_1a_int_fun   
   module procedure compose_3v_to_1a_int_fun
   module procedure compose_4v_to_1a_int_fun
   module procedure compose_2_1d_to_1_1d_int_fun_v1
   module procedure compose_2_1d_to_1_1d_int_fun_v2
   module procedure compose_2_1d_to_1_2d_int_fun

   module procedure compose_2v_to_1a_real_fun   
   module procedure compose_3v_to_1a_real_fun
   module procedure compose_4v_to_1a_real_fun
   module procedure compose_2_1d_to_1_1d_real_fun_v1
   module procedure compose_2_1d_to_1_1d_real_fun_v2
   module procedure compose_2_1d_to_1_2d_real_fun
end interface compose_funs

interface compose_subs
   module procedure compose_2v_to_1a_int_sub   
   module procedure compose_3v_to_1a_int_sub
   module procedure compose_4v_to_1a_int_sub
   module procedure compose_2_1d_to_1_1d_int_sub_v1
   module procedure compose_2_1d_to_1_1d_int_sub_v2
   module procedure compose_2_1d_to_1_1d_int_sub_v3
   module procedure compose_2_1d_to_1_2d_int_sub
   module procedure compose_1_1d_and_1_2d_to_1_2d_int_sub_v1
   module procedure compose_1_1d_and_1_2d_to_1_2d_int_sub_v2
   module procedure compose_2_2d_to_1_2d_int_sub
   
   module procedure compose_2v_to_1a_real_sub   
   module procedure compose_3v_to_1a_real_sub
   module procedure compose_4v_to_1a_real_sub
   module procedure compose_2_1d_to_1_1d_real_sub_v1
   module procedure compose_2_1d_to_1_1d_real_sub_v2
   module procedure compose_2_1d_to_1_1d_real_sub_v3
   module procedure compose_2_1d_to_1_2d_real_sub
   module procedure compose_1_1d_and_1_2d_to_1_2d_real_sub_v1
   module procedure compose_1_1d_and_1_2d_to_1_2d_real_sub_v2
   module procedure compose_2_2d_to_1_2d_real_sub
end interface compose_subs

contains

function compose_2v_to_1a_int_fun(d1,d2)   
integer :: compose_2v_to_1a_int_fun(2)
integer, intent(in) :: d1, d2
compose_2v_to_1a_int_fun(1)=d1
compose_2v_to_1a_int_fun(2)=d2
end function compose_2v_to_1a_int_fun

subroutine compose_2v_to_1a_int_sub(d1,d2,d3)   
integer, intent(in) :: d1, d2
integer, intent(out) :: d3(2)
d3(1)=d1
d3(2)=d2
end subroutine compose_2v_to_1a_int_sub

function compose_3v_to_1a_int_fun(d1,d2,d3)   
integer :: compose_3v_to_1a_int_fun(3)
integer, intent(in) :: d1, d2, d3
compose_3v_to_1a_int_fun(1)=d1
compose_3v_to_1a_int_fun(2)=d2
compose_3v_to_1a_int_fun(3)=d3
end function compose_3v_to_1a_int_fun

subroutine compose_3v_to_1a_int_sub(d1,d2,d3,d4)   
integer, intent(in) :: d1, d2, d3
integer, intent(out) :: d4(3)
d4(1)=d1
d4(2)=d2
d4(3)=d3
end subroutine compose_3v_to_1a_int_sub

function compose_4v_to_1a_int_fun(d1,d2,d3,d4)   
integer :: compose_4v_to_1a_int_fun(2,2)
integer, intent(in) :: d1, d2, d3, d4
compose_4v_to_1a_int_fun(:,1)=compose_funs(d1,d2)
compose_4v_to_1a_int_fun(:,2)=compose_funs(d3,d4)
end function compose_4v_to_1a_int_fun

subroutine compose_4v_to_1a_int_sub(d1,d2,d3,d4,d5)   
integer, intent(in) :: d1, d2, d3, d4
integer, intent(out) :: d5(2,2)
d5(:,1)=compose_funs(d1,d2)
d5(:,2)=compose_funs(d3,d4)
end subroutine compose_4v_to_1a_int_sub

function compose_2_1d_to_1_1d_int_fun_v1(d1,d2,n1,n2)
integer, intent(in) :: d1(n1), d2(n2), n1, n2
integer :: compose_2_1d_to_1_1d_int_fun_v1(n1+n2)
compose_2_1d_to_1_1d_int_fun_v1(1:n1)=d1
compose_2_1d_to_1_1d_int_fun_v1(n1+1:n1+n2)=d2
end function compose_2_1d_to_1_1d_int_fun_v1

function compose_2_1d_to_1_1d_int_fun_v2(d1,d2,n1,n2,n3)
integer, intent(in) :: d1(n1), d2(n2), n1, n2, n3
integer :: compose_2_1d_to_1_1d_int_fun_v2(n3)
compose_2_1d_to_1_1d_int_fun_v2(1:n1)=d1
compose_2_1d_to_1_1d_int_fun_v2(n1+1:n3)=d2
end function compose_2_1d_to_1_1d_int_fun_v2

subroutine compose_2_1d_to_1_1d_int_sub_v1(d1,d2,n1,n2,n3,d3)
integer, intent(in) :: d1(n1), d2(n2), n1, n2, n3
integer,intent(out) :: d3(n3)
d3(1:n1)=d1
d3(n1+1:n3)=d2
end subroutine compose_2_1d_to_1_1d_int_sub_v1

subroutine compose_2_1d_to_1_1d_int_sub_v2(d1,d2,n1,n2,d3)
integer, intent(in) :: d1(n1), d2(n2), n1, n2
integer,intent(out) :: d3(n1+n2)
d3(1:n1)=d1
d3(n1+1:n1+n2)=d2
end subroutine compose_2_1d_to_1_1d_int_sub_v2

subroutine compose_2_1d_to_1_1d_int_sub_v3(d1,d2,d3)
integer, intent(in) :: d1(:), d2(:)
integer,intent(out) :: d3(:)
integer :: n1, n2
n1=size(d1)
n2=size(d2)
d3(1:n1)=d1
d3(n1+1:n1+n2)=d2
end subroutine compose_2_1d_to_1_1d_int_sub_v3

function compose_2_1d_to_1_2d_int_fun(d1,d2,n1)
integer, intent(in) :: d1(n1), d2(n1), n1
integer :: compose_2_1d_to_1_2d_int_fun(n1,2)
compose_2_1d_to_1_2d_int_fun(:,1)=d1
compose_2_1d_to_1_2d_int_fun(:,2)=d2
end function compose_2_1d_to_1_2d_int_fun

subroutine compose_2_1d_to_1_2d_int_sub(d1,d2,d3)
integer, intent(in) :: d1(:), d2(:)
integer, intent(out):: d3(:,:)
d3(:,1)=d1
d3(:,2)=d2
end subroutine compose_2_1d_to_1_2d_int_sub

subroutine compose_1_1d_and_1_2d_to_1_2d_int_sub_v1(d1,d2,d3)
integer, intent(in) :: d1(:),d2(:,:)
integer, intent(out) :: d3(:,:)
integer :: n2
n2=size(d2,2)
d3(:,1)=d1
d3(:,2:n2)=d2
end subroutine compose_1_1d_and_1_2d_to_1_2d_int_sub_v1

subroutine compose_1_1d_and_1_2d_to_1_2d_int_sub_v2(d1,d2,d3)
integer, intent(in) :: d1(:,:),d2(:)
integer, intent(out) :: d3(:,:)
integer :: n1
n1=size(d1,2)
d3(:,1:n1)=d1
d3(:,n1+1)=d2
end subroutine compose_1_1d_and_1_2d_to_1_2d_int_sub_v2

subroutine compose_2_2d_to_1_2d_int_sub(d1,d2,d3)
integer, intent(in) :: d1(:,:),d2(:,:)
integer, intent(out) :: d3(:,:)
integer :: n1(2),n2(2)
n1=shape(d1)
n2=shape(d2)
d3(:,1:n1(2))=d1
d3(:,n1(2)+1:n1(2)+n2(2))=d2
end subroutine compose_2_2d_to_1_2d_int_sub

!=============================================================

function compose_2v_to_1a_real_fun(d1,d2)   
real :: compose_2v_to_1a_real_fun(2)
real, intent(in) :: d1, d2
compose_2v_to_1a_real_fun(1)=d1
compose_2v_to_1a_real_fun(2)=d2
end function compose_2v_to_1a_real_fun

subroutine compose_2v_to_1a_real_sub(d1,d2,d3)   
real, intent(in)  :: d1, d2
real, intent(out) :: d3(2)
d3(1)=d1
d3(2)=d2
end subroutine compose_2v_to_1a_real_sub

function compose_3v_to_1a_real_fun(d1,d2,d3)   
real :: compose_3v_to_1a_real_fun(3)
real, intent(in) :: d1, d2, d3
compose_3v_to_1a_real_fun(1)=d1
compose_3v_to_1a_real_fun(2)=d2
compose_3v_to_1a_real_fun(3)=d3
end function compose_3v_to_1a_real_fun

subroutine compose_3v_to_1a_real_sub(d1,d2,d3,d4)   
real, intent(in)  :: d1, d2, d3
real, intent(out) :: d4(3)
d4(1)=d1
d4(2)=d2
d4(3)=d3
end subroutine compose_3v_to_1a_real_sub

function compose_4v_to_1a_real_fun(d1,d2,d3,d4)   
real :: compose_4v_to_1a_real_fun(2,2)
real, intent(in) :: d1, d2, d3, d4
compose_4v_to_1a_real_fun(:,1)=compose_funs(d1,d2)
compose_4v_to_1a_real_fun(:,2)=compose_funs(d3,d4)
end function compose_4v_to_1a_real_fun

subroutine compose_4v_to_1a_real_sub(d1,d2,d3,d4,d5)   
real, intent(in)  :: d1, d2, d3, d4
real, intent(out) :: d5(2,2)
d5(:,1)=compose_funs(d1,d2)
d5(:,2)=compose_funs(d3,d4)
end subroutine compose_4v_to_1a_real_sub

function compose_2_1d_to_1_1d_real_fun_v1(d1,d2,n1,n2)
real,    intent(in) :: d1(n1), d2(n2)
integer, intent(in) :: n1, n2
real :: compose_2_1d_to_1_1d_real_fun_v1(n1+n2)
compose_2_1d_to_1_1d_real_fun_v1(1:n1)=d1
compose_2_1d_to_1_1d_real_fun_v1(n1+1:n1+n2)=d2
end function compose_2_1d_to_1_1d_real_fun_v1

function compose_2_1d_to_1_1d_real_fun_v2(d1,d2,n1,n2,n3)
real,    intent(in) :: d1(n1), d2(n2)
integer, intent(in) :: n1, n2, n3
integer :: compose_2_1d_to_1_1d_real_fun_v2(n3)
compose_2_1d_to_1_1d_real_fun_v2(1:n1)=d1
compose_2_1d_to_1_1d_real_fun_v2(n1+1:n3)=d2
end function compose_2_1d_to_1_1d_real_fun_v2

subroutine compose_2_1d_to_1_1d_real_sub_v1(d1,d2,n1,n2,n3,d3)
real,    intent(in)  :: d1(n1), d2(n2)
integer, intent(in)  :: n1, n2, n3
real,    intent(out) :: d3(n3)
d3(1:n1)=d1
d3(n1+1:n3)=d2
end subroutine compose_2_1d_to_1_1d_real_sub_v1

subroutine compose_2_1d_to_1_1d_real_sub_v2(d1,d2,n1,n2,d3)
real,    intent(in)  :: d1(n1), d2(n2)
integer, intent(in)  :: n1, n2
real,    intent(out) :: d3(n1+n2)
d3(1:n1)=d1
d3(n1+1:n1+n2)=d2
end subroutine compose_2_1d_to_1_1d_real_sub_v2

subroutine compose_2_1d_to_1_1d_real_sub_v3(d1,d2,d3)
real,   intent(in)  :: d1(:), d2(:)
real,   intent(out) :: d3(:)
integer :: n1, n2
n1=size(d1)
n2=size(d2)
d3(1:n1)=d1
d3(n1+1:n1+n2)=d2
end subroutine compose_2_1d_to_1_1d_real_sub_v3

function compose_2_1d_to_1_2d_real_fun(d1,d2,n1)
integer, intent(in) :: n1
real,    intent(in) :: d1(n1), d2(n1)
real ::  compose_2_1d_to_1_2d_real_fun(n1,2)
compose_2_1d_to_1_2d_real_fun(:,1)=d1
compose_2_1d_to_1_2d_real_fun(:,2)=d2
end function compose_2_1d_to_1_2d_real_fun

subroutine compose_2_1d_to_1_2d_real_sub(d1,d2,d3)
real,    intent(in) :: d1(:), d2(:)
real,    intent(out):: d3(:,:)
d3(:,1)=d1
d3(:,2)=d2
end subroutine compose_2_1d_to_1_2d_real_sub

subroutine compose_1_1d_and_1_2d_to_1_2d_real_sub_v1(d1,d2,d3)
real, intent(in) :: d1(:),d2(:,:)
real, intent(out) :: d3(:,:)
integer :: n2
n2=size(d2,2)
d3(:,1)=d1
d3(:,2:n2)=d2
end subroutine compose_1_1d_and_1_2d_to_1_2d_real_sub_v1

subroutine compose_1_1d_and_1_2d_to_1_2d_real_sub_v2(d1,d2,d3)
real, intent(in) :: d1(:,:),d2(:)
real, intent(out) :: d3(:,:)
integer :: n1
n1=size(d1,2)
d3(:,1:n1)=d1
d3(:,n1+1)=d2
end subroutine compose_1_1d_and_1_2d_to_1_2d_real_sub_v2

subroutine compose_2_2d_to_1_2d_real_sub(d1,d2,d3)
real, intent(in) :: d1(:,:),d2(:,:)
real, intent(out) :: d3(:,:)
integer :: n1(2),n2(2)
n1=shape(d1)
n2=shape(d2)
d3(:,1:n1(2))=d1
d3(:,n1(2)+1:n1(2)+n2(2))=d2
end subroutine compose_2_2d_to_1_2d_real_sub

end module module_array_compose
