module module_pweikfd3d_dir
use module_global, only: PI
use module_utility
use module_io
use module_string
use module_array_compose

implicit none
integer, private    :: nz, ny, nx ,n1, n2, n
real,    private    :: h,hh,tm
real,pointer,private:: s(:,:,:),t(:,:,:)
real,private,parameter :: c1=1.0/sqrt(2.0)
real, private, allocatable :: s1(:,:),s2(:,:),t1(:,:),t2(:,:),t1sort(:)
integer, private, allocatable :: indx(:)
!real,private,target,allocatable :: ss(:,:,:),tt(:,:,:)
real, private :: s11(2,2),s12(2,2),t11(2,2),t12(2,2)
real, private :: s21(3,2),s22(3,2),t21(3,2),t22(3,2)
real, private :: s31(3,3),s32(3,3),t31(3,3),t32(3,3)
logical, private :: is_t2near(4), is_t2far(4)
real, private, allocatable :: es1(:),es2(:),es3(:),es4(:)
real, private, allocatable :: et1(:),et2(:),et3(:),et4(:)

private :: surface_expand, eadge_expand, treat1, treat2, treat3, treat1_ur,&
   treat1_dr, treat1_dl, treat1_ul, treat2_u, treat2_r, treat2_d, treat2_l,&
   treat3_c, is_near, is_far, count_me

contains

!======================================================================
! the main subroutine for 3D FD solver of direct arrival of eikonal eqn
!======================================================================
subroutine pweikfd3d_dir(s0,t0,theta_x,xref,theta_y,yref,zp,nz0,ny0,nx0,h0)
real,target,intent(in)   :: s0(:,:,:)
real,intent(in)          :: theta_x,theta_y,xref,yref,h0,zp
integer, intent(in)      :: nz0,ny0,nx0
real,target,intent(out)  :: t0(:,:,:)
integer                  :: nmax,ix,iy,iz,i_outer,i_inner,izp
real :: sin_theta_x, sin_theta_y,temp_x,temp_y
call copy_variable(h0,h)
call copy_variable(nz0,nz,ny0,ny,nx0,nx)
tm=0.0001
s=>s0
t=>t0
hh=h*h
nmax=max(ny,nx)
t=0.0
! decide the initial value at the plane
izp=nint(zp/h)+1
sin_theta_x=sin(theta_x/180.0*PI)
sin_theta_y=sin(theta_y/180.0*PI)
do ix=1,nx
   temp_x=((ix-1)*h-xref)*sin_theta_x
   do iy=1,ny
      temp_y=((iy-1)*h-yref)*sin_theta_y
      t(izp,iy,ix)=(temp_x+temp_y)*s(izp,iy,ix)
   enddo
enddo
call allocate_and_initial(s1,s2,t1,t2,ny,nx)
call allocate_and_initial(es1,es2,es3,es4,nmax)
call allocate_and_initial(et1,et2,et3,et4,nmax)
call allocate_and_initial(t1sort,ny*nx)
call allocate_and_initial(indx,ny*nx)
!t=t+tm+0.00001
!goto 999
do iz=1,max(nz-izp,izp-1)
!  deal with the top face
   if (izp-iz>=1) then
      i_outer=izp-iz
      i_inner=i_outer+1
      s1=s(i_inner,:,:)
      t1=t(i_inner,:,:)
      s2=s(i_outer,:,:)
      t2=0.0
      call surface_expand
      t(i_outer,:,:)=t2
   endif
! the end of the dealing with the top face

! deal with the bottom face
   if (izp+iz<=nz) then
      i_outer=izp+iz
      i_inner=i_outer-1
      s1=s(i_inner,:,:)
      t1=t(i_inner,:,:)
      s2=s(i_outer,:,:)
      t2=0.0
      call surface_expand()
      t(i_outer,:,:)=t2
   endif
! end dealing with the bottom face
enddo

999 continue
call deallocate_and_free(s1,s2,t1,t2)
call deallocate_and_free(es1,es2,es3,es4)
call deallocate_and_free(et1,et2,et3,et4)
call deallocate_and_free(t1sort)
call deallocate_and_free(indx)
end subroutine pweikfd3d_dir

! total area has been divided into 9 parts
!                n2
!   o----o-----------------o----o   expand order: 1
!   | 9  |        2        | 6  |          -> 2,3,4 and 5
!   o----o-----------------o----o          -> 6,7,8 and 9
!   |    |                 |    |
!   |    |                 |    |
!n1 | 5  |        1        | 3  |
!   |    |                 |    |
!   |    |                 |    |
!   |    |                 |    |
!   o----o-----------------o----o
!   | 8  |        4        | 7  |
!   o----o-----------------o----o
!
subroutine surface_expand()
!integer, external :: count
integer :: nsort, i, i1, i2, count_is_t2near, count_is_t2far
! ----------------- dealing with area 1
nsort=(ny-2)*(nx-2)
t1sort(1:nsort)=arr2vec(t1(2:ny-1,2:nx-1),ny-2,nx-2)
call indexx(nsort,t1sort(1:nsort),indx(1:nsort))
do i=1,nsort
   call ind2sub(ny-2,indx(i),i1,i2)
   i1=i1+1
   i2=i2+1
   is_t2near=is_near(t2,i1,i2)
   count_is_t2near=count_me(is_t2near)
   if (count_is_t2near.eq.0) then ! none is know on t2
      t2(i1,i2)=treat3_c(i1,i2)
   elseif (count_is_t2near.eq.1) then
      if (is_t2near(1)) then 
         t2(i1,i2)=treat2_u(i1,i2)
      elseif (is_t2near(2)) then
         t2(i1,i2)=treat2_r(i1,i2)
      elseif (is_t2near(3)) then
         t2(i1,i2)=treat2_d(i1,i2)
      else   
         t2(i1,i2)=treat2_l(i1,i2)
      endif
   else
      is_t2far=is_far(t2,i1,i2)
      count_is_t2far=count_me(is_t2far)
      if (count_is_t2near.eq.2) then
         if (is_t2near(1)) then
            if (is_t2near(2)) then
               if (is_t2far(1)) then
                  t2(i1,i2)=treat1_ur(i1,i2)
               else
                  t2(i1,i2)=min(treat2_u(i1,i2),treat2_l(i1,i2))
               endif
            elseif (is_t2near(3)) then
                  t2(i1,i2)=min(treat2_u(i1,i2),treat2_d(i1,i2))
            elseif (is_t2near(4)) then
               if (is_t2far(4)) then
                  t2(i1,i2)=treat1_ul(i1,i2)
               else
                  t2(i1,i2)=min(treat2_u(i1,i2),treat2_r(i1,i2))
               endif
            endif
         else
            if (is_t2near(2)) then
               if (is_t2near(3)) then
                  if (is_t2far(2)) then
                     t2(i1,i2)=treat1_dr(i1,i2)
                  else
                     t2(i1,i2)=min(treat2_r(i1,i2),treat2_d(i1,i2))
                  endif
               elseif (is_t2near(4)) then
                  t2(i1,i2)=min(treat2_l(i1,i2),treat2_r(i1,i2))
               endif
            else
               if (is_t2far(3)) then
                  t2(i1,i2)=treat1_dl(i1,i2)
               else
                  t2(i1,i2)=min(treat2_d(i1,i2),treat2_r(i1,i2))
               endif
            endif
         endif
      elseif (count_is_t2near.eq.3) then
         if (is_t2near(1).eqv..false.) then ! is_t2near(2,3,4) known
            if (is_t2far(2)) then 
               if (is_t2far(3)) then ! is_t2far(2,3) known
                  t2(i1,i2)=min(treat1_dr(i1,i2),treat1_dl(i1,i2))
               else ! is_t2far(2) known
                  t2(i1,i2)=treat1_dr(i1,i2)
               endif
            else
               if (is_t2far(3)) then ! is_t2far(3) known
                  t2(i1,i2)=treat1_dl(i1,i2)
               else ! is_t2far all unknown
                  t2(i1,i2)=min(treat2_l(i1,i2),treat2_d(i1,i2),treat2_r(i1,i2))
               endif
            endif
         elseif (is_t2near(2).eqv..false.) then ! is_t2near(1,3,4) known
            if (is_t2far(3)) then
               if (is_t2far(4)) then ! is_t2far(3,4) known
                  t2(i1,i2)=min(treat1_dl(i1,i2),treat1_ul(i1,i2))
               else ! is_t2far(3) known
                  t2(i1,i2)=treat1_dl(i1,i2)
               endif
            else
               if (is_t2far(4)) then ! is_t2far(4) known
                  t2(i1,i2)=treat1_ul(i1,i2)
               else ! is_t2far all unknown
                  t2(i1,i2)=min(treat2_d(i1,i2),treat2_r(i1,i2),treat2_u(i1,i2))
               endif
            endif
         elseif (is_t2near(3).eqv..false.) then ! is_t2near(1,2,4) known
            if (is_t2far(4)) then
               if (is_t2far(1)) then ! is_t2far(1,4) known
                  t2(i1,i2)=min(treat1_ul(i1,i2),treat1_ur(i1,i2))
               else ! is_t2far(4) known
                  t2(i1,i2)=treat1_ul(i1,i2)
               endif
            else 
               if (is_t2far(1)) then ! is_t2far(1) known
                  t2(i1,i2)=treat1_ur(i1,i2)
               else ! is_t2far all unknown
                  t2(i1,i2)=min(treat2_l(i1,i2),treat2_u(i1,i2),treat2_r(i1,i2))
               endif
            endif
         else ! is_t2near(1,2,3) known 
            if (is_t2far(1)) then
               if (is_t2far(2)) then ! is_t2far(1,2) known
                  t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2))
               else  ! is_t2far(1) known
                  t2(i1,i2)=treat1_ur(i1,i2)
               endif
            else
               if (is_t2far(2)) then ! is_t2far(2) known
                  t2(i1,i2)=treat1_dr(i1,i2)
               else ! is_t2far all unknown
                  t2(i1,i2)=min(treat2_u(i1,i2),treat2_r(i1,i2),treat2_d(i1,i2))
               endif
            endif
         endif
      elseif (count_is_t2near.eq.4) then
         if (is_t2far(1)) then
            if (is_t2far(2)) then
               if (is_t2far(3)) then
                  if (is_t2far(4)) then ! is_t2far(1,2,3,4) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2),&
                                   treat1_dl(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(1,2,3) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2),treat1_dl(i1,i2))
                  endif
               else
                  if (is_t2far(4)) then ! is_t2far(1,2,4) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(1,2) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dr(i1,i2))
                  endif
               endif
            else
               if (is_t2far(3)) then 
                  if (is_t2far(4)) then ! is_t2far(1,3,4) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dl(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(1,3) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_dl(i1,i2))
                  endif
               else
                  if (is_t2far(4)) then ! is_t2far(1,4) known
                     t2(i1,i2)=min(treat1_ur(i1,i2),treat1_ul(i1,i2))
                  else  ! is_t2far(1) known
                     t2(i1,i2)=treat1_ur(i1,i2)
                  endif
               endif
            endif
         else
            if (is_t2far(2)) then
               if (is_t2far(3)) then
                  if (is_t2far(4)) then ! is_t2far(2,3,4) known
                     t2(i1,i2)=min(treat1_dr(i1,i2),treat1_dl(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(2,3) known
                     t2(i1,i2)=min(treat1_dr(i1,i2),treat1_dl(i1,i2))
                  endif
               else
                  if (is_t2far(4)) then ! is_t2far(2,4) known
                     t2(i1,i2)=min(treat1_dr(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(2) known
                     t2(i1,i2)=treat1_dr(i1,i2)
                  endif
               endif
            else
               if (is_t2far(3)) then
                  if (is_t2far(4)) then ! is_t2far(3,4) known
                     t2(i1,i2)=min(treat1_dl(i1,i2),treat1_ul(i1,i2))
                  else ! is_t2far(3) known
                     t2(i1,i2)=treat1_dl(i1,i2)
                  endif
               else 
                  if (is_t2far(4)) then ! is_t2far(4) known
                     t2(i1,i2)=treat1_ul(i1,i2)
                  else 
                     t2(i1,i2)=min(treat2_u(i1,i2),treat2_r(i1,i2),&
                                   treat2_d(i1,i2),treat2_l(i1,i2))
                  endif
               endif
            endif
         endif
      endif
   endif
enddo
! --------------end dealing with area 1

! ----------------- dealing with area 2
n=nx
es1(1:n)=s1(2,1:n)
es2(1:n)=s1(1,1:n)
es3(1:n)=s2(2,1:n)
es4(1:n)=s2(1,1:n)
et1(1:n)=t1(2,1:n)
et2(1:n)=t1(1,1:n)
et3(1:n)=t2(2,1:n)
et4(1:n)=0.0
call eadge_expand()
t2(1,2:n-1)=et4(2:n-1)
! --------------end dealing with area 2

! ----------------- dealing with area 3
n=ny
es1(1:n)=s1(1:n,nx-1)
es2(1:n)=s1(1:n,nx)
es3(1:n)=s2(1:n,nx-1)
es4(1:n)=s2(1:n,nx)
et1(1:n)=t1(1:n,nx-1)
et2(1:n)=t1(1:n,nx)
et3(1:n)=t2(1:n,nx-1)
et4(1:n)=0.0
call eadge_expand()
t2(2:n-1,nx)=et4(2:n-1)
! --------------end dealing with area 3

! ----------------- dealing with area 4
n=nx
es1(1:n)=s1(ny-1,1:n)
es2(1:n)=s1(ny,1:n)
es3(1:n)=s2(ny-1,1:n)
es4(1:n)=s2(ny,1:n)
et1(1:n)=t1(ny-1,1:n)
et2(1:n)=t1(ny,1:n)
et3(1:n)=t2(ny-1,1:n)
et4(1:n)=0.0
call eadge_expand()
t2(ny,2:n-1)=et4(2:n-1)
! --------------end dealing with area 4

! ----------------- dealing with area 5
n=ny
es1(1:n)=s1(1:n,2)
es2(1:n)=s1(1:n,1)
es3(1:n)=s2(1:n,2)
es4(1:n)=s2(1:n,1)
et1(1:n)=t1(1:n,2)
et2(1:n)=t1(1:n,1)
et3(1:n)=t2(1:n,2)
et4(1:n)=0.0
call eadge_expand()
t2(2:n-1,1)=et4(2:n-1)
! --------------end dealing with area 5

! ----------------- dealing with area 6
s11=s1(2:1:-1,nx-1:nx)
s12=s2(2:1:-1,nx-1:nx)
t11=t1(2:1:-1,nx-1:nx)
t12=t2(2:1:-1,nx-1:nx)
t2(1,nx)=treat1()
! --------------end dealing with area 6

! ----------------- dealing with area 7
s11=s1(ny-1:ny,nx-1:nx)
s12=s2(ny-1:ny,nx-1:nx)
t11=t1(ny-1:ny,nx-1:nx)
t12=t2(ny-1:ny,nx-1:nx)
t2(ny,nx)=treat1()
! --------------end dealing with area 7

! ----------------- dealing with area 8
s11=s1(ny-1:ny,2:1:-1)
s12=s2(ny-1:ny,2:1:-1)
t11=t1(ny-1:ny,2:1:-1)
t12=t2(ny-1:ny,2:1:-1)
t2(ny,1)=treat1()
! --------------end dealing with area 8

! ----------------- dealing with area 9
s11=s1(2:1:-1,2:1:-1)
s12=s2(2:1:-1,2:1:-1)
t11=t1(2:1:-1,2:1:-1)
t12=t2(2:1:-1,2:1:-1)
t2(1,1)=treat1()
! --------------end dealing with area 9
end subroutine surface_expand

! treat 1, most accurate type stencel
!     call ind2sub(n1-2,n2-2,indx(i),i1,i2)
!                      A (1,1) --
!       E--------F     B (1,2)   | t11 in plane
!      /|       /|     C (2,1)   |
!     / |      / |     D (2,2) --
!    /  |     /  |     
!   G--------H   |     E (1,1) --
!   |   A----|---B     F (1,2)   |
!   |  /     |  /      G (2,1)   | t12 out plane
!   | /      | /       H (2,2) --
!   |/       |/
!   C--------D         unknow is H, all other are known
! H  = A + 1/(2^0.5) * sq^0.5
! s_ave = the average of all 8 points
! sq = 6 * h^2 * s_ave^2 -(B-C)^2-(E-C)^2-(B-E)^2-(G-D)^2-(D-F)^2-(G-F)^2
real function treat1()
real :: s_ave, sq1, sq
s_ave=average(average(s11),average(s12))
sq1=6.0*hh*s_ave**2.0
sq=sq1-(t11(1,2)-t11(2,1))**2.0-(t12(1,1)-t11(2,1))**2.0&
      -(t11(1,2)-t12(1,1))**2.0-(t12(2,1)-t11(2,2))**2.0&
      -(t11(2,2)-t12(1,2))**2.0-(t12(2,1)-t12(1,2))**2.0
if (sq<0.0) then
   sq=sq1
endif
treat1=t11(1,1)+c1*sqrt(sq)
!if (sq<0.0) then
!   treat1=min(t11(2,2),t12(1,2),t12(2,1))+h*s_ave
!else
!   treat1=t11(1,1)+c1*sqrt(sq)
!endif
end function treat1

! treat2, second accurate type stencel
!
!             G---------H     A (1,1) ---
!            /|        /|     B (1,2)   |
!           / |       / |     C (2,1)   |_ t21 in plane
!          /  |      /  |     D (2,2)   |
!         I---------J   |     E (3,1)   |
!        /|  A/----/|---/B    F (3,2) ---  
!       / |  /    / |  /
!      /  | /    /  | /       G (1,1) ---
!     K---------L   |/        H (1,2)   |
!     |   C-----|---/D        I (2,1)   |_ t22 out plane
!     |  /      |  /          J (2,2)   |
!     | /       | /           K (3,1)   |
!     |/        |/            L (3,2) ---
!     E---------F         J is unknown    
! J = C + sq^0.5
! sq=2*h^2*s_ave^2-0.5*(A-E)^2-(I-D)^2
! s_ave= the average of all 12 points
real function treat2()
real :: s_ave, sq1, sq
!s_ave=average(average(s21),average(s22))
s_ave=average(average(s21(1:3,1)),s21(2,2),s22(2,1),s22(2,2))
sq1=2.0*hh*s_ave**2.0
sq=sq1-0.5*(t21(1,1)-t21(3,1))**2.0-(t22(2,1)-t21(2,2))**2.0
if (sq<0.0) then
   sq=sq1
endif
treat2=t21(2,1)+sqrt(sq)
!if (sq<0.0) then
!   treat2=min(t21(2,2),t22(2,1))+h*s_ave
!else
!   treat2=t21(2,1)+sqrt(sq)
!endif
end function treat2

! treat3, third type stencel            t1 in plane :
! 
!             J---------K---------L     A (1,1), B (1,2), C (1,3)   
!            /|        /|        /|     D (2,1), E (2,2), F (2,3)
!           / |       / |       / |     G (3,1), H (3,2), I (3,3)
!          /  |      /  |      /  |
!         M---A-----N---B-----O---C     t2 out plane : 
!        /|  /     /|  /     /|  /      J (1,1), K (1,2), L (1,3)
!       / | /     / | /     / | /       M (2,1), N (2,2), O (2,3)
!      /  |/     /  |/     /  |/        P (3,1), Q (3,2), R (3,3)
!     P---D-----Q---E-----R---F 
!     |  /      |  /      |  /          N is unknown
!     | /       | /       | /
!     |/        |/        |/
!     G---------H---------I   
! N = E + sq^0.5  
! sq=h^2*s_ave^2-0.25*[(D-F)^2+(B-H)^2]     
! s_ave= the average of 18 points    
real function treat3()
real :: s_ave, sq1, sq
!s_ave=average(average(s31),average(s32))
s_ave=average(average(s31(1,1),s31(2,1),s31(3,1),s31(1,2),s31(1,3)),s32(2,2))
sq1=hh*s_ave**2.0
sq=sq1-0.25*((t31(2,1)-t31(2,3))**2.0+(t31(1,2)-t31(3,2))**2.0)
if (sq<0.0) then
   sq=sq1
endif
treat3=t31(2,2)+sqrt(sq)
!if (sq<0.0) then
!   treat3=t31(2,2)+h*s_ave
!else
!   treat3=t31(2,2)+sqrt(sq)
!endif
end function treat3

! f4-----n1------f1    near(n1,n2,n3,n4)
!  |      |      |   
!  |      |      |     far(f1,f2,f3,f4)
! n4------o------n2
!  |      |      |
!  |      |      |
! f3-----n3------n4 
function is_near(t2,i1,i2)
logical :: is_near(4)
real,   intent(in) :: t2(:,:)
integer,intent(in) :: i1, i2
is_near=.false.
if (t2(i1-1,i2)>tm) is_near(1)=.true.
if (t2(i1,i2+1)>tm) is_near(2)=.true.
if (t2(i1+1,i2)>tm) is_near(3)=.true.
if (t2(i1,i2-1)>tm) is_near(4)=.true.
end function is_near

function is_far(t2,i1,i2)
logical :: is_far(4)
real,   intent(in) :: t2(:,:)
integer,intent(in) :: i1, i2
is_far=.false.
if (t2(i1-1,i2+1)>tm) is_far(1)=.true.
if (t2(i1+1,i2+1)>tm) is_far(2)=.true.
if (t2(i1+1,i2-1)>tm) is_far(3)=.true.
if (t2(i1-1,i2-1)>tm) is_far(4)=.true.
end function is_far

! f4-----n1------f1    near(n1,n2,n3,n4)
!  |      |      |   
!  |      |      |     far(f1,f2,f3,f4)
! n4------o------n2
!  |      |      |
!  |      |      |
! f3-----n3------f2 

real function treat1_ur(i1,i2)
integer, intent(in) :: i1, i2
s11=s1(i1-1:i1,i2+1:i2:-1) 
s12=s2(i1-1:i1,i2+1:i2:-1) 
t11=t1(i1-1:i1,i2+1:i2:-1) 
t12=t2(i1-1:i1,i2+1:i2:-1) 
treat1_ur=treat1()
end function treat1_ur

real function treat1_dr(i1,i2)
integer, intent(in) :: i1,i2
s11=s1(i1+1:i1:-1,i2+1:i2:-1)
s12=s2(i1+1:i1:-1,i2+1:i2:-1)
t11=t1(i1+1:i1:-1,i2+1:i2:-1)
t12=t2(i1+1:i1:-1,i2+1:i2:-1)
treat1_dr=treat1()
end function treat1_dr

real function treat1_dl(i1,i2)
integer, intent(in) :: i1,i2
s11=s1(i1+1:i1:-1,i2-1:i2)
s12=s2(i1+1:i1:-1,i2-1:i2)
t11=t1(i1+1:i1:-1,i2-1:i2)
t12=t2(i1+1:i1:-1,i2-1:i2)
treat1_dl=treat1()
end function treat1_dl

real function treat1_ul(i1,i2)
integer, intent(in) :: i1, i2
s11=s1(i1-1:i1,i2-1:i2)
s12=s2(i1-1:i1,i2-1:i2)
t11=t1(i1-1:i1,i2-1:i2)
t12=t2(i1-1:i1,i2-1:i2)
treat1_ul=treat1()
end function treat1_ul

real function treat2_u(i1,i2)
integer, intent(in) :: i1, i2
s21=transpose(s1(i1-1:i1,i2-1:i2+1))
s22=transpose(s2(i1-1:i1,i2-1:i2+1))
t21=transpose(t1(i1-1:i1,i2-1:i2+1))
t22=transpose(t2(i1-1:i1,i2-1:i2+1))
treat2_u=treat2()
end function treat2_u

real function treat2_r(i1,i2)
integer, intent(in) :: i1, i2
s21=s1(i1-1:i1+1,i2+1:i2:-1)
s22=s2(i1-1:i1+1,i2+1:i2:-1)
t21=t1(i1-1:i1+1,i2+1:i2:-1)
t22=t2(i1-1:i1+1,i2+1:i2:-1)
treat2_r=treat2()
end function treat2_r

real function treat2_d(i1,i2)
integer, intent(in) :: i1, i2
s21=transpose(s1(i1+1:i1:-1,i2-1:i2+1))
s22=transpose(s2(i1+1:i1:-1,i2-1:i2+1))
t21=transpose(t1(i1+1:i1:-1,i2-1:i2+1))
t22=transpose(t2(i1+1:i1:-1,i2-1:i2+1))
treat2_d=treat2()
end function treat2_d

real function treat2_l(i1,i2)
integer, intent(in) :: i1, i2
s21=s1(i1-1:i1+1,i2-1:i2)
s22=s2(i1-1:i1+1,i2-1:i2)
t21=t1(i1-1:i1+1,i2-1:i2)
t22=t2(i1-1:i1+1,i2-1:i2)
treat2_l=treat2()
end function treat2_l

real function treat3_c(i1,i2)
integer, intent(in) :: i1,i2
s31=s1(i1-1:i1+1,i2-1:i2+1)
s32=s2(i1-1:i1+1,i2-1:i2+1)
t31=t1(i1-1:i1+1,i2-1:i2+1)
treat3_c=treat3()
end function treat3_c

! eadge
!
!           K-1---------L-1      (I-1,I,I+1) -- et1,es1 
!            /|        /|        (J-1,J,J+1) -- et2,es2
!           / |       / |        (K-1,K,K+1) -- et3,es3
!          /  |      /  |        (L-1,L,L+1) -- et4,es4
!         K---------L   |     
!        /| I-1/---/|---/J-1     I(et1(i) is the sorting point  
!       / |  /    / |  /         L(et4(i) is unknown
!      /  | /    /  | /       
!   K+1---------L+1 |/        
!     |   I-----|---/J        
!     |  /      |  /          
!     | /       | /           
!     |/        |/            
!   I+1---------J+1             
subroutine eadge_expand()
integer :: nsort, ii, i
real :: et4i_tmp
nsort=n-2
t1sort(1:nsort)=et1(2:n-1)
call indexx(nsort,t1sort(1:nsort),indx(1:nsort))
do ii=1,nsort
   i=indx(ii)+1
   if (et4(i-1)>tm) then
      s11=compose_funs(es1(i-1),es2(i-1),es1(i),es2(i))
      s12=compose_funs(es3(i-1),es4(i-1),es3(i),es4(i))
      t11=compose_funs(et1(i-1),et2(i-1),et1(i),et2(i))
      t12=compose_funs(et3(i-1),et4(i-1),et3(i),et4(i))
      et4i_tmp=treat1()
      if (et4(i+1)>tm) then
         s11=compose_funs(es1(i+1),es2(i+1),es1(i),es2(i))
         s12=compose_funs(es3(i+1),es4(i+1),es3(i),es4(i))
         t11=compose_funs(et1(i+1),et2(i+1),et1(i),et2(i))
         t12=compose_funs(et3(i+1),et4(i+1),et3(i),et4(i))
         et4(i)=min(et4i_tmp,treat1())     
      else
         et4(i)=et4i_tmp
      endif
   else
      if (et4(i+1)>tm) then
         s11=compose_funs(es1(i+1),es2(i+1),es1(i),es2(i))
         s12=compose_funs(es3(i+1),es4(i+1),es3(i),es4(i))
         t11=compose_funs(et1(i+1),et2(i+1),et1(i),et2(i))
         t12=compose_funs(et3(i+1),et4(i+1),et3(i),et4(i))
         et4(i)=treat1()     
      else
         s21=compose_funs(es1(i-1:i+1),es2(i-1:i+1),3)
         s22=compose_funs(es3(i-1:i+1),es4(i-1:i+1),3)
         t21=compose_funs(et1(i-1:i+1),et2(i-1:i+1),3)
         t22=compose_funs(et3(i-1:i+1),et4(i-1:i+1),3)
         et4(i)=treat2()
      endif
   endif
enddo
end subroutine eadge_expand

integer function count_me(d_in)
logical, intent(in) :: d_in(4)
integer :: i
count_me=0
do i=1,4
   if (d_in(i)) then
      count_me=count_me+1
   endif
enddo
end function count_me

end module module_pweikfd3d_dir
