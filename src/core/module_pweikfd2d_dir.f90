module module_pweikfd2d_dir
use module_global , only : PI
use module_array
use module_utility
use module_math

implicit none
real, allocatable :: t_inner(:),s_inner(:),t_outer(:),s_outer(:)
integer,private   :: nx, nz
real,   private   :: h
private :: trat 

contains

!----------------------------------------------------
subroutine pweikfd2d_dir(s,t,theta_x,xref,zp,nz_in,nx_in,h_in)
real,    intent(in)  :: s(:,:),h_in,theta_x,xref,zp
integer, intent(in)  :: nz_in, nx_in
real,    intent(out) :: t(:,:)
integer  :: ix, iz, i_outer,i_inner,izp
real     :: sin_theta
call copy_variable(nz_in,nz,nx_in,nx)
call copy_variable(h_in,h)
t=0.0
! decide the initial time at the plane
izp=nint(zp/h)+1
sin_theta=sin(theta_x/180.0*PI)
do ix=1,nx
   t(izp,ix)=((ix-1)*h-xref)*sin_theta*s(izp,ix)
enddo
call allocate_and_initial(t_inner,t_outer,s_inner,s_outer,nx)
do iz=1,max(nz-izp,izp-1)
! deal with the top eadge
   if (izp-iz>=1) then
      i_outer=izp-iz
      i_inner=i_outer+1
      t_inner=t(i_inner,:)
      s_inner=s(i_inner,:)
      s_outer=s(i_outer,:)
      call trat ()
      t(i_outer,:)=t_outer
   endif

! deal with the bottom edge
   if (izp+iz<=nz) then
      i_outer=izp+iz
      i_inner=i_outer-1
      t_inner=t(i_inner,:)
      s_inner=s(i_inner,:)
      s_outer=s(i_outer,:)
      call trat ()
      t(i_outer,:)=t_outer
   endif
enddo
end subroutine pweikfd2d_dir

!----------------------------------------------------
subroutine trat()
integer :: min_number,max_number,min_loc(nx),max_loc(nx)
real    :: s_ave,sq,t_imax_temp
integer :: i,j,n_array,i_min,i_max 

n_array=nx
call local_min(t_inner,n_array,min_number,min_loc)
call local_max(t_inner,n_array,max_number,max_loc)

! calculate the traveltime just outside the local minimum points
do i=1,min_number
   i_min=min_loc(i)
   if (i_min==1 .or. i_min==n_array ) then
      !t_outer(i_min)=t_inner(i_min)+h*s_outer(i_min)
      t_outer(i_min)=t_inner(i_min)+h*average(s_outer(i_min),s_inner(i_min))
   else 
      !sq=(h*s_outer(i_min))**2-0.25*(t_inner(i_min-1)-t_inner(i_min+1))**2 
      s_ave=average(average(s_inner(i_min-1:i_min+1)),s_outer(i_min))
      sq=(h*s_ave)**2.0-0.25*(t_inner(i_min-1)-t_inner(i_min+1))**2 
      if (sq < 0.0) then
         !t_outer(i_min)=t_inner(i_min)+h*s_outer(i_min)
         t_outer(i_min)=t_inner(i_min)+h*s_ave
      else
         t_outer(i_min)=t_inner(i_min)+sqrt(sq)
      endif
   endif
enddo   

!sweep from left to right
do i=1,min_number
   i_min=min_loc(i)
   i_max=0
   loop1: do j=1,max_number
        if (max_loc(j) > i_min ) then
            i_max=max_loc(j)
            exit loop1
         endif
   enddo loop1
   do j=i_min+1,i_max
      !s_ave=0.25*(s_inner(j)+s_inner(j-1)+s_outer(j)+s_outer(j-1))
      s_ave=average(s_inner(j),s_inner(j-1),s_outer(j),s_outer(j-1))
      sq=2*(h*s_ave)**2-(t_inner(j)-t_outer(j-1))**2
      if (sq<0) then
         t_outer(j)=min(t_inner(j),t_outer(j-1))+h*s_ave
      else
         t_outer(j)=t_inner(j-1)+sqrt(sq)
      endif
   enddo 
enddo  

! sweep from right to left
do i=min_number,1,-1
   i_min=min_loc(i)
   i_max=2*nx
   loop2: do j=max_number,1,-1
      if ( max_loc(j) < i_min ) then
         i_max=max_loc(j)
         exit loop2
      endif
   enddo loop2
   if (i_max /= 2*nx) then
      do j=i_min-1,i_max,-1
         !s_ave=0.25*(s_inner(j)+s_inner(j+1)+s_outer(j)+s_outer(j+1))
         s_ave=average(s_inner(j),s_inner(j+1),s_outer(j),s_outer(j+1))
         sq=2*(h*s_ave)**2-(t_inner(j)-t_outer(j+1))**2
         if (j.eq.i_max) then
            if (sq < 0.0 ) then
               t_imax_temp=min(t_inner(j),t_outer(j+1))+h*s_ave 
            else
               t_imax_temp=t_inner(j+1)+sqrt(sq)
            endif
            if (t_outer(j) > 1e-6) then ! in case it is on the left end
               t_outer(j)=min(t_outer(j),t_imax_temp)
            else
               t_outer(j)=t_imax_temp
            endif
         else
            if (sq<0.0) then
               t_outer(j)=min(t_inner(j),t_outer(j+1))+h*s_ave 
            else
               t_outer(j)=t_inner(j+1)+sqrt(sq)
            endif
         endif
      enddo
   endif   
enddo 

end subroutine trat
!------------------

end module module_pweikfd2d_dir
