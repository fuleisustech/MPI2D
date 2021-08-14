module module_eikfd2d_dir
use module_array
use module_utility
use module_math

implicit none
integer, private    :: nz,nx
real,    private    :: sx,sz,h
integer, private    :: isx,isz,n_array,n_array_max,s_index,min_number,max_number
integer, allocatable, private :: min_loc(:),max_loc(:)
real,    allocatable, private :: s_inner(:),s_outer(:),t_inner(:),t_outer(:)

private :: aroundSource, squareExpand, trat
contains

!======================================================================
! the main subroutine for 2D FD solver of direct arrival of eikonal eqn
!======================================================================
subroutine eikfd2d_dir(s,t,sz0,sx0,nz0,nx0,h0)

real,    intent(in)  :: s(:,:),sx0,sz0,h0
integer, intent(in)  :: nz0,nx0
real,    intent(out) :: t(:,:)
!call date_and_time(date,time_now)
!write(*,*)time_now," eik start"
t=0.0
call copy_variable(sz0,sz,sx0,sx,h0,h)
call copy_variable(nz0,nz,nx0,nx)
isx=nint(sx/h)+1
isz=nint(sz/h)+1
n_array_max=max(nz,nx)
call allocate_and_initial(s_inner,s_outer,t_inner,t_outer,n_array_max)
call allocate_and_initial(min_loc,max_loc,n_array_max)
!call date_and_time(date,time_now)
!write(*,*)time_now," initial variables done"
call aroundSource(s,t)
!call date_and_time(date,time_now)
!write(*,*)time_now," around source done"
call squareExpand(s,t)
!call date_and_time(date,time_now)
!write(*,*)time_now," square done"
deallocate(s_inner,s_outer,t_inner,t_outer)
deallocate(min_loc,max_loc)
!call date_and_time(date,time_now)
!write(*,*)time_now," eik finished"
end subroutine eikfd2d_dir

!----------------------------------------------------------------------
! This subroutine is used to calculate the traveltime around the source
!----------------------------------------------------------------------
subroutine aroundSource(s,t)
real,  intent(in)     :: s(:,:)
real,  intent(inout)  :: t(:,:)
integer               :: ix1,ix2,iz1,iz2
real                  :: s_ave
integer               :: i1,i2
ix1=max(isx-2,1)
ix2=min(isx+2,nx)
iz1=max(isz-2,1)
iz2=min(isz+2,nz)

!write(*,*) "ix1===========",ix1
!write(*,*) "ix2===========",ix2
!write(*,*) "iz1===========",iz1
!write(*,*) "iz2===========",iz2
do i2=ix1,ix2
   do i1=iz1,iz2
      s_ave=sum(s(min(i1,isz):max(i1,isz),min(i2,isx):max(i2,isx))) &
           /real((max(i1,isz)-min(i1,isz)+1)*(max(i2,isx)-min(i2,isx)+1))
      !s1=0.0
      !ss1=0.0
      !do j=min(i2,isx),max(i2,isx)
      !   do k=min(i1,isz),max(i1,isz)
      !      s1=s1+s(k,j)
      !      ss1=ss1+1.0
      !   enddo
      !enddo
      !s_ave=s1/ss1
      t(i1,i2)=sqrt(((i2-1)*h-sx)**2+((i1-1)*h-sz)**2) * s_ave
   enddo
enddo
end subroutine aroundSource

!-----------------------------------------------------------------
! This subroutine is used to expand the traveltime table by square 
!-----------------------------------------------------------------
subroutine squareExpand(s,t)
real,   intent(in)    :: s(:,:)
real,   intent(inout) :: t(:,:)
integer :: i_outer,i_inner
integer :: j1,j2
integer :: i

do i=3,max(isx-1,isz-1,nx-isx,nz-isz)
!  The order will go through "Top" to "Right" to "Bottom" to "Left"        
!
!                        L---T---T---T---T---T---R
!                        |                       |
!                        L   *---*---*---*---*   R
!                        |   |               |   |
!                        L   *               *   R
!                        |   |               |   |
!                        L   *               *   R
!                        |   |               |   |
!                        L   *               *   R
!                        |   |               |   |
!                        L   *---*---*---*---*   R
!                        |                       |
!                        L---B---B---B---B---B---B

!  deal with the top edge
   if (isz-i>=1) then
      i_outer=isz-i
      i_inner=i_outer+1
      j1=max(1,isx-i+1)
      j2=min(isx+i-1,nx)
      n_array=j2-j1+1
      s_index=isx-j1+1
      s_inner(1:n_array)=s(i_inner,j1:j2)
      t_inner(1:n_array)=t(i_inner,j1:j2)
      s_outer(1:n_array)=s(i_outer,j1:j2)
      !t_outer(1:n_array)=0.0
      call trat()
      t(i_outer,j1:j2)=t_outer(1:n_array)
   endif    
! the end of the dealing with the top eadge

! deal with the right edge
   if (isx+i<=nx) then
      i_outer=isx+i
      i_inner=i_outer-1
      j1=max(1,isz-i)
      j2=min(isz+i-1,nz)
      n_array=j2-j1+1
      s_index=isz-j1+1
      s_inner(1:n_array)=s(j1:j2,i_inner)
      t_inner(1:n_array)=t(j1:j2,i_inner)
      s_outer(1:n_array)=s(j1:j2,i_outer)
      !t_outer(1:n_array)=0.0
      call trat()
      t(j1:j2,i_outer)=t_outer(1:n_array)
   endif    
! the end of the dealing with the right eadge

! deal with the bottom edge
   if (isz+i<=nz) then
      i_outer=isz+i
      i_inner=i_outer-1
      j1=max(1,isx-i+1)
      j2=min(isx+i,nx)
      n_array=j2-j1+1
      s_index=isx-j1+1
      s_inner(1:n_array)=s(i_inner,j1:j2)
      t_inner(1:n_array)=t(i_inner,j1:j2)
      s_outer(1:n_array)=s(i_outer,j1:j2)
      !t_outer(1:n_array)=0.0
!      write(*,*)"i_outer ",i_outer," j1 ",j1," j2 ",j2
!      call flush(6)
      call trat()
      t(i_outer,j1:j2)=t_outer(1:n_array)
   endif
! the end of the dealing with the bottom eadge

! deal with the left edge
   if (isx-i>=1) then
      i_outer=isx-i
      i_inner=i_outer+1
      j1=max(1,isz-i)
      j2=min(isz+i,nz)
      n_array=j2-j1+1
      s_index=isz-j1+1
      s_inner(1:n_array)=s(j1:j2,i_inner)
      t_inner(1:n_array)=t(j1:j2,i_inner)
      s_outer(1:n_array)=s(j1:j2,i_outer)
      !t_outer(1:n_array)=0.0
      call trat()      
      t(j1:j2,i_outer)=t_outer(1:n_array)
   endif
! the end of the dealing with the left eadge
enddo
end subroutine squareExpand

subroutine trat()
integer               :: i,j,i_min,i_max,i_close_min,diff 
real                  :: s_ave,sq
diff=nx+nz
call local_min(t_inner,n_array,min_number,min_loc)
call local_max(t_inner,n_array,max_number,max_loc)
! calculate the traveltime just outside the local minimum points
!call flush(6)
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
! find the closest local minimum point
i_close_min=0
do i=1,min_number    
   if (abs(min_loc(i)-s_index) < diff ) then
      diff=abs(min_loc(i)-s_index)
      i_close_min=i
   endif
enddo
! this i_close_min point will divide the eadge into 2 segments      

! deal with the left segment 
! deal with the incoming wave (sweep from left to right)
! no minimum left of the close point, no incoming wave
if (i_close_min==1) then 
   continue
else
   do i=1,i_close_min-1
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
endif

! deal with the outgoing wave (sweep from right to left)
! the close point is on left end, no outgoing
if (min_loc(i_close_min)==1) then 
   continue
else
   do i=i_close_min,1,-1
      i_min=min_loc(i)
      i_max=nx+nz
      loop2: do j=max_number,1,-1
         if ( max_loc(j) < i_min ) then
            i_max=max_loc(j)
            exit loop2
         endif
      enddo loop2
      if (i_max /= nx+nz) then
         do j=i_min-1,i_max,-1
            !s_ave=0.25*(s_inner(j)+s_inner(j+1)+s_outer(j)+s_outer(j+1))
            s_ave=average(s_inner(j),s_inner(j+1),s_outer(j),s_outer(j+1))
            sq=2*(h*s_ave)**2-(t_inner(j)-t_outer(j+1))**2
            if (sq < 0.0 ) then
               t_outer(j)=min(t_inner(j),t_outer(j+1))+h*s_ave 
            else
               t_outer(j)=t_inner(j+1)+sqrt(sq)
            endif
         enddo
      endif   
   enddo 
endif
! deal with the right segment
! deal with the incoming wave, sweep from right to left
! no local min right to the close point,no incoming wave
if (i_close_min==min_number) then 
   continue
else 
   do i=min_number,i_close_min+1,-1
      i_min=min_loc(i)
      i_max=nx+nz
      loop3: do j=max_number,1,-1
         if (max_loc(j) < i_min) then
            i_max=max_loc(j)
            exit loop3
         endif
      enddo loop3
      do j=i_min-1,i_max,-1
         !s_ave=0.25*(s_inner(j)+s_inner(j+1)+s_outer(j)+s_outer(j+1))
         s_ave=average(s_inner(j),s_inner(j+1),s_outer(j),s_outer(j+1))
         sq=2*(h*s_ave)**2-(t_inner(j)-t_outer(j+1))**2
         if (sq < 0.0 ) then
            t_outer(j)=min(t_inner(j),t_outer(j+1))+h*s_ave
         else
            t_outer(j)=t_inner(j+1)+sqrt(sq)
         endif
      enddo
   enddo  
endif

! deal with the outgoing wave, sweep from left to right
! the close min on the right end of array, no outgoing wave
if (min_loc(i_close_min)==n_array) then 
   continue
else 
   do i=i_close_min,min_number
      i_min=min_loc(i)
      i_max=0
      loop4: do j=1,max_number
         if (max_loc(j) > i_min ) then
            i_max=max_loc(j)
            exit loop4
         endif
      enddo loop4  
      if (i_max /= 0) then
         do j=i_min+1,i_max
            !s_ave=0.25*(s_inner(j)+s_inner(j-1)+s_outer(j)+s_outer(j-1))
            s_ave=average(s_inner(j),s_inner(j-1),s_outer(j),s_outer(j-1))
            sq=2*(h*s_ave)**2-(t_inner(j)-t_outer(j-1))**2
            if (sq < 0.0) then
               t_outer(j)=min(t_inner(j),t_outer(j-1))+h*s_ave
            else
               t_outer(j)=t_inner(j-1)+sqrt(sq)
            endif
         enddo
      endif
   enddo 
endif   

end subroutine trat

end module module_eikfd2d_dir
