module module_pweikfd2d_fir
use module_array; use module_utility
use module_global, only : PI
implicit none
integer,private :: nz,nx,izp
real,   private :: dz,dx,h,hh,dt,ddx,theta_x,xref,zp,sin_theta,t00
real,allocatable,private :: s(:,:),t(:,:)
private :: around_source_fir,eband,trat1,trat2,trat3,trat4,trat5,&
   trat6,trat7,trat8,trat9,trat10,trat11,trat12,trat13,trat14,&
   trat_plane,coners,arr1

contains

subroutine pweikfd2d_fir(s0,t0,theta_x0,xref0,zp0,nz0,nx0,dz0,dx0)
real,intent(in)     :: s0(:,:),theta_x0,xref0,zp0,dz0,dx0
integer,intent(in)  :: nz0,nx0
real,intent(out)    :: t0(:,:)
integer,allocatable :: mx(:),mz(:),mxwk(:),mzwk(:) 
integer             :: ntwk
real                :: tmin,ssmin
integer             :: i,j
integer             :: nt,k1,n1

call copy_variable(nz0,nz,nx0,nx);   call copy_variable(dz0,dz,dx0,dx)
call copy_variable(theta_x0,theta_x,xref0,xref,zp0,zp)
!pad the velocity 2 points
nz=nz+2; nx=nx+2; xref=xref+dx; zp=zp+dz
call allocate_and_initial(mx,mz,mxwk,mzwk,20*(nx+nz))
call allocate_and_initial(t,s,nz,nx)
call padarray(s0,nz-2,1,1,nx-2,1,1,s,"replicate")
!s=s0

h=dx;   hh=h*h;   ssmin=minval(s)
!vmax=1./ssmin

!Try to let all the arrival times of the cauculated points be nonzero.
t=0.0

! decide the initial time at the plane
izp=nint(zp/h)+1;   sin_theta=sin(theta_x/180.0*PI)
write(*,*)"theta_x0,xref,izp"
write(*,*)theta_x0,xref,izp
do j=1,nx
   t(izp,j)=((j-1)*h-xref)*sin_theta*s(izp,j)
enddo
do j=1,nx
   write(60,*)t(izp,j)
enddo
t00=0.001
t=t+t00

ddx=1.0;   dt=ddx*h*ssmin

!mx(1)=ix0
!mz(1)=iz0
!tmin=t(iz0,ix0)
!ntwk=0
      
call around_source_fir(mxwk,mzwk,ntwk)
!do j=1,nx
!   write(61,*)t(izp-1,j)
!   write(62,*)t(izp+1,j)
!   write(63,*)t(izp+2,j)
!enddo

nt=ntwk;   mx(1:nt)=mxwk(1:nt);   mz(1:nt)=mzwk(1:nt)

101 ntwk=0

mx(1:nt)=mxwk(1:nt);   mz(1:nt)=mzwk(1:nt)      
      
call arr1(nt,t,mx,mz)      
tmin=t(mz(1),mx(1))
!do i=1,nt
!   write(79,*)mz(i)      
!   write(80,*)mx(i) 
!   write(81,*)t(mz(i),mx(i))
!enddo
do i=1,nt
   n1=mx(i);   k1=mz(i)
   !c*To check whether the traveltime of this point is too big. 
   if (t(k1,n1).gt.tmin+dt) then
      do j=i,nt
         n1=mx(j); k1=mz(j); ntwk=ntwk+1; mxwk(ntwk)=n1; mzwk(ntwk)=k1
      enddo
      goto 11
   endif      
   call eband(i,mx,mz,mxwk,mzwk,ntwk)
enddo

11 continue 
      
nt=ntwk      
if (nt.ge.1) goto 101

!c****
!c* Final treatment of the four coners.
!c****
100 call coners() 

t0=t(2:nz-1,2:nx-1)
!t0=t-t00
call deallocate_and_free(t,s)
call deallocate_and_free(mx,mz,mxwk,mzwk)
 
return

end subroutine pweikfd2d_fir

subroutine around_source_fir(mxwk,mzwk,ntwk)
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
integer               :: up, down, up_to_down,i,j!,nt
real :: t_inner(nx-2),t_outer(nx-2),s_inner(nx-2),s_outer(nx-2)
!real,   allocatable :: t_sort(:)
!integer,allocatable :: indx(:)
!c***********************************************************
!c This subroutine deals with the source
!c***********************************************************
up=izp;   down=izp
do i=1,2
   ! deal with top
   if (izp-i>=1) then
      up=up-1;   t_inner=t(izp-i+1,2:nx-1)
      s_inner=s(izp-i+1,2:nx-1);   s_outer=s(izp-i,2:nx-1)
      call trat_plane(s_inner,s_outer,t_inner,t_outer)
      t(izp-i,2:nx-1)=t_outer
   endif
   ! deal with bottom
   if (izp+i<=nz) then
      down=down+1;   t_inner=t(izp+i-1,2:nx-1)
      s_inner=s(izp+i-1,2:nx-1);   s_outer=s(izp+i,2:nx-1)
      call trat_plane(s_inner,s_outer,t_inner,t_outer)
      t(izp+i,2:nx-1)=t_outer
   endif
enddo
up_to_down=down-up+1!;  nt=nx*up_to_down
!call allocata_and_initial(t_sort,nt)
!call allocata_and_initial(indx,nt)
!t_sort=arr2vec(t(up:down,:),down-up+1,nx)
!call indexx(nt,t_sort,indx)
!do i=1,nt
!   call ind2sub(up_to_down,indx(i),i1,i2)
!   mzwk(i)=i1; mxwk(i)=i2
!enddo
!ntwk=nt
!call deallocate_and_free(t_sort); call deallocate_and_free(indx)
ntwk=0
do j=2,nx-1;   do i=max(up,2),min(down,nz-1)
   ntwk=ntwk+1;   mxwk(ntwk)=j;   mzwk(ntwk)=i
enddo;   enddo
end subroutine around_source_fir


!c                        *************
!c                        * Sub eband *
!c                        *************

subroutine eband(i,mx,mz,mxwk,mzwk,ntwk)

!c**********************************************************************
!c* This part wants to expand the wavefront boundary from the minimum  *
!c* arriving time point.                                               * 
!c**********************************************************************
integer,intent(in)  :: i
integer,intent(in)    :: mx(:),mz(:)
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
integer             :: n1,k1

n1=mx(i);   k1=mz(i)
!write(*,*)"n1",n1,"k1",k1
!c#1a
if (t(k1+1,n1).gt.t00) then
!c#2a
   if (t(k1-1,n1).gt.t00) then
!c#3a
      if (t(k1,n1+1).gt.t00) then
!c#4a      
         if (t(k1,n1-1).gt.t00) then
            return
!c#4b
         else
!c****
!c* find the value of t(k1,n1-1) and other treatment.
!c****
            call trat1(n1,k1,mxwk,mzwk,ntwk)
!c**
!c#4c
         endif
!c#3b
      else
!c#4a
         if (t(k1,n1-1).gt.t00) then
!c****
!c* To find the value of t(k1,n1+1) and etc.
!c****
            call trat2(n1,k1,mxwk,mzwk,ntwk)
!c#4b
         else
!c****
!c* find the value of t(k1,n1+1) & t(k1,n1-1) and etc.
!c****
            call trat3(n1,k1,mxwk,mzwk,ntwk)
!c#4c
         end if
!c#3c
      end if
!c#2b
   else
!c#3a
      if (t(k1,n1+1).gt.t00) then
!c#4a
         if (t(k1,n1-1).gt.t00) then
!c****
!c* find t(k1-1,n1) & etc.
!c****
            call trat4(n1,k1,mxwk,mzwk,ntwk)
!c#4b
         else
!c****
!c* find t(k1,n1-1),t(k1-1,n1) & etc.
!c****
            call trat5(n1,k1,mxwk,mzwk,ntwk)
!c#4c      
         endif
!c#3b      
      else
!c#4a      
         if (t(k1,n1-1).gt.t00) then
!c****
!c* find t(k1-1,n1),t(k1,n1+1) & etc.
!c****
            call trat6(n1,k1,mxwk,mzwk,ntwk)
!c#4b
         else
!c****
!c* find t(k1-1,n1),t(k1,n1+1),t(k1,n1-1) & etc.
!c****
            call trat7(n1,k1,mxwk,mzwk,ntwk)
!c#4c
         end if
!c#3c      
      endif
!c#2c      
   endif
!c#1b      
else
!c#2a      
   if (t(k1-1,n1).gt.t00) then
!c#3a      
      if (t(k1,n1+1).gt.t00) then
!c#4a      
         if (t(k1,n1-1).gt.t00) then
!c****
!c* find t(k1+1,n1) & etc.
!c****
            call trat8(n1,k1,mxwk,mzwk,ntwk)
!c#4b
         else
!c****
!c* find t(k1,n1-1) & t(k1+1,n1) & etc.
!c****
            call trat9(n1,k1,mxwk,mzwk,ntwk)
!c#4c
         endif
!c#3b
      else
         if (t(k1,n1-1).gt.t00) then
!c****
!c* find t(k1,n1+1) t(k1+1,n1) & etc.
!c****
            call trat10(n1,k1,mxwk,mzwk,ntwk)
!c#4b      
         else
!c****
!c* find t(k1,n1+1),t(k1,n1-1) t(k1+1,n1) & etc.
!c****
            call trat11(n1,k1,mxwk,mzwk,ntwk)
!c#4c
         endif
!c#3c
      endif
!c#2b
   else
!c#3a      
      if (t(k1,n1+1).gt.t00) then
!c#4a      
         if (t(k1,n1-1).gt.t00) then
!c****
!c* find t(k1-1,n1) t(k1+1,n1) & etc.
!c****
            call trat12(n1,k1,mxwk,mzwk,ntwk)
!c#4b      
         else
!c****
!c* find t(k1,n1-1) t(k1-1,n1),t(k1+1,n1) & etc.
!c****
            call trat13(n1,k1,mxwk,mzwk,ntwk)
!c#4c
         endif
!c#3b      
      else
!c#4a
         if (t(k1,n1-1).gt.t00) then
!c****
!c* find t(k1-1,n1) t(k1,n1+1),t(k1+1,n1) & etc.
!c****
            call trat14(n1,k1,mxwk,mzwk,ntwk)
!c#4b
         else
!c****
!c* Source point treatment
!c****
!c Impossible
!c
!c#4c      
         end if
!c#3c    
      end if
!c#2c      
   end if
!c#1c      
end if
      
return

end subroutine eband

!c                        *************
!c                        * Sub trat1 *
!c                        *************
       
subroutine trat1(n1,k1,mxwk,mzwk,ntwk)
!c***********************************************************
!c* This sub is to find the value of t(k1,n1-1)             *
!c***********************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
if (t(k1-1,n1-1).ne.0.0) then
   s1=(s(k1-1,n1)+s(k1-1,n1-1)+s(k1,n1)+s(k1,n1-1))/4
   sq=hh*s1*s1*2-(t(k1-1,n1-1)-t(k1,n1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1,n1-1)=t(k1-1,n1)+sqrt(sq)
   t(k1,n1-1)=max(t(k1,n1-1),t(k1,n1))
   goto 1
endif

if (t(k1+1,n1-1).ne.0.0) then
   s1=(s(k1+1,n1)+s(k1,n1-1)+s(k1+1,n1-1)+s(k1,n1))/4
   sq=hh*s1*s1*2-(t(k1+1,n1-1)-t(k1,n1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1,n1-1)=t(k1+1,n1)+sqrt(sq)
   t(k1,n1-1)=max(t(k1,n1-1),t(k1,n1))
   goto 1
endif
    
s1=(s(k1,n1-1)+(s(k1,n1)+s(k1+1,n1)+s(k1-1,n1))/3)/2
sq=hh*s1*s1-.25*(t(k1+1,n1)-t(k1-1,n1))**2

if (sq.lt.0.0) sq=s1*s1*h*h
t(k1,n1-1)=t(k1,n1)+sqrt(sq)
!c
!c* The (k1,n1-1) point is on the boundary.
!c

1 continue

if (n1-1.eq.1) goto 8
ntwk=ntwk+1      
mxwk(ntwk)=n1-1      
mzwk(ntwk)=k1

8 continue
      
return
      
end subroutine trat1

!c                        *************
!c                        * Sub trat2 *
!c                        *************
       
subroutine trat2(n1,k1,mxwk,mzwk,ntwk)
!c************************************************************
!c* This sub is to find the value of t(k1,n1+1)              *
!c************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
if (t(k1-1,n1+1).ne.0.0) then
   s1=(s(k1-1,n1)+s(k1,n1+1)+s(k1-1,n1+1)+s(k1,n1))/4
   sq=hh*s1*s1*2-(t(k1-1,n1+1)-t(k1,n1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1,n1+1)=t(k1-1,n1)+sqrt(sq)
   t(k1,n1+1)=max(t(k1,n1+1),t(k1,n1))
   goto 1
endif
   
if (t(k1+1,n1+1).ne.0.0) then
   s1=(s(k1+1,n1)+s(k1,n1+1)+s(k1,n1)+s(k1+1,n1+1))/4
   sq=hh*s1*s1*2-(t(k1,n1)-t(k1+1,n1+1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1,n1+1)=t(k1+1,n1)+sqrt(sq)
   t(k1,n1+1)=max(t(k1,n1+1),t(k1,n1))
   goto 1
endif
      
s1=(s(k1,n1+1)+(s(k1,n1)+s(k1+1,n1)+s(k1-1,n1))/3)/2      
sq=h*s1*h*s1-.25*(t(k1+1,n1)-t(k1-1,n1))**2
      
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1,n1+1)=t(k1,n1)+sqrt(sq)

!c* The (k1,n1+1) is on boundary?
1 continue
 
if (n1+1.eq.nx) goto 10
ntwk=ntwk+1 
mxwk(ntwk)=n1+1 
mzwk(ntwk)=k1 

10  continue

return

end subroutine trat2

!c                        *************
!c                        * Sub trat3 *
!c                        *************
       
subroutine trat3(n1,k1,mxwk,mzwk,ntwk)
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
!write(*,*)'sub=3'

!c**************************************************************
!c* This part is to find the values of t(k1,n1+1) & t(k1,n1-1) *
!c**************************************************************
       
s1=(s(k1,n1-1)+(s(k1,n1)+s(k1+1,n1)+s(k1-1,n1))/3)*.5
sq=h*s1*s1*h-.25*(t(k1+1,n1)-t(k1-1,n1))**2
      
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1,n1-1)=t(k1,n1)+sqrt(sq)

!c* The (k1,n1-1) point is on the boundary?
if (n1-1.eq.1) goto 18
ntwk=ntwk+1  
mxwk(ntwk)=n1-1  
mzwk(ntwk)=k1  
18 continue
      
s1=(s(k1,n1+1)+(s(k1,n1)+s(k1+1,n1)+s(k1-1,n1))/3)*.5
sq=hh*s1*s1-.25*(t(k1+1,n1)-t(k1-1,n1))**2
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1,n1+1)=t(k1,n1)+sqrt(sq)

if (n1+1.eq.nx) goto 20
ntwk=ntwk+1   
mxwk(ntwk)=n1+1   
mzwk(ntwk)=k1   
20 continue

return
      
end subroutine trat3

!c                        *************
!c                        * Sub trat4 *
!c                        *************
subroutine trat4(n1,k1,mxwk,mzwk,ntwk)
!c******************************************************************
!c* This sub is to find the values of t(k1-1,n1)                   *
!c******************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
if (t(k1-1,n1-1).ne.0.0) then
   s1=(s(k1-1,n1)+s(k1,n1-1)+s(k1,n1)+s(k1-1,n1-1))*.25
   sq=hh*s1*s1*2-(t(k1,n1)-t(k1-1,n1-1))**2
   if (sq.lt.0.0)sq=2*s1*s1*h*h
   t(k1-1,n1)=t(k1,n1-1)+sqrt(sq)
   t(k1-1,n1)=max(t(k1-1,n1),t(k1,n1))
   goto 1
endif
      
if (t(k1-1,n1+1).ne.0.0) then
   s1=(s(k1-1,n1)+s(k1,n1+1)+s(k1,n1)+s(k1-1,n1+1))*.25
   sq=hh*s1*s1*2.-(t(k1,n1)-t(k1-1,n1+1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1-1,n1)=t(k1,n1+1)+sqrt(sq)
   t(k1-1,n1)=max(t(k1-1,n1),t(k1,n1))
   goto 1
endif
      
s1=(s(k1-1,n1)+(s(k1,n1)+s(k1,n1-1)+s(k1,n1+1))/3)*.5      
sq=hh*s1*s1-.25*(t(k1,n1-1)-t(k1,n1+1))**2
      
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1-1,n1)=t(k1,n1)+sqrt(sq)
!c**

1 continue 

if (k1-1.eq.1) return
ntwk=ntwk+1      
mxwk(ntwk)=n1      
mzwk(ntwk)=k1-1
      
return

end subroutine trat4

!c                        *************
!c                        * Sub trat5 *
!c                        *************

subroutine trat5(n1,k1,mxwk,mzwk,ntwk)
!c******************************************************************
!c* This sub is to find the values of t(k1,n1-1),t(k1-1,n1)        *
!c******************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
integer             :: m1,m2
m1=0
m2=0
1 continue 
if (t(k1-1,n1+1).ne.0) then
     s1=(s(k1,n1+1)+s(k1-1,n1)+s(k1-1,n1+1)+s(k1,n1))*.25
     sq=hh*s1*s1*2-(t(k1-1,n1+1)-t(k1,n1))**2
     if (sq.lt.0.0)sq=2*s1*s1*h*h
     t(k1-1,n1)=t(k1,n1+1)+sqrt(sq)
     t(k1-1,n1)=max(t(k1-1,n1),t(k1,n1))
     m1=1
endif
       
if (t(k1+1,n1-1).ne.0.0)then
   s1=(s(k1,n1-1)+s(k1+1,n1)+s(k1,n1)+s(k1+1,n1-1))*.25
   sq=hh*s1*s1*2-(t(k1,n1)-t(k1+1,n1-1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1,n1-1)=t(k1+1,n1)+sqrt(sq)
   t(k1,n1-1)=max(t(k1,n1-1),t(k1,n1))
   m2=1
endif
      
if (m1.eq.0.and.m2.eq.1) then
   s1=((s(k1,n1)+s(k1,n1+1)+s(k1,n1-1))/3+s(k1-1,n1))*.5
   sq=hh*s1*s1-.25*(t(k1,n1+1)-t(k1,n1-1))**2
   if (sq.lt.0.0) sq=s1*s1*h*h
   t(k1-1,n1)=t(k1,n1)+sqrt(sq)
endif

if (m1.eq.1.and.m2.eq.0) then
   s1=(s(k1,n1-1)+(s(k1,n1)+s(k1+1,n1)+s(k1-1,n1))/3)*.5
   sq=hh*s1*s1-.25*(t(k1+1,n1)-t(k1-1,n1))**2
   if (sq.lt.0.0) sq=s1*s1*h*h
   t(k1,n1-1)=t(k1,n1)+sqrt(sq)
endif
  
if (m1+m2.eq.0) then
   if (t(k1,n1+2).ne.0.0) then
      s1=(s(k1-1,n1+1)+(s(k1,n1+1)+s(k1,n1)+s(k1,n1+2))/3)*.5
      sq=hh*s1*s1-.25*(t(k1,n1)-t(k1,n1+2))**2
      if (sq.lt.0.0) sq=s1*s1*h*h
      t(k1-1,n1+1)=t(k1,n1+1)+sqrt(sq)
      t(k1-1,n1+1)=max(t(k1,n1+1),t(k1,n1))

      if (n1+1.eq.nx.or.k1-1.eq.1) goto 11
      ntwk=ntwk+1
      mxwk(ntwk)=n1+1
      mzwk(ntwk)=k1-1
   11 continue
   endif

   if (t(k1+2,n1).ne.0.0) then
      s1=((s(k1+1,n1)+s(k1,n1)+s(k1+2,n1))/3+s(k1+1,n1-1))*.5
      sq=hh*s1*s1-.25*(t(k1,n1)-t(k1+2,n1))**2
      if (sq.lt.0.0) sq=s1*s1*h*h
      t(k1+1,n1-1)=t(k1+1,n1)+sqrt(sq)
      t(k1+1,n1-1)=max(t(k1+1,n1-1),t(k1,n1))

      if (n1-1.eq.1.or.k1+1.eq.nz) goto 13
      ntwk=ntwk+1
      mxwk(ntwk)=n1-1
      mzwk(ntwk)=k1+1
   13 continue
   endif

   goto 1
endif
      
if (n1-1.eq.1) goto 20
ntwk=ntwk+1
mxwk(ntwk)=n1-1
mzwk(ntwk)=k1

20 continue 
if (k1-1.eq.1) goto 21

ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1-1

21    continue
      
return
      
end subroutine trat5

!c                        *************
!c                        * Sub trat6 *
!c                        *************

subroutine trat6(n1,k1,mxwk,mzwk,ntwk)
!c***************************************************************
!c* This sub is to find the values of t(k1-1,n1) & t(k1,n1+1)   *
!c***************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
integer             :: m1,m2      

m1=0      
m2=0

1 continue 

if (t(k1+1,n1+1).ne.0.0) then
   m1=1
   s1=(s(k1+1,n1)+s(k1,n1+1)+s(k1,n1)+s(k1+1,n1+1))*.25
   sq=hh*s1*s1*2-(t(k1+1,n1+1)-t(k1,n1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1,n1+1)=t(k1+1,n1)+sqrt(sq)
   t(k1,n1+1)=max(t(k1,n1+1),t(k1,n1))
endif

if (t(k1-1,n1-1).ne.0.0) then
   m2=1
   s1=(s(k1,n1-1)+s(k1-1,n1)+s(k1,n1)+s(k1-1,n1-1))*.25
   sq=hh*s1*s1*2-(t(k1,n1)-t(k1-1,n1-1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1-1,n1)=t(k1,n1-1)+sqrt(sq)
   t(k1-1,n1)=max(t(k1-1,n1),t(k1,n1))
endif

if (m1.eq.0.and.m2.eq.1) then
   s1=(s(k1,n1+1)+(s(k1,n1)+s(k1-1,n1)+s(k1+1,n1))/3)*.5
   sq=hh*s1*s1-.25*(t(k1-1,n1)-t(k1+1,n1))**2
   if (sq.lt.0.0) sq=s1*s1*h*h
   t(k1,n1+1)=t(k1,n1)+sqrt(sq)
endif
      
if (m1.eq.1.and.m2.eq.0) then
   s1=(s(k1-1,n1)+(s(k1,n1)+s(k1,n1-1)+s(k1,n1+1))/3)*.5
   sq=hh*s1*s1-.25*(t(k1,n1-1)-t(k1,n1+1))**2
   if (sq.lt.0.0) sq=s1*s1*h*h
   t(k1-1,n1)=t(k1,n1)+sqrt(sq)
endif
      
if (m1+m2.eq.0) then
   if (t(k1+2,n1).ne.0.0) then
      s1=(s(k1+1,n1+1)+(s(k1+1,n1)+s(k1,n1)+s(k1+2,n1))/3)*.5
      sq=hh*s1*s1-.25*(t(k1,n1)-t(k1+2,n1))**2
      if (sq.lt.0.0) sq=s1*s1*h*h
      t(k1+1,n1+1)=t(k1+1,n1)+sqrt(sq)
      t(k1+1,n1+1)=max(t(k1+1,n1+1),t(k1,n1))

      if (n1+1.eq.nx.or.k1+1.eq.nz) goto 11
      ntwk=ntwk+1
      mxwk(ntwk)=n1+1
      mzwk(ntwk)=k1+1
   11 continue
   endif
      
   if (t(k1,n1-2).ne.0.0) then
      s1=((s(k1,n1-1)+s(k1,n1)+s(k1,n1-2))/3+s(k1-1,n1-1))*.5
      sq=hh*s1*s1-.25*(t(k1,n1)-t(k1,n1-2))**2
      if (sq.lt.0.0) sq=s1*s1*h*h
      t(k1-1,n1-1)=t(k1,n1-1)+sqrt(sq)
      t(k1-1,n1-1)=max(t(k1-1,n1-1),t(k1,n1))

      if (n1-1.eq.1.or.k1-1.eq.1) goto 21
      ntwk=ntwk+1
      mxwk(ntwk)=n1-1
      mzwk(ntwk)=k1-1
   21 continue
   endif

   goto 1
endif

if (k1-1.eq.1) goto 25
ntwk=ntwk+1      
mxwk(ntwk)=n1      
mzwk(ntwk)=k1-1
25 continue 
if (n1+1.eq.nx) goto 26      
ntwk=ntwk+1      
mxwk(ntwk)=n1+1      
mzwk(ntwk)=k1
26 continue
return
      
end subroutine trat6

!c                        *************
!c                        * Sub trat7 *
!c                        *************
       
subroutine trat7(n1,k1,mxwk,mzwk,ntwk)
!c*************************************************************
!c* This sub is to find the values of t(k1-1,n1) & t(k1,n1+1) *
!c*           & t(k1,n1-1)                                    *
!c*************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
s1=(s(k1,n1-1)+s(k1+1,n1)+s(k1,n1)+s(k1+1,n1-1))*.25      
sq=2*hh*s1*s1-(t(k1,n1)-t(k1+1,n1-1))**2      
if (sq.lt.0.0) sq=2*s1*s1*h*h
t(k1,n1-1)=t(k1+1,n1)+sqrt(sq)
t(k1,n1-1)=max(t(k1,n1-1),t(k1,n1))

s1=(s(k1+1,n1)+s(k1,n1+1)+s(k1,n1)+s(k1+1,n1+1))*.25
sq=hh*s1*s1*2-(t(k1,n1)-t(k1+1,n1+1))**2
if (sq.lt.0.0) sq=2*s1*s1*h*h
t(k1,n1+1)=t(k1+1,n1)+sqrt(sq)      
t(k1,n1+1)=max(t(k1,n1+1),t(k1,n1))
      
s1=((s(k1,n1)+s(k1,n1+1)+s(k1,n1-1))/3+s(k1-1,n1))*.5      
sq=hh*s1*s1-.25*(t(k1,n1+1)-t(k1,n1-1))**2      
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1-1,n1)=t(k1,n1)+sqrt(sq)

if (k1-1.eq.1) goto 60
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1-1

60 continue
if (n1+1.eq.nx) goto 62
ntwk=ntwk+1
mxwk(ntwk)=n1+1
mzwk(ntwk)=k1

62 continue
if (n1-1.eq.1) goto 159
ntwk=ntwk+1
mxwk(ntwk)=n1-1
mzwk(ntwk)=k1

159 continue
return
    
end subroutine trat7

!c                        *************
!c                        * Sub trat8 *
!c                        *************

subroutine trat8(n1,k1,mxwk,mzwk,ntwk)
!c**************************************************************
!c* This sub is to find the values of t(k1+1,n1)               *
!c**************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq

if (t(k1+1,n1-1).ne.0.0) then
   s1=(s(k1+1,n1)+s(k1,n1-1)+s(k1,n1)+s(k1+1,n1-1))*.25
   sq=hh*s1*s1*2-(t(k1,n1)-t(k1+1,n1-1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1+1,n1)=t(k1,n1-1)+sqrt(sq)
   t(k1+1,n1)=max(t(k1+1,n1),t(k1,n1))
   goto 1
endif
      
if (t(k1+1,n1+1).ne.0.0) then
   s1=(s(k1,n1+1)+s(k1+1,n1)+s(k1,n1)+s(k1+1,n1+1))*.25
   sq=2*hh*s1*s1-(t(k1,n1)-t(k1+1,n1+1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1+1,n1)=t(k1,n1+1)+sqrt(sq)
   t(k1+1,n1)=max(t(k1+1,n1),t(k1,n1))
   goto 1
endif
      
s1=(s(k1+1,n1)+(s(k1,n1)+s(k1,n1-1)+s(k1,n1+1))/3)*.5      
sq=hh*s1*s1-.25*(t(k1,n1-1)-t(k1,n1+1))**2
      
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1+1,n1)=t(k1,n1)+sqrt(sq)

1 continue
if (k1+1.eq.nz) goto 161
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1+1

161 continue
      
return
      
end subroutine trat8


!c                        *************
!c                        * Sub trat9 *
!c                        *************
       
subroutine trat9(n1,k1,mxwk,mzwk,ntwk)
!c**************************************************************
!c* This sub is to find the values of t(k1,n1-1) & t(k1+1,n1)  *
!c**************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
integer             :: m1,m2      
m1=0      
m2=0
1 continue
if (t(k1+1,n1+1).ne.0.0) then
   m1=1
   s1=(s(k1+1,n1)+s(k1,n1+1)+s(k1+1,n1+1)+s(k1,n1))*.25
   sq=2*hh*s1*s1-(t(k1+1,n1+1)-t(k1,n1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1+1,n1)=t(k1,n1+1)+sqrt(sq)
   t(k1+1,n1)=max(t(k1+1,n1),t(k1,n1))
endif
      
if (t(k1-1,n1-1).ne.0.0) then
   m2=1
   s1=(s(k1,n1-1)+s(k1-1,n1)+s(k1,n1)+s(k1-1,n1-1))*.25
   sq=2*hh*s1*s1-(t(k1-1,n1-1)-t(k1,n1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1,n1-1)=t(k1-1,n1)+sqrt(sq)
   t(k1,n1-1)=max(t(k1,n1-1),t(k1,n1))
endif
      
if (m1.eq.0.and.m2.eq.1) then
   s1=(s(k1+1,n1)+(s(k1,n1)+s(k1,n1+1)+s(k1,n1-1))/3)*.5
   sq=hh*s1*s1-.25*(t(k1,n1+1)-t(k1,n1-1))**2
   if (sq.lt.0.0) sq=s1*s1*h*h
   t(k1+1,n1)=t(k1,n1)+sqrt(sq)
endif
      
if (m1.eq.1.and.m2.eq.0) then
   s1=(s(k1,n1-1)+(s(k1,n1)+s(k1-1,n1)+s(k1+1,n1))/3)*.5
   sq=hh*s1*s1-.25*(t(k1-1,n1)-t(k1+1,n1))**2
   if (sq.lt.0.0) sq=s1*s1*h*h
   t(k1,n1-1)=t(k1,n1)+sqrt(sq)
endif
      
if (m1+m2.eq.0) then
   if (t(k1,n1+2).ne.0.0) then
      s1=((s(k1,n1+1)+s(k1,n1+2)+s(k1,n1))/3+s(k1+1,n1+1))*.5
      sq=hh*s1*s1-.25*(t(k1,n1+2)-t(k1,n1))**2
      if (sq.lt.0.0) sq=s1*s1*h*h
      t(k1+1,n1+1)=t(k1,n1+1)+sqrt(sq)
      t(k1+1,n1+1)=max(t(k1+1,n1+1),t(k1,n1))

      if (n1+1.eq.nx.or.k1+1.eq.nz) goto 6
      ntwk=ntwk+1
      mxwk(ntwk)=n1+1
      mzwk(ntwk)=k1+1
   6  continue
   endif

   if (t(k1-2,n1).ne.0.0) then
      s1=((s(k1-1,n1)+s(k1-2,n1)+s(k1,n1))/3+s(k1-1,n1-1))*.5
      sq=hh*s1*s1-.25*(t(k1-2,n1)-t(k1,n1))**2
      if (sq.lt.0.0) sq=s1*s1*h*h
      t(k1-1,n1-1)=t(k1-1,n1)+sqrt(sq)
      t(k1-1,n1-1)=max(t(k1-1,n1-1),t(k1,n1))

      if (n1-1.eq.1.or.k1-1.eq.1) goto 8
      ntwk=ntwk+1
      mxwk(ntwk)=n1-1
      mzwk(ntwk)=k1-1
   8  continue
   endif

   goto 1
endif
      
if (n1-1.eq.1) goto 100
   ntwk=ntwk+1
   mxwk(ntwk)=n1-1
   mzwk(ntwk)=k1

100 continue

if (k1+1.eq.nz) goto 110
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1+1

110   continue
      
return
      
end subroutine trat9

!c                        *************
!c                        * Sub trat10 *
!c                        *************

subroutine trat10(n1,k1,mxwk,mzwk,ntwk)
!c*************************************************************
!c* This sub is to find the values of t(k1,n1+1) & t(k1+1,n1) *
!c*************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
integer             :: m1,m2      
      
m1=0      
m2=0
1 continue
if (t(k1-1,n1+1).ne.0.0) then
   m1=1
   s1=(s(k1-1,n1)+s(k1,n1+1)+s(k1,n1)+s(k1-1,n1+1))*.25
   sq=2*hh*s1*s1-(t(k1,n1)-t(k1-1,n1+1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1,n1+1)=t(k1-1,n1)+sqrt(sq)
   t(k1,n1+1)=max(t(k1,n1+1),t(k1,n1))
endif
      
if (t(k1+1,n1-1).ne.0.0) then
   m2=1
   s1=(s(k1,n1-1)+s(k1+1,n1)+s(k1,n1)+s(k1+1,n1-1))*.25
   sq=2*hh*s1*s1-(t(k1,n1)-t(k1+1,n1-1))**2
   if (sq.lt.0.0) sq=2*s1*s1*h*h
   t(k1+1,n1)=t(k1,n1-1)+sqrt(sq)
   t(k1+1,n1)=max(t(k1+1,n1),t(k1,n1))
endif
   
if (m1.eq.0.and.m2.eq.1) then
   s1=(s(k1,n1+1)+(s(k1,n1)+s(k1-1,n1)+s(k1+1,n1))/3)*.5
   sq=hh*s1*s1-.25*(t(k1-1,n1)-t(k1+1,n1))**2
   if (sq.lt.0.0) sq=s1*s1*h*h
   t(k1,n1+1)=t(k1,n1)+sqrt(sq)
endif
  
if (m1.eq.1.and.m2.eq.0) then
   s1=((s(k1,n1)+s(k1,n1+1)+s(k1,n1-1))/3+s(k1+1,n1))*.5
   sq=hh*s1*s1-.25*(t(k1,n1+1)-t(k1,n1-1))**2
   if (sq.lt.0.0) sq=s1*s1*h*h
   t(k1+1,n1)=t(k1,n1)+sqrt(sq)
endif
      
if (m1+m2.eq.0.0) then
   if (t(k1-2,n1).ne.0.0) then
      s1=(s(k1-1,n1+1)+(s(k1-1,n1)+s(k1,n1)+s(k1-2,n1))/3)*.5
      sq=hh*s1*s1-.25*(t(k1-2,n1)-t(k1,n1))**2
      if (sq.lt.0.0) sq=s1*s1*h*h
      t(k1-1,n1+1)=t(k1-1,n1)+sqrt(sq)
      t(k1-1,n1+1)=max(t(k1-1,n1+1),t(k1,n1))

      if (n1+1.eq.nx.or.k1-1.eq.1) goto 11
      ntwk=ntwk+1
      mxwk(ntwk)=n1+1
      mzwk(ntwk)=k1-1
   11 continue
   endif

   if (t(k1,n1-2).ne.0.0) then
      s1=((s(k1,n1-1)+s(k1,n1)+s(k1,n1-2))/3+s(k1+1,n1-1))*.5
      sq=hh*s1*s1-.25*(t(k1,n1)-t(k1,n1-2))**2
      if (sq.lt.0.0) sq=s1*s1*h*h
      t(k1+1,n1-1)=t(k1,n1-1)+sqrt(sq)
      t(k1+1,n1-1)=max(t(k1+1,n1-1),t(k1,n1))

      if (n1-1.eq.1.or.k1+1.eq.nz) goto 21
      ntwk=ntwk+1
      mxwk(ntwk)=n1-1
      mzwk(ntwk)=k1+1
   21 continue
   endif

   goto 1
endif 

if (n1+1.eq.nx) goto 100
ntwk=ntwk+1
mxwk(ntwk)=n1+1
mzwk(ntwk)=k1

100 continue
if (k1+1.eq.nz) goto 110
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1+1

110  continue
      
return
      
end subroutine trat10

!c                        *************
!c                        * Sub trat11 *
!c                        *************
       
subroutine trat11(n1,k1,mxwk,mzwk,ntwk)
!c*************************************************************
!c* This sub is to find the values of t(k1,n1+1) & t(k1,n1-1) *
!c*          & t(k1+1,n1)                                     *
!c*************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
      
s1=(s(k1,n1+1)+s(k1-1,n1)+s(k1,n1)+s(k1-1,n1+1))*.25
sq=2*hh*s1*s1-(t(k1,n1)-t(k1-1,n1+1))**2
if (sq.lt.0.0) sq=2*s1*s1*h*h
t(k1,n1+1)=t(k1-1,n1)+sqrt(sq)
t(k1,n1+1)=max(t(k1,n1+1),t(k1,n1))

s1=(s(k1,n1-1)+s(k1-1,n1)+s(k1,n1)+s(k1-1,n1-1))*.25
sq=2*hh*s1*s1-(t(k1,n1)-t(k1-1,n1-1))**2
if (sq.lt.0.0) sq=2*s1*s1*h*h
t(k1,n1-1)=t(k1-1,n1)+sqrt(sq)
t(k1,n1-1)=max(t(k1,n1-1),t(k1,n1))

s1=(s(k1+1,n1)+(s(k1,n1)+s(k1,n1+1)+s(k1,n1-1))/3)*.5
sq=hh*s1*s1-.25*(t(k1,n1+1)-t(k1,n1-1))**2
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1+1,n1)=t(k1,n1)+sqrt(sq)

if (n1+1.eq.nx) goto 10
ntwk=ntwk+1      
mxwk(ntwk)=n1+1      
mzwk(ntwk)=k1

10 continue
if (n1-1.eq.1) goto 20
ntwk=ntwk+1
mxwk(ntwk)=n1-1      
mzwk(ntwk)=k1

20 continue
if (k1+1.eq.nz) goto 30
ntwk=ntwk+1
mxwk(ntwk)=n1      
mzwk(ntwk)=k1+1

30 continue
      
return
      
end subroutine trat11

!c                        *************
!c                        * Sub trat12 *
!c                        *************
       
subroutine trat12(n1,k1,mxwk,mzwk,ntwk)
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
!write(*,*)'sub=12'

!c***************************************************************
!c* This sub is to find the values of t(k1-1,n1) & t(k1+1,n1)   *
!c***************************************************************
      
s1=(s(k1-1,n1)+(s(k1,n1)+s(k1,n1+1)+s(k1,n1-1))/3)*.5      
sq=hh*s1*s1-.25*(t(k1,n1-1)-t(k1,n1+1))**2      
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1-1,n1)=t(k1,n1)+sqrt(sq)
      
s1=(s(k1+1,n1)+(s(k1,n1)+s(k1,n1-1)+s(k1,n1+1))/3)*.5      
sq=hh*s1*s1-.25*(t(k1,n1-1)-t(k1,n1+1))**2     
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1+1,n1)=t(k1,n1)+sqrt(sq)

if (k1-1.eq.1) goto 10
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1-1

10 continue
if (k1+1.eq.nz) goto 20
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1+1

20 continue
      
return
      
end subroutine trat12

!c                        *************
!c                        * Sub trat13 *
!c                        *************
       
subroutine trat13(n1,k1,mxwk,mzwk,ntwk)
!c*************************************************************
!c* This sub is to find the values of t(k1,n1-1) & t(k1-1,n1) *
!c*                   & t(k1+1,n1)                            *
!c*************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq
     
s1=(s(k1,n1+1)+s(k1-1,n1)+s(k1,n1)+s(k1-1,n1+1))*.25      
sq=2*hh*s1*s1-(t(k1,n1)-t(k1-1,n1+1))**2
if (sq.lt.0.0) sq=2*s1*s1*h*h
t(k1-1,n1)=t(k1,n1+1)+sqrt(sq)
t(k1-1,n1)=max(t(k1-1,n1),t(k1,n1))
      
s1=(s(k1+1,n1)+s(k1,n1+1)+s(k1+1,n1+1)+s(k1,n1))*.25      
sq=2*hh*s1*s1-(t(k1+1,n1+1)-t(k1,n1))**2       
if (sq.lt.0.0) sq=2*s1*s1*h*h
t(k1+1,n1)=t(k1,n1+1)+sqrt(sq)
t(k1+1,n1)=max(t(k1+1,n1),t(k1,n1))
      
s1=((s(k1,n1)+s(k1-1,n1)+s(k1+1,n1))/3+s(k1,n1-1))*.5
sq=hh*s1*s1-.25*(t(k1-1,n1)-t(k1+1,n1))**2
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1,n1-1)=t(k1,n1)+sqrt(sq)
      
if (n1-1.eq.1) goto 10
ntwk=ntwk+1
mxwk(ntwk)=n1-1
mzwk(ntwk)=k1

10 continue

if (k1-1.eq.1) goto 20
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1-1
 
20 continue
if (k1+1.eq.nz) goto 30
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1+1

30 continue

return

end subroutine trat13

!c                        *************
!c                        * Sub trat14 *
!c                        *************

subroutine trat14(n1,k1,mxwk,mzwk,ntwk)
!c**************************************************************
!c* This sub is to find the values of t(k1-1,n1) & t(k1+1,n1)  *
!c*                   & t(k1,n1+1)                             *
!c**************************************************************
integer, intent(in) :: n1,k1
integer,intent(inout) :: mxwk(:),mzwk(:),ntwk
real                :: s1,sq

s1=(s(k1,n1-1)+s(k1+1,n1)+s(k1+1,n1-1)+s(k1,n1))*.25
sq=2*hh*s1*s1-(t(k1+1,n1-1)-t(k1,n1))**2
if (sq.lt.0.0) sq=2*s1*s1*h*h
t(k1+1,n1)=t(k1,n1-1)+sqrt(sq)
t(k1+1,n1)=max(t(k1+1,n1),t(k1,n1))

s1=(s(k1,n1-1)+s(k1-1,n1)+s(k1,n1)+s(k1-1,n1-1))*.25
sq=2*hh*s1*s1-(t(k1-1,n1-1)-t(k1,n1))**2
if (sq.lt.0.0) sq=2*s1*s1*h*h
t(k1-1,n1)=t(k1,n1-1)+sqrt(sq)
t(k1-1,n1)=max(t(k1-1,n1),t(k1,n1))

s1=((s(k1,n1)+s(k1-1,n1)+s(k1+1,n1))/3+s(k1,n1+1))*.5
sq=hh*s1*s1-.25*(t(k1-1,n1)-t(k1+1,n1))**2
if (sq.lt.0.0) sq=s1*s1*h*h
t(k1,n1+1)=t(k1,n1)+sqrt(sq)

if (k1-1.eq.1) goto 10
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1-1

10 continue
if (k1+1.eq.nz) goto 20
ntwk=ntwk+1
mxwk(ntwk)=n1
mzwk(ntwk)=k1+1

20 continue
if (n1+1.eq.nx) goto 30
ntwk=ntwk+1
mxwk(ntwk)=n1+1
mzwk(ntwk)=k1

30 continue
      
return
      
end subroutine trat14

!c                        *************
!c                        * Sub coners *
!c                        *************
      
subroutine coners()

!c*****************************************************************
!c* This sub is to consider the four coners of the model          *
!c*****************************************************************
real                :: s1,sq

s1=(s(1,1)+s(2,2))*.5
s1=(s(1,1)+s(2,2)+s(2,1)+s(1,2))*.25
sq=2*hh*s1*s1-(t(1,2)-t(2,1))**2
t(1,1)=t(2,2)+sqrt(sq)
s1=(s(1,nx)+s(2,nx-1))*.5
s1=(s(1,nx)+s(2,nx-1)+s(1,nx-1)+s(2,nx))*.25
sq=2*h*s1*s1*h-(t(1,nx-1)-t(2,nx))**2
t(1,nx)=t(2,nx-1)+sqrt(sq)
s1=(s(nz,nx)+s(nz-1,nx-1)+s(nz-1,nx)+s(nz,nx-1))*.25
sq=2*h*s1*s1*h-(t(nz-1,nx)-t(nz,nx-1))**2
t(nz,nx)=t(nz-1,nx-1)+sqrt(sq)
s1=(s(nz,1)+s(nz-1,2)+s(nz-1,1)+s(nz,2))*.25
sq=2*hh*s1*s1-(t(nz-1,1)-t(nz,2))**2
t(nz,1)=t(nz-1,2)+sqrt(sq)
      
return
      
end subroutine coners

subroutine arr1(n,arr,mx,mz)
real,parameter :: aln2i=1./0.69314718,tiny=1.e-5
integer, intent(in) :: n
real   , intent(in) :: arr(:,:)
integer, intent(inout):: mx(:),mz(:)
integer        :: lognb2, m, nn, k, j, i, l, n1, n2
       
lognb2=int(alog(float(n))*aln2i+tiny)       
m=n
       
do 12 nn=1,lognb2
   m=m/2
   k=n-m
   do 11 j=1,k
      i=j
3     continue
      l=i+m
      if (arr(mz(l),mx(l)).lt.arr(mz(i),mx(i))) then
         n1=mx(i)
         n2=mz(i)

         mx(i)=mx(l)
         mz(i)=mz(l)

         mx(l)=n1
         mz(l)=n2
         i=i-m
         if (i.ge.1) goto 3
      endif
11    continue
12    continue

return

end subroutine arr1

!----------------------------------------------------
subroutine trat_plane(s_inner,s_outer,t_inner,t_outer)
real,   intent(in)  :: s_inner(:),s_outer(:),t_inner(:)
real,   intent(out) :: t_outer(:)
integer :: min_number,max_number,min_loc(nx),max_loc(nx)
real    :: s_ave,sq,t_imax_temp
integer :: i,j,n_array,i_min,i_max
n_array=nx-2
call local_min(t_inner,n_array,min_number,min_loc)
call local_max(t_inner,n_array,max_number,max_loc)
! calculate the traveltime just outside the local minimum points
do i=1,min_number
   i_min=min_loc(i)
   if (i_min==1 .or. i_min==n_array ) then
      t_outer(i_min)=t_inner(i_min)+h*average(s_outer(i_min),s_inner(i_min))
   else
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
   i_min=min_loc(i);   i_max=0
   loop1: do j=1,max_number
      if (max_loc(j) > i_min ) then
         i_max=max_loc(j);   exit loop1
      endif
   enddo loop1
   do j=i_min+1,i_max
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
   i_min=min_loc(i);   i_max=2*nx
   loop2: do j=max_number,1,-1
      if ( max_loc(j) < i_min ) then
         i_max=max_loc(j);   exit loop2
      endif
   enddo loop2
   if (i_max /= 2*nx) then
      do j=i_min-1,i_max,-1
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
end subroutine trat_plane

end module module_pweikfd2d_fir
