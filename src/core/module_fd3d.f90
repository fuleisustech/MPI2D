module module_fd3d
use module_global
use module_boundary3d
implicit none

integer, private :: iz,iy,ix
real,    private :: alpha,temp1,temp2,temp3
integer, private :: sp   ! Used in elastic 3D modeling Starting point in the z direction


interface a3d_abc_kernal
   module procedure a3d_abc_kernal_whole
   module procedure a3d_abc_kernal_working
   module procedure a2d_abc_kernal_whole_Lax_Wendroff     ! modified by Bowen Guo  
end interface a3d_abc_kernal

interface e3d_abc_kernal_vel
   module procedure e3d_abc_kernal_u_whole
   module procedure e3d_abc_kernal_v_whole
   module procedure e3d_abc_kernal_w_whole
end interface e3d_abc_kernal_vel                   ! modified by Bowen Guo

interface e3d_abc_kernal_stress
   module procedure e3d_abc_kernal_totalstress_whole
   module procedure e3d_abc_kernal_tauxx_whole
   module procedure e3d_abc_kernal_tauyy_whole
   module procedure e3d_abc_kernal_tauzz_whole
   module procedure e3d_abc_kernal_tauxy_whole  
   module procedure e3d_abc_kernal_tauxz_whole
   module procedure e3d_abc_kernal_tauyz_whole
end interface e3d_abc_kernal_stress                ! modifidd by Bowen Guo














contains

!=====================================================================

subroutine a3d_abc_kernal_whole(p,p0,p1,v,kappa,dtdx,c1t3,&
   nz_pml,ny_pml,nx_pml,fd_order)
integer, intent(in) :: nz_pml,ny_pml,nx_pml,fd_order
real, intent(inout) :: p(:,:,:)
real, intent(in)    :: dtdx,c1t3,p0(:,:,:),p1(:,:,:),v(:,:,:),kappa(:,:,:)
if (fd_order==22) then
   !forall (iz=2:nz_pml-1,iy=2:ny_pml-1,ix=2:nx_pml-1)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha,temp1,temp2,temp3)
   do ix=2,nx_pml-1
      do iy=2,ny_pml-1
         do iz=2,nz_pml-1
            alpha=(v(iz,iy,ix)*dtdx)**2.0
            temp1=2.0+c1t3*alpha-kappa(iz,iy,ix)*kappa(iz,iy,ix)
            temp2=1.0 - kappa(iz,iy,ix)
            temp3=1.0 + kappa(iz,iy,ix)
            p(iz,iy,ix)=(temp1*p1(iz,iy,ix)-temp2*p0(iz,iy,ix)&
               +alpha*(C22*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                           +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                           +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))))/temp3
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
elseif (fd_order==24) then
   !forall (iz=3:nz_pml-2,iy=3:ny_pml-2,ix=3:nx_pml-2)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha,temp1,temp2,temp3)
   do ix=3,nx_pml-2
      do iy=3,ny_pml-2
         do iz=3,nz_pml-2
            alpha=(v(iz,iy,ix)*dtdx)**2.0
            temp1=2.0+c1t3*alpha-kappa(iz,iy,ix)*kappa(iz,iy,ix)
            temp2=1.0 - kappa(iz,iy,ix)
            temp3=1.0 + kappa(iz,iy,ix)
            p(iz,iy,ix)=(temp1*p1(iz,iy,ix)-temp2*p0(iz,iy,ix) &
               +alpha*(C42*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                           +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                           +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                      +C43*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                           +p1(iz,iy+2,ix)+p1(iz,iy-2,ix) &
                           +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))))/temp3
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
elseif (fd_order==26) then
   !forall (iz=4:nz_pml-3,iy=4:ny_pml-3,ix=4:nx_pml-3)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha,temp1,temp2,temp3)
   do ix=4,nx_pml-3
      do iy=4,ny_pml-3
         do iz=4,nz_pml-3
            alpha=(v(iz,iy,ix)*dtdx)**2.0
            temp1=2.0+c1t3*alpha-kappa(iz,iy,ix)*kappa(iz,iy,ix)
            temp2=1.0 - kappa(iz,iy,ix)
            temp3=1.0 + kappa(iz,iy,ix)
            p(iz,iy,ix)=(temp1*p1(iz,iy,ix) - temp2*p0(iz,iy,ix) &
               +alpha*(C62*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                           +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                           +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                      +C63*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                           +p1(iz,iy+2,ix)+p1(iz,iy-2,ix) &
                           +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))&
                      +C64*(p1(iz,iy,ix+3)+p1(iz,iy,ix-3) &
                           +p1(iz,iy+3,ix)+p1(iz,iy-3,ix) &
                           +p1(iz+3,iy,ix)+p1(iz-3,iy,ix))))/temp3
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
elseif (fd_order==28) then
   !forall (iz=5:nz_pml-4,iy=5:ny_pml-4,ix=5:nx_pml-4)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha,temp1,temp2,temp3)
   do ix=5,nx_pml-4
      do iy=5,ny_pml-4
         do iz=5,nz_pml-4
            alpha=(v(iz,iy,ix)*dtdx)**2.0
            temp1=2.0+c1t3*alpha-kappa(iz,iy,ix)*kappa(iz,iy,ix)
            temp2=1.0 - kappa(iz,iy,ix)
            temp3=1.0 + kappa(iz,iy,ix)
            p(iz,iy,ix)=(temp1*p1(iz,iy,ix) - temp2*p0(iz,iy,ix) &
               +alpha*(C82*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                           +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                           +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                      +C83*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                           +p1(iz,iy+2,ix)+p1(iz,iy-2,ix) &
                           +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))&
                      +C84*(p1(iz,iy,ix+3)+p1(iz,iy,ix-3) &
                           +p1(iz,iy+3,ix)+p1(iz,iy-3,ix) &
                           +p1(iz+3,iy,ix)+p1(iz-3,iy,ix))&
                      +C85*(p1(iz,iy,ix+4)+p1(iz,iy,ix-4) &
                           +p1(iz,iy+4,ix)+p1(iz,iy-4,ix) &
                           +p1(iz+4,iy,ix)+p1(iz-4,iy,ix))))/temp3
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
endif
end subroutine a3d_abc_kernal_whole

!=====================================================================

subroutine a3d_abc_kernal_working(p,p0,p1,v,kappa,dtdx,c1t3,&
   npml,nzp,nyp,nxp,fd_order)
integer, intent(in) :: npml,nzp,nyp,nxp,fd_order
real, intent(inout) :: p(:,:,:)
real, intent(in)    :: p0(:,:,:),p1(:,:,:),v(:,:,:),kappa(:,:,:),dtdx,c1t3
if (fd_order==22) then
   !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha,temp1,temp2,temp3)
   do ix=npml+1,nxp
      do iy=npml+1,nyp
         do iz=npml+1,nzp
            alpha=(v(iz,iy,ix)*dtdx)**2.0
            temp1=2.0+c1t3*alpha-kappa(iz,iy,ix)*kappa(iz,iy,ix)
            temp2=1.0 - kappa(iz,iy,ix)
            temp3=1.0 + kappa(iz,iy,ix)
            p(iz,iy,ix)=(temp1*p1(iz,iy,ix) - temp2*p0(iz,iy,ix)&
               +alpha*(C22*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                           +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                           +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))))/temp3
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
elseif (fd_order==24) then
   !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha,temp1,temp2,temp3)
   do ix=npml+1,nxp
      do iy=npml+1,nyp
         do iz=npml+1,nzp
            alpha=(v(iz,iy,ix)*dtdx)**2.0
            temp1=2.0+c1t3*alpha-kappa(iz,iy,ix)*kappa(iz,iy,ix)
            temp2=1.0 - kappa(iz,iy,ix)
            temp3=1.0 + kappa(iz,iy,ix)
            p(iz,iy,ix)=(temp1*p1(iz,iy,ix) - temp2*p0(iz,iy,ix) &
               +alpha*(C42*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                           +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                           +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                      +C43*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                           +p1(iz,iy+2,ix)+p1(iz,iy+2,ix) &
                           +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))))/temp3
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
elseif (fd_order==26) then
   !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha,temp1,temp2,temp3)
   do ix=npml+1,nxp
      do iy=npml+1,nyp
         do iz=npml+1,nzp
            alpha=(v(iz,iy,ix)*dtdx)**2.0
            temp1=2.0+c1t3*alpha-kappa(iz,iy,ix)*kappa(iz,iy,ix)
            temp2=1.0 - kappa(iz,iy,ix)
            temp3=1.0 + kappa(iz,iy,ix)
            p(iz,iy,ix)=(temp1*p1(iz,iy,ix) - temp2*p0(iz,iy,ix) &
               +alpha*(C62*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                           +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                           +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                      +C63*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                           +p1(iz,iy+2,ix)+p1(iz,iy-2,ix) &
                           +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))&
                      +C64*(p1(iz,iy,ix+3)+p1(iz,iy,ix-3) &
                           +p1(iz,iy+3,ix)+p1(iz,iy-3,ix) &
                           +p1(iz+3,iy,ix)+p1(iz-3,iy,ix))))/temp3
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
elseif (fd_order==28) then
   !forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha,temp1,temp2,temp3)
   do ix=npml+1,nxp
      do iy=npml+1,nyp
         do iz=npml+1,nzp 
            alpha=(v(iz,iy,ix)*dtdx)**2.0
            temp1=2.0+c1t3*alpha-kappa(iz,iy,ix)*kappa(iz,iy,ix)
            temp2=1.0 - kappa(iz,iy,ix)
            temp3=1.0 + kappa(iz,iy,ix)
            p(iz,iy,ix)=(temp1*p1(iz,iy,ix) - temp2*p0(iz,iy,ix) &
               +alpha*(C82*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                           +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                           +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                      +C83*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                           +p1(iz,iy+2,ix)+p1(iz,iy-2,ix) &
                           +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))&
                      +C84*(p1(iz,iy,ix+3)+p1(iz,iy,ix-3) &
                           +p1(iz,iy+3,ix)+p1(iz,iy-3,ix) &
                           +p1(iz+3,iy,ix)+p1(iz-3,iy,ix))&
                      +C85*(p1(iz,iy,ix+4)+p1(iz,iy,ix-4) &
                           +p1(iz,iy+4,ix)+p1(iz,iy-4,ix) &
                           +p1(iz+4,iy,ix)+p1(iz-4,iy,ix))))/temp3
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
   !end forall
endif
end subroutine a3d_abc_kernal_working
!===================================================================================
! Written by Bowen Guo to calculate wave field by Lax-Wendroff method
subroutine a2d_abc_kernal_whole_Lax_Wendroff(p,p0,p1,temp1,temp2,temp3,temp4,c,dx,dy,dz,nz_pml,ny_pml,nx_pml,fd_order)
integer,intent(in)    :: nz_pml,ny_pml,nx_pml,fd_order
real,intent(in)       :: dx,dy,dz
real,intent(in)       :: c(:,:,:),temp1(:,:,:),temp2(:,:,:),temp3(:,:,:),temp4(:,:,:),p0(:,:,:),&
p1(:,:,:)
real,intent(out)      :: p(:,:,:)
real,allocatable      :: d2p(:,:,:),d4p(:,:,:)
call allocate_and_initial(d2p,d4p,nz_pml,ny_pml,nx_pml)
call cal_space_laplace_3d_whole(p1,d2p,dx,dy,dz,fd_order,nz_pml,ny_pml,nx_pml)
! Calculate the 4th time derivative
d2p=d2p*(c**2)
call cal_space_laplace_3d_whole(d2p,d4p,dx,dy,dz,fd_order,nz_pml,ny_pml,nx_pml)
p=temp1*d2p+temp2*p0+temp3*p1+temp4*d4p
call deallocate_and_free(d2p,d4p)
end subroutine a2d_abc_kernal_whole_Lax_Wendroff
!============================================================================================
subroutine cal_space_laplace_3d_whole(p,d2p,dx,dy,dz,fd_order,nz_pml,ny_pml,nx_pml)
integer,intent(in)    :: nz_pml,ny_pml,nx_pml,fd_order
real,intent(in)       :: dx,dz,dy
real,intent(in)       :: p(:,:,:)
real,intent(out)      :: d2p(:,:,:)
real,allocatable      :: d2p_dx2(:,:,:),d2p_dy2(:,:,:),d2p_dz2(:,:,:)

call allocate_and_initial(d2p_dx2,d2p_dy2,d2p_dz2,nz_pml,ny_pml,nx_pml)
call cal_space_2nddx_3d_whole(p,d2p_dx2,dx,fd_order,nz_pml,ny_pml,nx_pml)
call cal_space_2nddy_3d_whole(p,d2p_dy2,dy,fd_order,nz_pml,ny_pml,nx_pml)
call cal_space_2nddz_3d_whole(p,d2p_dz2,dz,fd_order,nz_pml,ny_pml,nx_pml)
d2p=d2p_dx2+d2p_dz2+d2p_dy2
call deallocate_and_free(d2p_dx2,d2p_dy2,d2p_dz2)
end subroutine cal_space_laplace_3d_whole

!=====================================================================
subroutine cal_space_2nddx_3d_whole(p,d2p_dx2,dx,fd_order,nz_pml,ny_pml,nx_pml)
integer,intent(in)    :: fd_order,nz_pml,ny_pml,nx_pml
real,intent(in)       :: dx
real,intent(in)       :: p(:,:,:)
real,intent(out)      :: d2p_dx2(:,:,:)
real                  :: C21_temp,C22_temp,C41_temp,C42_temp,C43_temp,C61_temp,C62_temp,C63_temp,&
C64_temp,C65_temp,C81_temp,C82_temp,C83_temp,C84_temp,C85_temp


if (fd_order==42) then
   C21_temp=C21/(dx**2);C22_temp=C22/(dx**2)
   forall (iz=1:nz_pml,iy=1:ny_pml,ix=2:nx_pml-1)
      d2p_dx2(iz,iy,ix)=C21_temp*p(iz,iy,ix)+C22_temp*(p(iz,iy,ix+1)+p(iz,iy,ix-1))
   end forall
elseif (fd_order==44) then
   C41_temp=C41/(dx**2);C42_temp=C42/(dx**2);C43_temp=C43/(dx**2)
   forall (iz=1:nz_pml,iy=1:ny_pml,ix=3:nx_pml-2)
      d2p_dx2(iz,iy,ix)=C41_temp*p(iz,iy,ix)+C42_temp*(p(iz,iy,ix+1)+p(iz,iy,ix-1))+C43_temp*(p(iz,iy,ix+2)+p(iz,iy,ix-2))
   end forall
elseif (fd_order==46) then
   C61_temp=C61/(dx**2);C62_temp=C62/(dx**2);C63_temp=C63/(dx**2);C64_temp=C64/(dx**2)
   forall (iz=1:nz_pml,iy=1:ny_pml,ix=4:nx_pml-3)
      d2p_dx2(iz,iy,ix)=C61_temp*p(iz,iy,ix)+C62_temp*(p(iz,iy,ix+1)+p(iz,iy,ix-1))+C63_temp*(p(iz,iy,ix+2)+p(iz,iy,ix-2))+&
      C64_temp*(p(iz,iy,ix+3)+p(iz,iy,ix-3))
   end forall
elseif (fd_order==48) then
   C81_temp=C81/(dx**2);C82_temp=C82/(dx**2);C83_temp=C83/(dx**2);C84_temp=C84/(dx**2);C85_temp=C85/&
   (dx**2)
   forall (iz=1:nz_pml,iy=1:ny_pml,ix=5:nx_pml-4)
      d2p_dx2(iz,iy,ix)=C81_temp*p(iz,iy,ix)+C82_temp*(p(iz,iy,ix+1)+p(iz,iy,ix-1))+C83_temp*(p(iz,iy,ix+2)+p(iz,iy,ix-2))+&
      C84_temp*(p(iz,iy,ix+3)+p(iz,iy,ix-3))+C85_temp*(p(iz,iy,ix+4)+p(iz,iy,ix-4))
   end forall
endif
end subroutine cal_space_2nddx_3d_whole
!=============================================================================================
subroutine cal_space_2nddy_3d_whole(p,d2p_dx2,dy,fd_order,nz_pml,ny_pml,nx_pml)
integer,intent(in)    :: fd_order,nz_pml,ny_pml,nx_pml
real,intent(in)       :: dy
real,intent(in)       :: p(:,:,:)
real,intent(out)      :: d2p_dx2(:,:,:)
real                  :: C21_temp,C22_temp,C41_temp,C42_temp,C43_temp,C61_temp,C62_temp,C63_temp,&
C64_temp,C65_temp,C81_temp,C82_temp,C83_temp,C84_temp,C85_temp


if (fd_order==42) then
   C21_temp=C21/(dy**2);C22_temp=C22/(dy**2)
   forall (iz=1:nz_pml,iy=2:ny_pml-1,ix=1:nx_pml)
      d2p_dx2(iz,iy,ix)=C21_temp*p(iz,iy,ix)+C22_temp*(p(iz,iy+1,ix)+p(iz,iy-1,ix))
   end forall
elseif (fd_order==44) then
   C41_temp=C41/(dy**2);C42_temp=C42/(dy**2);C43_temp=C43/(dy**2)
   forall (iz=1:nz_pml,iy=3:ny_pml-2,ix=1:nx_pml)
      d2p_dx2(iz,iy,ix)=C41_temp*p(iz,iy,ix)+C42_temp*(p(iz,iy+1,ix)+p(iz,iy-1,ix))+C43_temp*(p(iz,iy+2,ix)+p(iz,iy-2,ix))
   end forall
elseif (fd_order==46) then
   C61_temp=C61/(dy**2);C62_temp=C62/(dy**2);C63_temp=C63/(dy**2);C64_temp=C64/(dy**2)
   forall (iz=1:nz_pml,iy=4:ny_pml-3,ix=1:nx_pml)
d2p_dx2(iz,iy,ix)=C61_temp*p(iz,iy,ix)+C62_temp*(p(iz,iy+1,ix)+p(iz,iy-1,ix))+C63_temp*(p(iz,iy+2,ix)+p(iz,iy-2,ix))+&
      C64_temp*(p(iz,iy+3,ix)+p(iz,iy-3,ix))
   end forall
elseif (fd_order==48) then
   C81_temp=C81/(dy**2);C82_temp=C82/(dy**2);C83_temp=C83/(dy**2);C84_temp=C84/(dy**2);C85_temp=C85/&
   (dy**2)
   forall (iz=1:nz_pml,iy=5:ny_pml-4,ix=1:nx_pml)
      d2p_dx2(iz,iy,ix)=C81_temp*p(iz,iy,ix)+C82_temp*(p(iz,iy+1,ix)+p(iz,iy-1,ix))+C83_temp*(p(iz,iy+2,ix)+p(iz,iy-2,ix))+&
      C84_temp*(p(iz,iy+3,ix)+p(iz,iy-3,ix))+C85_temp*(p(iz,iy+4,ix)+p(iz,iy-4,ix))
   end forall
endif
end subroutine cal_space_2nddy_3d_whole
!=================================================================================================
subroutine cal_space_2nddz_3d_whole(p,d2p_dz2,dz,fd_order,nz_pml,ny_pml,nx_pml)
integer,intent(in)    :: fd_order,nz_pml,ny_pml,nx_pml
real,intent(in)       :: dz
real,intent(in)       :: p(:,:,:)
real,intent(out)      :: d2p_dz2(:,:,:)
real                  :: C21_temp,C22_temp,C41_temp,C42_temp,C43_temp,C61_temp,C62_temp,C63_temp,&
C64_temp,C65_temp,C81_temp,C82_temp,C83_temp,C84_temp,C85_temp
if (fd_order==42) then
   C21_temp=C21/(dz**2);C22_temp=C22/(dz**2)
   forall (iz=2:nz_pml-1,iy=1:ny_pml,ix=1:nx_pml)
      d2p_dz2(iz,iy,ix)=C21_temp*p(iz,iy,ix)+C22_temp*(p(iz+1,iy,ix)+p(iz-1,iy,ix))
   end forall
elseif (fd_order==44) then
   C41_temp=C41/(dz**2);C42_temp=C42/(dz**2);C43_temp=C43/(dz**2)
   forall (iz=3:nz_pml-2,iy=1:ny_pml,ix=1:nx_pml)
      d2p_dz2(iz,iy,ix)=C41_temp*p(iz,iy,ix)+C42_temp*(p(iz+1,iy,ix)+p(iz-1,iy,ix))+C43_temp*(p(iz+2,iy,ix)+p(iz-2,iy,ix))
   end forall
elseif (fd_order==46) then
   C61_temp=C61/(dz**2);C62_temp=C62/(dz**2);C63_temp=C63/(dz**2);C64_temp=C64/(dz**2)
   forall (iz=4:nz_pml-3,iy=1:ny_pml,ix=1:nx_pml)
      d2p_dz2(iz,iy,ix)=C61_temp*p(iz,iy,ix)+C62_temp*(p(iz+1,iy,ix)+p(iz-1,iy,ix))+C63_temp*(p(iz+2,iy,ix)+p(iz-2,iy,ix))+&
      C64_temp*(p(iz+3,iy,ix)+p(iz-3,iy,ix))
   end forall
elseif (fd_order==48) then
C81_temp=C81/(dz**2);C82_temp=C82/(dz**2);C83_temp=C83/(dz**2);C84_temp=C84/(dz**2);C85_temp=C85/&
   (dz**2)
   forall (iz=5:nz_pml-4,iy=1:ny_pml,ix=1:nx_pml)
      d2p_dz2(iz,iy,ix)=C81_temp*p(iz,iy,ix)+C82_temp*(p(iz+1,iy,ix)+p(iz-1,iy,ix))+C83_temp*(p(iz+2,iy,ix)+p(iz-2,iy,ix))+&
      C84_temp*(p(iz+3,iy,ix)+p(iz-3,iy,ix))+C85_temp*(p(iz+4,iy,ix)+p(iz-4,iy,ix))
   end forall
endif
end subroutine cal_space_2nddz_3d_whole



!================================================================================
subroutine e3d_abc_kernal_u_whole(u,tau_x,tau_y,tau_z,temp,den,dtdx,nx_pml,ny_pml,nz_pml,fd_order,isFS,npml)

integer,intent(in)   :: nx_pml,ny_pml,nz_pml,fd_order,npml
logical,intent(in)   :: isFS
real,intent(inout)   :: u(:,:,:)
real,intent(in)      :: tau_x(:,:,:),tau_y(:,:,:),tau_z(:,:,:),den(:,:,:),temp(:,:,:)
real,intent(in)      :: dtdx

if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2+1
endif


if (fd_order==22) then
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha)
   do ix=1,nx_pml-1
      do iy=2,ny_pml
          do iz=sp,nz_pml
             alpha=dtdx/(den(iz,iy,ix)+den(iz,iy,ix+1))*2
             u(iz,iy,ix)=temp(iz,iy,ix)*u(iz,iy,ix)+alpha*S21*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix)+tau_y(iz,iy,ix)-tau_y(iz,iy-1,ix)+tau_z(iz,iy,ix)-tau_z(iz-1,iy,ix))
          enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24) then 
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha)
   do ix=2,nx_pml-2
      do iy=3,ny_pml-1
         do iz=sp,nz_pml-1
             alpha=dtdx/(den(iz,iy,ix)+den(iz,iy,ix+1))*2
             u(iz,iy,ix)=temp(iz,iy,ix)*u(iz,iy,ix)+alpha*(S41*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix)+tau_y(iz,iy,ix)-tau_y(iz,iy-1,ix)+tau_z(iz,iy,ix)-tau_z(iz-1,iy,ix))+S42*(tau_x(iz,iy,ix+2)-tau_x(iz,iy,ix-1)+tau_y(iz,iy+1,ix)-tau_y(iz,iy-2,ix)+tau_z(iz+1,iy,ix)-tau_z(iz-2,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26) then
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha)
   do ix=3,nx_pml-3
      do iy=4,ny_pml-2
         do iz=sp,nz_pml-2
             alpha=dtdx/(den(iz,iy,ix)+den(iz,iy,ix+1))*2 
             u(iz,iy,ix)=temp(iz,iy,ix)*u(iz,iy,ix)+alpha*(S61*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix)+tau_y(iz,iy,ix)-tau_y(iz,iy-1,ix)+tau_z(iz,iy,ix)-tau_z(iz-1,iy,ix))+S62*(tau_x(iz,iy,ix+2)-tau_x(iz,iy,ix-1)+tau_y(iz,iy+1,ix)-tau_y(iz,iy-2,ix)+tau_z(iz+1,iy,ix)-tau_z(iz-2,iy,ix))+S63*(tau_x(iz,iy,ix+3)-tau_x(iz,iy,ix-2)+tau_y(iz,iy+2,ix)-tau_y(iz,iy-3,ix)+tau_z(iz+2,iy,ix)-tau_z(iz-3,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28) then
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha)
   do ix=4,nx_pml-4
      do iy=5,ny_pml-3
         do iz=sp,nz_pml-3
             alpha=dtdx/(den(iz,iy,ix)+den(iz,iy,ix+1))*2 
             u(iz,iy,ix)=temp(iz,iy,ix)*u(iz,iy,ix)+alpha*(S81*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix)+tau_y(iz,iy,ix)-tau_y(iz,iy-1,ix)+tau_z(iz,iy,ix)-tau_z(iz-1,iy,ix))+S82*(tau_x(iz,iy,ix+2)-tau_x(iz,iy,ix-1)+tau_y(iz,iy+1,ix)-tau_y(iz,iy-2,ix)+tau_z(iz+1,iy,ix)-tau_z(iz-2,iy,ix))+S83*(tau_x(iz,iy,ix+3)-tau_x(iz,iy,ix-2)+tau_y(iz,iy+2,ix)-tau_y(iz,iy-3,ix)+tau_z(iz+2,iy,ix)-tau_z(iz-3,iy,ix))+S84*(tau_x(iz,iy,ix+4)-tau_x(iz,iy,ix-3)+tau_y(iz,iy+3,ix)-tau_y(iz,iy-4,ix)+tau_z(iz+3,iy,ix)-tau_z(iz-4,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif 
end subroutine e3d_abc_kernal_u_whole
!========================================================================================
subroutine e3d_abc_kernal_v_whole(fd_order,v,tau_x,tau_y,tau_z,temp,den,dtdx,nx_pml,ny_pml,nz_pml,isFS,npml)

integer,intent(in)   :: nx_pml,ny_pml,nz_pml,fd_order,npml
logical,intent(in)   :: isFS
real,intent(inout)   :: v(:,:,:)
real,intent(in)      :: tau_x(:,:,:),tau_y(:,:,:),tau_z(:,:,:),den(:,:,:),temp(:,:,:)
real,intent(in)      :: dtdx
if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2+1
endif

if (fd_order==22) then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=2,nx_pml
      do iy=1,ny_pml-1
          do iz=sp,nz_pml
             alpha=dtdx/(den(iz,iy,ix)+den(iz,iy+1,ix))*2
             v(iz,iy,ix)=temp(iz,iy,ix)*v(iz,iy,ix)+alpha*S21*(tau_x(iz,iy,ix)-tau_x(iz,iy,ix-1)+tau_y(iz,iy+1,ix)-tau_y(iz,iy,ix)+tau_z(iz,iy,ix)-tau_z(iz-1,iy,ix))
          enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24) then 
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=3,nx_pml-1
      do iy=2,ny_pml-2
         do iz=sp,nz_pml-1
             alpha=dtdx/(den(iz,iy,ix)+den(iz,iy+1,ix))*2
             v(iz,iy,ix)=temp(iz,iy,ix)*v(iz,iy,ix)+alpha*(S41*(tau_x(iz,iy,ix)-tau_x(iz,iy,ix-1)+tau_y(iz,iy+1,ix)-tau_y(iz,iy,ix)+tau_z(iz,iy,ix)-tau_z(iz-1,iy,ix))+S42*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix-2)+tau_y(iz,iy+2,ix)-tau_y(iz,iy-1,ix)+tau_z(iz+1,iy,ix)-tau_z(iz-2,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26) then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=4,nx_pml-2
      do iy=3,ny_pml-3
         do iz=sp,nz_pml-2
             alpha=dtdx/(den(iz,iy,ix)+den(iz,iy+1,ix))*2
             v(iz,iy,ix)=temp(iz,iy,ix)*v(iz,iy,ix)+alpha*(S61*(tau_x(iz,iy,ix)-tau_x(iz,iy,ix-1)+tau_y(iz,iy+1,ix)-tau_y(iz,iy,ix)+tau_z(iz,iy,ix)-tau_z(iz-1,iy,ix))+S62*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix-2)+tau_y(iz,iy+2,ix)-tau_y(iz,iy-1,ix)+tau_z(iz+1,iy,ix)-tau_z(iz-2,iy,ix))+S63*(tau_x(iz,iy,ix+2)-tau_x(iz,iy,ix-3)+tau_y(iz,iy+3,ix)-tau_y(iz,iy-2,ix)+tau_z(iz+2,iy,ix)-tau_z(iz-3,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28) then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=5,nx_pml-3
      do iy=4,ny_pml-4
         do iz=sp,nz_pml-3
             alpha=dtdx/(den(iz,iy,ix)+den(iz,iy+1,ix))*2
             v(iz,iy,ix)=temp(iz,iy,ix)*v(iz,iy,ix)+alpha*(S81*(tau_x(iz,iy,ix)-tau_x(iz,iy,ix-1)+tau_y(iz,iy+1,ix)-tau_y(iz,iy,ix)+tau_z(iz,iy,ix)-tau_z(iz-1,iy,ix))+S82*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix-2)+tau_y(iz,iy+2,ix)-tau_y(iz,iy-1,ix)+tau_z(iz+1,iy,ix)-tau_z(iz-2,iy,ix))+S83*(tau_x(iz,iy,ix+2)-tau_x(iz,iy,ix-3)+tau_y(iz,iy+3,ix)-tau_y(iz,iy-2,ix)+tau_z(iz+2,iy,ix)-tau_z(iz-3,iy,ix))+S84*(tau_x(iz,iy,ix+3)-tau_x(iz,iy,ix-4)+tau_y(iz,iy+4,ix)-tau_y(iz,iy-3,ix)+tau_z(iz+3,iy,ix)-tau_z(iz-4,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif 

end subroutine e3d_abc_kernal_v_whole

!========================================================================================
subroutine e3d_abc_kernal_w_whole(nx_pml,ny_pml,nz_pml,w,tau_x,tau_y,tau_z,temp,den,dtdx,fd_order,isFS,npml)

integer,intent(in)   :: nx_pml,ny_pml,nz_pml,fd_order,npml
logical,intent(in)   :: isFS
real,intent(inout)   :: w(:,:,:)
real,intent(in)      :: tau_x(:,:,:),tau_y(:,:,:),tau_z(:,:,:),den(:,:,:),temp(:,:,:)
real,intent(in)      :: dtdx
if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2
endif

if (fd_order==22) then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=2,nx_pml
      do iy=2,ny_pml
          do iz=sp,nz_pml-1
             alpha=dtdx/(den(iz,iy,ix)+den(iz+1,iy,ix))*2
             w(iz,iy,ix)=temp(iz,iy,ix)*w(iz,iy,ix)+alpha*S21*(tau_x(iz,iy,ix)-tau_x(iz,iy,ix-1)+tau_y(iz,iy,ix)-tau_y(iz,iy-1,ix)+tau_z(iz+1,iy,ix)-tau_z(iz,iy,ix))
          enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24) then 
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=3,nx_pml-1
      do iy=3,ny_pml-1
         do iz=sp,nz_pml-2
             alpha=dtdx/(den(iz,iy,ix)+den(iz+1,iy,ix))*2
             w(iz,iy,ix)=temp(iz,iy,ix)*w(iz,iy,ix)+alpha*(S41*(tau_x(iz,iy,ix)-tau_x(iz,iy,ix-1)+tau_y(iz,iy,ix)-tau_y(iz,iy-1,ix)+tau_z(iz+1,iy,ix)-tau_z(iz,iy,ix))+S42*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix-2)+tau_y(iz,iy+1,ix)-tau_y(iz,iy-2,ix)+tau_z(iz+2,iy,ix)-tau_z(iz-1,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26) then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=4,nx_pml-2
      do iy=4,ny_pml-2
         do iz=sp,nz_pml-3
             alpha=dtdx/(den(iz,iy,ix)+den(iz+1,iy,ix))*2
             w(iz,iy,ix)=temp(iz,iy,ix)*w(iz,iy,ix)+alpha*(S61*(tau_x(iz,iy,ix)-tau_x(iz,iy,ix-1)+tau_y(iz,iy,ix)-tau_y(iz,iy-1,ix)+tau_z(iz+1,iy,ix)-tau_z(iz,iy,ix))+S62*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix-2)+tau_y(iz,iy+1,ix)-tau_y(iz,iy-2,ix)+tau_z(iz+2,iy,ix)-tau_z(iz-1,iy,ix))+S63*(tau_x(iz,iy,ix+2)-tau_x(iz,iy,ix-3)+tau_y(iz,iy+2,ix)-tau_y(iz,iy-3,ix)+tau_z(iz+3,iy,ix)-tau_z(iz-2,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28) then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=5,nx_pml-3
      do iy=5,ny_pml-3
         do iz=sp,nz_pml-4
             alpha=dtdx/(den(iz,iy,ix)+den(iz+1,iy,ix))*2
             w(iz,iy,ix)=temp(iz,iy,ix)*w(iz,iy,ix)+alpha*(S81*(tau_x(iz,iy,ix)-tau_x(iz,iy,ix-1)+tau_y(iz,iy,ix)-tau_y(iz,iy-1,ix)+tau_z(iz+1,iy,ix)-tau_z(iz,iy,ix))+S82*(tau_x(iz,iy,ix+1)-tau_x(iz,iy,ix-2)+tau_y(iz,iy+1,ix)-tau_y(iz,iy-2,ix)+tau_z(iz+2,iy,ix)-tau_z(iz-1,iy,ix))+S83*(tau_x(iz,iy,ix+2)-tau_x(iz,iy,ix-3)+tau_y(iz,iy+2,ix)-tau_y(iz,iy-3,ix)+tau_z(iz+3,iy,ix)-tau_z(iz-2,iy,ix))+S84*(tau_x(iz,iy,ix+3)-tau_x(iz,iy,ix-4)+tau_y(iz,iy+3,ix)-tau_y(iz,iy-4,ix)+tau_z(iz+4,iy,ix)-tau_z(iz-3,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif 
end subroutine e3d_abc_kernal_w_whole

!======================================================================================
subroutine e3d_abc_kernal_totalstress_whole(tau,u,v,w,lambda,dtdx,nx_pml,ny_pml,nz_pml,fd_order,isFS,npml)
integer,intent(in)        :: nx_pml,ny_pml,nz_pml,fd_order,npml
logical,intent(in)        :: isFS
real,intent(out)          :: tau(:,:,:)
real,intent(in)           :: u(:,:,:),v(:,:,:),w(:,:,:),lambda(:,:,:)
real,intent(in)           :: dtdx

if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2+1
endif


if (fd_order==22)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=2,nx_pml
      do iy=2,ny_pml
         do iz=sp,nz_pml
            alpha=lambda(iz,iy,ix)*dtdx
            tau(iz,iy,ix)=alpha*S21*(u(iz,iy,ix)-u(iz,iy,ix-1)+v(iz,iy,ix)-v(iz,iy-1,ix)+w(iz,iy,ix)-w(iz-1,iy,ix))
        enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=3,nx_pml-1
      do iy=3,ny_pml-1
         do iz=sp,nz_pml-1
            alpha=lambda(iz,iy,ix)*dtdx
            tau(iz,iy,ix)=alpha*(S41*(u(iz,iy,ix)-u(iz,iy,ix-1)+v(iz,iy,ix)-v(iz,iy-1,ix)+w(iz,iy,ix)-w(iz-1,iy,ix))+S42*(u(iz,iy,ix+1)-u(iz,iy,ix-2)+v(iz,iy+1,ix)-v(iz,iy-2,ix)+w(iz+1,iy,ix)-w(iz-2,iy,ix)))         
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=4,nx_pml-2
      do iy=4,ny_pml-2
         do iz=sp,nz_pml-2
            alpha=lambda(iz,iy,ix)*dtdx
            tau(iz,iy,ix)=alpha*(S61*(u(iz,iy,ix)-u(iz,iy,ix-1)+v(iz,iy,ix)-v(iz,iy-1,ix)+w(iz,iy,ix)-w(iz-1,iy,ix))+S62*(u(iz,iy,ix+1)-u(iz,iy,ix-2)+v(iz,iy+1,ix)-v(iz,iy-2,ix)+w(iz+1,iy,ix)-w(iz-2,iy,ix))+S63*(u(iz,iy,ix+2)-u(iz,iy,ix-3)+v(iz,iy+2,ix)-v(iz,iy-3,ix)+w(iz+2,iy,ix)-w(iz-3,iy,ix)))         
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=5,nx_pml-3
      do iy=5,ny_pml-3
         do iz=sp,nz_pml-3
            alpha=lambda(iz,iy,ix)*dtdx
            tau(iz,iy,ix)=alpha*(S81*(u(iz,iy,ix)-u(iz,iy,ix-1)+v(iz,iy,ix)-v(iz,iy-1,ix)+w(iz,iy,ix)-w(iz-1,iy,ix))+S82*(u(iz,iy,ix+1)-u(iz,iy,ix-2)+v(iz,iy+1,ix)-v(iz,iy-2,ix)+w(iz+1,iy,ix)-w(iz-2,iy,ix))+S83*(u(iz,iy,ix+2)-u(iz,iy,ix-3)+v(iz,iy+2,ix)-v(iz,iy-3,ix)+w(iz+2,iy,ix)-w(iz-3,iy,ix))+S84*(u(iz,iy,ix+3)-u(iz,iy,ix-4)+v(iz,iy+3,ix)-v(iz,iy-4,ix)+w(iz+3,iy,ix)-w(iz-4,iy,ix)))         
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif       
end subroutine e3d_abc_kernal_totalstress_whole

!=================================================================================
subroutine e3d_abc_kernal_tauxx_whole(total,fd_order,u,tau_xx,temp,mu,dtdx,nx_pml,ny_pml,nz_pml,isFS,npml)
integer,intent(in)    :: fd_order,nx_pml,ny_pml,nz_pml,npml
logical,intent(in)    :: isFS
real,intent(in)       :: dtdx
real,intent(in)       :: total(:,:,:),temp(:,:,:),u(:,:,:),mu(:,:,:)
real,intent(inout)    :: tau_xx(:,:,:)
if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2+1
endif

if (fd_order==22)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=2,nx_pml
      do iy=2,ny_pml
         do iz=sp,nz_pml
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_xx(iz,iy,ix)=temp(iz,iy,ix)*tau_xx(iz,iy,ix)+total(iz,iy,ix)+alpha*S21*(u(iz,iy,ix)-u(iz,iy,ix-1))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=3,nx_pml-1
      do iy=3,ny_pml-1
         do iz=sp,nz_pml-1
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_xx(iz,iy,ix)=temp(iz,iy,ix)*tau_xx(iz,iy,ix)+total(iz,iy,ix)+alpha*(S41*(u(iz,iy,ix)-u(iz,iy,ix-1))+S42*(u(iz,iy,ix+1)-u(iz,iy,ix-2)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=4,nx_pml-2
      do iy=4,ny_pml-2
         do iz=sp,nz_pml-2
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_xx(iz,iy,ix)=temp(iz,iy,ix)*tau_xx(iz,iy,ix)+total(iz,iy,ix)+alpha*(S61*(u(iz,iy,ix)-u(iz,iy,ix-1))+S62*(u(iz,iy,ix+1)-u(iz,iy,ix-2))+S63*(u(iz,iy,ix+2)-u(iz,iy,ix-3)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=5,nx_pml-3
      do iy=5,ny_pml-3
         do iz=sp,nz_pml-3
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_xx(iz,iy,ix)=temp(iz,iy,ix)*tau_xx(iz,iy,ix)+total(iz,iy,ix)+alpha*(S81*(u(iz,iy,ix)-u(iz,iy,ix-1))+S82*(u(iz,iy,ix+1)-u(iz,iy,ix-2))+S83*(u(iz,iy,ix+2)-u(iz,iy,ix-3))+S84*(u(iz,iy,ix+3)-u(iz,iy,ix-4)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif   
end subroutine e3d_abc_kernal_tauxx_whole
!=================================================================================   
subroutine e3d_abc_kernal_tauyy_whole(fd_order,total,v,tau_yy,temp,mu,dtdx,nx_pml,ny_pml,nz_pml,isFS,npml)
integer,intent(in)    :: fd_order,nx_pml,ny_pml,nz_pml,npml
logical,intent(in)    :: isFS
real,intent(in)       :: dtdx
real,intent(in)       :: total(:,:,:),temp(:,:,:),v(:,:,:),mu(:,:,:)
real,intent(inout)    :: tau_yy(:,:,:)

if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2+1
endif
if (fd_order==22)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=2,nx_pml
      do iy=2,ny_pml
         do iz=sp,nz_pml
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_yy(iz,iy,ix)=temp(iz,iy,ix)*tau_yy(iz,iy,ix)+total(iz,iy,ix)+alpha*S21*(v(iz,iy,ix)-v(iz,iy-1,ix))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=3,nx_pml-1
      do iy=3,ny_pml-1
         do iz=sp,nz_pml-1
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_yy(iz,iy,ix)=temp(iz,iy,ix)*tau_yy(iz,iy,ix)+total(iz,iy,ix)+alpha*(S41*(v(iz,iy,ix)-v(iz,iy-1,ix))+S42*(v(iz,iy+1,ix)-v(iz,iy-2,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=4,nx_pml-2
      do iy=4,ny_pml-2
         do iz=sp,nz_pml-2
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_yy(iz,iy,ix)=temp(iz,iy,ix)*tau_yy(iz,iy,ix)+total(iz,iy,ix)+alpha*(S61*(v(iz,iy,ix)-v(iz,iy-1,ix))+S62*(v(iz,iy+1,ix)-v(iz,iy-2,ix))+S63*(v(iz,iy+2,ix)-v(iz,iy-3,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
   do ix=5,nx_pml-3
      do iy=5,ny_pml-3
         do iz=sp,nz_pml-3
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_yy(iz,iy,ix)=temp(iz,iy,ix)*tau_yy(iz,iy,ix)+total(iz,iy,ix)+alpha*(S81*(v(iz,iy,ix)-v(iz,iy-1,ix))+S82*(v(iz,iy+1,ix)-v(iz,iy-2,ix))+S83*(v(iz,iy+2,ix)-v(iz,iy-3,ix))+S84*(v(iz,iy+3,ix)-v(iz,iy-4,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif   
end subroutine e3d_abc_kernal_tauyy_whole
!=================================================================================
subroutine e3d_abc_kernal_tauzz_whole(nx_pml,ny_pml,nz_pml,total,w,tau_zz,temp,mu,dtdx,fd_order,isFS,npml)
integer,intent(in)    :: fd_order,nx_pml,ny_pml,nz_pml,npml
logical,intent(in)    :: isFS
real,intent(in)       :: dtdx
real,intent(in)       :: total(:,:,:),temp(:,:,:),w(:,:,:),mu(:,:,:)
real,intent(inout)    :: tau_zz(:,:,:)

if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2+1
endif

if (fd_order==22)   then
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha)
   do ix=2,nx_pml
      do iy=2,ny_pml
         do iz=sp,nz_pml
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_zz(iz,iy,ix)=temp(iz,iy,ix)*tau_zz(iz,iy,ix)+total(iz,iy,ix)+alpha*S21*(w(iz,iy,ix)-w(iz-1,iy,ix))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24)   then
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha)
   do ix=3,nx_pml-1
      do iy=3,ny_pml-1
         do iz=sp,nz_pml-1
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_zz(iz,iy,ix)=temp(iz,iy,ix)*tau_zz(iz,iy,ix)+total(iz,iy,ix)+alpha*(S41*(w(iz,iy,ix)-w(iz-1,iy,ix))+S42*(w(iz+1,iy,ix)-w(iz-2,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26)   then
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha)
   do ix=4,nx_pml-2
      do iy=4,ny_pml-2
         do iz=sp,nz_pml-2
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_zz(iz,iy,ix)=temp(iz,iy,ix)*tau_zz(iz,iy,ix)+total(iz,iy,ix)+alpha*(S61*(w(iz,iy,ix)-w(iz-1,iy,ix))+S62*(w(iz+1,iy,ix)-w(iz-2,iy,ix))+S63*(w(iz+2,iy,ix)-w(iz-3,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28)   then
   !$OMP PARALLEL DO PRIVATE(iz,iy,ix,alpha)
   do ix=5,nx_pml-3
      do iy=5,ny_pml-3
         do iz=sp,nz_pml-3
             alpha=2*mu(iz,iy,ix)*dtdx
             tau_zz(iz,iy,ix)=temp(iz,iy,ix)*tau_zz(iz,iy,ix)+total(iz,iy,ix)+alpha*(S81*(w(iz,iy,ix)-w(iz-1,iy,ix))+S82*(w(iz+1,iy,ix)-w(iz-2,iy,ix))+S83*(w(iz+2,iy,ix)-w(iz-3,iy,ix))+S84*(w(iz+3,iy,ix)-w(iz-4,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif   
end subroutine e3d_abc_kernal_tauzz_whole
!=================================================================================
subroutine e3d_abc_kernal_tauxy_whole(fd_order,nx_pml,tau_xy,u,v,temp,mu,dtdx,ny_pml,nz_pml,isFS,npml)

integer,intent(in)     :: fd_order,nx_pml,ny_pml,nz_pml,npml
logical,intent(in)     :: isFS
real,intent(in)        :: temp(:,:,:),u(:,:,:),v(:,:,:),mu(:,:,:)
real,intent(inout)     :: tau_xy(:,:,:)
real,intent(in)        :: dtdx
real                   :: mu_temp


if (isFS) then
   sp=npml
else
   sp=1
endif
if (fd_order==22)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,mu_temp,alpha)
   do ix=1,nx_pml-1
      do iy=1,ny_pml-1
         do iz=sp,nz_pml
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy,ix+1)+1/mu(iz,iy+1,ix)+1/mu(iz,iy+1,ix+1))/4)
            alpha=dtdx*mu_temp
            tau_xy(iz,iy,ix)=temp(iz,iy,ix)*tau_xy(iz,iy,ix)+alpha*S21*(u(iz,iy+1,ix)-u(iz,iy,ix)+v(iz,iy,ix+1)-v(iz,iy,ix))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,mu_temp,alpha)
   do ix=2,nx_pml-2
      do iy=2,ny_pml-2
         do iz=sp,nz_pml
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy,ix+1)+1/mu(iz,iy+1,ix)+1/mu(iz,iy+1,ix+1))/4)
            alpha=dtdx*mu_temp
            tau_xy(iz,iy,ix)=temp(iz,iy,ix)*tau_xy(iz,iy,ix)+alpha*(S41*(u(iz,iy+1,ix)-u(iz,iy,ix)+v(iz,iy,ix+1)-v(iz,iy,ix))+S42*(u(iz,iy+2,ix)-u(iz,iy-1,ix)+v(iz,iy,ix+2)-v(iz,iy,ix-1)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,mu_temp,alpha)
   do ix=3,nx_pml-3
      do iy=3,ny_pml-3
         do iz=sp,nz_pml
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy,ix+1)+1/mu(iz,iy+1,ix)+1/mu(iz,iy+1,ix+1))/4)
            alpha=dtdx*mu_temp
            tau_xy(iz,iy,ix)=temp(iz,iy,ix)*tau_xy(iz,iy,ix)+alpha*(S61*(u(iz,iy+1,ix)-u(iz,iy,ix)+v(iz,iy,ix+1)-v(iz,iy,ix))+S62*(u(iz,iy+2,ix)-u(iz,iy-1,ix)+v(iz,iy,ix+2)-v(iz,iy,ix-1))+S63*(u(iz,iy+3,ix)-u(iz,iy-2,ix)+v(iz,iy,ix+3)-v(iz,iy,ix-2)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,mu_temp,alpha)
   do ix=4,nx_pml-4
      do iy=4,ny_pml-4
         do iz=sp,nz_pml
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy,ix+1)+1/mu(iz,iy+1,ix)+1/mu(iz,iy+1,ix+1))/4)
            alpha=dtdx*mu_temp
            tau_xy(iz,iy,ix)=temp(iz,iy,ix)*tau_xy(iz,iy,ix)+alpha*(S81*(u(iz,iy+1,ix)-u(iz,iy,ix)+v(iz,iy,ix+1)-v(iz,iy,ix))+S82*(u(iz,iy+2,ix)-u(iz,iy-1,ix)+v(iz,iy,ix+2)-v(iz,iy,ix-1))+S83*(u(iz,iy+3,ix)-u(iz,iy-2,ix)+v(iz,iy,ix+3)-v(iz,iy,ix-2))+S84*(u(iz,iy+4,ix)-u(iz,iy-3,ix)+v(iz,iy,ix+4)-v(iz,iy,ix-3)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif 

end subroutine e3d_abc_kernal_tauxy_whole
!==================================================================================
subroutine e3d_abc_kernal_tauxz_whole(fd_order,nx_pml,ny_pml,nz_pml,tau_xz,u,w,temp,mu,dtdx,isFS,npml)
integer,intent(in)     :: fd_order,nx_pml,ny_pml,nz_pml,npml
logical,intent(in)     :: isFS
real,intent(in)        :: temp(:,:,:),u(:,:,:),w(:,:,:),mu(:,:,:)
real,intent(inout)     :: tau_xz(:,:,:)
real,intent(in)        :: dtdx
real                   :: mu_temp

integer                :: OMP_GET_NUM_THREADS
integer                :: yaya

if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2
endif
if (fd_order==22)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha,mu_temp)
   do ix=1,nx_pml-1
      do iy=1,ny_pml
         do iz=sp,nz_pml-1
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy,ix+1)+1/mu(iz+1,iy,ix)+1/mu(iz+1,iy,ix+1))/4)
            alpha=dtdx*mu_temp
            tau_xz(iz,iy,ix)=temp(iz,iy,ix)*tau_xz(iz,iy,ix)+alpha*S21*(u(iz+1,iy,ix)-u(iz,iy,ix)+w(iz,iy,ix+1)-w(iz,iy,ix))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO 
elseif (fd_order==24)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha,mu_temp)
   do ix=2,nx_pml-2
      do iy=1,ny_pml
         do iz=sp,nz_pml-2
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy,ix+1)+1/mu(iz+1,iy,ix)+1/mu(iz+1,iy,ix+1))/4)
            alpha=dtdx*mu_temp
            tau_xz(iz,iy,ix)=temp(iz,iy,ix)*tau_xz(iz,iy,ix)+alpha*(S41*(u(iz+1,iy,ix)-u(iz,iy,ix)+w(iz,iy,ix+1)-w(iz,iy,ix))+S42*(u(iz+2,iy,ix)-u(iz-1,iy,ix)+w(iz,iy,ix+2)-w(iz,iy,ix-1)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO 
elseif (fd_order==26)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha,mu_temp)
   do ix=3,nx_pml-3
      do iy=1,ny_pml
         do iz=sp,nz_pml-3
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy,ix+1)+1/mu(iz+1,iy,ix)+1/mu(iz+1,iy,ix+1))/4)
            alpha=dtdx*mu_temp
            tau_xz(iz,iy,ix)=temp(iz,iy,ix)*tau_xz(iz,iy,ix)+alpha*(S61*(u(iz+1,iy,ix)-u(iz,iy,ix)+w(iz,iy,ix+1)-w(iz,iy,ix))+S62*(u(iz+2,iy,ix)-u(iz-1,iy,ix)+w(iz,iy,ix+2)-w(iz,iy,ix-1))+S63*(u(iz+3,iy,ix)-u(iz-2,iy,ix)+w(iz,iy,ix+3)-w(iz,iy,ix-2)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO 
elseif (fd_order==28)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha,mu_temp) 
 
   do ix=4,nx_pml-4
      do iy=1,ny_pml
         do iz=sp,nz_pml-4
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy,ix+1)+1/mu(iz+1,iy,ix)+1/mu(iz+1,iy,ix+1))/4)
            alpha=dtdx*mu_temp
            tau_xz(iz,iy,ix)=temp(iz,iy,ix)*tau_xz(iz,iy,ix)+alpha*(S81*(u(iz+1,iy,ix)-u(iz,iy,ix)+w(iz,iy,ix+1)-w(iz,iy,ix))+S82*(u(iz+2,iy,ix)-u(iz-1,iy,ix)+w(iz,iy,ix+2)-w(iz,iy,ix-1))+S83*(u(iz+3,iy,ix)-u(iz-2,iy,ix)+w(iz,iy,ix+3)-w(iz,iy,ix-2))+S84*(u(iz+4,iy,ix)-u(iz-3,iy,ix)+w(iz,iy,ix+4)-w(iz,iy,ix-3)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO 
endif 

end subroutine e3d_abc_kernal_tauxz_whole
!==================================================================================

subroutine e3d_abc_kernal_tauyz_whole(fd_order,nx_pml,ny_pml,tau_yz,nz_pml,v,w,temp,mu,dtdx,isFS,npml)
integer,intent(in)     :: fd_order,nx_pml,ny_pml,nz_pml,npml
logical,intent(in)     :: isFS
real,intent(in)        :: temp(:,:,:),w(:,:,:),v(:,:,:),mu(:,:,:)
real,intent(inout)     :: tau_yz(:,:,:)
real,intent(in)        :: dtdx
real                   :: mu_temp

if (isFS) then
   sp=npml
else
   sp=(fd_order-20)/2
endif
if (fd_order==22)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha,mu_temp)
   do ix=1,nx_pml
      do iy=1,ny_pml-1
         do iz=sp,nz_pml-1
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy+1,ix)+1/mu(iz+1,iy,ix)+1/mu(iz+1,iy+1,ix))/4)
            alpha=dtdx*mu_temp
            tau_yz(iz,iy,ix)=temp(iz,iy,ix)*tau_yz(iz,iy,ix)+alpha*S21*(w(iz,iy+1,ix)-w(iz,iy,ix)+v(iz+1,iy,ix)-v(iz,iy,ix))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==24)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha,mu_temp)
   do ix=1,nx_pml
      do iy=2,ny_pml-2
         do iz=sp,nz_pml-2
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy+1,ix)+1/mu(iz+1,iy,ix)+1/mu(iz+1,iy+1,ix))/4)
            alpha=dtdx*mu_temp
            tau_yz(iz,iy,ix)=temp(iz,iy,ix)*tau_yz(iz,iy,ix)+alpha*(S41*(w(iz,iy+1,ix)-w(iz,iy,ix)+v(iz+1,iy,ix)-v(iz,iy,ix))+S42*(w(iz,iy+2,ix)-w(iz,iy-1,ix)+v(iz+2,iy,ix)-v(iz-1,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==26)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha,mu_temp)
   do ix=1,nx_pml
      do iy=3,ny_pml-3
         do iz=sp,nz_pml-3
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy+1,ix)+1/mu(iz+1,iy,ix)+1/mu(iz+1,iy+1,ix))/4)
            alpha=dtdx*mu_temp
            tau_yz(iz,iy,ix)=temp(iz,iy,ix)*tau_yz(iz,iy,ix)+alpha*(S61*(w(iz,iy+1,ix)-w(iz,iy,ix)+v(iz+1,iy,ix)-v(iz,iy,ix))+S62*(w(iz,iy+2,ix)-w(iz,iy-1,ix)+v(iz+2,iy,ix)-v(iz-1,iy,ix))+S63*(w(iz,iy+3,ix)-w(iz,iy-2,ix)+v(iz+3,iy,ix)-v(iz-2,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
elseif (fd_order==28)   then
   !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha,mu_temp)
   do ix=1,nx_pml
      do iy=4,ny_pml-4
         do iz=sp,nz_pml-4
            mu_temp=1/((1/mu(iz,iy,ix)+1/mu(iz,iy+1,ix)+1/mu(iz+1,iy,ix)+1/mu(iz+1,iy+1,ix))/4)
            alpha=dtdx*mu_temp
            tau_yz(iz,iy,ix)=temp(iz,iy,ix)*tau_yz(iz,iy,ix)+alpha*(S81*(w(iz,iy+1,ix)-w(iz,iy,ix)+v(iz+1,iy,ix)-v(iz,iy,ix))+S82*(w(iz,iy+2,ix)-w(iz,iy-1,ix)+v(iz+2,iy,ix)-v(iz-1,iy,ix))+S83*(w(iz,iy+3,ix)-w(iz,iy-2,ix)+v(iz+3,iy,ix)-v(iz-2,iy,ix))+S84*(w(iz,iy+4,ix)-w(iz,iy-3,ix)+v(iz+4,iy,ix)-v(iz-3,iy,ix)))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
endif 

end subroutine e3d_abc_kernal_tauyz_whole
!================================================================================




























end module module_fd3d
