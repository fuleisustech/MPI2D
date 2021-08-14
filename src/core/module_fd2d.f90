module module_fd2d

use module_global
use module_boundary2d
use module_array   ! modified by Bowen Guo
implicit none

integer, private :: iz,iy,ix,sp ! sp is used in the elastic 3D modeling starting point in the z direction

interface a2d_abc_kernal
   module procedure a2d_abc_kernal_whole
   module procedure a2d_abc_kernal_working
   module procedure a2d_abc_kernal_whole_Lax_Wendroff     ! modified by BOwen Guo  
end interface a2d_abc_kernal

interface a2d_abc_kernal_SG_p                             ! modified by Bowen Guo
   module procedure a2d_abc_kernal_SG_p_working           ! modified by Bowen Guo
end interface a2d_abc_kernal_SG_p                         ! modified by Bowen Guo

interface a2d_abc_kernal_SG_vel                           ! modified by Bowen Guo
   module procedure a2d_abc_kernal_SG_vel_working         ! modified by Bowen Guo
end interface a2d_abc_kernal_SG_vel                       ! modified by Bowen Guo

interface e2d_abc_kernal_vel                              ! modified by Bowen Guo
   module procedure e2d_abc_kernal_vel_working            ! modified by Bowen Guo
end interface e2d_abc_kernal_vel                          ! modified by Bowen Guo

interface e2d_abc_kernal_stress                           ! modified by Bowen Guo
   module procedure e2d_abc_kernal_stress_working         ! modified by Bowen Guo
end interface e2d_abc_kernal_stress                       ! modified by Bowen Guo

interface a3d_abc_kernal
   module procedure a3d_abc_kernal_whole
   module procedure a3d_abc_kernal_working
end interface a3d_abc_kernal

interface a2d_sponge_kernal
   module procedure a2d_sponge_kernal_whole
   module procedure a2d_sponge_kernal_working
end interface a2d_sponge_kernal

contains

!=====================================================================

subroutine a2d_abc_kernal_whole(p,p0,p1,temp1,temp2,&
   alpha,nz_pml,nx_pml,fd_order)
integer, intent(in) :: nz_pml,nx_pml,fd_order
real, intent(inout) :: p(:,:)
real, intent(in)    :: p0(:,:),p1(:,:),temp1(:,:),temp2(:,:),alpha(:,:)

!write(*,*) "size of temp1 = ", size(temp1,1), size(temp1,2)
!write(*,*) "size of temp2 = ", size(temp2,1), size(temp2,2)
!write(*,*) "size of alpha = ", size(alpha,1), size(alpha,2)
!write(*,*) "nz_pml = ",nz_pml
!write(*,*) "nx_pml = ",nx_pml
if (fd_order==22) then
   forall (iz=2:nz_pml-1,ix=2:nx_pml-1)
      p(iz,ix)=temp1(iz,ix)*p1(iz,ix)-temp2(iz,ix)*p0(iz,ix)+alpha(iz,ix)*&
                (C22*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix)))
   end forall
elseif (fd_order==24) then
   forall (iz=3:nz_pml-2,ix=3:nx_pml-2)
      p(iz,ix)=temp1(iz,ix)*p1(iz,ix)-temp2(iz,ix)*p0(iz,ix)+alpha(iz,ix)*&
                 (C42*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
                  C43*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix)))
   end forall
elseif (fd_order==26) then
   forall (iz=4:nz_pml-3,ix=4:nx_pml-3)
      p(iz,ix)=temp1(iz,ix)*p1(iz,ix)-temp2(iz,ix)*p0(iz,ix)+alpha(iz,ix)*&
                (C62*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
                 C63*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix))+&
                 C64*(p1(iz,ix+3)+p1(iz,ix-3)+p1(iz+3,ix)+p1(iz-3,ix)))
   end forall
elseif (fd_order==28) then
   forall (iz=5:nz_pml-4,ix=5:nx_pml-4)
      p(iz,ix)=temp1(iz,ix)*p1(iz,ix)-temp2(iz,ix)*p0(iz,ix)+alpha(iz,ix)*&
         (C82*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+ &
          C83*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix))+ &
          C84*(p1(iz,ix+3)+p1(iz,ix-3)+p1(iz+3,ix)+p1(iz-3,ix))+ &
          C85*(p1(iz,ix+4)+p1(iz,ix-4)+p1(iz+4,ix)+p1(iz-4,ix)))
   end forall
endif
end subroutine a2d_abc_kernal_whole

!=====================================================================

subroutine a2d_abc_kernal_working(p,p0,p1,temp1,temp2,&
   alpha,npml,nzp,nxp,fd_order)
integer, intent(in) :: npml,nzp,nxp,fd_order
real, intent(inout) :: p(:,:)
real, intent(in)    :: p0(:,:),p1(:,:),temp1(:,:),temp2(:,:),alpha(:,:)
if (fd_order==22) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      p(iz,ix)=temp1(iz,ix)*p1(iz,ix)-temp2(iz,ix)*p0(iz,ix)+alpha(iz,ix)*&
                (C22*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix)))
   end forall
elseif (fd_order==24) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      p(iz,ix)=temp1(iz,ix)*p1(iz,ix)-temp2(iz,ix)*p0(iz,ix)+alpha(iz,ix)*&
                 (C42*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
                  C43*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix)))
   end forall
elseif (fd_order==26) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      p(iz,ix)=temp1(iz,ix)*p1(iz,ix)-temp2(iz,ix)*p0(iz,ix)+alpha(iz,ix)*&
                (C62*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
                 C63*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix))+&
                 C64*(p1(iz,ix+3)+p1(iz,ix-3)+p1(iz+3,ix)+p1(iz-3,ix)))
   end forall
elseif (fd_order==28) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      p(iz,ix)=temp1(iz,ix)*p1(iz,ix)-temp2(iz,ix)*p0(iz,ix)+alpha(iz,ix)*&
         (C82*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+ &
          C83*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix))+ &
          C84*(p1(iz,ix+3)+p1(iz,ix-3)+p1(iz+3,ix)+p1(iz-3,ix))+ &
          C85*(p1(iz,ix+4)+p1(iz,ix-4)+p1(iz+4,ix)+p1(iz-4,ix)))
   end forall
endif
end subroutine a2d_abc_kernal_working

!=====================================================================
! modified by Bowen Guo
subroutine a2d_abc_kernal_SG_p_working(u,w,p,temp,alpha1,nz_pml,nx_pml,fd_order)
integer,intent(in)   :: nz_pml,nx_pml,fd_order
real,intent(inout)   :: u(:,:),w(:,:),p(:,:)
real,intent(in)      :: temp(:,:),alpha1(:,:)

! update pressure 
if (fd_order==22) then
   forall (iz=2:nz_pml,ix=2:nx_pml)
   p(iz,ix)=temp(iz,ix)*p(iz,ix)-alpha1(iz,ix)*(S21*(u(iz,ix)-u(iz,ix-1)+w(iz,ix)&
            -w(iz-1,ix)))
   end forall

elseif (fd_order==24) then
   forall (iz=3:nz_pml-1,ix=3:nx_pml-1)
   p(iz,ix)=temp(iz,ix)*p(iz,ix)-alpha1(iz,ix)*(S41*(u(iz,ix)-u(iz,ix-1))+S42* &
            (u(iz,ix+1)-u(iz,ix-2))+S41*(w(iz,ix)-w(iz-1,ix))+S42*(w(iz+1 &
            ,ix)-w(iz-2,ix)))
   end forall
elseif (fd_order==26) then
   forall (iz=4:nz_pml-2,ix=4:nx_pml-2)
   p(iz,ix)=temp(iz,ix)*p(iz,ix)-alpha1(iz,ix)*(S61*(u(iz,ix)-u(iz,ix-1))+S62*(u( &
           iz,ix+1)-u(iz,ix-2))+S63*(u(iz,ix+2)-u(iz,ix-3))+S61*(w(iz,ix) &
            -w(iz-1,ix))+S62*(w(iz+1,ix)-w(iz-2,ix))+S63*(w(iz+2,ix)-w(iz- &
            3,ix)))
   end forall
elseif (fd_order==28) then
   forall (iz=5:nz_pml-3,ix=5:nx_pml-3)
   p(iz,ix)=temp(iz,ix)*p(iz,ix)-alpha1(iz,ix)*(S81*(u(iz,ix)-u(iz,ix-1))+S82*(u( &
            iz,ix+1)-u(iz,ix-2))+S83*(u(iz,ix+2)-u(iz,ix-3))+S84*(u(iz,ix+ &
            3)-u(iz,ix-4))+S81*(w(iz,ix)-w(iz-1,ix))+S82*(w(iz+1,ix)-w(iz- &
            2,ix))+S83*(w(iz+2,ix)-w(iz-3,ix))+S84*(w(iz+3,ix)-w(iz-4,ix)))
   end forall
endif
end subroutine a2d_abc_kernal_SG_p_working
!=====================================================================================  
! modified by Bowen Guo
subroutine a2d_abc_kernal_SG_vel_working(u,w,p,temp,alpha2,nz_pml,nx_pml,fd_order)
integer,intent(in)   :: nz_pml,nx_pml,fd_order
real,intent(inout)   :: u(:,:),w(:,:),p(:,:)
real,intent(in)      :: temp(:,:),alpha2(:,:)
! update particle velocity u and w
if (fd_order==22) then
   forall (iz=1:nz_pml-1,ix=1:nx_pml-1)
   u(iz,ix)=temp(iz,ix)*u(iz,ix)-alpha2(iz,ix)*(S21*(p(iz,ix+1)-p(iz,ix)))
   w(iz,ix)=temp(iz,ix)*u(iz,ix)-alpha2(iz,ix)*(S21*(p(iz+1,ix)-p(iz,ix)))
   end forall
elseif (fd_order==24) then
   forall  (iz=2:nz_pml-2,ix=2:nx_pml-2)
   u(iz,ix)=temp(iz,ix)*u(iz,ix)-alpha2(iz,ix)*(S41*(p(iz,ix+1)-p(iz,ix))+S42*(p(iz&
            ,ix+2)-p(iz,ix-1)))
   w(iz,ix)=temp(iz,ix)*w(iz,ix)-alpha2(iz,ix)*(S41*(p(iz+1,ix)-p(iz,ix))+S42*(p(iz&
            +2,ix)-p(iz-1,ix)))
   end forall
elseif (fd_order==26) then
   forall (iz=3:nz_pml-3,ix=3:nx_pml-3)
   u(iz,ix)=temp(iz,ix)*u(iz,ix)-alpha2(iz,ix)*(S61*(p(iz,ix+1)-p(iz,ix))+S62*(p(iz&
            ,ix+2)-p(iz,ix-1))+S63*(p(iz,ix+3)-p(iz,ix-2)))
   w(iz,ix)=temp(iz,ix)*w(iz,ix)-alpha2(iz,ix)*(S61*(p(iz+1,ix)-p(iz,ix))+S62*(p(iz&
            +2,ix)-p(iz-1,ix))+S63*(p(iz+3,ix)-p(iz-2,ix)))
   end forall
elseif (fd_order==28) then
   forall (iz=4:nz_pml-4,ix=4:nx_pml-4)
   u(iz,ix)=temp(iz,ix)*u(iz,ix)-alpha2(iz,ix)*(S81*(p(iz,ix+1)-p(iz,ix))+S82*(p(iz&
            ,ix+2)-p(iz,ix-1))+S83*(p(iz,ix+3)-p(iz,ix-2))+S84*(p(iz,ix+4)-p&
            (iz,ix-3)))
   w(iz,ix)=temp(iz,ix)*w(iz,ix)-alpha2(iz,ix)*(S81*(p(iz+1,ix)-p(iz,ix))+S82*(p(iz&
            +2,ix)-p(iz-1,ix))+S83*(p(iz+3,ix)-p(iz-2,ix))+S84*(p(iz+4,ix)-p&
            (iz-3,ix)))
   end forall
endif

end subroutine a2d_abc_kernal_SG_vel_working
!========================================================================================
subroutine a2d_sponge_kernal_whole(p,p0,p1,u,c1t2,beta_dtdx,&
   spgx,spgz,spcx1,spcx2,spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
integer, intent(in) :: npml,nz_pml,nx_pml,fd_order
real, intent(inout) :: p(:,:),p1(:,:),p0(:,:),u(:,:)
real, intent(in)    :: c1t2,beta_dtdx(:,:),spgx(:),spgz(:),&
   spcx1(:),spcx2(:),spcz1(:),spcz2(:)
u=0.0
if (fd_order==22) then
   forall (iz=2:nz_pml-1,ix=2:nx_pml-1)
      u(iz,ix)=beta_dtdx(iz,ix) * (c1t2 * p1(iz,ix) + &
                C22*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix)))
   end forall
elseif (fd_order==24) then
   forall (iz=3:nz_pml-2,ix=3:nx_pml-2)
      u(iz,ix)=beta_dtdx(iz,ix) * (c1t2 * p1(iz,ix) + &
                C42*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
                C43*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix)))
   end forall
elseif (fd_order==26) then
   forall (iz=4:nz_pml-3,ix=4:nx_pml-3)
      u(iz,ix)=beta_dtdx(iz,ix) * (c1t2 * p1(iz,ix) + &
                C62*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
                C63*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix))+&
                C64*(p1(iz,ix+3)+p1(iz,ix-3)+p1(iz+3,ix)+p1(iz-3,ix)))
   end forall
elseif (fd_order==28) then
   forall (iz=5:nz_pml-4,ix=5:nx_pml-4)
      u(iz,ix)=beta_dtdx(iz,ix) * (c1t2 * p1(iz,ix) + &
                C82*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
                C83*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix))+&
                C84*(p1(iz,ix+3)+p1(iz,ix-3)+p1(iz+3,ix)+p1(iz-3,ix))+&
                C85*(p1(iz,ix+4)+p1(iz,ix-4)+p1(iz+4,ix)+p1(iz-4,ix)))
   end forall
endif
call sponge2D(p1,spgx,spgz,spcx1,spcx2,spcz1,spcz2,nx_pml,nz_pml,npml)
call sponge2D(p0,spgx,spgz,spcx1,spcx2,spcz1,spcz2,nx_pml,nz_pml,npml)
forall (iz=1:nz_pml,ix=1:nx_pml)
   p(iz,ix)=2.0*p1(iz,ix)-p0(iz,ix)+u(iz,ix)
end forall
end subroutine a2d_sponge_kernal_whole

!=====================================================================

subroutine a2d_sponge_kernal_working(p,p0,p1,u,c1t2,beta_dtdx,npml,nzp,nxp,fd_order)
integer, intent(in) :: npml,nzp,nxp,fd_order
real, intent(inout) :: p(:,:),p0(:,:),p1(:,:),u(:,:)
real, intent(in)    :: c1t2,beta_dtdx(:,:)
if (fd_order==22) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      u(iz,ix)=beta_dtdx(iz,ix) * (c1t2 * p1(iz,ix) + &
               C22*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix)))
   end forall
elseif (fd_order==24) then
   forall (iz=npml+1:nzp,ix=npml:nxp)
      u(iz,ix)=beta_dtdx(iz,ix) * (c1t2 * p1(iz,ix) + &
               C42*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
               C43*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix)))
   end forall
elseif (fd_order==26) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      u(iz,ix)=beta_dtdx(iz,ix) * (c1t2 * p1(iz,ix) + &
               C62*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
               C63*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix))+&
               C64*(p1(iz,ix+3)+p1(iz,ix-3)+p1(iz+3,ix)+p1(iz-3,ix)))
   end forall
elseif (fd_order==28) then
   forall (iz=npml+1:nzp,ix=npml+1:nxp)
      u(iz,ix)=beta_dtdx(iz,ix) * (c1t2 * p1(iz,ix) + &
               C82*(p1(iz,ix+1)+p1(iz,ix-1)+p1(iz+1,ix)+p1(iz-1,ix))+&
               C83*(p1(iz,ix+2)+p1(iz,ix-2)+p1(iz+2,ix)+p1(iz-2,ix))+&
               C84*(p1(iz,ix+3)+p1(iz,ix-3)+p1(iz+3,ix)+p1(iz-3,ix))+&
               C85*(p1(iz,ix+4)+p1(iz,ix-4)+p1(iz+4,ix)+p1(iz-4,ix)))
   end forall
endif
forall (iz=npml+1:nzp,ix=npml+1:nxp)
   p(iz,ix)=2.0*p1(iz,ix)-p0(iz,ix)+u(iz,ix)
end forall
end subroutine a2d_sponge_kernal_working

!=====================================================================

subroutine a3d_abc_kernal_whole(p,p0,p1,temp1,temp2,temp3,&
   alpha,nz_pml,ny_pml,nx_pml,fd_order)
integer, intent(in) :: nz_pml,ny_pml,nx_pml,fd_order
real, intent(inout) :: p(:,:,:)
real, intent(in)    :: p0(:,:,:),p1(:,:,:),temp1(:,:,:),temp2(:,:,:),&
   temp3(:,:,:),alpha(:,:,:)
if (fd_order==22) then
   forall (iz=2:nz_pml-1,iy=2:ny_pml-1,ix=2:nx_pml-1)
      p(iz,iy,ix)=(temp1(iz,iy,ix)*p1(iz,iy,ix) &
         -temp2(iz,iy,ix)*p0(iz,iy,ix)&
         +alpha(iz,iy,ix)*(C22*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                               +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                               +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))))/temp3(iz,iy,ix)
   end forall
elseif (fd_order==24) then
   forall (iz=3:nz_pml-2,iy=3:ny_pml-2,ix=3:nx_pml-2)
      p(iz,iy,ix)=(temp1(iz,iy,ix)*p1(iz,iy,ix) &
         -temp2(iz,iy,ix)*p0(iz,iy,ix) &
         +alpha(iz,iy,ix)*(C42*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                               +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                               +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                          +C43*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                               +p1(iz,iy+2,ix)+p1(iz,iy+2,ix) &
                               +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))))/temp3(iz,iy,ix)
   end forall
elseif (fd_order==26) then
   forall (iz=4:nz_pml-3,iy=4:ny_pml-3,ix=4:nx_pml-3)
      p(iz,iy,ix)=(temp1(iz,iy,ix)*p1(iz,iy,ix) &
         -temp2(iz,iy,ix)*p0(iz,iy,ix) &
         +alpha(iz,iy,ix)*(C62*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                               +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                               +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                          +C63*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                               +p1(iz,iy+2,ix)+p1(iz,iy-2,ix) &
                               +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))&
                          +C64*(p1(iz,iy,ix+3)+p1(iz,iy,ix-3) &
                               +p1(iz,iy+3,ix)+p1(iz,iy-3,ix) &
                               +p1(iz+3,iy,ix)+p1(iz-3,iy,ix))))/temp3(iz,iy,ix)
   end forall
elseif (fd_order==28) then
   forall (iz=5:nz_pml-4,iy=5:ny_pml-4,ix=5:nx_pml-4)
      p(iz,iy,ix)=(temp1(iz,iy,ix)*p1(iz,iy,ix) &
         -temp2(iz,iy,ix)*p0(iz,iy,ix) &
         +alpha(iz,iy,ix)*(C82*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
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
                               +p1(iz+4,iy,ix)+p1(iz-4,iy,ix))))/temp3(iz,iy,ix)
   end forall
endif
end subroutine a3d_abc_kernal_whole

!=====================================================================

subroutine a3d_abc_kernal_working(p,p0,p1,temp1,temp2,temp3,&
   alpha,npml,nzp,nyp,nxp,fd_order)
integer, intent(in) :: npml,nzp,nyp,nxp,fd_order
real, intent(inout) :: p(:,:,:)
real, intent(in)    :: p0(:,:,:),p1(:,:,:),temp1(:,:,:),temp2(:,:,:),&
   temp3(:,:,:),alpha(:,:,:)
if (fd_order==22) then
   forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
      p(iz,iy,ix)=(temp1(iz,iy,ix)*p1(iz,iy,ix) &
         -temp2(iz,iy,ix)*p0(iz,iy,ix)&
         +alpha(iz,iy,ix)*(C22*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                               +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                               +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))))/temp3(iz,iy,ix)
   end forall
elseif (fd_order==24) then
   forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
      p(iz,iy,ix)=(temp1(iz,iy,ix)*p1(iz,iy,ix) &
         -temp2(iz,iy,ix)*p0(iz,iy,ix) &
         +alpha(iz,iy,ix)*(C42*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                               +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                               +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                          +C43*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                               +p1(iz,iy+2,ix)+p1(iz,iy+2,ix) &
                               +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))))/temp3(iz,iy,ix)
   end forall
elseif (fd_order==26) then
   forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
      p(iz,iy,ix)=(temp1(iz,iy,ix)*p1(iz,iy,ix) &
         -temp2(iz,iy,ix)*p0(iz,iy,ix) &
         +alpha(iz,iy,ix)*(C62*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
                               +p1(iz,iy+1,ix)+p1(iz,iy-1,ix) &
                               +p1(iz+1,iy,ix)+p1(iz-1,iy,ix))&
                          +C63*(p1(iz,iy,ix+2)+p1(iz,iy,ix-2) &
                               +p1(iz,iy+2,ix)+p1(iz,iy-2,ix) &
                               +p1(iz+2,iy,ix)+p1(iz-2,iy,ix))&
                          +C64*(p1(iz,iy,ix+3)+p1(iz,iy,ix-3) &
                               +p1(iz,iy+3,ix)+p1(iz,iy-3,ix) &
                               +p1(iz+3,iy,ix)+p1(iz-3,iy,ix))))/temp3(iz,iy,ix)
   end forall
elseif (fd_order==28) then
   forall (iz=npml+1:nzp,iy=npml+1:nyp,ix=npml+1:nxp)
      p(iz,iy,ix)=(temp1(iz,iy,ix)*p1(iz,iy,ix) &
         -temp2(iz,iy,ix)*p0(iz,iy,ix) &
         +alpha(iz,iy,ix)*(C82*(p1(iz,iy,ix+1)+p1(iz,iy,ix-1) &
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
                               +p1(iz+4,iy,ix)+p1(iz-4,iy,ix))))/temp3(iz,iy,ix)
   end forall
endif
end subroutine a3d_abc_kernal_working

!=====================================================================
! modified by Bowen Guo
subroutine e2d_abc_kernal_vel_working(u,w,tau_xx,tau_xz,tau_zz,temp,alpha1,nz_pml,nx_pml,fd_order,isFS,npml)
logical,intent(in)         :: isFS
integer,intent(in)         :: nz_pml,nx_pml,npml,fd_order
real,intent(out)           :: u(:,:),w(:,:)
real,intent(in)            :: tau_xx(:,:),tau_xz(:,:),tau_zz(:,:),temp(:,:),alpha1(:,:)
! update particle velocity u and w 
if (isFS) then
   sp=npml
else 
   sp=(fd_order-20)/2+1
endif


if (fd_order==22) then
   forall (iz=sp:nz_pml,ix=2:nx_pml)
      u(iz,ix)=temp(iz,ix)*u(iz,ix)+alpha1(iz,ix)*(S21*(tau_xx(iz,ix)-tau_xx(iz,ix-1)+tau_xz(iz,ix)&
            -tau_xz(iz-1,ix)))
   endforall
   forall (iz=sp:nz_pml-1,ix=1:nx_pml-1)
      w(iz,ix)=temp(iz,ix)*w(iz,ix)+alpha1(iz,ix)*(S21*(tau_xz(iz,ix+1)-tau_xz(iz,ix)+tau_zz(iz+1,ix)&
            -tau_zz(iz,ix)))
   endforall
elseif (fd_order==24) then
   forall (iz=sp:nz_pml-1,ix=3:nx_pml-1)
      u(iz,ix)=temp(iz,ix)*u(iz,ix)+alpha1(iz,ix)*(S41*(tau_xx(iz,ix)-tau_xx(iz,ix-1))+S42* &
            (tau_xx(iz,ix+1)-tau_xx(iz,ix-2))+S41*(tau_xz(iz,ix)-tau_xz(iz-1,ix))+S42*(tau_xz(iz+1 &
            ,ix)-tau_xz(iz-2,ix)))
   endforall    
   forall (iz=sp:nz_pml-2,ix=2:nx_pml-2)  
      w(iz,ix)=temp(iz,ix)*w(iz,ix)+alpha1(iz,ix)*(S41*(tau_xz(iz,ix+1)-tau_xz(iz,ix))+S42* &
            (tau_xz(iz,ix+2)-tau_xz(iz,ix-1))+S41*(tau_zz(iz+1,ix)-tau_zz(iz,ix))+S42*(tau_zz(iz+2 &
            ,ix)-tau_zz(iz-1,ix)))
   end forall
elseif (fd_order==26) then
   forall (iz=sp:nz_pml-2,ix=4:nx_pml-2)
      u(iz,ix)=temp(iz,ix)*u(iz,ix)+alpha1(iz,ix)*(S61*(tau_xx(iz,ix)-tau_xx(iz,ix-1))+S62*(tau_xx( &
           iz,ix+1)-tau_xx(iz,ix-2))+S63*(tau_xx(iz,ix+2)-tau_xx(iz,ix-3))+S61*(tau_xz(iz,ix) &
            -tau_xz(iz-1,ix))+S62*(tau_xz(iz+1,ix)-tau_xz(iz-2,ix))+S63*(tau_xz(iz+2,ix)-tau_xz(iz- &
            3,ix)))
   endforall   
   forall (iz=sp:nz_pml-3,ix=3:nx_pml-3)
      w(iz,ix)=temp(iz,ix)*w(iz,ix)+alpha1(iz,ix)*(S61*(tau_xz(iz,ix+1)-tau_xz(iz,ix))+S62*(tau_xz( &
           iz,ix+2)-tau_xz(iz,ix-1))+S63*(tau_xz(iz,ix+3)-tau_xz(iz,ix-2))+S61*(tau_zz(iz+1,ix) &
            -tau_zz(iz,ix))+S62*(tau_zz(iz+2,ix)-tau_zz(iz-1,ix))+S63*(tau_zz(iz+3,ix)-tau_zz(iz- &
           2,ix)))
   end forall
elseif (fd_order==28) then
   forall (iz=sp:nz_pml-3,ix=5:nx_pml-3)
      u(iz,ix)=temp(iz,ix)*u(iz,ix)+alpha1(iz,ix)*(S81*(tau_xx(iz,ix)-tau_xx(iz,ix-1))+S82*(tau_xx( &
            iz,ix+1)-tau_xx(iz,ix-2))+S83*(tau_xx(iz,ix+2)-tau_xx(iz,ix-3))+S84*(tau_xx(iz,ix+ &
            3)-tau_xx(iz,ix-4))+S81*(tau_xz(iz,ix)-tau_xz(iz-1,ix))+S82*(tau_xz(iz+1,ix)-tau_xz(iz- &
            2,ix))+S83*(tau_xz(iz+2,ix)-tau_xz(iz-3,ix))+S84*(tau_xz(iz+3,ix)-tau_xz(iz-4,ix)))
   endforall   
   forall (iz=sp:nz_pml-4,ix=4:nx_pml-4)
      w(iz,ix)=temp(iz,ix)*w(iz,ix)+alpha1(iz,ix)*(S81*(tau_xz(iz,ix+1)-tau_xz(iz,ix))+S82*(tau_xz( &
            iz,ix+2)-tau_xz(iz,ix-1))+S83*(tau_xz(iz,ix+3)-tau_xz(iz,ix-2))+S84*(tau_xz(iz,ix+ &
            4)-tau_xz(iz,ix-3))+S81*(tau_zz(iz+1,ix)-tau_zz(iz,ix))+S82*(tau_zz(iz+2,ix)-tau_zz(iz- &
            1,ix))+S83*(tau_zz(iz+3,ix)-tau_zz(iz-2,ix))+S84*(tau_zz(iz+4,ix)-tau_zz(iz-3,ix)))
   end forall
endif
end subroutine e2d_abc_kernal_vel_working

!===================================================================================================
!modified by Bowen Guo
subroutine e2d_abc_kernal_stress_working(u,w,tau_xx,tau_xz,tau_zz,temp,alpha_lambda,alpha_mu,nz_pml,nx_pml,fd_order,isFS,npml)
logical,intent(in)      :: isFS
integer,intent(in)      :: nz_pml,nx_pml,npml,fd_order
real,intent(in)         :: u(:,:),w(:,:),temp(:,:),alpha_lambda(:,:),alpha_mu(:,:)
real,intent(inout)      :: tau_xx(:,:),tau_xz(:,:),tau_zz(:,:)
real,allocatable        :: TAU(:,:)

if (isFS) then
   sp = npml
else 
   sp = (fd_order-20)/2+1
endif
call allocate_and_initial(TAU,nz_pml,nx_pml)
! update normal stress 
if (fd_order==22) then
   forall (iz=sp:nz_pml,ix=1:nx_pml-1)
      TAU(iz,ix)=alpha_lambda(iz,ix)*(S21*(u(iz,ix+1)-u(iz,ix)+w(iz,ix)-w(iz-1,ix)))
      tau_xx(iz,ix)=temp(iz,ix)*tau_xx(iz,ix)+(TAU(iz,ix)+2*alpha_mu(iz,ix)*(S21*(u(iz,ix+1)-u(iz,ix))))
      tau_zz(iz,ix)=temp(iz,ix)*tau_zz(iz,ix)+(TAU(iz,ix)+2*alpha_mu(iz,ix)*(S21*(w(iz,ix)-w(iz-1,ix))))
   end forall 
elseif (fd_order==24) then
   forall (iz=sp:nz_pml-1,ix=2:nx_pml-2)
      TAU(iz,ix)=alpha_lambda(iz,ix)*(S41*(u(iz,ix+1)-u(iz,ix))+S42* &
            (u(iz,ix+2)-u(iz,ix-1))+S41*(w(iz,ix)-w(iz-1,ix))+S42*(w(iz+1 &
            ,ix)-w(iz-2,ix)))
      tau_xx(iz,ix)=temp(iz,ix)*tau_xx(iz,ix)+(TAU(iz,ix)+2*alpha_mu(iz,ix)*(S41*(u(iz,ix+1)-u(iz,ix))+S42*(u(iz,ix+2)-u(iz,ix-1))))
      tau_zz(iz,ix)=temp(iz,ix)*tau_zz(iz,ix)+(TAU(iz,ix)+2*alpha_mu(iz,ix)*(S41*(w(iz,ix)-w(iz-1,ix))+S42*(w(iz+1,ix)-w(iz-2,ix))))
   end forall
elseif (fd_order==26) then
   forall (iz=sp:nz_pml-2,ix=3:nx_pml-3)
      TAU(iz,ix)=alpha_lambda(iz,ix)*(S61*(u(iz,ix+1)-u(iz,ix))+S62*(u( &
           iz,ix+2)-u(iz,ix-1))+S63*(u(iz,ix+3)-u(iz,ix-2))+S61*(w(iz,ix) &
            -w(iz-1,ix))+S62*(w(iz+1,ix)-w(iz-2,ix))+S63*(w(iz+2,ix)-w(iz- &
            3,ix)))
      tau_xx(iz,ix)=temp(iz,ix)*tau_xx(iz,ix)+(TAU(iz,ix)+2*alpha_mu(iz,ix)*(S61*(u(iz,ix+1)-u(iz,ix))+S62*(u(iz,ix+2)-u(iz,ix-1))+&
           S63*(u(iz,ix+3)-u(iz,ix-2))))
      tau_zz(iz,ix)=temp(iz,ix)*tau_zz(iz,ix)+(TAU(iz,ix)+2*alpha_mu(iz,ix)*(S61*(w(iz,ix)-w(iz-1,ix))+S62*(w(iz+1,ix)-w(iz-2,ix))+&
           S63*(w(iz+2,ix)-w(iz-3,ix))))
    end forall
elseif (fd_order==28) then
   forall (iz=sp:nz_pml-3,ix=4:nx_pml-4)
      TAU(iz,ix)=alpha_lambda(iz,ix)*(S81*(u(iz,ix+1)-u(iz,ix))+S82*(u( &
            iz,ix+2)-u(iz,ix-1))+S83*(u(iz,ix+3)-u(iz,ix-2))+S84*(u(iz,ix+ &
            4)-u(iz,ix-3))+S81*(w(iz,ix)-w(iz-1,ix))+S82*(w(iz+1,ix)-w(iz- &
            2,ix))+S83*(w(iz+2,ix)-w(iz-3,ix))+S84*(w(iz+3,ix)-w(iz-4,ix)))
      tau_xx(iz,ix)=temp(iz,ix)*tau_xx(iz,ix)+(TAU(iz,ix)+2*alpha_mu(iz,ix)*(S81*(u(iz,ix+1)-u(iz,ix))+S82*(u(iz,ix+2)-u(iz,ix-1))+&
           S83*(u(iz,ix+3)-u(iz,ix-2))+S84*(u(iz,ix+4)-u(iz,ix-3))))
      tau_zz(iz,ix)=temp(iz,ix)*tau_zz(iz,ix)+(TAU(iz,ix)+2*alpha_mu(iz,ix)*(S81*(w(iz,ix)-w(iz-1,ix))+S82*(w(iz+1,ix)-w(iz-2,ix))+&
           S83*(w(iz+2,ix)-w(iz-3,ix))+S84*(w(iz+3,ix)-w(iz-4,ix))))
   end forall
endif
call deallocate_and_free(TAU)
! update shear stress
if (fd_order==22) then
   forall (iz=sp:nz_pml-1,ix=2:nx_pml)
      tau_xz(iz,ix)=temp(iz,ix)*tau_xz(iz,ix)+alpha_mu(iz,ix)*(S21*(u(iz+1,ix)-u(iz,ix)+w(iz,ix)-w(iz,ix-1)))
   end forall
elseif (fd_order==24) then
   forall (iz=sp:nz_pml-2,ix=3:nx_pml-1)
      tau_xz(iz,ix)=temp(iz,ix)*tau_xz(iz,ix)+alpha_mu(iz,ix)*(S41*(u(iz+1,ix)-u(iz,ix)+w(iz,ix)-w(iz,ix-1))+S42*(u(iz+2,ix)&
                    -u(iz-1,ix)+w(iz,ix+1)-w(iz,ix-2)))
   end forall
elseif (fd_order==26) then
   forall (iz=sp:nz_pml-3,ix=4:nx_pml-2)
      tau_xz(iz,ix)=temp(iz,ix)*tau_xz(iz,ix)+alpha_mu(iz,ix)*(S61*(u(iz+1,ix)-u(iz,ix)+w(iz,ix)-w(iz,ix-1))+S62*(u(iz+2,ix)&
                    -u(iz-1,ix)+w(iz,ix+1)-w(iz,ix-2))+S63*(u(iz+3,ix)-u(iz-2,ix)+w(iz,ix+2)-w(iz,ix-3)))
   end forall
elseif (fd_order==28) then
   forall (iz=sp:nz_pml-4,ix=5:nx_pml-3)
      tau_xz(iz,ix)=temp(iz,ix)*tau_xz(iz,ix)+alpha_mu(iz,ix)*(S81*(u(iz+1,ix)-u(iz,ix)+w(iz,ix)-w(iz,ix-1))+S82*(u(iz+2,ix)&
                    -u(iz-1,ix)+w(iz,ix+1)-w(iz,ix-2))+S83*(u(iz+3,ix)-u(iz-2,ix)+w(iz,ix+2)-w(iz,ix-3))+S84&
                     *(u(iz+4,ix)-u(iz-3,ix)+w(iz,ix+3)-w(iz,ix-4)))
   end forall
endif
      
end subroutine e2d_abc_kernal_stress_working
!============================================================================================\
! Written by Bowen Guo to calculate wave field by Lax-Wendroff method
subroutine a2d_abc_kernal_whole_Lax_Wendroff(p,p0,p1,temp1,temp2,temp3,temp4,c,dx,dz,nz_pml,nx_pml,fd_order)            
integer,intent(in)    :: nz_pml,nx_pml,fd_order
real,intent(in)       :: dx,dz
real,intent(in)       :: c(:,:),temp1(:,:),temp2(:,:),temp3(:,:),temp4(:,:),p0(:,:),p1(:,:)
real,intent(out)      :: p(:,:)
real,allocatable      :: d2p(:,:),d4p(:,:)
call allocate_and_initial(d2p,d4p,nz_pml,nx_pml)
call cal_space_laplace_2d_whole(p1,d2p,dx,dz,fd_order,nz_pml,nx_pml)
! Calculate the 4th time derivative
d2p=d2p*(c**2)
call cal_space_laplace_2d_whole(d2p,d4p,dx,dz,fd_order,nz_pml,nx_pml)
p=temp1*d2p+temp2*p0+temp3*p1+temp4*d4p
call deallocate_and_free(d2p,d4p)
end subroutine a2d_abc_kernal_whole_Lax_Wendroff


!=====================================================================================
subroutine cal_space_laplace_2d_whole(p,d2p,dx,dz,fd_order,nz_pml,nx_pml)
integer,intent(in)    :: nz_pml,nx_pml,fd_order
real,intent(in)       :: dx,dz
real,intent(in)       :: p(:,:)
real,intent(out)      :: d2p(:,:)
real,allocatable      :: d2p_dx2(:,:),d2p_dz2(:,:)

call allocate_and_initial(d2p_dx2,d2p_dz2,nz_pml,nx_pml)
call cal_space_2nddx_2d_whole(p,d2p_dx2,dx,fd_order,nz_pml,nx_pml)
call cal_space_2nddz_2d_whole(p,d2p_dz2,dz,fd_order,nz_pml,nx_pml)
d2p=d2p_dx2+d2p_dz2
call deallocate_and_free(d2p_dx2,d2p_dz2)
end subroutine cal_space_laplace_2d_whole

!======================================================================================================
subroutine cal_space_2nddx_2d_whole(p,d2p_dx2,dx,fd_order,nz_pml,nx_pml)
integer,intent(in)    :: fd_order,nz_pml,nx_pml
real,intent(in)       :: dx
real,intent(in)       :: p(:,:)
real,intent(out)      :: d2p_dx2(:,:)
real                  :: C21_temp,C22_temp,C41_temp,C42_temp,C43_temp,C61_temp,C62_temp,C63_temp,&
C64_temp,C65_temp,C81_temp,C82_temp,C83_temp,C84_temp,C85_temp


if (fd_order==42) then
   C21_temp=C21/(dx**2);C22_temp=C22/(dx**2)
   forall (iz=1:nz_pml,ix=2:nx_pml-1)
      d2p_dx2(iz,ix)=C21_temp*p(iz,ix)+C22_temp*(p(iz,ix+1)+p(iz,ix-1))
   end forall
elseif (fd_order==44) then 
   C41_temp=C41/(dx**2);C42_temp=C42/(dx**2);C43_temp=C43/(dx**2)
   forall (iz=1:nz_pml,ix=3:nx_pml-2)
      d2p_dx2(iz,ix)=C41_temp*p(iz,ix)+C42_temp*(p(iz,ix+1)+p(iz,ix-1))+C43_temp*(p(iz,ix+2)+p(iz,ix-2))
   end forall
elseif (fd_order==46) then 
   C61_temp=C61/(dx**2);C62_temp=C62/(dx**2);C63_temp=C63/(dx**2);C64_temp=C64/(dx**2)
   forall (iz=1:nz_pml,ix=4:nx_pml-3)
      d2p_dx2(iz,ix)=C61_temp*p(iz,ix)+C62_temp*(p(iz,ix+1)+p(iz,ix-1))+C63_temp*(p(iz,ix+2)+p(iz,ix-2))+&
      C64_temp*(p(iz,ix+3)+p(iz,ix-3))
   end forall
elseif (fd_order==48) then 
   C81_temp=C81/(dx**2);C82_temp=C82/(dx**2);C83_temp=C83/(dx**2);C84_temp=C84/(dx**2);C85_temp=C85/&
   (dx**2)
   forall (iz=1:nz_pml,ix=5:nx_pml-4)
      d2p_dx2(iz,ix)=C81_temp*p(iz,ix)+C82_temp*(p(iz,ix+1)+p(iz,ix-1))+C83_temp*(p(iz,ix+2)+p(iz,ix-2))+&
      C84_temp*(p(iz,ix+3)+p(iz,ix-3))+C85_temp*(p(iz,ix+4)+p(iz,ix-4))
   end forall
endif
end subroutine cal_space_2nddx_2d_whole

!=============================================================================================
subroutine cal_space_2nddz_2d_whole(p,d2p_dz2,dz,fd_order,nz_pml,nx_pml)
integer,intent(in)    :: fd_order,nz_pml,nx_pml
real,intent(in)       :: dz
real,intent(in)       :: p(:,:)
real,intent(out)      :: d2p_dz2(:,:)
real                  :: C21_temp,C22_temp,C41_temp,C42_temp,C43_temp,C61_temp,C62_temp,C63_temp,&
C64_temp,C65_temp,C81_temp,C82_temp,C83_temp,C84_temp,C85_temp
if (fd_order==42) then
   C21_temp=C21/(dz**2);C22_temp=C22/(dz**2)
   forall (iz=2:nz_pml-1,ix=1:nx_pml)
      d2p_dz2(iz,ix)=C21_temp*p(iz,ix)+C22_temp*(p(iz+1,ix)+p(iz-1,ix))
   end forall
elseif (fd_order==44) then
   C41_temp=C41/(dz**2);C42_temp=C42/(dz**2);C43_temp=C43/(dz**2)
   forall (iz=3:nz_pml-2,ix=1:nx_pml)
      d2p_dz2(iz,ix)=C41_temp*p(iz,ix)+C42_temp*(p(iz+1,ix)+p(iz-1,ix))+C43_temp*(p(iz+2,ix)+p(iz-2,ix))
   end forall
elseif (fd_order==46) then
   C61_temp=C61/(dz**2);C62_temp=C62/(dz**2);C63_temp=C63/(dz**2);C64_temp=C64/(dz**2)
   forall (iz=4:nz_pml-3,ix=1:nx_pml)
      d2p_dz2(iz,ix)=C61_temp*p(iz,ix)+C62_temp*(p(iz+1,ix)+p(iz-1,ix))+C63_temp*(p(iz+2,ix)+p(iz-2,ix))+&
      C64_temp*(p(iz+3,ix)+p(iz-3,ix))
   end forall
elseif (fd_order==48) then
   C81_temp=C81/(dz**2);C82_temp=C82/(dz**2);C83_temp=C83/(dz**2);C84_temp=C84/(dz**2);C85_temp=C85/&
   (dz**2)
   forall (iz=5:nz_pml-4,ix=1:nx_pml)
      d2p_dz2(iz,ix)=C81_temp*p(iz,ix)+C82_temp*(p(iz+1,ix)+p(iz-1,ix))+C83_temp*(p(iz+2,ix)+p(iz-2,ix))+&
      C84_temp*(p(iz+3,ix)+p(iz-3,ix))+C85_temp*(p(iz+4,ix)+p(iz-4,ix))
   end forall
endif
end subroutine cal_space_2nddz_2d_whole
end module module_fd2d
