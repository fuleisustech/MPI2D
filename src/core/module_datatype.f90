module module_datatype
use module_utility
implicit none

! Basic datatype
type fname
   character(len=200) :: vel_file, vels_file, den_file, den_homo_file, &   
      vel_homo_file, vel_water_file,vel_whomo_file, &   
      coord_in_file,coord_out_file,coord_all_file,&
      pw1_coord_file,pw2_coord_file,log_file,&
      pw1_par_file,pw2_par_file,pw1_coord_all_file,pw2_coord_all_file,&
      csg_out,csg_homo,csg_recons,csg_recons_ps,seis_in,seis_out,wf_file,gf_file,&
      ttt_in,ttt_out,sttt_in,gttt_in,sttt_s_in,gttt_s_in,&   !sttt_s_in,gttt_s_in are used in LSMF
      pw1ttt_in,pw2sttt_in,pw2gttt_in,&
      source_file,topo_file,mig0_file,mig0_pp_file,mig0_ps_file,& ! mig0_pp_file/mig0_ps_file are used in lsmf_PpPs
      refl_file,refl_ob_file,refl_fs_file,&
      img_file,img_ps_file,img_illum_file,pre_illum_out,pre_img_out,pre_img_out_ps,& ! LSMF
      pre_illum_img_out,z0_mig_file,res_file,res_tot_file,smwin_file,&
      vel0_file,gk1_file,prec_gk1_file,dk1_file,&
      lambda_file,mu_file,vp_file,vs_file,&   ! LSMF and elastic
      vp_homo_file,vs_homo_file,&             ! elastic muting
      csg_out_u,csg_out_v,csg_out_w,csg_tmp,csg_out_p,csg_dir,csg_mask,mask_csg,&         ! Used in multicomponent data output
      green_out,f0m_out,green_shot,csg_cal,&          ! Marchenko
      repath, &                                  ! result path
      fin,fout,&                                ! divi_csg
      green_total,green_d,green_bs,nm_image,&     ! added by Lei Fu ! Green_T: Direct Rayleigh wave, Green_BS: Back Scattered R Wave  
      vp0_file,vs0_file,mask_file,mask_shot_file,tt_obs_file,& ! Lei   
      csg_w,csg0_w,mask_salt_file
end type fname


type nsnr
  integer :: nshot, nrec 

end type nsnr


type inv
   character(len=200) :: inv_type,cg_type
   integer :: nit, reset,nearoffset_cut,nf,nmax , ewi_type,no_ofset ! nf = num. of high frequency, nmax= maxium integrations num.
   real    :: stl_coe,vmin,vmax,fl1,fl2,fl3,fl4,fl5,fl6,fh1,fh2,fh3,fh4,fh5,fh6,ofset_step,max_ofset ! fh1,fh2,fh3,fh4 4 different high frequencies
   logical :: isMig0, isNMLS, isReset, isSeisRecons, isPreCondition,&
              isBackTrack,isAdjointTest
end type inv

type trm
  character(len=200) :: fn_csg_w,fn_csg0_w
  integer :: nt,ng,ns,nmute,dmute
end type trm


type other
   character(len=200) :: pw1_side,sort_in, sort_out, eik_type, coordfile_type,&
      dfor
   logical :: isMask,isIllum, isIllumReceiver, isSaveIllum, isPreStackImg,isEikSep,&
              isShotConti, isSortData, isPreCondition,isCut,isDepthPC,isMask_Salt,&
              isSmooth_gk,isRemove_dir,isTaper_dir,isMask_csg,isWindow_csg,isCut_far_trace,&
              is_csg_filt,grad_flip!,is_output_pre,is_output_obs,is_output_residual
   real    :: xs_taper,ys_taper,illum_order,sort_tol,cut_distance,near_cut,far_cut
   integer :: nt, nt_in, nt_out, ns_taper,np1,np2,ncut,ncut2,hw,hw_m,hw_t
   real    :: dt, dt_in, dt_out,dfor_coe
   integer :: win,nh !window size for muteing first arrival by ::Lei
   integer :: lagmax,nlag,wt_res !by Lei
   integer :: Tmax,nit0,n_taper,nt0,n0,T_p !Max samples for tau, Tmax = tau_max/dt 
   real    :: taumax,toler,T_window !tau max in seconds 
   real    :: A_err,N_period ! A_err = 0.0001, N_period = 1.5 ~ 2
   real    :: fpk,fl,fh
   integer :: is_conv_sou,use_amp 
end type other

type aperture
   logical :: isAper
   real    :: x0, y0, z0, aper_x,aper_y,aper_z
end type aperture

type source
   real,allocatable :: sou(:)
   logical :: isSource
   integer          :: nw,source_shift,period_int    ! source_shift,period_int is added by Bowen Guo
   real             :: dt, freq ,period              ! period is added by Bowen Guo       
end type source

type fd2d
   integer, allocatable :: fs(:)
   character(len=200) :: fd_type,source_type,data_type
   integer :: nt, npml, ic, fd_order, bc_type
   logical :: isFS, isDipole, isSaveBC, isSaveWF, &
      isWriteData, isReadData
   logical :: sourcestress           ! modified by Bowen Guo
   real    :: dt
   real,allocatable     :: z0_mig(:)  ! modified by Bowen Guo
end type fd2d

type fd3d
   integer, allocatable :: fs(:,:)
   character(len=200) :: fd_type,source_type,data_type ! Modified by Bowen Guo
                                                       ! Used in module_e3d
   integer :: nt, npml, ic, fd_order, bc_type
   logical :: isFS, isDipole, isSaveBC, isSaveWF, & 
      isWriteData, isReadData
   real    :: dt
   
end type fd3d

!type admva2d
!   in 

type km2d
   complex,allocatable :: wfft(:)
   real,   allocatable :: z0_mig(:)
   character(len=200) :: eik_type,sort_in,dfor
   integer :: nt, nt_conv
   logical :: isSaveTTT,isSaveBothTTT,isTheta,isDip,isWriteData,isReadData,issvs,isgvs
   logical :: isWriteGreen,isWritef0m,isgttt,isSaveGreen,isSaveGreenShot   ! Used in marchenko 
   real          :: dt,theta, dip, dfor_coe
end type km2d

type km3d
   complex,allocatable :: wfft(:)
   real,   allocatable :: z0_mig(:,:)
   character(len=200) :: eik_type,sort_in,dfor
   integer :: nt, nt_conv
   logical :: isSaveTTT,isSaveBothTTT,isTheta,isDip,isWriteData,isReadData
   real          :: dt,theta, dip, dfor_coe
end type km3d

type gdm2d
   complex,allocatable :: wfft(:)
   real,   allocatable :: z0_mig(:)
   integer             :: nt,ntw_conv,nt_conv
   real                :: dt
   character(len=200)  :: sort_in
end type gdm2d

type gdm3d
   complex,allocatable :: wfft(:)
   real,   allocatable :: z0_mig(:,:)
   integer             :: nt,ntw_conv,nt_conv
   real                :: dt
   character(len=200)  :: sort_in
end type gdm3d

type coord2d
   integer, allocatable :: proid(:),ng(:),sid(:,:),gid(:,:)
   real, allocatable    :: xs(:,:), zs(:,:), xg(:,:), zg(:,:), t(:,:)
   integer              :: npro, ngmax, ntrace
end type coord2d

type coord3d
   integer, allocatable :: proid(:),ng(:),sid(:,:),gid(:,:)
   real, allocatable    :: xs(:,:),ys(:,:),zs(:,:),xg(:,:),yg(:,:),zg(:,:),t(:,:)
   integer              :: npro, ngmax, ntrace
end type coord3d

type eik_coord2d
   integer, allocatable :: proid(:),sid(:)
   real, allocatable    :: xs(:), zs(:)
   integer              :: npro
end type eik_coord2d

type eik_coord3d
   integer, allocatable :: proid(:),sid(:)
   real, allocatable    :: xs(:),ys(:),zs(:)
   integer              :: npro
end type eik_coord3d

type ms_coord2d
   integer, allocatable :: ns(:), ng(:), proid(:),sid(:,:),gid(:,:)
   real,    allocatable :: xs(:,:),zs(:,:),xg(:,:),zg(:,:),d(:,:),p(:,:),t(:,:)
   integer              :: npro, nsmax, ngmax, ntrace
end type ms_coord2d

type ms_coord3d
   integer, allocatable :: ns(:), ng(:),proid(:),sid(:,:),gid(:,:)
   real, allocatable    :: xs(:,:), ys(:,:), zs(:,:), xg(:,:), yg(:,:),&
      zg(:,:), d(:,:), p(:,:),t(:,:)
   integer              :: npro, nsmax, ngmax, ntrace
end type ms_coord3d

type pw1_coord2d
   integer, allocatable :: proid(:),ng(:),gid(:,:)
   real,    allocatable :: zp1(:),theta_x1(:),xref_1(:),xg(:,:),zg(:,:),t(:,:)
   integer              :: npro,ngmax,ntrace
end type pw1_coord2d

type pw1_shot2d
   integer, allocatable :: gid(:)
   real,allocatable :: xg(:), zg(:), seis(:,:)
   integer :: nt, ng, proid
   real    :: dt, zp1, theta_x1, xref_1
end type pw1_shot2d

type pw1_coord3d
   integer, allocatable :: proid(:),ng(:),gid(:,:)
   real,    allocatable :: zp1(:),theta_x1(:),theta_y1(:),xref_1(:),yref_1(:),&
                           xg(:,:),yg(:,:),zg(:,:),t(:,:)
   integer              :: npro,ngmax,ntrace
end type pw1_coord3d

type pw1_shot3d
   integer, allocatable :: gid(:)
   real,allocatable :: xg(:), yg(:),zg(:), seis(:,:)
   integer :: nt, ng, proid
   real    :: dt, zp1, theta_x1, xref_1, theta_y1, yref_1
end type pw1_shot3d

type pw2_coord2d
   integer, allocatable :: proid(:),ng(:),gid(:,:)
   real,    allocatable :: zp1(:),theta_x1(:),xref_1(:),&
      zp2(:,:),theta_x2(:,:),xref_2(:,:),t(:,:)
   integer              :: npro,ngmax,ntrace
end type pw2_coord2d

type pw2_shot2d
   integer, allocatable :: gid(:)
   real,allocatable :: zp2(:),theta_x2(:),xref_2(:),seis(:,:)
   integer :: nt, ng, proid
   real    :: dt, zp1, theta_x1, xref_1
end type pw2_shot2d

type pw2_coord3d
   integer              :: npro,ngmax,ntrace
   integer, allocatable :: proid(:),ng(:),gid(:,:)
   real,    allocatable :: zp1(:),theta_x1(:),theta_y1(:),xref_1(:),yref_1(:),&
      zp2(:,:),theta_x2(:,:),theta_y2(:,:),xref_2(:,:),yref_2(:,:),t(:,:)
end type pw2_coord3d

type pw2_shot3d
   integer, allocatable :: gid(:)
   real,allocatable :: seis(:,:),zp2(:),theta_x2(:),theta_y2(:),xref_2(:),yref_2(:)
   integer :: nt, ng, proid
   real    :: dt,zp1,theta_x1,theta_y1,xref_1,yref_1
end type pw2_shot3d

type model2d
   real,allocatable :: v(:,:), refl(:,:),den(:,:),vp(:,:),vs(:,:),vs0(:,:),lambda(:,:),mu(:,:),refl_pp(:,:),refl_ps(:,:),mask(:,:),mask_salt(:,:),mask_shot(:,:)! modified by Bowen Guo  
   integer :: nx, nz,nx_s,nz_s,nx_e,nz_e   ! modified by Bowen Guo (nz_sp,nx_sp, are the part of the model) 
   real    :: dx, dz, vmin, vmax
end type model2d

type model3d
   real,allocatable :: v(:,:,:), refl(:,:,:),den(:,:,:),vs(:,:,:),lambda(:,:,:),mu(:,:,:),vp(:,:,:)  ! modified by Bowen Guo
   integer :: nx, ny, nz
   real    :: dx, dy, dz, vmin, vmax
end type model3d

type water_model2d
   real, allocatable :: v_w(:,:),v_h(:,:),refl_fs(:,:),refl_ob(:,:)
   type(model2d) :: m
   integer :: nzw
end type water_model2d

type water_model3d
   real, allocatable :: v_w(:,:,:),v_h(:,:,:),refl_fs(:,:,:),refl_ob(:,:,:)
   type(model3d) :: m
   integer :: nzw
end type water_model3d

type shot2d
   integer, allocatable :: sid(:),gid(:)
   real,allocatable :: xs(:),zs(:),xg(:), zg(:), seis(:,:),seis_u(:,:),seis_w(:,:),seis_p(:,:),seis_all(:,:,:)
   real,allocatable :: seis_ps(:,:)  ! Used in LSMF
   real,allocatable :: seis_dir(:,:)  ! Used in LSMF
   ! seis_all is added by Bowen Guo
   integer :: nt, ng, proid, ns
   real    :: dt
   integer :: nps
end type shot2d

type shot3d
   integer, allocatable :: sid(:), gid(:)
   real,allocatable :: xs(:),ys(:),zs(:),xg(:),yg(:),zg(:),seis(:,:),seis_u(:,:),seis_v(:,:),seis_w(:,:)
   integer :: nt, ng, proid
   real    :: dt
end type shot3d

type ms_shot2d
   integer, allocatable :: sid(:), gid(:)
   real,allocatable :: xs(:), zs(:), xg(:), zg(:), d(:), p(:),seis(:,:)
   integer :: nt, ns, ng
   real    :: dt
end type ms_shot2d

type ms_shot3d
   integer, allocatable :: sid(:), gid(:)
   real,allocatable :: xs(:), ys(:), zs(:), xg(:), yg(:), zg(:), &
      d(:), p(:), seis(:,:)
   integer :: nt, ns, ng
   real    :: dt
end type ms_shot3d

type image1d
   real,allocatable    :: img(:),illum(:)
 !  logical :: isIllum, isIllumReceiver, isPreStackImg
 !  real    :: illum_order
   integer             :: nx,npml,ix1,ix2,iz1,iz2
end type image1d
type image2d
   real,allocatable    :: img(:,:),img_h(:,:,:),illum(:,:),img_t(:,:),img_m(:,:)
   logical :: isIllum, isIllumReceiver, isPreStackImg
   real    :: illum_order
   integer             :: nx,nz,npml,ix1,ix2,iz1,iz2
end type image2d

type image3d
   real,allocatable    :: img(:,:,:),illum(:,:,:)
   logical :: isIllum, isIllumReceiver, isPreStackImg
   real    :: illum_order
   integer             :: nx,ny,nz,npml,ix1,ix2,iy1,iy2,iz1,iz2
end type image3d

! Advanced datatype, can be composed by basic datatypes
type eik2d
   type(model2d) :: m
   real,allocatable :: ttt(:,:)
   real          :: xs, zs
end type eik2d

type eik3d
   type(model3d) :: m
   real,allocatable :: ttt(:,:,:)
   real          :: xs, ys, zs
end type eik3d

type pweik2d
   type(model2d) :: m
   real,allocatable :: ttt(:,:), t_init(:)
   real          :: xref, zp, theta
end type pweik2d

type pweik3d
   type(model3d) :: m
   real,allocatable :: ttt(:,:,:), t_init(:,:)
   real          :: xref, yref, zp, theta_x, theta_y
end type pweik3d

type fdmod2d
   type(model2d) :: m
   type(fd2d)    :: f
   type(source)  :: s
end type fdmod2d

type fdmod3d
   type(model3d) :: m
   type(fd3d)    :: f
   type(source)  :: s
end type fdmod3d

type water_fdmod2d
   type(water_model2d) :: m
   type(fd2d)    :: f
   type(source)  :: s
end type water_fdmod2d

type water_fdmod3d
   type(water_model3d) :: m
   type(fd3d)    :: f
   type(source)  :: s
end type water_fdmod3d

type kmmod2d
   type(model2d) :: m
   type(km2d)    :: k
   type(source)  :: s
end type kmmod2d

type kmmod3d
   type(model3d) :: m
   type(km3d)    :: k
   type(source)  :: s
end type kmmod3d

type ref2d
   integer,allocatable :: fs(:)
   real,allocatable :: ttt(:,:),v(:,:),xs(:),zs(:),xg(:),zg(:)
   integer       :: nz, nx, ng,ix1,ix2,iz1,iz2
   real          :: dz, dx
end type ref2d

type ref3d
   integer,allocatable :: fs(:,:)
   real,allocatable :: ttt(:,:,:),v(:,:,:),xs(:),ys(:),zs(:),xg(:),yg(:),zg(:)
   integer       :: nz, ny, nx, ng, ix1,ix2,iy1,iy2,iz1,iz2
   real          :: dz, dy,dx
end type ref3d

type gdmod2d
   type(gdm2d)   :: g
   type(model2d) :: m
   type(fd2d)    :: f
   type(source)  :: s
end type gdmod2d

type gdmod3d
   type(gdm3d)   :: g
   type(model3d) :: m
   type(fd3d)    :: f
   type(source)  :: s
end type gdmod3d

type wf2d
   integer           :: nx_2, nz_2, nt
   real, allocatable :: pwf(:,:,:),uwf(:,:,:),wwf(:,:,:) ! modified by Bowen Guo
end type wf2d

type wf3d
   integer           :: nx_2, ny_2,nz_2, nt
   real, allocatable :: pwf(:,:,:,:)
end type wf3d

interface operator(+)
   module procedure shot2d_add
   module procedure image2d_add
end interface

interface operator(-)
   module procedure shot2d_shot_sub_shot
   module procedure shot2d_real_sub_shot
   module procedure shot2d_shot_sub_real
   module procedure image2d_img_sub_img
end interface

interface operator(*)
   module procedure image2d_img_times_img
   module procedure image2d_real_times_img
end interface

interface operator(/)
   module procedure image2d_img_divide_img
end interface

interface operator(**)
   module procedure shot2d_exp
   module procedure image2d_exp
end interface

interface copytype
   module procedure copytype_source
   module procedure copytype_fd2d
   module procedure copytype_fd3d
   module procedure copytype_km2d
   module procedure copytype_km3d
   module procedure copytype_gdm2d
   module procedure copytype_gdm3d
   module procedure copytype_model2d
   module procedure copytype_model3d
   module procedure copytype_water_model2d
   module procedure copytype_water_model3d
   module procedure copytype_shot2d
   module procedure copytype_shot3d
   module procedure copytype_ms_shot2d
   module procedure copytype_ms_shot3d
   module procedure copytype_pw1_shot2d
   module procedure copytype_pw1_shot3d
   module procedure copytype_pw2_shot2d
   module procedure copytype_pw2_shot3d
   module procedure copytype_image2d
   module procedure copytype_image3d
end interface copytype

interface deallocate_type
   module procedure deallocate_km2d_m2d_s
   module procedure deallocate_km3d_m3d_s
   module procedure deallocate_fd2d_m2d_s
   module procedure deallocate_fd3d_m3d_s
   module procedure deallocate_2d_coord
   module procedure deallocate_3d_coord
   module procedure deallocate_eik_2d_coord
   module procedure deallocate_eik_3d_coord
   module procedure deallocate_ms2d_coord
   module procedure deallocate_ms3d_coord
   module procedure deallocate_pw1_2d_coord
   module procedure deallocate_pw1_2d_shot
   module procedure deallocate_pw2_2d_coord
   module procedure deallocate_pw2_2d_shot
   module procedure deallocate_pw1_3d_coord
   module procedure deallocate_pw1_3d_shot
   module procedure deallocate_pw2_3d_coord
   module procedure deallocate_pw2_3d_shot
   module procedure deallocate_source
   module procedure deallocate_km2d
   module procedure deallocate_km3d
   module procedure deallocate_gdm2d
   module procedure deallocate_gdm3d
   module procedure deallocate_shot2d
   module procedure deallocate_shot3d
   module procedure deallocate_model2d
   module procedure deallocate_model3d
   module procedure deallocate_ref2d
   module procedure deallocate_ref3d
   module procedure deallocate_gdmod2d
   module procedure deallocate_gdmod3d
   module procedure deallocate_gdm2d_m2_s
end interface deallocate_type

contains

type(shot2d) function shot2d_add(shot1,shot2)
type(shot2d), intent(in) :: shot1, shot2
shot2d_add%seis = shot1%seis+shot2%seis
end function shot2d_add

type(image2d) function image2d_add(img1,img2)
type(image2d), intent(in) :: img1, img2
image2d_add%img = img1%img + img2%img
end function image2d_add

type(shot2d) function shot2d_shot_sub_shot(shot1,shot2)
type(shot2d), intent(in) :: shot1, shot2
shot2d_shot_sub_shot%seis = shot1%seis - shot2%seis
end function shot2d_shot_sub_shot

type(shot2d) function shot2d_real_sub_shot(num,shot2)
real,         intent(in) :: num
type(shot2d), intent(in) :: shot2
shot2d_real_sub_shot%seis = num-shot2%seis
end function shot2d_real_sub_shot

type(shot2d) function shot2d_shot_sub_real(shot1,num)
real,         intent(in) :: num
type(shot2d), intent(in) :: shot1
shot2d_shot_sub_real%seis = shot1%seis-num
end function shot2d_shot_sub_real

type(image2d) function image2d_img_sub_img(img1,img2)
type(image2d), intent(in) :: img1, img2
image2d_img_sub_img%img = img1%img - img2%img
end function image2d_img_sub_img

type(image2d) function image2d_img_times_img(img1,img2)
type(image2d), intent(in) :: img1, img2
image2d_img_times_img%img = img1%img * img2%img
end function image2d_img_times_img

type(image2d) function image2d_real_times_img(num,img)
type(image2d), intent(in) :: img
real,          intent(in) :: num
image2d_real_times_img%img = num * img%img
end function image2d_real_times_img

type(image2d) function image2d_img_divide_img(img1,img2)
type(image2d), intent(in) :: img1, img2
image2d_img_divide_img%img = img1%img / img2%img
end function image2d_img_divide_img

type(image2d) function image2d_img_divide_array(img,array)
type(image2d), intent(in) :: img
real,          intent(in) :: array(:,:)
image2d_img_divide_array%img = img%img / array
end function image2d_img_divide_array

type(shot2d) function shot2d_exp(shot,num)
type(shot2d), intent(in) :: shot
real,         intent(in) :: num
shot2d_exp%seis = shot%seis**num
end function shot2d_exp

type(image2d) function image2d_exp(img,num)
type(image2d), intent(in) :: img
real,          intent(in) :: num
image2d_exp%img = img%img**num
end function image2d_exp

subroutine copytype_fd2d(f_in,f_out)
type(fd2d),  intent(in)  :: f_in
type(fd2d),  intent(out) :: f_out
call copy_variable(f_in%nt,f_out%nt,f_in%npml,f_out%npml,f_in%ic,f_out%ic)
call copy_variable(f_in%fd_type,f_out%fd_type)
call copy_variable(f_in%fd_order,f_out%fd_order,f_in%bc_type,f_out%bc_type)
call copy_variable(f_in%isFS,f_out%isFS,f_in%isDipole,f_out%isDipole)
call copy_variable(f_in%isWriteData,f_out%isWriteData,f_in%isReadData,f_out%isReadData)
call copy_variable(f_in%isSaveBC,f_out%isSaveBC,f_in%isSaveWF,f_out%isSaveWF)
call copy_variable(f_in%dt,f_out%dt)
call copy_variable(f_in%source_type,f_out%source_type)        ! modified by Bowen Guo
call copy_variable(f_in%data_type,f_out%data_type)      ! modified by Bowen Guo
if (allocated(f_in%fs).and.allocated(f_out%fs)) then
   if (size(f_in%fs).eq.size(f_out%fs)) then
      call copy_array(f_in%fs,f_out%fs)
   endif
endif
end subroutine copytype_fd2d

subroutine copytype_fd3d(f_in,f_out)
type(fd3d),  intent(in)  :: f_in
type(fd3d),  intent(out) :: f_out
call copy_variable(f_in%nt,f_out%nt,f_in%npml,f_out%npml,f_in%ic,f_out%ic)
call copy_variable(f_in%fd_type,f_out%fd_type)
call copy_variable(f_in%data_type,f_out%data_type) ! Modified by Bowen Guo
call copy_variable(f_in%source_type,f_out%source_type) ! Modified by Bowen Guo
call copy_variable(f_in%fd_order,f_out%fd_order,f_in%bc_type,f_out%bc_type)
call copy_variable(f_in%isFS,f_out%isFS,f_in%isDipole,f_out%isDipole)
call copy_variable(f_in%isWriteData,f_out%isWriteData,f_in%isReadData,f_out%isReadData)
call copy_variable(f_in%isSaveBC,f_out%isSaveBC,f_in%isSaveWF,f_out%isSaveWF)
call copy_variable(f_in%dt,f_out%dt)
if (allocated(f_in%fs).and.allocated(f_out%fs)) then
   if ((size(f_in%fs,1).eq.size(f_out%fs,1)) .and. &
       (size(f_in%fs,2).eq.size(f_out%fs,2)) ) then
      call copy_array(f_in%fs,f_out%fs)
   endif
endif
end subroutine copytype_fd3d

subroutine copytype_source(s_in,s_out)
type(source),intent(in)  :: s_in
type(source),intent(out) :: s_out
s_out%nw=s_in%nw;s_out%dt=s_in%dt;s_out%freq=s_in%freq;s_out%isSource=s_in%isSource
if (allocated(s_in%sou).and.allocated(s_out%sou)) then
   if (size(s_in%sou).eq.size(s_out%sou)) then
      call copy_array(s_in%sou,s_out%sou)
   endif
endif
end subroutine copytype_source

subroutine copytype_model2d(m_in,m_out)
type(model2d),   intent(in)  :: m_in
type(model2d),   intent(out) :: m_out
call copy_variable(m_in%nx,m_out%nx,m_in%nz,m_out%nz)
call copy_variable(m_in%dx,m_out%dx,m_in%dz,m_out%dz)
call copy_variable(m_in%vmin,m_out%vmin,m_in%vmax,m_out%vmax)
if (allocated(m_in%v).and.allocated(m_out%v)) then
   if ((size(m_in%v,1).eq.size(m_out%v,1)) .and. &
       (size(m_in%v,2).eq.size(m_out%v,2)) ) then
      call copy_array(m_in%v,m_out%v)
   endif
endif
if (allocated(m_in%vs).and.allocated(m_out%vs)) then
   if ((size(m_in%vs,1).eq.size(m_out%vs,1)) .and. &
       (size(m_in%vs,2).eq.size(m_out%vs,2)) ) then
      call copy_array(m_in%vs,m_out%vs)
   endif
endif
if (allocated(m_in%v).and.allocated(m_out%v)) then
   if ((size(m_in%v,1).eq.size(m_out%v,1)) .and. &
       (size(m_in%v,2).eq.size(m_out%v,2)) ) then
      call copy_array(m_in%v,m_out%v)
   endif
endif
if (allocated(m_in%refl).and.allocated(m_out%refl)) then
   if ((size(m_in%refl,1).eq.size(m_out%refl,1)) .and. &
       (size(m_in%refl,2).eq.size(m_out%refl,2)) ) then
      call copy_array(m_in%refl,m_out%refl)
   endif
endif
if (allocated(m_in%refl_ps).and.allocated(m_out%refl_ps)) then
   if ((size(m_in%refl_ps,1).eq.size(m_out%refl_ps,1)) .and. &
       (size(m_in%refl_ps,2).eq.size(m_out%refl_ps,2)) ) then
      call copy_array(m_in%refl_ps,m_out%refl_ps)
   endif
endif
if (allocated(m_in%den).and.allocated(m_out%den)) then
   if ((size(m_in%den,1).eq.size(m_out%den,1)) .and. &
       (size(m_in%den,2).eq.size(m_out%den,2)) ) then
      call copy_array(m_in%den,m_out%den)
   endif
endif 
if (allocated(m_in%lambda).and.allocated(m_out%lambda)) then      ! modified by Bowen Guo
   if ((size(m_in%lambda,1).eq.size(m_out%lambda,1)) .and. &     ! modified by Bowen Guo
       (size(m_in%lambda,2).eq.size(m_out%lambda,2)) ) then      ! modified by Bowen Guo
      call copy_array(m_in%lambda,m_out%lambda)                  ! modified by Bowen Guo
   endif                                                         ! modified by Bowen Guo
endif                                                            ! modified by Bowen Guo
if (allocated(m_in%mu).and.allocated(m_out%mu)) then              ! modified by Bowen Guo
   if ((size(m_in%mu,1).eq.size(m_out%mu,1)) .and. &             ! modified by Bowen Guo
       (size(m_in%mu,2).eq.size(m_out%mu,2)) ) then              ! modified by Bowen Guo
      call copy_array(m_in%mu,m_out%mu)                          ! modified by Bowen Guo
   endif                                                         ! modified by Bowen Guo
endif                                                            ! modified by Bowen Guo
if (allocated(m_in%vp).and.allocated(m_out%vp)) then              ! modified by Bowen Guo
   if ((size(m_in%vp,1).eq.size(m_out%vp,1)) .and. &             ! modified by Bowen Guo
       (size(m_in%vp,2).eq.size(m_out%vp,2)) ) then              ! modified by Bowen Guo
      call copy_array(m_in%vp,m_out%vp)                          ! modified by Bowen Guo
   endif                                                         ! modified by Bowen Guo
endif                                                            ! modified by Bowen Guo
if (allocated(m_in%vs).and.allocated(m_out%vs)) then
   if ((size(m_in%vs,1).eq.size(m_out%vs,1)) .and. &
       (size(m_in%vs,2).eq.size(m_out%vs,2)) ) then
      call copy_array(m_in%vs,m_out%vs)
   endif
endif
end subroutine copytype_model2d

!==========================================================================================
subroutine copytype_model3d(m_in,m_out)
type(model3d),   intent(in)  :: m_in
type(model3d),   intent(out) :: m_out
call copy_variable(m_in%nx,m_out%nx,m_in%ny,m_out%ny,m_in%nz,m_out%nz)
call copy_variable(m_in%dx,m_out%dx,m_in%dy,m_out%dy,m_in%dz,m_out%dz)
call copy_variable(m_in%vmin,m_out%vmin,m_in%vmax,m_out%vmax)
if (allocated(m_in%v).and.allocated(m_out%v)) then
   if ((size(m_in%v,1).eq.size(m_out%v,1)) .and. &
       (size(m_in%v,2).eq.size(m_out%v,2)) .and. &
       (size(m_in%v,3).eq.size(m_out%v,3)) ) then
      call copy_array(m_in%v,m_out%v)
   endif
endif
if (allocated(m_in%refl).and.allocated(m_out%refl)) then
   if ((size(m_in%refl,1).eq.size(m_out%refl,1)) .and. & 
       (size(m_in%refl,2).eq.size(m_out%refl,2)) .and. &
       (size(m_in%refl,3).eq.size(m_out%refl,3)) ) then
      call copy_array(m_in%refl,m_out%refl)
   endif
endif
if (allocated(m_in%den).and.allocated(m_out%den)) then
   if ((size(m_in%den,1).eq.size(m_out%den,1)) .and. &
       (size(m_in%den,2).eq.size(m_out%den,2)) .and. &
       (size(m_in%den,3).eq.size(m_out%den,3)) ) then
      call copy_array(m_in%den,m_out%den)
   endif
endif
if (allocated(m_in%vs).and.allocated(m_out%vs)) then
   if ((size(m_in%vs,1).eq.size(m_out%vs,1)) .and. & 
       (size(m_in%vs,2).eq.size(m_out%vs,2)) .and. &
       (size(m_in%vs,3).eq.size(m_out%vs,3)) ) then
      call copy_array(m_in%vs,m_out%vs)
   endif
endif
if (allocated(m_in%vp).and.allocated(m_out%vp)) then
   if ((size(m_in%vp,1).eq.size(m_out%vp,1)) .and. & 
       (size(m_in%vp,2).eq.size(m_out%vp,2)) .and. &
       (size(m_in%vp,3).eq.size(m_out%vp,3))) then
      call copy_array(m_in%vp,m_out%vp)
   endif
endif                 ! Modified by Bowen Guo
if (allocated(m_in%lambda).and.allocated(m_out%lambda)) then
   if ((size(m_in%lambda,1).eq.size(m_out%lambda,1)) .and. & 
       (size(m_in%lambda,2).eq.size(m_out%lambda,2)) .and. &
       (size(m_in%lambda,3).eq.size(m_out%lambda,3)) ) then
      call copy_array(m_in%lambda,m_out%lambda)
   endif
endif                 ! Modified by Bowen Guo
if (allocated(m_in%mu).and.allocated(m_out%mu)) then
   if ((size(m_in%mu,1).eq.size(m_out%mu,1)) .and. & 
       (size(m_in%mu,2).eq.size(m_out%mu,2)) .and. &
       (size(m_in%mu,3).eq.size(m_out%mu,3)) ) then
      call copy_array(m_in%mu,m_out%mu)
   endif
endif                 ! Modified by Bowen Guo

end subroutine copytype_model3d

subroutine copytype_water_model2d(m_in,m_out)
type(water_model2d),   intent(in)  :: m_in
type(water_model2d),   intent(out) :: m_out
call copytype_model2d(m_in%m,m_out%m)
call copy_variable(m_in%nzw,m_out%nzw)
end subroutine copytype_water_model2d

subroutine copytype_water_model3d(m_in,m_out)
type(water_model3d),   intent(in)  :: m_in
type(water_model3d),   intent(out) :: m_out
call copytype_model3d(m_in%m,m_out%m)
call copy_variable(m_in%nzw,m_out%nzw)
end subroutine copytype_water_model3d

subroutine copytype_km2d(k_in,k_out)
type(km2d),     intent(in)  :: k_in
type(km2d),     intent(out) :: k_out
call copy_variable(k_in%isSaveTTT,k_out%isSaveTTT,k_in%isSaveBothTTT,k_out%isSaveBothTTT)
call copy_variable(k_in%eik_type,k_out%eik_type,k_in%sort_in,k_out%sort_in)
call copy_variable(k_in%dfor,k_out%dfor);
call copy_variable(k_in%isTheta,k_out%isTheta,k_in%isDip,k_out%isDip)
call copy_variable(k_in%isWriteData,k_out%isWriteData,k_in%isReadData,k_out%isReadData)
call copy_variable(k_in%nt,k_out%nt,k_in%nt_conv,k_out%nt_conv)
call copy_variable(k_in%dt,k_out%dt,k_in%theta,k_out%theta,k_in%dip,k_out%dip)
call copy_variable(k_in%issvs,k_out%issvs) !  Modified by Bowen Guo
call copy_variable(k_in%isgvs,k_out%isgvs) !  Modified by Bowen Guo
if (allocated(k_in%wfft).and.allocated(k_out%wfft)) then
   if (size(k_in%wfft).eq.size(k_out%wfft)) then
      call copy_array(k_in%wfft,k_out%wfft)
   endif
endif
if (allocated(k_in%z0_mig).and.allocated(k_out%z0_mig)) then
   if (size(k_in%z0_mig).eq.size(k_out%z0_mig)) then
      call copy_array(k_in%z0_mig,k_out%z0_mig)
   endif
endif
call copy_variable(k_in%dfor,k_out%dfor);
call copy_variable(k_in%dfor_coe,k_out%dfor_coe)
end subroutine copytype_km2d

subroutine copytype_km3d(k_in,k_out)
type(km3d),     intent(in)  :: k_in
type(km3d),     intent(out) :: k_out
call copy_variable(k_in%isSaveTTT,k_out%isSaveTTT,k_in%isSaveBothTTT,k_out%isSaveBothTTT)
call copy_variable(k_in%eik_type,k_out%eik_type,k_in%sort_in,k_out%sort_in)
call copy_variable(k_in%isTheta,k_out%isTheta,k_in%isDip,k_out%isDip)
call copy_variable(k_in%isWriteData,k_out%isWriteData,k_in%isReadData,k_out%isReadData)
call copy_variable(k_in%nt,k_out%nt,k_in%nt_conv,k_out%nt_conv)
call copy_variable(K_in%dt,k_out%dt,k_in%theta,k_out%theta,k_in%dip,k_out%dip)
if (allocated(k_in%wfft).and.allocated(k_out%wfft)) then
   if (size(k_in%wfft).eq.size(k_out%wfft)) then
      call copy_array(k_in%wfft,k_out%wfft)
   endif
endif
if (allocated(k_in%z0_mig).and.allocated(k_out%z0_mig)) then
   if ((size(k_in%z0_mig,1).eq.size(k_out%z0_mig,1)) .and. &
       (size(k_in%z0_mig,2).eq.size(k_out%z0_mig,2)) ) then
      call copy_array(k_in%z0_mig,k_out%z0_mig)
   endif
endif
call copy_variable(k_in%dfor,k_out%dfor);
call copy_variable(k_in%dfor_coe,k_out%dfor_coe)
end subroutine copytype_km3d

subroutine copytype_gdm2d(g_in,g_out)
type(gdm2d),     intent(in)  :: g_in
type(gdm2d),     intent(out) :: g_out
call copy_variable(g_in%nt,g_out%nt,g_in%ntw_conv,g_out%ntw_conv,g_in%nt_conv,g_out%nt_conv)
call copy_variable(g_in%dt,g_out%dt)
call copy_variable(g_in%sort_in,g_out%sort_in)
if (allocated(g_in%wfft).and.allocated(g_out%wfft)) then
   if (size(g_in%wfft).eq.size(g_out%wfft)) then
      call copy_array(g_in%wfft,g_out%wfft)
   endif
endif
if (allocated(g_in%z0_mig).and.allocated(g_out%z0_mig)) then
   if (size(g_in%z0_mig).eq.size(g_out%z0_mig)) then
      call copy_array(g_in%z0_mig,g_out%z0_mig)
   endif
endif
end subroutine copytype_gdm2d

subroutine copytype_gdm3d(g_in,g_out)
type(gdm3d),     intent(in)  :: g_in
type(gdm3d),     intent(out) :: g_out
call copy_variable(g_in%nt,g_out%nt,g_in%ntw_conv,g_out%ntw_conv,g_in%nt_conv,g_out%nt_conv)
call copy_variable(g_in%dt,g_out%dt)
call copy_variable(g_in%sort_in,g_out%sort_in)
if (allocated(g_in%wfft).and.allocated(g_out%wfft)) then
   if (size(g_in%wfft).eq.size(g_out%wfft)) then
      call copy_array(g_in%wfft,g_out%wfft)
   endif
endif
if (allocated(g_in%z0_mig).and.allocated(g_out%z0_mig)) then
   if (size(g_in%z0_mig).eq.size(g_out%z0_mig)) then
      call copy_array(g_in%z0_mig,g_out%z0_mig)
   endif
endif
end subroutine copytype_gdm3d

subroutine copytype_shot2d(sh_in,sh_out)
type(shot2d),     intent(in)  :: sh_in
type(shot2d),     intent(out) :: sh_out
call copy_variable(sh_in%nt,sh_out%nt,sh_in%ng,sh_out%ng)
call copy_variable(sh_in%dt,sh_out%dt)
if (allocated(sh_in%sid).and.allocated(sh_out%sid)) then
   call copy_array(sh_in%sid,sh_out%sid,sh_in%gid,sh_out%gid)
   call copy_array(sh_in%xs,sh_out%xs,sh_in%zs,sh_out%zs)
   call copy_array(sh_in%xg,sh_out%xg,sh_in%zg,sh_out%zg)
   call copy_array(sh_in%seis,sh_out%seis)
endif
end subroutine copytype_shot2d

subroutine copytype_shot3d(sh_in,sh_out)
type(shot3d),     intent(in)  :: sh_in
type(shot3d),     intent(out) :: sh_out
call copy_variable(sh_in%nt,sh_out%nt,sh_in%ng,sh_out%ng)
call copy_variable(sh_in%dt,sh_out%dt)
if (allocated(sh_in%sid).and.allocated(sh_out%sid)) then
   call copy_array(sh_in%sid,sh_out%sid,sh_in%gid,sh_out%gid)
   call copy_array(sh_in%xs,sh_out%xs,sh_in%zs,sh_out%zs)
   call copy_array(sh_in%xg,sh_out%xg,sh_in%zg,sh_out%zg)
   call copy_array(sh_in%ys,sh_out%ys,sh_in%yg,sh_out%yg)
   call copy_array(sh_in%seis,sh_out%seis)
endif
end subroutine copytype_shot3d

subroutine copytype_ms_shot2d(sh_in,sh_out)
type(ms_shot2d),     intent(in)  :: sh_in
type(ms_shot2d),     intent(out) :: sh_out
call copy_variable(sh_in%nt,sh_out%nt,sh_in%ns,sh_out%ns,sh_in%ng,sh_out%ng)
call copy_variable(sh_in%dt,sh_out%dt)
if (allocated(sh_in%xs).and.allocated(sh_out%xs)) then
   call copy_array(sh_in%sid,sh_out%sid,sh_in%gid,sh_out%gid)
   call copy_array(sh_in%d,sh_out%d,sh_in%p,sh_out%p)
   call copy_array(sh_in%xs,sh_out%xs,sh_in%zs,sh_out%zs)
   call copy_array(sh_in%xg,sh_out%xg,sh_in%zg,sh_out%zg)
   call copy_array(sh_in%seis,sh_out%seis)
endif
end subroutine copytype_ms_shot2d

subroutine copytype_ms_shot3d(sh_in,sh_out)
type(ms_shot3d),     intent(in)  :: sh_in
type(ms_shot3d),     intent(out) :: sh_out
call copy_variable(sh_in%nt,sh_out%nt,sh_in%ns,sh_out%ns,sh_in%ng,sh_out%ng)
call copy_variable(sh_in%dt,sh_out%dt)
if (allocated(sh_in%xs).and.allocated(sh_out%xs)) then
   call copy_array(sh_in%sid,sh_out%sid,sh_in%gid,sh_out%gid)
   call copy_array(sh_in%d,sh_out%d,sh_in%p,sh_out%p)
   call copy_array(sh_in%xs,sh_out%xs,sh_in%zs,sh_out%zs)
   call copy_array(sh_in%xg,sh_out%xg,sh_in%zg,sh_out%zg)
   call copy_array(sh_in%ys,sh_out%ys,sh_in%yg,sh_out%yg)
   call copy_array(sh_in%seis,sh_out%seis)
endif
end subroutine copytype_ms_shot3d

subroutine copytype_pw1_shot2d(sh_in,sh_out)
type(pw1_shot2d),     intent(in)  :: sh_in
type(pw1_shot2d),     intent(out) :: sh_out
call copy_variable(sh_in%nt,sh_out%nt,sh_in%proid,sh_out%proid,&
   sh_in%ng,sh_out%ng)
call copy_variable(sh_in%dt,sh_out%dt,sh_in%zp1,sh_out%zp1)
if (allocated(sh_in%xg).and.allocated(sh_out%xg)) then
   call copy_array(sh_in%gid,sh_out%gid)
   call copy_array(sh_in%xg,sh_out%xg,sh_in%zg,sh_out%zg)
endif
if (allocated(sh_in%seis).and.allocated(sh_out%seis)) then
   call copy_array(sh_in%seis,sh_out%seis)
endif
end subroutine copytype_pw1_shot2d

subroutine copytype_pw1_shot3d(sh_in,sh_out)
type(pw1_shot3d),     intent(in)  :: sh_in
type(pw1_shot3d),     intent(out) :: sh_out
call copy_variable(sh_in%nt,sh_out%nt,sh_in%proid,sh_out%proid,&
   sh_in%ng,sh_out%ng)
call copy_variable(sh_in%dt,sh_out%dt,sh_in%zp1,sh_out%zp1)
if (allocated(sh_in%xg).and.allocated(sh_out%xg)) then
   call copy_array(sh_in%gid,sh_out%gid)
   call copy_array(sh_in%xg,sh_out%xg,sh_in%yg,sh_out%yg,&
        sh_in%zg,sh_out%zg)
endif
if (allocated(sh_in%seis).and.allocated(sh_out%seis)) then
   call copy_array(sh_in%seis,sh_out%seis)
endif
end subroutine copytype_pw1_shot3d

subroutine copytype_pw2_shot2d(sh_in,sh_out)
type(pw2_shot2d),     intent(in)  :: sh_in
type(pw2_shot2d),     intent(out) :: sh_out
call copy_variable(sh_in%nt,sh_out%nt,sh_in%proid,sh_out%proid,&
   sh_in%ng,sh_out%ng)
call copy_variable(sh_in%dt,sh_out%dt,sh_in%zp1,sh_out%zp1,&
       sh_in%theta_x1,sh_out%theta_x1,sh_in%xref_1,sh_out%xref_1)
if (allocated(sh_in%theta_x2).and.allocated(sh_out%theta_x2)) then
   call copy_array(sh_in%gid,sh_out%gid)
   call copy_array(sh_in%theta_x2,sh_out%theta_x2)
   call copy_array(sh_in%xref_2,sh_out%xref_2,sh_in%zp2,sh_out%zp2)
endif
if (allocated(sh_in%seis).and.allocated(sh_out%seis)) then
   call copy_array(sh_in%seis,sh_out%seis)
endif
end subroutine copytype_pw2_shot2d

subroutine copytype_pw2_shot3d(sh_in,sh_out)
type(pw2_shot3d),     intent(in)  :: sh_in
type(pw2_shot3d),     intent(out) :: sh_out
call copy_variable(sh_in%nt,sh_out%nt,sh_in%proid,sh_out%proid,&
   sh_in%ng,sh_out%ng)
call copy_variable(sh_in%dt,sh_out%dt,sh_in%zp1,sh_out%zp1,&
        sh_in%theta_x1,sh_out%theta_x1,sh_in%xref_1,sh_out%xref_1,&
        sh_in%theta_y1,sh_out%theta_y1,sh_in%yref_1,sh_out%yref_1)
if (allocated(sh_in%theta_x2).and.allocated(sh_out%theta_x2)) then
   call copy_array(sh_in%gid,sh_out%gid)
   call copy_array(sh_in%theta_x2,sh_out%theta_x2)
   call copy_array(sh_in%theta_y2,sh_out%theta_y2)
   call copy_array(sh_in%yref_2,sh_out%yref_2)
   call copy_array(sh_in%xref_2,sh_out%xref_2,sh_in%zp2,sh_out%zp2)
endif
if (allocated(sh_in%seis).and.allocated(sh_out%seis)) then
   call copy_array(sh_in%seis,sh_out%seis)
endif
end subroutine copytype_pw2_shot3d

subroutine copytype_image2d(img_in,img_out)
type(image2d), intent(in) :: img_in
type(image2d), intent(out):: img_out
call copy_variable(img_in%isIllum,img_out%isIllum)
call copy_variable(img_in%isIllumReceiver,img_out%isIllumReceiver)
call copy_variable(img_in%isPreStackImg,img_out%isPreStackImg)
call copy_variable(img_in%illum_order,img_out%illum_order)
call copy_variable(img_in%nx,img_out%nx,img_in%nz,img_out%nz,&
   img_in%npml,img_out%npml)
call copy_variable(img_in%ix1,img_out%ix1,img_in%ix2,img_out%ix2,&
                   img_in%iz1,img_out%iz1,img_in%iz2,img_out%iz2)
if (allocated(img_in%img).and.allocated(img_out%img)) then
   call copy_array(img_in%img,img_out%img)
endif
if (allocated(img_in%illum).and.allocated(img_out%illum)) then
   call copy_array(img_in%illum,img_out%illum)
endif
end subroutine copytype_image2d

subroutine copytype_image3d(img_in,img_out)
type(image3d), intent(in) :: img_in
type(image3d), intent(out):: img_out
call copy_variable(img_in%isIllum,img_out%isIllum)
call copy_variable(img_in%isIllumReceiver,img_out%isIllumReceiver)
call copy_variable(img_in%isPreStackImg,img_out%isPreStackImg)
call copy_variable(img_in%illum_order,img_out%illum_order)
call copy_variable(img_in%nx,img_out%nx,img_in%ny,img_out%ny,&
        img_in%nz,img_out%nz,img_in%npml,img_out%npml)
call copy_variable(img_in%ix1,img_out%ix1,img_in%ix2,img_out%ix2)
call copy_variable(img_in%iy1,img_out%iy1,img_in%iy2,img_out%iy2)
call copy_variable(img_in%iz1,img_out%iz1,img_in%iz2,img_out%iz2)
if (allocated(img_in%img).and.allocated(img_out%img)) then
   call copy_array(img_in%img,img_out%img)
endif
if (allocated(img_in%illum).and.allocated(img_out%illum)) then
   call copy_array(img_in%illum,img_out%illum)
endif
end subroutine copytype_image3d

subroutine deallocate_2d_coord(c)
type(coord2d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid,c%ng)
if (allocated(c%sid)) call deallocate_and_free(c%sid,c%gid)
if (allocated(c%xs)) call deallocate_and_free(c%xs,c%zs,c%xg,c%zg,c%t)
end subroutine deallocate_2d_coord

subroutine deallocate_3d_coord(c)
type(coord3d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid,c%ng)
if (allocated(c%sid)) call deallocate_and_free(c%sid,c%gid)
if (allocated(c%xs)) call deallocate_and_free(c%xs,c%ys,c%zs,c%xg,c%zg,c%t)
end subroutine deallocate_3d_coord

subroutine deallocate_eik_2d_coord(c)
type(eik_coord2d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid)
if (allocated(c%sid)) call deallocate_and_free(c%sid)
if (allocated(c%xs)) call deallocate_and_free(c%xs,c%zs)
end subroutine deallocate_eik_2d_coord

subroutine deallocate_eik_3d_coord(c)
type(eik_coord3d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid)
if (allocated(c%sid)) call deallocate_and_free(c%sid)
if (allocated(c%xs)) call deallocate_and_free(c%xs,c%ys,c%zs)
end subroutine deallocate_eik_3d_coord

subroutine deallocate_ms2d_coord(c)
type(ms_coord2d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid,c%ns,c%ng)
if (allocated(c%sid)) call deallocate_and_free(c%sid,c%gid)
if (allocated(c%xs)) call deallocate_and_free(c%xs,c%zs,c%xg,c%zg,c%d,c%p,c%t)
end subroutine deallocate_ms2d_coord

subroutine deallocate_ms3d_coord(c)
type(ms_coord3d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid,c%ns,c%ng)
if (allocated(c%sid)) call deallocate_and_free(c%sid,c%gid)
if (allocated(c%xs)) call deallocate_and_free(c%xs,c%ys,c%zs,c%xg,c%yg,c%zg,c%d,c%p,c%t)
end subroutine deallocate_ms3d_coord

subroutine deallocate_pw1_2d_coord(c)
type(pw1_coord2d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid)
if (allocated(c%ng)) call deallocate_and_free(c%ng)
if (allocated(c%gid)) call deallocate_and_free(c%gid)
if (allocated(c%zp1)) call deallocate_and_free(c%zp1,c%theta_x1,c%xref_1)
if (allocated(c%xg)) call deallocate_and_free(c%xg,c%zg,c%t)
end subroutine deallocate_pw1_2d_coord

subroutine deallocate_pw1_3d_coord(c)
type(pw1_coord3d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid)
if (allocated(c%ng)) call deallocate_and_free(c%ng)
if (allocated(c%gid)) call deallocate_and_free(c%gid)
if (allocated(c%zp1)) call deallocate_and_free(c%zp1,c%theta_x1,c%xref_1,c%theta_y1,c%yref_1)
if (allocated(c%xg)) call deallocate_and_free(c%xg,c%yg,c%zg,c%t)
end subroutine deallocate_pw1_3d_coord

subroutine deallocate_pw2_2d_coord(c)
type(pw2_coord2d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid)
if (allocated(c%ng)) call deallocate_and_free(c%ng)
if (allocated(c%gid)) call deallocate_and_free(c%gid)
if (allocated(c%zp1)) call deallocate_and_free(c%zp1,c%theta_x1,c%xref_1)
if (allocated(c%zp2)) call deallocate_and_free(c%zp2,c%theta_x2,c%xref_2,c%t)
end subroutine deallocate_pw2_2d_coord

subroutine deallocate_pw2_3d_coord(c)
type(pw2_coord3d), intent(inout) :: c
if (allocated(c%proid)) call deallocate_and_free(c%proid)
if (allocated(c%ng)) call deallocate_and_free(c%ng)
if (allocated(c%gid)) call deallocate_and_free(c%gid)
if (allocated(c%zp1)) call deallocate_and_free(c%zp1,c%theta_x1,c%xref_1,c%theta_y1,c%yref_1)
if (allocated(c%zp2)) call deallocate_and_free(c%zp2,c%theta_x2,c%xref_2,c%theta_y2,c%yref_2,c%t)
end subroutine deallocate_pw2_3d_coord

subroutine deallocate_pw1_2d_shot(sh)
type(pw1_shot2d), intent(inout) :: sh
if (allocated(sh%xg)) call deallocate_and_free(sh%xg,sh%zg)
if (allocated(sh%gid)) call deallocate_and_free(sh%gid)
if (allocated(sh%seis)) deallocate(sh%seis)
end subroutine deallocate_pw1_2d_shot

subroutine deallocate_pw1_3d_shot(sh)
type(pw1_shot3d), intent(inout) :: sh
if (allocated(sh%xg)) call deallocate_and_free(sh%xg,sh%yg,sh%zg)
if (allocated(sh%gid)) call deallocate_and_free(sh%gid)
if (allocated(sh%seis)) deallocate(sh%seis)
end subroutine deallocate_pw1_3d_shot

subroutine deallocate_pw2_2d_shot(sh)
type(pw2_shot2d), intent(inout) :: sh
if (allocated(sh%gid)) call deallocate_and_free(sh%gid)
if (allocated(sh%zp2)) call deallocate_and_free(sh%theta_x2,sh%xref_2,sh%zp2)
if (allocated(sh%seis)) deallocate(sh%seis)
end subroutine deallocate_pw2_2d_shot

subroutine deallocate_pw2_3d_shot(sh)
type(pw2_shot3d), intent(inout) :: sh
if (allocated(sh%gid)) call deallocate_and_free(sh%gid)
if (allocated(sh%zp2)) then
   call deallocate_and_free(sh%theta_x2,sh%xref_2,sh%theta_y2,sh%yref_2,sh%zp2)
endif
if (allocated(sh%seis)) deallocate(sh%seis)
end subroutine deallocate_pw2_3d_shot

subroutine deallocate_source(s)
type(source), intent(inout) :: s
deallocate(s%sou)
end subroutine deallocate_source

subroutine deallocate_km2d(k)
type(km2d), intent(inout) :: k
if (allocated(k%wfft)) then 
   deallocate(k%wfft)
endif
if (allocated(k%z0_mig)) then
   deallocate(k%z0_mig)
endif
end subroutine deallocate_km2d

subroutine deallocate_km3d(k)
type(km3d), intent(inout) :: k
if (allocated(k%wfft)) then 
   deallocate(k%wfft)
endif
if (allocated(k%z0_mig)) then
   deallocate(k%z0_mig)
endif
end subroutine deallocate_km3d

subroutine deallocate_shot2d(sh)
type(shot2d),  intent(inout) :: sh
call deallocate_and_free(sh%xs,sh%zs,sh%xg,sh%zg)
!if (allocated(sh%t)) then
!   call deallocate_and_free(sh%t)
!endif
if (allocated(sh%seis)) then
   call deallocate_and_free(sh%seis)
endif
end subroutine deallocate_shot2d
 
subroutine deallocate_shot3d(sh)
type(shot3d),  intent(inout) :: sh
call deallocate_and_free(sh%xs,sh%ys,sh%zs,sh%xg,sh%yg,sh%zg)
!if (allocated(sh%t)) then
!   call deallocate_and_free(sh%t)
!endif
if (allocated(sh%seis)) then
   call deallocate_and_free(sh%seis)
endif
end subroutine deallocate_shot3d

subroutine deallocate_model2d(m)
type(model2d), intent(inout) :: m
if (allocated(m%v)) then
   deallocate(m%v)
endif
if (allocated(m%refl)) then
   deallocate(m%refl)
endif
if (allocated(m%lambda)) then  ! modified by Bowen Guo 
   deallocate(m%lambda)        ! modified by Bowen Guo
endif                          ! modified by Bowen Guo
if (allocated(m%mu)) then      ! modified by Bowen Guo
   deallocate(m%mu)            ! modified by Bowen Guo
endif                          ! modified by Bowen Guo
if (allocated(m%vp)) then      ! modified by Bowen Guo
   deallocate(m%vp)            ! modified by Bowen Guo
endif                          ! modified by Bowen Guo
if (allocated(m%den)) then     
   deallocate(m%den)
endif
if (allocated(m%vs)) then
   deallocate(m%vs)
endif
end subroutine deallocate_model2d

subroutine deallocate_model3d(m)
type(model3d), intent(inout) :: m
if (allocated(m%v)) then
   deallocate(m%v)
endif
if (allocated(m%refl)) then
   deallocate(m%refl)
endif
if (allocated(m%den)) then
   deallocate(m%den)
endif
if (allocated(m%vs)) then
   deallocate(m%vs)
endif
if (allocated(m%vp)) then
   deallocate(m%vp)
endif
if (allocated(m%den)) then
   deallocate(m%den)
endif
if (allocated(m%lambda)) then
   deallocate(m%lambda)
endif
if (allocated(m%mu)) then
   deallocate(m%mu)
endif
end subroutine deallocate_model3d

subroutine deallocate_fd2d(f)
type(fd2d), intent(inout) :: f
if (allocated(f%fs)) then
   call deallocate_and_free(f%fs)
endif
end subroutine deallocate_fd2d

subroutine deallocate_fd3d(f)
type(fd3d), intent(inout) :: f
if (allocated(f%fs)) then
   call deallocate_and_free(f%fs)
endif
end subroutine deallocate_fd3d

subroutine deallocate_gdm2d(gdm)
type(gdm2d), intent(inout) :: gdm
if (allocated(gdm%wfft)) then
   deallocate(gdm%wfft)
endif
if (allocated(gdm%z0_mig)) then
   deallocate(gdm%z0_mig)
endif
end subroutine deallocate_gdm2d

subroutine deallocate_gdm3d(gdm)
type(gdm3d), intent(inout) :: gdm
if (allocated(gdm%wfft)) then
   deallocate(gdm%wfft)
endif
if (allocated(gdm%z0_mig)) then
   deallocate(gdm%z0_mig)
endif
end subroutine deallocate_gdm3d

subroutine deallocate_gdmod2d(gdm)
type(gdmod2d), intent(inout) :: gdm
call deallocate_gdm2d(gdm%g)
call deallocate_model2d(gdm%m)
call deallocate_fd2d(gdm%f)
call deallocate_source(gdm%s)
end subroutine deallocate_gdmod2d

subroutine deallocate_gdmod3d(gdm)
type(gdmod3d), intent(inout) :: gdm
call deallocate_gdm3d(gdm%g)
call deallocate_model3d(gdm%m)
call deallocate_fd3d(gdm%f)
call deallocate_source(gdm%s)
end subroutine deallocate_gdmod3d

subroutine deallocate_ref2d(r)
type(ref2d), intent(inout) :: r
call deallocate_and_free(r%fs)
call deallocate_and_free(r%xs,r%zs,r%xg,r%zg)
call deallocate_and_free(r%ttt)
call deallocate_and_free(r%v)
end subroutine deallocate_ref2d

subroutine deallocate_ref3d(r)
type(ref3d), intent(inout) :: r
call deallocate_and_free(r%fs)
call deallocate_and_free(r%xs,r%ys,r%zs,r%xg,r%yg,r%zg)
call deallocate_and_free(r%ttt)
call deallocate_and_free(r%v)
end subroutine deallocate_ref3d

subroutine deallocate_km2d_m2d_s(k,m,s)
type(km2d),      intent(inout) :: k
type(model2d),   intent(inout) :: m
type(source),    intent(inout) :: s
call deallocate_km2d(k)
call deallocate_model2d(m)
call deallocate_source(s)
end subroutine deallocate_km2d_m2d_s

subroutine deallocate_km3d_m3d_s(k,m,s)
type(km3d),      intent(inout) :: k
type(model3d),   intent(inout) :: m
type(source),    intent(inout) :: s
call deallocate_km3d(k)
call deallocate_model3d(m)
call deallocate_source(s)
end subroutine deallocate_km3d_m3d_s

subroutine deallocate_fd2d_m2d_s(f,m,s)
type(fd2d),      intent(inout) :: f
type(model2d),   intent(inout) :: m
type(source),    intent(inout) :: s
call deallocate_fd2d(f)
call deallocate_model2d(m)
call deallocate_source(s)
end subroutine deallocate_fd2d_m2d_s

subroutine deallocate_fd3d_m3d_s(f,m,s)
type(fd3d),      intent(inout) :: f
type(model3d),   intent(inout) :: m
type(source),    intent(inout) :: s
call deallocate_fd3d(f)
call deallocate_model3d(m)
call deallocate_source(s)
end subroutine deallocate_fd3d_m3d_s

subroutine deallocate_gdm2d_m2_s(g,m2,s)
type(gdm2d),     intent(inout) :: g
type(model2d),   intent(inout) :: m2
type(source),    intent(inout) :: s
call deallocate_gdm2d(g)
call deallocate_model2d(m2)
call deallocate_source(s)
end subroutine deallocate_gdm2d_m2_s

end module module_datatype
