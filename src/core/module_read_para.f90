module module_read_para
use module_parser
use module_datatype
use module_sf90_mpi
use module_utility
implicit none

interface readparamfile_subs
  
   module procedure readparamfile_sub_fname
   module procedure readparamfile_sub_inv
   module procedure readparamfile_sub_other
!   module procedure readparamfile_sub_trm
   module procedure readparamfile_sub_nsnr   ! Added by Lei
   module procedure readparamfile_sub_aperture
   module procedure readparamfile_sub_source
   module procedure readparamfile_sub_fd2d
   module procedure readparamfile_sub_fd3d
   module procedure readparamfile_sub_model2d
   module procedure readparamfile_sub_model3d
   module procedure readparamfile_sub_water_model2d
   module procedure readparamfile_sub_water_model3d
   module procedure readparamfile_sub_km2d
   module procedure readparamfile_sub_km3d
   module procedure readparamfile_sub_gdm2d
   module procedure readparamfile_sub_gdm3d
end interface readparamfile_subs

interface readparamfile
   module procedure readparamfile_gen       ! Added by Bowen Guo
   module procedure readparamfile_nm2d      ! Added by Lei Fu 2015
   module procedure readparamfile_a2d
   module procedure readparamfile_a3d
   module procedure readparamfile_eik2d
   module procedure readparamfile_eik3d
   module procedure readparamfile_km2d
   module procedure readparamfile_km3d
   module procedure readparamfile_lsrtm2d     ! modified by Bowen Guo
   module procedure readparamfile_gdm2d
   module procedure readparamfile_gdm3d
   module procedure readparamfile_lskm2d
   module procedure readparamfile_lskm3d
   module procedure readparamfile_fn_o
   module procedure readparamfile_fn_o_m
   module procedure readparamfile_fn_o_i_a_m2
   module procedure readparamfile_fn_o_i_a_m3
   module procedure readparamfile_fn_o_a ! Natural Mig by Lei Fu
end interface readparamfile

private readparamfile_subs

contains
!=======================================================
! Used for apps of csg2multi 1         ! Added by Bowen Guo
subroutine readparamfile_gen(parfile,fn,o,ap,s,m2)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(model2d),      intent(out) :: m2
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,m2)
end subroutine readparamfile_gen

!=====================================
! use for apps of nm2d: natural migration  ! Added by Lei Fu
subroutine readparamfile_nm2d(parfile,fn,o,n)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(nsnr),         intent(out)  :: n
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,n)

end subroutine readparamfile_nm2d

!=====================================
subroutine readparamfile_fn_o_a(parfile,fn,o,ap)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
end subroutine readparamfile_fn_o_a
!==================================================================
! use for apps of a2d_modeling/modeling_mute/born_modeling/_rtm
subroutine readparamfile_a2d(parfile,fn,o,ap,s,f2,m2)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(fd2d),         intent(out) :: f2
type(model2d),      intent(out) :: m2
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,f2)
call readparamfile_subs(parfile,m2)
end subroutine readparamfile_a2d

!=====================================
! use for apps of a3d_modeling/modeling_mute/born_modeling/_rtm
subroutine readparamfile_a3d(parfile,fn,o,ap,s,f3,m3)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(fd3d),         intent(out) :: f3
type(model3d),      intent(out) :: m3
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,f3)
call readparamfile_subs(parfile,m3)
end subroutine readparamfile_a3d
!=====================================
! use for apps of eik2d
subroutine readparamfile_eik2d(parfile,fn,o,ap,m2)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(model2d),      intent(out) :: m2
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,m2)
end subroutine readparamfile_eik2d

!=====================================
! use for apps of eik3d
subroutine readparamfile_eik3d(parfile,fn,o,ap,m3)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(model3d),      intent(out) :: m3
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,m3)
end subroutine readparamfile_eik3d

!=====================================
! use for apps of km2d_modeling/km2d_migration
subroutine readparamfile_km2d(parfile,fn,o,ap,s,k2,m2)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(km2d),         intent(out) :: k2
type(model2d),      intent(out) :: m2
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,k2)
call readparamfile_subs(parfile,m2)
end subroutine readparamfile_km2d

!=====================================
! use for apps of km3d_modeling/km3d_migration
subroutine readparamfile_km3d(parfile,fn,o,ap,s,k3,m3)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(km3d),         intent(out) :: k3
type(model3d),      intent(out) :: m3
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,k3)
call readparamfile_subs(parfile,m3)
end subroutine readparamfile_km3d

!=====================================
! use for apps of gdm2d_modeling/gdm2d_migration
subroutine readparamfile_gdm2d(parfile,fn,o,ap,s,g2,m2)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(gdm2d),        intent(out) :: g2
type(model2d),      intent(out) :: m2
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,g2)
call readparamfile_subs(parfile,m2)
end subroutine readparamfile_gdm2d

!=====================================
! use for apps of gdm3d_modeling/gdm3d_migration
subroutine readparamfile_gdm3d(parfile,fn,o,ap,s,g3,m3)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(gdm3d),        intent(out) :: g3
type(model3d),      intent(out) :: m3
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,g3)
call readparamfile_subs(parfile,m3)
end subroutine readparamfile_gdm3d

!=====================================
! use for apps of lskm2d
subroutine readparamfile_lskm2d(parfile,fn,o,i,ap,s,k2,m2)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(inv),          intent(out) :: i
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(km2d),         intent(out) :: k2
type(model2d),      intent(out) :: m2
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,i)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,k2)
call readparamfile_subs(parfile,m2)
end subroutine readparamfile_lskm2d

!=====================================
! use for apps of km3d_modeling/km3d_migration
subroutine readparamfile_lskm3d(parfile,fn,o,i,ap,s,k3,m3)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(inv),          intent(out) :: i
type(aperture),     intent(out) :: ap
type(source),       intent(out) :: s
type(km3d),         intent(out) :: k3
type(model3d),      intent(out) :: m3
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,i)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,k3)
call readparamfile_subs(parfile,m3)
end subroutine readparamfile_lskm3d

!=====================================
! modified by Bowen Guo
subroutine readparamfile_lsrtm2d(parfile,fn,o,i,ap,s,f2,m2)
character(len=*), intent(in)     :: parfile
type(fname),      intent(out)    :: fn
type(other),      intent(out)    :: o
type(inv),        intent(out)    :: i
type(aperture),   intent(out)    :: ap
type(source),     intent(out)    :: s
type(fd2d),       intent(out)    :: f2
type(model2d),    intent(out)    :: m2
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,i)
call readparamfile_subs(parfile,ap)
call readparamfile_subs(parfile,s)
call readparamfile_subs(parfile,f2)
call readparamfile_subs(parfile,m2)
end subroutine readparamfile_lsrtm2d
!==============================================
subroutine readparamfile_fn_o(parfile,fn,o)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
end subroutine readparamfile_fn_o

!=====================================
! use for apps of pw1_generate_2d
subroutine readparamfile_fn_o_m(parfile,fn,o,m)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(model2d),      intent(out) :: m
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,m)
end subroutine readparamfile_fn_o_m

!=====================================
! use for 2d ref traveltime tomography
subroutine readparamfile_fn_o_i_a_m2(parfile,fn,o,i,a,m)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(inv),          intent(out) :: i
type(aperture),     intent(out) :: a
type(model2d),      intent(out) :: m
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,i)
call readparamfile_subs(parfile,a)
call readparamfile_subs(parfile,m)
end subroutine readparamfile_fn_o_i_a_m2
!=====================================

!=====================================
! use for 3d ref traveltime tomography
subroutine readparamfile_fn_o_i_a_m3(parfile,fn,o,i,a,m)
character(len=*),   intent(in) :: parfile
type(fname),        intent(out) :: fn
type(other),        intent(out) :: o
type(inv),          intent(out) :: i
type(aperture),     intent(out) :: a
type(model3d),      intent(out) :: m
call readparamfile_subs(parfile,fn)
call readparamfile_subs(parfile,o)
call readparamfile_subs(parfile,i)
call readparamfile_subs(parfile,a)
call readparamfile_subs(parfile,m)
end subroutine readparamfile_fn_o_i_a_m3
!=====================================

subroutine readparamfile_sub_fname(parfile,fn)
character(len=*),    intent(in)  :: parfile
type(fname),         intent(out) :: fn
if (rank.eq.0) then
   ! Only for Master processor
   call readParFile(parfile, "CSG_W",           fn%csg_w,          'n/a')
   call readParFile(parfile, "CSG0_W",           fn%csg0_w,          'n/a')
   call readParFile(parfile, "VELS_FILE",           fn%vels_file,          'n/a')
   call readParFile(parfile, "VEL_WATER_FILE_",     fn%vel_water_file,     'n/a')
   call readParFile(parfile, "VEL_WHOME_FILE_",     fn%vel_whomo_file,     'n/a')
   call readParFile(parfile, "SOURCE_FILE",         fn%source_file, 'source.bin')
   call readParFile(parfile, "TOPO_FILE",           fn%topo_file,          'n/a') 
   call readParFile(parfile, "REFL_OB_FILE",        fn%refl_ob_file,       'n/a')
   call readParFile(parfile, "REFL_FS_FILE",        fn%refl_fs_file,       'n/a')
   call readParFile(parfile, "RES_FILE",            fn%res_file,       'res.bin')
   call readParFile(parfile, "RES_TOT_FILE",        fn%res_tot_file,   'res_tot.bin')! prelsm
   call readParFile(parfile, "PW1_PAR_FILE",        fn%pw1_par_file,       'n/a')
   call readParFile(parfile, "PW2_PAR_FILE",        fn%pw2_par_file,       'n/a')
   call readParFile(parfile, "SMWIN_FILE",          fn%smwin_file,         'n/a')

   ! For all processors
   ! Optional String parameters
   call readParFile(parfile, "VEL_FILE",          fn%vel_file,           'n/a')
   call readParFile(parfile, "TT_OBS_FILE",        fn%tt_obs_file,      'n/a') ! Lei
   call readParFile(parfile, "MASK_FILE",          fn%mask_file,           'n/a') ! Lei
   call readParFile(parfile, "MASK_SALT_FILE",     fn%mask_salt_file,      'n/a') ! Lei
   call readParFile(parfile, "MASK_SHOT_FILE",     fn%mask_shot_file,      'n/a') ! Lei
   call readParFile(parfile, "REFL_FILE",         fn%refl_file,          'n/a')
   call readParFile(parfile, "DEN_FILE",          fn%den_file,           'n/a')
   call readParFile(parfile, "VEL_HOMO_FILE",     fn%vel_homo_file,      'n/a')
   call readParFile(parfile, "LAMBDA_FILE",       fn%lambda_file,         'n/a')
   call readParFile(parfile, "MU_FILE",           fn%mu_file,            'n/a')
   call readParFile(parfile, "VP_FILE",           fn%vp_file,            'n/a')
   call readParFile(parfile, "VS_FILE",           fn%vs_file,            'n/a')
   call readParFile(parfile, "VP_HOMO_FILE",      fn%vp_homo_file,            'n/a')
   call readParFile(parfile, "VS_HOMO_FILE",      fn%vs_homo_file,            'n/a')
   call readParFile(parfile, "GREEN_TOTAL",      fn%green_total,            'n/a')
   call readParFile(parfile, "GREEN_D",      fn%green_d,            'n/a')
   call readParFile(parfile, "GREEN_BS",      fn%green_bs,            'n/a')
   call readParFile(parfile, "NM_IMAGE",      fn%nm_image,            'n/a')
   ! modified by Bowen Guo
   call readParfile(parfile, "DEN_HOMO_FILE",     fn%den_homo_file,      'n/a')
   call readParFile(parfile, "MIG0_FILE",         fn%mig0_file,          'n/a') 
   call readParFile(parfile, "MIG0_PP_FILE",      fn%mig0_pp_file,    'n/a') ! Added by Bowen Guo 
   call readParFile(parfile, "MIG0_PS_FILE",      fn%mig0_ps_file,     'n/a') ! Added by Bowen Guo
   call readParFile(parfile, "Z0_MIG_FILE",       fn%z0_mig_file,        'n/a')
   call readParFile(parfile, "CSG_OUT",           fn%csg_out,           "csg_")
   call readParFile(parfile, "CSG_TMP",           fn%csg_tmp,           "csgtmp_")
   call readParFile(parfile, "CSG_CAL",           fn%csg_cal,           "csg_")
   call readParFile(parfile, "CSG_OUT_U",         fn%csg_out_u,          'n/a')! modified by Bowen Guo
   call readParFile(parfile, "CSG_OUT_V",         fn%csg_out_v,          'n/a')! modified by Bowen Guo
   call readParFile(parfile, "CSG_OUT_W",         fn%csg_out_w,          'n/a')! modified by Bowen Guo
   call readParFile(parfile, "CSG_OUT_P",         fn%csg_out_p,          'n/a')! modified by Bowen Guo
   call readParFile(parfile, "CSG_DIR",         fn%csg_dir,          'n/a')! modified by Bowen Guo
   call readParFile(parfile, "CSG_MASK",         fn%csg_mask,          'n/a')! modified by Bowen Guo
   call readParFile(parfile, "MASK_CSG",         fn%mask_csg,          'n/a')! modified by Bowen Guo

   call readParFile(parfile, "CSG_HOMO",          fn%csg_homo,     "csg_homo_")
   call readParFile(parfile, "CSG_RECONS",        fn%csg_recons, "csg_recons_")
   call readParFile(parfile, "CSG_RECONS_PS",     fn%csg_recons_ps, "csg_recons_ps")!LSMF
   call readParFile(parfile, "SEIS_IN",           fn%seis_in,           "csg_")
   call readParFile(parfile, "SEIS_OUT",          fn%seis_out,          "csg_")
   call readParFile(parfile, "PRE_IMAGE_OUT",     fn%pre_img_out,       "img_")
   call readParFile(parfile, "PRE_IMAGE_OUT_PS",  fn%pre_img_out_ps,     "img_")! LSMF
   call readParFile(parfile, "PRE_ILLUM_OUT",     fn%pre_illum_out,   "illum_")
   call readParFile(parfile, "PRE_ILLUM_IMG_OUT", fn%pre_illum_img_out,"img_illum_")
   call readParFile(parfile, "TTT_IN",            fn%ttt_in,            "ttt_")
   call readParFile(parfile, "TTT_OUT",           fn%ttt_out,           "ttt_")
   call readParFile(parfile, "STTT_IN",           fn%sttt_in,          "sttt_")
   call readParFile(parfile, "GTTT_IN",           fn%gttt_in,          "gttt_")
   call readParFile(parfile, "STTT_S_IN",         fn%sttt_s_in,        "sttt_s_")! Used in LSMF
   call readParFile(parfile, "GTTT_S_IN",         fn%gttt_s_in,        "gttt_s_")! Used in LSMF

   call readParFile(parfile, "PW1TTT_IN",         fn%pw1ttt_in,      "pw1ttt_")
   call readParFile(parfile, "PW2STTT_IN",        fn%pw2sttt_in,    "pw2sttt_")
   call readParFile(parfile, "PW2GTTT_IN",        fn%pw2gttt_in,   "pw2g:ttt_")
   call readParFile(parfile, "IMAGE_FILE",        fn%img_file,       "img.bin")
   call readParFile(parfile, "IMAGE_PS_FILE",     fn%img_ps_file,    "img_ps.bin") ! LSMF
   call readParFile(parfile, "IMAGE_ILLUM_FILE",  fn%img_illum_file,"img_illum.bin")
   call readParFile(parfile, "COORD_IN_FILE",     fn%coord_in_file,      "n/a")
   call readParFile(parfile, "COORD_OUT_FILE",    fn%coord_out_file,     "n/a")
   call readParFile(parfile, "COORD_ALL_FILE",    fn%coord_all_file,     "n/a")
   call readParFile(parfile, "PW1_COORD_ALL_FILE",fn%pw1_coord_all_file, "n/a")
   call readParFile(parfile, "PW2_COORD_ALL_FILE",fn%pw2_coord_all_file, "n/a")
   call readParFile(parfile, "PW1_COORD_FILE",    fn%pw1_coord_file,     "n/a")
   call readParFile(parfile, "PW2_COORD_FILE",    fn%pw2_coord_file,     "n/a")
   call readParFile(parfile, "VEL0_FILE",         fn%vel0_file,          "n/a")
   call readParFile(parfile, "GK1_FILE",          fn%gk1_file,          "gk1_")
   call readParFile(parfile, "PREC_GK1_FILE",     fn%prec_gk1_file,"prec_gk1_")
   call readParFile(parfile, "DK1_FILE",          fn%Dk1_file,          "dk1_")
   call readParFile(parfile, "WF_FILE",           fn%wf_file,            "wf_")
   call readParFile(parfile, "GF_FILE",           fn%gf_file,            "gf_")
   call readParFile(parfile, "LOG_FILE",          fn%log_file,   "working.log")
   call readParFile(parfile, "GREEN_OUT",         fn%green_out,          'n/a') ! Marchenko 
   call readParFile(parfile, "F0M_OUT",           fn%f0m_out,            'n/a') ! Marchenko
   call readParFile(parfile, "GREEN_SHOT",        fn%green_shot,         'n/a') ! Marchenko
   call readParFile(parfile, "REPATH",            fn%repath,              "n/a") ! csg2multi1
endif
call MPI_BCAST(fn%source_file,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%green_total,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! NM
call MPI_BCAST(fn%green_d,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)     ! NM
call MPI_BCAST(fn%green_bs,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)    ! NM
call MPI_BCAST(fn%vel_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%tt_obs_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%mask_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) !Lei
call MPI_BCAST(fn%mask_shot_file,100, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) !Lei
call MPI_BCAST(fn%topo_file,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%refl_file,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%den_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%vel_homo_file,100, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%mig0_file,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%z0_mig_file,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%csg_out_u,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! added by Bowen Guo
call MPI_BCAST(fn%csg_out_v,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! added by Bowen Guo
call MPI_BCAST(fn%csg_out_w,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! added by Bowen Guo
call MPI_BCAST(fn%csg_out_p,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! added by Bowen Guo
call MPI_BCAST(fn%csg_out,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%csg_dir,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%csg_mask,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%mask_csg,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%csg_tmp,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%csg_cal,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%csg_homo,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%csg_recons,100,    MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%csg_recons_ps,100, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! LSMF
call MPI_BCAST(fn%seis_in,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%seis_out,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%ttt_in,100,        MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%ttt_out,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%sttt_in,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%gttt_in,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%sttt_s_in,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! added by Bowen Guo
call MPI_BCAST(fn%gttt_s_in,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! added by Bowen Guo
call MPI_BCAST(fn%pw1ttt_in,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pw2sttt_in,100,    MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pw2gttt_in,100,    MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%img_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%img_ps_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! LSMF
call MPI_BCAST(fn%img_illum_file,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pre_img_out,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pre_img_out_ps,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) !LSMF
call MPI_BCAST(fn%pre_illum_out,100, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pre_illum_img_out,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%coord_in_file,100, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%coord_out_file,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%coord_all_file,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pw1_coord_all_file,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pw2_coord_all_file,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pw1_coord_file,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%pw2_coord_file,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%vel0_file,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%gk1_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%prec_gk1_file,100, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%dk1_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%wf_file,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%gf_file,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%log_file,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(fn%green_out,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! Marchenko
call MPI_BCAST(fn%f0m_out,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! Marchenko
call MPI_BCAST(fn%green_shot,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! Marchenko
call MPI_BCAST(fn%repath,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! csg2multi1
end subroutine readparamfile_sub_fname

!=====================================

subroutine readparamfile_sub_inv(parfile, i)
character(len=*),      intent(in)  :: parfile
type(inv),             intent(out) :: i
! Logical parameters
if (rank.eq.0) then
   call readParFile(parfile, "IsMig0",                i%isMig0,         .false.)
   call readParFile(parfile, "IsNumericalLineSearch", i%isNMLS,         .false.)
   call readParFile(parfile, "IsResetCGDirection",    i%isReset,        .false.)
   call readParFile(parfile, "IsSeisRecons",          i%isSeisRecons,   .false.)
   call readParFile(parfile, "IsPreCondition",        i%isPreCondition, .false.)
   call readParFile(parfile, "IsBackTrack",           i%isBackTrack,    .false.)
   call readParFile(parfile, "IsAdjointTest",         i%isAdjointTest,    .false.)  ! Added by Bowen Guo

endif
call MPI_BCAST(i%isMig0,1,          MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%isNMLS,1,          MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%isReset,1,         MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%isSeisRecons,1,    MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%isPreCondition,1,  MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%isBackTrack,1,     MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%isAdjointTest,1,     MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)   ! Added by Bowen Guo
! Character parameters
if (rank.eq.0) then
   call readParFile(parfile, "INV_TYPE",              i%inv_type,       "CG")
   call readParFile(parfile, "CG_TYPE",               i%cg_type,"Polak-Rebiere")
endif
call MPI_BCAST(i%inv_type,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%cg_type,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
! Integer parameters
if (rank.eq.0) then
   call readParFile(parfile, "NF",                   i%nf,                4)
   call readParFile(parfile, "NMAX",                 i%nmax,               2)
   call readParFile(parfile, "NIT",                   i%nit,                -1)
   call readParFile(parfile, "RESET",                 i%reset,              -1)
   call readParFile(parfile, "NEAROFFSET_CUT",                 i%nearoffset_cut,              10)
   call readParFile(parfile, "EWI_TYPE",                   i%ewi_type,                1)
   call readParFile(parfile, "NO_OFSET",                   i%no_ofset,                10)
endif
call MPI_BCAST(i%nf,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%nmax,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%nit,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%reset,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%nearoffset_cut,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%ewi_type,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%no_ofset,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! Real parameters
if (rank.eq.0) then
   call readParFile(parfile, "STL_COE",               i%stl_coe,          0.05)
   call readParFile(parfile, "VMIN",                  i%vmin,             0.05)
   call readParFile(parfile, "VMAX",                  i%vmax,             0.05)
   call readParFile(parfile, "FL1",                  i%fl1,             1.0)
   call readParFile(parfile, "FL2",                  i%fl2,             1.6)
   call readParFile(parfile, "FL3",                  i%fl3,             2.5)
   call readParFile(parfile, "FL4",                  i%fl4,             4.0)
   call readParFile(parfile, "FL5",                  i%fl5,             6.3)
   call readParFile(parfile, "FL6",                  i%fl6,             9.8)
   call readParFile(parfile, "FH1",                  i%fh1,             2.1)
   call readParFile(parfile, "FH2",                  i%fh2,             3.3)
   call readParFile(parfile, "FH3",                  i%fh3,             5.4)
   call readParFile(parfile, "FH4",                  i%fh4,             8.5)
   call readParFile(parfile, "FH5",                  i%fh5,             13.1)
   call readParFile(parfile, "FH6",                  i%fh6,             25.1)
   call readParFile(parfile, "OFSET_STEP",           i%ofset_step,      250.0)
   call readParFile(parfile, "MAX_OFSET",            i%max_ofset,      5000.0)
endif
call MPI_BCAST(i%stl_coe,1,         MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%vmin,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%vmax,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fl1,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fl2,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fl3,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fl4,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fl5,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fl6,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fh1,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fh2,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fh3,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fh4,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fh5,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%fh6,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%ofset_step,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(i%max_ofset,1,            MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_inv

!============================================
subroutine readparamfile_sub_nsnr(parfile,n)
character(len=*),      intent(in)  :: parfile
type(nsnr),           intent(out) :: n
! Integer Parameter
if (rank.eq.0) then
   call readParFile(parfile, "NSHOT",    n%nshot,          -1)
   call readParFile(parfile, "NREC",     n%nrec,          -1)

endif
call MPI_BCAST(n%nshot,     1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(n%nrec,     1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_nsnr

!========================================================================
subroutine readparamfile_sub_other(parfile,o)
character(len=*),      intent(in)  :: parfile
type(other),           intent(out) :: o
! Character parameter
if (rank.eq.0) then
   call readParFile(parfile, "PW1_SIDE",        o%pw1_side,          "SOURCE")
   call readParFile(parfile, "SORT_IN",         o%sort_in,           "CSG")
   call readParFile(parfile, "SORT_OUT",        o%sort_out,          "CRG")
   call readParFile(parfile, "EIK_TYPE",        o%eik_type,          "DIR")
   call readParFile(parfile, "COORDFILE_TYPE",  o%coordfile_type,    "ASCII")
   call readParFile(parfile, "DATA_FORMAT",     o%dfor,              "DEFAULT")
endif
call MPI_BCAST(o%pw1_side,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%sort_in,100,       MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%sort_out,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%eik_type,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%coordfile_type,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%dfor,100,          MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
! Logical parameter
if (rank.eq.0) then
   call readParFile(parfile, "IsIllum",               o%isIllum,        .false.)
   call readParFile(parfile, "IsMask",                o%isMask,        .false.)
   call readParFile(parfile, "IsMask_Salt",           o%isMask_Salt,    .false.)
   call readParFile(parfile, "IsIllumReceiver",       o%isIllumReceiver,.false.)
   call readParFile(parfile, "IsSaveIllum",           o%isSaveIllum,    .false.)
   call readParFile(parfile, "IsPreStackImg",         o%isPreStackImg,  .false.)
   call readParFile(parfile, "IsEikSep",              o%isEikSep,       .true.)
   call readParFile(parfile, "IsShotConti",           o%isShotConti,    .true.)
   call readParFile(parfile, "IsSortData",            o%isSortData,     .true.)
   call readParFile(parfile, "IsPreCondition",        o%isPreCondition, .false.)
   call readParFile(parfile, "IsDepthPC",                 o%isDepthPC, .false.)
   call readParFile(parfile, "IsSmooth_gk",              o%isSmooth_gk, .false.)
   call readParFile(parfile, "IsCut",                          o%isCut, .false.)
   call readParFile(parfile, "IsRemove_dir",            o%isRemove_dir, .false.)
   call readParFile(parfile, "IsTaper_dir",            o%isTaper_dir, .false.)
   call readParFile(parfile, "IsMask_csg",                o%isMask_csg, .false.)
   call readParFile(parfile, "IsWindow_csg",            o%isWindow_csg, .false.)
   call readParFile(parfile, "IsCut_far_trace",        o%isCut_far_trace, .false.)
   call readParFile(parfile, "IS_CSG_FILT",              o%is_csg_filt, .false.)
   call readParFile(parfile, "GRAD_FLIP",                  o%grad_flip, .false.)
 !  call readParFile(parfile, "IS_OUTPUT_PRE",          o%is_output_pre, .false.)
 !  call readParFile(parfile, "IS_OUTPUT_OBS",          o%is_output_obs, .false.)
 !  call readParFile(parfile, "IS_OUTPUT_RESIDUAL", o%is_output_residual, .false.)
endif
call MPI_BCAST(o%isIllum,1,          MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isMask,1,           MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isMask_Salt,1,      MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isIllumReceiver,1,  MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isSaveIllum,1,      MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isPreStackImg,1,    MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isEikSep,1,         MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isShotConti,1,      MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isSortData,1,       MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isPreCondition,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isDepthPC,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isSmooth_gk,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isCut,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isRemove_dir,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isTaper_dir,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isMask_csg,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isWindow_csg,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%isCut_far_trace,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%is_csg_filt,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%grad_flip,1,     MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!call MPI_BCAST(o%is_output_pre,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!call MPI_BCAST(o%is_output_obs,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!call MPI_BCAST(o%is_output_residual,1,   MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
! Integer Parameter
if (rank.eq.0) then
   call readParFile(parfile, "NT",     o%nt,          -1)
   call readParFile(parfile, "NIT0",     o%nit0,      10)
   call readParFile(parfile, "NT0",     o%nt0,      2)
   call readParFile(parfile, "N0",     o%n0,      20)
   call readParFile(parfile, "T_P",     o%T_p,      300)
   call readParFile(parfile, "NT_IN",  o%nt_in,     o%nt)
   call readParFile(parfile, "NT_OUT", o%nt_out,    o%nt)
   call readParFile(parfile, "NS_TAPER",o%ns_taper,    0)
   call readParFile(parfile, "TMAX",o%Tmax,   200) ! Lei 
   call readParFile(parfile, "WIN",o%win,   200) ! Lei 
   call readParFile(parfile, "NH",o%nh,   51) ! Lei 
   call readParFile(parfile, "LAGMAX",o%lagmax,   100) ! Lei 
   call readParFile(parfile, "NLAG",o%nlag,   21) ! Lei 
   call readParFile(parfile, "WT_RES",o%wt_res,   2) ! Lei 
   call readParFile(parfile, "NCUT",o%ncut,   0) ! Lei 
   call readParFile(parfile, "NCUT2",o%ncut2,   0) ! Lei 
   call readParFile(parfile, "HW",o%hw,   5) ! Lei 
   call readParFile(parfile, "HW_M",o%hw_m,   7) ! Lei 
   call readParFile(parfile, "HW_T",o%hw_t,   19) ! Lei 
   call readParFile(parfile, "N_TAPER",o%n_taper,   100) ! Lei 
   call readParFile(parfile, "IS_CONV_SOU",o%is_conv_sou,   0) ! Lei 
   call readParFile(parfile, "USE_AMP",o%use_amp,   0) ! Lei 
endif
call MPI_BCAST(o%nt,     1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%nit0,   1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%nt0,   1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%n0,   1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%T_p,   1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%nt_in,  1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%nt_out, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%ns_taper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%win,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%nh,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%Tmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%lagmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%nlag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%wt_res,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%ncut,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%ncut2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%hw,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%hw_m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%hw_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%n_taper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%is_conv_sou,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%use_amp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! Real Parameter
if (rank.eq.0) then
   call readParFile(parfile, "FLOW",       o%fl,             0.2)
   call readParFile(parfile, "FHIGH",       o%fh,             20.0)
   call readParFile(parfile, "FPK",       o%fpk,             1.5)
   call readParFile(parfile, "N_PERIOD",  o%N_period,        1.5)
   call readParFile(parfile, "A_ERR",     o%A_err,        0.0001)
   call readParFile(parfile, "DT",     o%dt,        -1.0)
   call readParFile(parfile, "TAUMAX",     o%taumax,        0.5)
   call readParFile(parfile, "TOLER",     o%toler,        0.02)
   call readParFile(parfile, "DT_IN",  o%dt_in,     o%dt)
   call readParFile(parfile, "DT_OUT", o%dt_out,    o%dt)
   call readParFile(parfile, "Illum_Order",   o%illum_order,        1.0)
   call readParFile(parfile, "SORT_TOLERANCE",o%sort_tol,           5.0)
   call readParFile(parfile, "XS_TAPER",      o%xs_taper,           0.0)
   call readParFile(parfile, "YS_TAPER",      o%ys_taper,           0.0)
   call readParFile(parfile, "DATA_FORMAT_COE",o%dfor_coe,          0.0)
   call readParFile(parfile, "CUT_DISTANCE",o%cut_distance,         0.0)
   call readParFile(parfile, "T_WINDOW",        o%T_window,         1.5)
   call readParFile(parfile, "NEAR_CUT",        o%near_cut,         0.0)
   call readParFile(parfile, "FAR_CUT",         o%far_cut,         0.0)
endif

call MPI_BCAST(o%fl,               1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%fh,               1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%fpk,               1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%N_period,          1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%A_err,          1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%dt,          1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%taumax,       1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%toler,       1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%dt_in,       1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%dt_out,      1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%illum_order, 1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%sort_tol,    1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%xs_taper,    1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%ys_taper,    1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%dfor_coe,    1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%cut_distance,1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%T_window,    1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%near_cut,    1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(o%far_cut,     1, MPI_REAL,   0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_other
   
!============================================

subroutine readparamfile_sub_aperture(parfile,a)
character(len=*),      intent(in)  :: parfile
type(aperture),        intent(out) :: a
if (rank.eq.0) then
   call readParFile(parfile, "IsAper",a%isAper,     .false.)
   call readParFile(parfile, "X0",    a%x0,         -1.0)
   call readParFile(parfile, "Y0",    a%y0,         -1.0)
   call readParFile(parfile, "Z0",    a%z0,         -1.0)
   call readParFile(parfile, "APER_X",a%aper_x,      0.0)
   call readParFile(parfile, "APER_Y",a%aper_y, a%aper_x)
   call readParFile(parfile, "APER_Z",a%aper_z,      0.0)
endif
call MPI_BCAST(a%isAper, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(a%x0,     1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(a%y0,     1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(a%z0,     1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(a%aper_x, 1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(a%aper_y, 1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(a%aper_z, 1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_aperture   

!============================================

subroutine readparamfile_sub_source(parfile,s)
character(len=*),   intent(in)  :: parfile
type(source),       intent(out) :: s
if (rank.eq.0) then
   call readParFile(parfile, "IsSource",    s%isSource,   .false.)
   call readParFile(parfile, "NW",          s%nw,         -1)
   call readParFile(parfile, "DTW",         s%dt,         -1.0)
   call readParFile(parfile, "FREQ",        s%freq,       -1.0)
endif   
call MPI_BCAST(s%nw,1,        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(s%dt,1,        MPI_REAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(s%isSource,1,  MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(s%freq,1,      MPI_REAL,   0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_source   

!============================================

subroutine readparamfile_sub_fd2d(parfile,f)
character(len=*),   intent(in)  :: parfile
type(fd2d),         intent(out) :: f
! String parameters
if (rank.eq.0) then
   call readParFile(parfile, "FD_TYPE",        f%fd_type,          "2ND")
endif
call MPI_BCAST(f%fd_type,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
! Logical parameters
if (rank.eq.0) then
   call readParFile(parfile, "IsFS",                  f%isFS,           .false.)
   call readParFile(parfile, "IsDipole",              f%isDipole,       .false.)
   call readParFile(parfile, "IsSaveBC",              f%isSaveBC,       .true.)
   call readParFile(parfile, "IsSaveWF",              f%isSaveWF,       .false.)
   call readParFile(parfile, "IsWriteData",           f%isWriteData,    .false.)
   call readParFile(parfile, "IsReadData",            f%isReadData,     .false.)
   if (f%isSaveBC) then
      if (f%isSaveWF) then
         call message(0,"IsSaveBC and IsSaveWF cannot be true at the same time")
         stop
      endif
   endif
endif
call MPI_BCAST(f%isFS,1,            MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isDipole,1,        MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isSaveWF,1,        MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isSaveBC,1,        MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isWriteData,1,     MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isReadData,1,      MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%sourcestress,1,    MPI_LOGIcal,0,MPI_COMM_WORLD,ierr)   ! modified by Bowen Guo
! Integer/Character parameters
if (rank.eq.0) then
   call readParFile(parfile, 'NT',               f%nt,                 -1)
   call readParFile(parfile, 'NPML',             f%npml,               40)
  ! call readParFile(parfile, 'IMAGE_CONDITION',  f%ic,                  2)
   call readParFile(parfile, 'IMAGE_CONDITION',  f%ic,                  1) !modified by Bowen Guo
   call readParFile(parfile, 'FD_ORDER',         f%fd_order,           28)
   call readParFile(parfile, 'BC_TYPE',          f%bc_type,             1)
   call readParFile(parfile, "SOURCE_TYPE",      f%source_type,         "w" ) ! modified by Bowen Guo
   call readParFile(parfile, "DATA_TYPE",       f%data_type,         "w") ! modified by Bowen Guo
endif
call MPI_BCAST(f%nt,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%npml,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%ic,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%fd_order,1,       MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%bc_type,1,        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%source_type,1,      MPI_INTEGER,0,MPI_COMM_WORLD,ierr) ! modified by Bowen Guo
call MPI_BCAST(f%data_type,1,     MPI_INTEGER,0,MPI_COMM_WORLD,ierr) ! modified by Bowen Guo
! Real parameters
if (rank.eq.0) then
   call readParFile(parfile, 'DT',               f%dt,             -1.0)
endif
call MPI_BCAST(f%data_type,100,        MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%source_type,100,        MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%dt,1,        MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_fd2d

!============================================

subroutine readparamfile_sub_fd3d(parfile,f)
character(len=*),   intent(in)  :: parfile
type(fd3d),         intent(out) :: f
! String parameters
if (rank.eq.0) then
   call readParFile(parfile, "FD_TYPE",        f%fd_type,          "2ND")
   call readParFile(parfile, "SOURCE_TYPE",    f%source_type,      "w")
   call readParFile(parfile, "DATA_TYPE",      f%data_type,        "w")
endif
call MPI_BCAST(f%fd_type,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%source_type,100,  MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%data_type,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
! Logical parameters
if (rank.eq.0) then
   call readParFile(parfile, "IsFS",                  f%isFS,           .false.)
   call readParFile(parfile, "IsDipole",              f%isDipole,       .false.)
   call readParFile(parfile, "IsSaveBC",              f%isSaveBC,       .false.)
   call readParFile(parfile, "IsSaveWF",              f%isSaveWF,       .false.)
   call readParFile(parfile, "IsWriteData",           f%isWriteData,    .true.)
   call readParFile(parfile, "IsReadData",            f%isReadData,     .true.)
   if (f%isSaveBC) then
      if (f%isSaveWF) then
         call message(0,"IsSaveBC and IsSaveWF cannot be true at the same time")
         stop
      endif
   endif
endif
call MPI_BCAST(f%isFS,1,            MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isDipole,1,        MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isSaveWF,1,        MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isSaveBC,1,        MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isWriteData,1,     MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%isReadData,1,      MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
! Integer parameters
if (rank.eq.0) then
   call readParFile(parfile, 'NT',               f%nt,                 -1)
   call readParFile(parfile, 'NPML',             f%npml,               40)
   call readParFile(parfile, 'IMAGE_CONDITION',  f%ic,                  2)
   call readParFile(parfile, 'FD_ORDER',         f%fd_order,           28)
   call readParFile(parfile, 'BC_TYPE',          f%bc_type,             1)
endif
call MPI_BCAST(f%nt,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%npml,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%ic,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%fd_order,1,       MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(f%bc_type,1,        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! Real parameters
if (rank.eq.0) then
   call readParFile(parfile, 'DT',               f%dt,             -1.0)
endif
call MPI_BCAST(f%dt,1,        MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_fd3d

!============================================

subroutine readparamfile_sub_model2d(parfile,m)
character(len=*),      intent(in)  :: parfile
type(model2d),         intent(out) :: m
if (rank.eq.0) then
   ! Integer parameters
   call readParFile(parfile, 'NX',               m%nx,                 -1)
   call readParFile(parfile, 'NZ',               m%nz,                 -1)
   call readParfile(parfile, 'NX_S',             m%nx_s,                1 ) ! modified by Bowen Guo
   call readParfile(parfile, 'NZ_S',             m%nz_s,                1 ) ! modified by Bowen Guo
   call readParfile(parfile, 'NX_E',             m%nx_e,                m%nx)! modified by Bowen Guo
   call readParfile(parfile, 'NZ_E',             m%nz_e,                m%nz)! modified by Bowen Guo
   ! real parameters
   call readParFile(parfile, 'DX',               m%dx,                 -1.0)
   call readParFile(parfile, 'DZ',               m%dz,                 m%dx)
   call readParFile(parfile, 'VMIN',             m%vmin,               -1.0)
   call readParFile(parfile, 'VMAX',             m%vmax,               -1.0)
endif
call MPI_BCAST(m%nx,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%nz,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%nx_s,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%nx_e,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%nz_s,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%nz_e,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%dx,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%dz,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%vmin,1,         MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%vmax,1,         MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_model2d

!============================================

subroutine readparamfile_sub_model3d(parfile,m)
character(len=*),      intent(in)  :: parfile
type(model3d),         intent(out) :: m
if (rank.eq.0) then
   ! Integer parameters
   call readParFile(parfile, 'NX',               m%nx,                 -1)
   call readParFile(parfile, 'NY',               m%ny,                 -1)
   call readParFile(parfile, 'NZ',               m%nz,                 -1)
   ! real parameters
   call readParFile(parfile, 'DX',               m%dx,                 -1.0)
   call readParFile(parfile, 'DX',               m%dy,                 m%dx)
   call readParFile(parfile, 'DZ',               m%dz,                 m%dx)
   call readParFile(parfile, 'VMIN',             m%vmin,               -1.0)
   call readParFile(parfile, 'VMAX',             m%vmax,               -1.0)
endif
call MPI_BCAST(m%nx,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%ny,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%nz,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%dx,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%dy,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%dz,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%vmin,1,         MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(m%vmax,1,         MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_model3d

!=======================================================

subroutine readparamfile_sub_water_model2d(parfile,wm2d)
character(len=*),    intent(in)          :: parfile
type(water_model2d), intent(out)         :: wm2d
call readparamfile_sub_model2d(parfile,wm2d%m)
if (rank.eq.0) then
   call readParFile(parfile, 'NZW',              wm2d%nzw,                -1)
endif
call MPI_BCAST(wm2d%nzw,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_water_model2d

!=======================================================

subroutine readparamfile_sub_water_model3d(parfile,wm3d)
character(len=*),    intent(in)          :: parfile
type(water_model3d), intent(out)         :: wm3d
call readparamfile_sub_model3d(parfile,wm3d%m)
if (rank.eq.0) then
   ! Integer parameters
   call readParFile(parfile, 'NZW',              wm3d%nzw,                -1)
endif
call MPI_BCAST(wm3d%nzw,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_water_model3d

!=======================================================

subroutine readparamfile_sub_km2d(parfile, k)
character(len=*),  intent(in)  :: parfile
type(km2d),          intent(out) :: k
if (rank.eq.0) then
   call readParFile(parfile, "NT",        k%nt,        -1)
   call readParFile(parfile, "NT_CONV",   k%nt_conv,    2*k%nt-1)
!   call readParFile(parfile, "NW",        k%nw,        -1)
   call readParFile(parfile, "IsSaveTTT", k%isSaveTTT, .false.)
   call readParFile(parfile, "IsSaveBothTTT", k%isSaveBothTTT, .false.)
   call readParFile(parfile, "IsTheta",   k%isTheta,   .false.)
   call readParFile(parfile, "IsDip",     k%isDip,     .false.)
   call readParFile(parfile, "IsWriteData",k%isWriteData,.false.)
   call readParFile(parfile, "DT",        k%dt,        -1.0)
   call readParFile(parfile, "THETA",     k%theta,     -1.0)
   call readParFile(parfile, "DIP",       k%dip,       -1.0)
   call readParFile(parfile, "EIK_TYPE",  k%eik_type,  "DIR")
   call readParFile(parfile, "SORT_IN",   k%sort_in,   "CSG")
   call readParFile(parfile, "DATA_FORMAT",k%dfor,  "DEFAULT")
   call readParFile(parfile, "DATA_FORMAT_COE",k%dfor_coe,0.0)
   call readParFile(parfile, "IsWriteGreen",k%isWriteGreen,  .false.) ! Used in marchenko
   call readParFile(parfile, "IsWritef0m",k%isWritef0m,  .false.) ! Used in marchenko
   call readParFile(parfile, "Isgttt",k%isgttt,  .true.) ! Used in marchenko
   call readParFile(parfile, "IsSaveGreen",k%isSaveGreen,  .false.) ! Used in marchenko
   call readParFile(parfile, "IsSaveGreenShot",k%isSaveGreenShot,  .false.) ! Used in marchenko

endif
call MPI_BCAST(k%nt,          1, MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%nt_conv,          1, MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
!call MPI_BCAST(k%nw,          1, MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isSaveTTT,   1, MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isSaveBothTTT,1,MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isTheta,     1, MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isDip,       1, MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isWriteData,1,  MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isWriteGreen,1,  MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr) ! Used in marchenko
call MPI_BCAST(k%isWritef0m,  1,  MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr) ! Used in marchenko
call MPI_BCAST(k%isgttt,      1,  MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr) ! Used in marchenko
call MPI_BCAST(k%isSaveGreen,      1,  MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr) ! Used in marchenko
call MPI_BCAST(k%isSaveGreenShot,      1,  MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr) ! Used in marchenko

call MPI_BCAST(k%dt,          1, MPI_REAL,      0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%theta,       1, MPI_REAL,      0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%dip,         1, MPI_REAL,      0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%eik_type,  100, MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%sort_in,   100, MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%dfor,      100, MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%dfor_coe,    1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_km2d

!=======================================================

subroutine readparamfile_sub_km3d(parfile, k)
character(len=*),  intent(in)  :: parfile
type(km3d),          intent(out) :: k
if (rank.eq.0) then
   call readParFile(parfile, "NT",        k%nt,        -1)
!   call readParFile(parfile, "NW",        k%nw,        -1)
   call readParFile(parfile, "IsSaveTTT", k%isSaveTTT, .false.)
   call readParFile(parfile, "IsSaveBothTTT", k%isSaveBothTTT, .false.)
   call readParFile(parfile, "IsTheta",   k%isTheta,   .false.)
   call readParFile(parfile, "IsDip",     k%isDip,     .false.)
   call readParFile(parfile, "IsWriteData",k%isWriteData,.false.)
   call readParFile(parfile, "DT",        k%dt,        -1.0)
   call readParFile(parfile, "THETA",     k%theta,     -1.0)
   call readParFile(parfile, "DIP",       k%dip,       -1.0)
   call readParFile(parfile, "EIK_TYPE",  k%eik_type,  "DIR")
   call readParFile(parfile, "SORT_IN",   k%sort_in,   "CSG")
   call readParFile(parfile, "DATA_FORMAT",k%dfor,  "DEFAULT")
   call readParFile(parfile, "DATA_FORMAT_COE",k%dfor_coe,0.0)
endif
call MPI_BCAST(k%nt,          1, MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
!call MPI_BCAST(k%nw,          1, MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isSaveTTT,   1, MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isSaveBothTTT,1,MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isTheta,     1, MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isDip,       1, MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%isWriteData, 1, MPI_LOGICAL,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%dt,          1, MPI_REAL,      0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%theta,       1, MPI_REAL,      0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%dip,         1, MPI_REAL,      0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%eik_type,  100, MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%sort_in,   100, MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%dfor,      100, MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(k%dfor_coe,    1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_km3d

!=======================================================

subroutine readparamfile_sub_gdm2d(parfile, g)
character(len=*),  intent(in)  :: parfile
type(gdm2d),       intent(out) :: g
if (rank.eq.0) then
   call readParFile(parfile, "NT",        g%nt,        -1)
   call readParFile(parfile, "DT",        g%dt,        -1.0)
endif
call MPI_BCAST(g%nt,          1, MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(g%dt,          1, MPI_REAL,      0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_gdm2d

!=======================================================

subroutine readparamfile_sub_gdm3d(parfile, g)
character(len=*),  intent(in)  :: parfile
type(gdm3d),       intent(out) :: g
if (rank.eq.0) then
   call readParFile(parfile, "NT",        g%nt,        -1)
   call readParFile(parfile, "DT",        g%dt,        -1.0)
endif
call MPI_BCAST(g%nt,          1, MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(g%dt,          1, MPI_REAL,      0,MPI_COMM_WORLD,ierr)
end subroutine readparamfile_sub_gdm3d

!==============================================================================

end module module_read_para

