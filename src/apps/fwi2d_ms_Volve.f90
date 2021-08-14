! 2D frequency multiscale full waveform inversion based on constant density acoustic equation
! Written by Lei Fu based on Xin Wang's matlab version
! @2016-5-13
! for Volve cable2  data, cut traces near sources, and mute the direct waves

program fwi2d_ms 
use module_global, only : parfile
use module_sf90_mpi
use module_parser
use module_datatype
use module_inversion
use module_acquisition2d
use module_read_para
use module_string
use module_source
use module_io
use module_aperture2d
use module_a2d
use module_utility
use module_math
use module_array
implicit none
type(fname)          :: fn
type(other)          :: o
type(inv)            :: i
type(aperture)       :: a
type(source)         :: s
type(fd2d)           :: f
type(bc2d_general)   :: bc
type(wf2d)           :: wf
type(model2d)        :: m,m_temp
type(coord2d)        :: c_all,c
type(shot2d)         :: sh
type(image2d)        :: img
type(fdmod2d)        :: fdm,fdm_h
integer              :: ipro,npro_me,it,res_nit,ig,itt,nt,ng,Ncut,fr_Ind
real                 :: t1, t2,pur=0.02,f1,f2,f3,fff ! fff is the true step length ratio
real                 :: f11,f22,f33,f44,ftmp ! frequencies used for bandpass filtering
logical              :: isError
integer,allocatable  :: ipro_me(:)
real,    allocatable :: data_true(:,:),gk(:,:),gk1(:,:),gk1_s(:,:),gk1_temp(:,:)&
   ,prec_gk1(:,:),prec_gk1_s(:,:),prec_gk1_temp(:,:),prec_gk(:,:),prec_gk_temp(:,:),dk(:,:),dk1(:,:),mk(:,:),temp1(:,:),ss(:,:),ss1(:,:),vtemp(:,:),mask_csg(:,:)! ss is the slowness 
real,allocatable     :: data_temp(:,:),res_data(:,:),temp_out(:,:)
real, allocatable    :: res(:),res_dt(:),fL(:),fH(:),xs(:,:),xg(:,:) ! FH(:) , store the high frequncy
double precision     :: misfit_temp,misfit,gTLTLg_temp,gTLTLg,alpha,misfit1,misfit2,misfit3,err,res0
character(len=slen)   :: output, output_temp

call cpu_time(t1)
call readCmd("par",parfile,"parfile")
call init_mpi
! Read in parameters
isError=.false.
call readparamfile(parfile,fn,o,i,a,s,f,m)
f%isReadData=.false.;f%isWriteData=.false.
! Read acquisition geometry data
call read_coordfile_mpi(fn%coord_in_file,o%coordfile_type,c_all)
! Set up free surface
call allocate_and_initial(f%fs,m%nx)
call read_topofile_mpi(fn%topo_file,m%nx,f%npml,m%dx,a%z0,f%fs)
! Read velocity model
call allocate_and_initial(m%v,m%nz,m%nx)
call read_binfile_mpi(fn%vel_file, m%v, m%nz,m%nx)
!call read_binfile_mpi(fn%vel_homo_file,m_h%v,m%nz,m%nx)
! read mask file for gradient
if (o%isMask) then
call allocate_and_initial(m%mask,m%nz,m%nx)
call read_binfile_mpi(fn%mask_file,m%mask,m%nz,m%nx)
endif

call allocate_and_initial(temp_out,f%nt,c%ngmax)

m%vmin=minval(m%v)
m%vmax=maxval(m%v)

if (rank.eq.0) then
   write(*,*) " Frequency Multiscale FWI Based on Constant Density Acoustic Equation "
   call display_parameters(m,f,o,i,c_all,a,fn,s) 
   call log_parameters(m,f,o,i,c_all,a,s,fn)
endif

call allocate_and_initial(fL,i%nf)
call allocate_and_initial(fH,i%nf)
fL(1)=i%fl1
fL(2)=i%fl2
fL(3)=i%fl3
fL(4)=i%fl4
fL(5)=i%fl5
fL(6)=i%fl6

fH(1)=i%fh1
fH(2)=i%fh2
fH(3)=i%fh3
fH(4)=i%fh4
fH(5)=i%fh5
fH(6)=i%fh6

! Setup Source Wavelet
call allocate_and_initial(s%sou,s%nw)
if (s%isSource) then
   call ricker(s)
   if (f%fd_type.eq."SG") call integrate(s%sou,s%nw)
   if (rank.eq.0) call write_binfile(fn%source_file,s%sou,s%nw)
else
   call read_sourcefile_mpi(fn%source_file,s%sou,s%nw)
endif

! Static load balancing
call allocate_and_initial(ipro_me,c_all%npro)
if (o%isShotConti) then
   call get_assigned_continue(c_all%npro, npro_me, ipro_me)
else
   call get_assigned_evenly(c_all%npro, npro_me, ipro_me)
endif
call assign_coord(c_all,c,npro_me,ipro_me)
call deallocate_type(c_all)
if (i%isSeisRecons) then
   res_nit=i%nit+1
else
   res_nit=i%nit
endif
 
if (rank.eq.0) then
   call allocate_and_initial(gk1,dk1,m%nz,m%nx)
   call allocate_and_initial(gk1_s,m%nz,m%nx)
   call allocate_and_initial(res,res_nit)
   if (o%isIllum) then
      call allocate_and_initial(prec_gk,prec_gk1,m%nz,m%nx)
      call allocate_and_initial(prec_gk1_s,m%nz,m%nx)
   endif
   if (i%inv_type.eq."CG") then
      call allocate_and_initial(gk,dk,m%nz,m%nx)
   endif
endif

!Ncut=0
call copytype(m,m_temp)
! ######### Inversion Loop ###############
! ########################################
it=1
fr_Ind=1

do while (it<=i%nit)
   f11 = fL(fr_Ind);
   f22 = f11+0.01*fH(fr_Ind);
   f33 = fH(fr_Ind);
   f44 = f33+0.05*f33; 

   o%fpk = 0.5*(f22+f33)
 
   if (rank.eq.0)     write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 
   if (rank.eq.0)     write(*,*) "=============================================================================" 
   if (rank.eq.0)     write(*,*) "===========  Frequency Multiscale (==>",fr_Ind,"        <==) FWI " 
   if (rank.eq.0)     write(*,*) "(" ,f11,f22,f33,f44,") Hz"
   if (rank.eq.0)     write(*,*) "!   PEAK FREQ " ,o%fpk, " Hz !" 
   if (rank.eq.0)     write(99,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 
   if (rank.eq.0)     write(99,*) "=============================================================================" 
   if (rank.eq.0)     write(99,*) "===========  Frequency Multiscale (==>",fr_Ind,"        <==) FWI " 
   if (rank.eq.0)     write(99,*) "(" ,f11,f22,f33,f44,") Hz" 
   if (rank.eq.0)     write(99,*) "!   PEAK FREQ " ,o%fpk, " Hz !" 
   if ((rank.eq.0).and.(o%grad_flip))     write(99,*) "Flip the gradient " 
   
   ! calculate the new data residual res = L m_k - d_obs
   misfit_temp=0.0
   call allocate_and_initial(gk1_temp,m%nz,m%nx)
   if (o%isIllum) then                                        ! Preconditioning 
      call allocate_and_initial(prec_gk1_temp,m%nz,m%nx)
   endif
   
   call allocate_and_initial(res_data,f%nt,c%ngmax)
   call allocate_and_initial(res_dt,c%ngmax)
   
   do ipro=1,npro_me
      call aperture_init(c,a,m,f,s,o,ipro,sh,img,fdm)
      call allocate_and_initial(data_true,sh%seis,sh%nt,sh%ng)
      call allocate_and_initial(mask_csg,sh%nt,sh%ng)
      call allocate_and_initial(sh%seis_dir,sh%nt,sh%ng)
      call allocate_and_initial(xs,sh%ng,1)
      call allocate_and_initial(xg,sh%ng,1)
      ! new data residual
      call a2d_mod(fdm,bc,wf,sh,output_temp)
      ! read the direct wave, and removed from the predicted data 
      if (o%isRemove_dir) then
      call filename(output, fn%csg_dir, c%proid(ipro), ".bin")
      call read_binfile(output,sh%seis_dir,sh%nt,sh%ng)
      sh%seis=sh%seis-sh%seis_dir
      endif

      if (o%isMask_csg) then
      call filename(output, fn%csg_mask, c%proid(ipro), ".bin")
      call read_binfile(output,mask_csg,sh%nt,sh%ng)
      sh%seis = sh%seis*mask_csg
      endif

      call filename(output, fn%csg_out, c%proid(ipro), ".bin")
      call read_binfile(output,data_true,sh%nt,sh%ng)
   
      xs=c%xs
      xg=c%xg
      if (o%isCut) call csg_cut(sh%seis,data_true,sh%nt,sh%ng, xs,xg, o%near_cut,o%far_cut) 
      if (o%isWindow_csg)  call csg_window(sh%seis,data_true,sh%nt,sh%ng, sh%dt,o%T_window,o%n_taper)
      if (o%is_conv_sou.eq.1) call csg_conv_sou(sh%seis,data_true,sh%nt,sh%dt,sh%ng,o%fpk) 
      !call filename(output, fn%csg_out, c%proid(ipro), "_pre.bin")
      !call write_binfile(output,sh%seis,sh%nt, sh%ng)
      !call filename(output, fn%csg_out, c%proid(ipro), "_obs.bin")
      !call write_binfile(output,data_true,sh%nt, sh%ng)

      call res_fwi(sh%seis,data_true,sh%nt,sh%ng,res_data,res_dt,sh%dt,f11,f22,f33,f44)
      sh%seis=res_data
      !call filename(output, fn%csg_out, c%proid(ipro), "_res.bin")
      !call write_binfile(output,sh%seis ,sh%nt, sh%ng)
  
      call deallocate_and_free(data_true)
      call deallocate_and_free(mask_csg)
      call deallocate_and_free(xs,xg)

      misfit_temp=misfit_temp+0.5*sum(res_dt)
      ! new gradient
      call a2d_mig(sh,bc,wf,fdm,img,output_temp)
      call output_illum(o%isIllum,o%isSaveIllum,fn%pre_illum_out,&  ! Preconditioning
           c%proid(ipro),img%illum,img%ix1,img%ix2,img%nz,img%nx)
      call output_prestack_image(o%isPreStackImg,fn%pre_illum_img_out,&
           c%proid(ipro),img%img,img%ix1,img%ix2,img%nz,img%nx)
      call allocate_and_initial(temp1,m%nz,m%nx)
      call aperture_image_l2g(a%isAper,img%img,img%iz1,img%iz2,img%ix1,img%ix2,temp1,m%nz,m%nx)
      gk1_temp=gk1_temp+temp1
      call deallocate_and_free(temp1)
      if (o%isIllum) then
         img%img=img%img/(img%illum+0.001*maxval(img%illum))
         call allocate_and_initial(temp1,m%nz,m%nx)
         call aperture_image_l2g(a%isAper,img%img,img%iz1,img%iz2,img%ix1,img%ix2,temp1,m%nz,m%nx)
         prec_gk1_temp=prec_gk1_temp+temp1
         call deallocate_and_free(temp1)
      endif
   enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   misfit=0.0
   call MPI_REDUCE(misfit_temp,misfit,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call sf90_bcast(misfit)
   call message(0,"it=",it,"		 misfit=",misfit)
   if (rank.eq.0)     write(99,*) "it = " ,it," misfit = ",misfit 
   
   res0 = misfit  !initial misfit

   call MPI_REDUCE(gk1_temp,gk1,m%nz*m%nx,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   if (o%isIllum) then
      call MPI_REDUCE(prec_gk1_temp,prec_gk1,m%nz*m%nx,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   endif
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 ! smoothe the gradient
 if (rank.eq.0) then
    if (o%isSmooth_gk) then
        call smooth(gk1,m%nz,m%nx,gk1_s,o%hw)
        if (o%isIllum) then
        call smooth(prec_gk1,m%nz,m%nx,prec_gk1_s,o%hw)
        endif
     else
         gk1_s=gk1
         if (o%isIllum) then
           prec_gk1_s = prec_gk1
         endif
    endif
 endif

! mask the gradient
  if (rank.eq.0) then
   if (o%isMask) then
    gk1_s=gk1_s*m%mask
    if (o%isIllum) then
      prec_gk1_s=prec_gk1_s*m%mask
    endif
   endif
  endif

   call deallocate_and_free(gk1_temp)
   if (o%isIllum) then
      call deallocate_and_free(prec_gk1_temp)
   endif
  
   call allocate_and_initial(m%refl,m%nz,m%nx)
   if (rank.eq.0) then
      res(it)=misfit
   !   call gradient(i%inv_type,i%cg_type,o%isIllum,it,gk,gk1,prec_gk1,dk,dk1)
      m%refl=gk1_s   ! for conv data
       if (o%isIllum) then
        m%refl=prec_gk1_s
       endif
   endif
   ! write out the gradient for each iteration
   if (rank.eq.0) then 
     call filename(output, fn%gk1_file,it, ".bin")
     call write_binfile(output,m%refl,m%nz,m%nx)
   endif
   
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_BCAST(m%refl,m%nz*m%nx,MPI_REAL,0,MPI_COMM_WORLD,ierr)
   
   ! ############ Calculate the step length ############### 
   call allocate_and_initial(ss,ss1,m%nz,m%nx)
   ss=1/m%v

   alpha=pur*sqrt(sum(ss**2))/(0.001*maxval(m%refl)+sqrt(sum(m%refl**2)))
   if (o%grad_flip)  alpha=-1.0*alpha
   if (rank.eq.0) write(*,*)"alpha = ",alpha
   if (rank.eq.0) write(99,*)"alpha = ",alpha

   ! Back tracking to find the numerical step length
   !====================================================================================================
   f1=0.5 
   ss1=ss+alpha*f1*m%refl
   call allocate_and_initial(vtemp,m%nz,m%nx)
   if (.NOT.(allocated(m_temp%v))) call allocate_and_initial(m_temp%v,m%nz,m%nx)
   m_temp%v=1/ss1
   call control_value(m_temp%v,i%vmax,i%vmin)         
   misfit_temp=0.0

   do ipro=1,npro_me
      call aperture_init(c,a,m_temp,f,s,o,ipro,sh,img,fdm)
      call allocate_and_initial(data_true,sh%seis,sh%nt,sh%ng)
      call allocate_and_initial(mask_csg,sh%nt,sh%ng)
      call allocate_and_initial(sh%seis_dir,sh%nt,sh%ng)
      call allocate_and_initial(xs,sh%ng,1)
      call allocate_and_initial(xg,sh%ng,1)
      call a2d_mod(fdm,bc,wf,sh,output_temp)
      
      ! read the direct wave, and removed from the predicted data 
      if (o%isRemove_dir) then
      call filename(output, fn%csg_dir, c%proid(ipro), ".bin")
      call read_binfile(output,sh%seis_dir,sh%nt,sh%ng)
      sh%seis=sh%seis-sh%seis_dir
      endif
 
      if (o%isMask_csg) then
      call filename(output, fn%csg_mask, c%proid(ipro), ".bin")
      call read_binfile(output,mask_csg,sh%nt,sh%ng)
      sh%seis = sh%seis*mask_csg
      endif
      call filename(output, fn%csg_out, c%proid(ipro), ".bin")
      call read_binfile(output,data_true,sh%nt,sh%ng)
      
      xs=c%xs
      xg=c%xg
      if (o%isCut) call csg_cut(sh%seis,data_true,sh%nt,sh%ng, xs,xg, o%near_cut,o%far_cut) 
      if (o%isWindow_csg)  call csg_window(sh%seis,data_true,sh%nt,sh%ng, sh%dt,o%T_window,o%n_taper)
      if (o%is_conv_sou.eq.1) call csg_conv_sou(sh%seis,data_true,sh%nt,sh%dt,sh%ng,o%fpk) 
      call res_fwi(sh%seis,data_true,sh%nt,sh%ng,res_data,res_dt,sh%dt,f11,f22,f33,f44)
      call deallocate_and_free(data_true)
      call deallocate_and_free(mask_csg)
      call deallocate_and_free(xs)
      call deallocate_and_free(xg)
      misfit_temp=misfit_temp+0.5*sum(res_dt)
   enddo
   !call aperture_finalize(sh,fdm)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   misfit1=0.0 
   call MPI_REDUCE(misfit_temp,misfit1,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call sf90_bcast(misfit1)
   if (rank.eq.0)  write(*,*) "f1 = ",f1, " residual 1 = ", misfit1
   if (rank.eq.0)  write(99,*) "f1 = ",f1, " residual 1 = ", misfit1

   if (misfit1>misfit)   then
      !do while ((misfit1>misfit) .and. (f1>0.001))
      do while ((misfit1>misfit) .and. (f1>0.015))
          f2=f1;misfit2=misfit1
          f1=f1/2.0
          ss1=ss+alpha*f1*m%refl
          m_temp%v=1/ss1
          call control_value(m_temp%v,i%vmax,i%vmin)
          misfit_temp=0.0
          do ipro=1,npro_me
             call aperture_init(c,a,m_temp,f,s,o,ipro,sh,img,fdm)
             call allocate_and_initial(data_true,sh%seis,sh%nt,sh%ng)
             call allocate_and_initial(mask_csg,sh%nt,sh%ng)
             call allocate_and_initial(sh%seis_dir,sh%nt,sh%ng)
             call allocate_and_initial(xs,sh%ng,1)
             call allocate_and_initial(xg,sh%ng,1)
             call a2d_mod(fdm,bc,wf,sh,output_temp)
      
             ! read the direct wave, and removed from the predicted data 
             if (o%isRemove_dir) then
             call filename(output, fn%csg_dir, c%proid(ipro), ".bin")
             call read_binfile(output,sh%seis_dir,sh%nt,sh%ng)
             sh%seis=sh%seis-sh%seis_dir
             endif

             if (o%isMask_csg) then
             call filename(output, fn%csg_mask, c%proid(ipro), ".bin")
             call read_binfile(output,mask_csg,sh%nt,sh%ng)
             sh%seis = sh%seis*mask_csg
             endif

             call filename(output, fn%csg_out, c%proid(ipro), ".bin")
             call read_binfile(output,data_true,sh%nt,sh%ng)
            
             xs=c%xs
             xg=c%xg
             if (o%isCut) call csg_cut(sh%seis,data_true,sh%nt,sh%ng, xs,xg, o%near_cut,o%far_cut) 
             if (o%isWindow_csg)  call csg_window(sh%seis,data_true,sh%nt,sh%ng, sh%dt,o%T_window,o%n_taper)
             if (o%is_conv_sou.eq.1) call csg_conv_sou(sh%seis,data_true,sh%nt,sh%dt,sh%ng,o%fpk) 
             call res_fwi(sh%seis,data_true,sh%nt,sh%ng,res_data,res_dt,sh%dt,f11,f22,f33,f44)
             call deallocate_and_free(data_true)
             call deallocate_and_free(mask_csg)
             call deallocate_and_free(xs)
             call deallocate_and_free(xg)
             misfit_temp=misfit_temp+0.5*sum(res_dt)
          enddo
          !call aperture_finalize(sh,fdm)
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          misfit1=0.0 
          call MPI_REDUCE(misfit_temp,misfit1,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          call sf90_bcast(misfit1)
          if (rank.eq.0)  write(*,*) "f1 = ",f1, " residual 1 = ", misfit1
          if (rank.eq.0)  write(99,*) "f1 = ",f1, " residual 1 = ", misfit1
      enddo
   else
      
      f2=f1*2  
      ss1=ss+alpha*f2*m%refl
      m_temp%v=1/ss1
      call control_value(m_temp%v,i%vmax,i%vmin)
      misfit_temp=0.0
      do ipro=1,npro_me
         call aperture_init(c,a,m_temp,f,s,o,ipro,sh,img,fdm)
         call allocate_and_initial(data_true,sh%seis,sh%nt,sh%ng)
         call allocate_and_initial(mask_csg,sh%nt,sh%ng)
         call allocate_and_initial(sh%seis_dir,sh%nt,sh%ng)
         call allocate_and_initial(xs,sh%ng,1)
         call allocate_and_initial(xg,sh%ng,1)
         call a2d_mod(fdm,bc,wf,sh,output_temp)
             
         ! read the direct wave, and removed from the predicted data 
         if (o%isRemove_dir) then
         call filename(output, fn%csg_dir, c%proid(ipro), ".bin")
         call read_binfile(output,sh%seis_dir,sh%nt,sh%ng)
         sh%seis=sh%seis-sh%seis_dir
         endif
 
         if (o%isMask_csg) then
         call filename(output, fn%csg_mask, c%proid(ipro), ".bin")
         call read_binfile(output,mask_csg,sh%nt,sh%ng)
         sh%seis = sh%seis*mask_csg
         endif

         call filename(output, fn%csg_out, c%proid(ipro), ".bin")
         call read_binfile(output,data_true,sh%nt,sh%ng)
         xs=c%xs
         xg=c%xg
         if (o%isCut) call csg_cut(sh%seis,data_true,sh%nt,sh%ng, xs,xg, o%near_cut,o%far_cut) 
         if (o%isWindow_csg)  call csg_window(sh%seis,data_true,sh%nt,sh%ng, sh%dt,o%T_window,o%n_taper)
         if (o%is_conv_sou.eq.1) call csg_conv_sou(sh%seis,data_true,sh%nt,sh%dt,sh%ng,o%fpk) 
         call res_fwi(sh%seis,data_true,sh%nt,sh%ng,res_data,res_dt,sh%dt,f11,f22,f33,f44)
         call deallocate_and_free(data_true)
         call deallocate_and_free(mask_csg)
         call deallocate_and_free(xs)
         call deallocate_and_free(xg)
         misfit_temp=misfit_temp+0.5*sum(res_dt)
      enddo
      !call aperture_finalize(sh,fdm)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      misfit2=0.0 
      call MPI_REDUCE(misfit_temp,misfit2,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call sf90_bcast(misfit2)
      if (rank.eq.0)  write(*,*) "f2 = ",f2, " residual 2 = ", misfit2
      if (rank.eq.0)  write(99,*) "f2 = ",f2, " residual 2 = ", misfit2
   endif
   
   fff=((f1**2)*(misfit-misfit2)+(f2**2)*(misfit1-misfit))/(2*misfit*(f1-f2)+2*misfit1*f2-2*misfit2*f1)
   ss1=ss+alpha*fff*m%refl
   m_temp%v=1/ss1
   call control_value(m_temp%v,i%vmax,i%vmin)
   misfit_temp=0.0
   do ipro=1,npro_me
      call aperture_init(c,a,m_temp,f,s,o,ipro,sh,img,fdm)
      call allocate_and_initial(data_true,sh%seis,sh%nt,sh%ng)
      call allocate_and_initial(mask_csg,sh%nt,sh%ng)
      call allocate_and_initial(sh%seis_dir,sh%nt,sh%ng)
      call allocate_and_initial(xs,sh%ng,1)
      call allocate_and_initial(xg,sh%ng,1)
      call a2d_mod(fdm,bc,wf,sh,output_temp)
         
      ! read the direct wave, and removed from the predicted data 
      if (o%isRemove_dir) then
      call filename(output, fn%csg_dir, c%proid(ipro), ".bin")
      call read_binfile(output,sh%seis_dir,sh%nt,sh%ng)
      sh%seis=sh%seis-sh%seis_dir
      endif 

      if (o%isMask_csg) then
      call filename(output, fn%csg_mask, c%proid(ipro), ".bin")
      call read_binfile(output,mask_csg,sh%nt,sh%ng)
      sh%seis = sh%seis*mask_csg
      endif

      call filename(output, fn%csg_out, c%proid(ipro), ".bin")
      call read_binfile(output,data_true,sh%nt,sh%ng)
      xs=c%xs
      xg=c%xg
      if (o%isCut) call csg_cut(sh%seis,data_true,sh%nt,sh%ng, xs,xg, o%near_cut,o%far_cut) 
      if (o%isWindow_csg)  call csg_window(sh%seis,data_true,sh%nt,sh%ng, sh%dt,o%T_window,o%n_taper)
      if (o%is_conv_sou.eq.1) call csg_conv_sou(sh%seis,data_true,sh%nt,sh%dt,sh%ng,o%fpk) 
      call res_fwi(sh%seis,data_true,sh%nt,sh%ng,res_data,res_dt,sh%dt,f11,f22,f33,f44)
      call deallocate_and_free(data_true)
      call deallocate_and_free(mask_csg)
      call deallocate_and_free(xs)
      call deallocate_and_free(xg)
      misfit_temp=misfit_temp+0.5*sum(res_dt)
   enddo
   !call aperture_finalize(sh,fdm)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   misfit3=0.0 
   call MPI_REDUCE(misfit_temp,misfit3,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call sf90_bcast(misfit3)
   if (rank.eq.0)  write(*,*) "fff = ",fff, " residual 3 = ", misfit3
   if (rank.eq.0)  write(99,*) "fff = ",fff, " residual 3 = ", misfit3

   if (rank.eq.0) then
      if ((misfit3>misfit1) .or. (misfit3>misfit2)) then
         if (misfit2>misfit1) then
            misfit=misfit1;fff=f1
         else
            misfit=misfit2;fff=f2
         endif
      else
         misfit=misfit3
      endif
      write(*,*) "The misfit = ", misfit
      write(99,*) "The misfit = ", misfit
   endif

   call sf90_bcast(fff)
   call sf90_bcast(misfit)
   
   ss1=ss+alpha*fff*m%refl
   m%v=1/ss1
   call control_value(m%v,i%vmax,i%vmin)
   m_temp%v=1/ss1
   call control_value(m_temp%v,i%vmax,i%vmin)
   
   ! write down the velocity    
   if (rank.eq.0) then
   call filename(output, fn%img_file, it, ".bin")
   call write_binfile(output,m%v,m%nz,m%nx)
   endif

   ! Refresh gradient
   if (rank.eq.0) then
     if (i%inv_type.eq."CG") then
        gk=gk1;dk=dk1
     endif
     if (o%isIllum) then
       prec_gk=prec_gk1
     endif
   endif
  
   err=1.0*abs(misfit-res0)/res0
   if (rank.eq.0) then
     write(*,*),"Misfit Change = "   , err 
     write(99,*),"Misfit Change = "   , err 
     write(*,*) "  "
     write(99,*) "  "
   end if
   if (err<o%toler) then
      fr_Ind = fr_Ind+1
      if (fr_Ind>i%nf) then
          exit
      endif         
   endif

  it=it+1
enddo ! Inversion loop end 

   ! release memory
   if (rank.eq.0) then
      call deallocate_and_free(gk1,dk1)
      if (i%inv_type.eq."CG") then
         call deallocate_and_free(gk,dk)
       endif
    endif

   ! reconstruct shot gather
   if (i%isSeisRecons) then
      misfit_temp=0.0
      do ipro=1,npro_me
         call aperture_init(c,a,m,f,s,o,ipro,sh,img,fdm)
         call allocate_and_initial(data_true,sh%seis,sh%nt,sh%ng)
         call allocate_and_initial(mask_csg,sh%nt,sh%ng)
         call allocate_and_initial(sh%seis_dir,sh%nt,sh%ng)
         call allocate_and_initial(xs,sh%ng,1)
         call allocate_and_initial(xg,sh%ng,1)
         call a2d_mod(fdm,bc,wf,sh,output_temp)
      
         ! read the direct wave, and removed from the predicted data 
         if (o%isRemove_dir) then
         call filename(output, fn%csg_dir, c%proid(ipro), ".bin")
         call read_binfile(output,sh%seis_dir,sh%nt,sh%ng)
         sh%seis=sh%seis-sh%seis_dir
         endif

         if (o%isMask_csg) then
         call filename(output, fn%csg_mask, c%proid(ipro), ".bin")
         call read_binfile(output,mask_csg,sh%nt,sh%ng)
         sh%seis = sh%seis*mask_csg
         endif

         call filename(output, fn%csg_recons, c%proid(ipro), ".bin")
         call write_binfile(output, sh%seis, sh%nt, sh%ng)
         call allocate_and_initial(data_true,sh%nt,sh%ng)
         call filename(output, fn%csg_out, c%proid(ipro), ".bin")
         call read_binfile(output,data_true,sh%nt, sh%ng)
         xs=c%xs
         xg=c%xg
         if (o%isCut) call csg_cut(sh%seis,data_true,sh%nt,sh%ng, xs,xg, o%near_cut,o%far_cut) 
         if (o%isWindow_csg)  call csg_window(sh%seis,data_true,sh%nt,sh%ng, sh%dt,o%T_window,o%n_taper)
         if (o%is_conv_sou.eq.1) call csg_conv_sou(sh%seis,data_true,sh%nt,sh%dt,sh%ng,o%fpk) 
         call res_fwi(sh%seis,data_true,sh%nt,sh%ng,res_data,res_dt,sh%dt,f11,f22,f33,f44)
         misfit_temp=misfit_temp+0.5*sum(res_dt)
         call deallocate_and_free(data_true)
         call deallocate_and_free(mask_csg)
         call deallocate_and_free(xs)
         call deallocate_and_free(xg)
       enddo   
       !call aperture_finalize(sh,fdm)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       misfit=0.0
       call MPI_REDUCE(misfit_temp,misfit,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       call message(0,"final misfit=",misfit)
      if (rank.eq.0) then
        res(i%nit+1)=misfit
      endif  
   endif
   ! write the misfit
   if (rank.eq.0) then
     !write(*,*) "res",res
     call write_binfile(fn%res_file,res,i%nit+1)
     call deallocate_and_free(res)
   endif

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call deallocate_type(f,m,s);call deallocate_type(c);
   call deallocate_type(m_temp); deallocate(ipro_me);

   if (rank.eq.0) call log_end(99)
   call stop_mpi
   call cpu_time(t2)
   if (rank == 0) write(*,*) 'Runtime = ', t2-t1
contains

subroutine display_parameters(m,f,o,i,c,a,fn,s)
implicit none
type(model2d),       intent(in) :: m
type(fd2d),          intent(in) :: f
type(other),         intent(in) :: o
type(inv),           intent(in) :: i
type(coord2d),       intent(in) :: c
type(aperture),      intent(in) :: a
type(fname),         intent(in) :: fn
type(source),        intent(in) :: s
write(*,*)"nx=",m%nx, " nz=",m%nz, " dx=",m%dx
write(*,*)"nt=",f%nt," dt=",f%dt
if (s%isSource) then
   write(*,*)"Make source wavelet, frequency=",s%freq
else
   write(*,*)"Read in wavelt file, source file:",fn%source_file
endif
write(*,*)"Sort in:",o%sort_in
write(*,*)"Coordinats file type:",o%coordfile_type
write(*,*)"Coordinats file:",fn%coord_in_file
write(*,*)"ns=",c%npro," ngmax=",c%ngmax
write(*,*)
write(*,*)"x0=",a%x0," z0=",a%z0
write(*,*)"IsAper=",a%isAper
if (a%isAper) then
   write(*,*) "Aper_x=",a%aper_x
   write(*,*) "Aper_z=",a%aper_z
endif
write(*,*)
write(*,*)"Velocity file:",fn%vel_file
write(*,*)"ISMASK:",o%isMask
write(*,*)"Mask file:",fn%mask_file
write(*,*)"Z0 file:",fn%z0_mig_file
write(*,*)
write(*,*)"VMIN=",m%vmin," ,VMAX=",m%vmax
write(*,*)
if (o%isCut) write(*,*)"Cut Near-off Traces=",o%ncut,'Cut Far-offset Traces=',o%ncut2
write(*,*)
if (o%isRemove_dir) write(*,*)"Remove the direct wave from predicted data!"

write(*,*)
!if (o%wt_res.eq.1) write(*,*) "Residual Method: Hilbert Envelope "
!if (o%wt_res.eq.2) write(*,*) "Residual Method: Shihang's Method "
if (o%isIllum) then
   write(*,*)"Use the Illumination Compensation!"
write(*,*)
   write(*,*)"Use the Prestack Illumination!"
write(*,*)
else
   write(*,*)"No Use the Illumination Compensation!"
write(*,*)
endif
if (o%isSmooth_gk) then
   write(*,*)"Smoothing Gradent!"
write(*,*)
endif
!call stop_mpi
call cpu_time(t2)
!if (rank == 0) write(*,*) 'Runtime = ', t2-t1
!if (rank == 0) write(*,*) 'Runtime = ', (t2-t1)/nsize
if (rank == 0) write(*,*)"Inversion type:",i%inv_type
if (i%inv_type.eq."CG") then
   write(*,*)"CG type=",i%cg_type
endif
write(*,*)"IsMig0:",i%isMig0
if (i%isMig0) then
   write(*,*)"Mig0_file:",fn%mig0_file
endif
write(*,*)"Iteration Number:",i%nit
write(*,*)
end subroutine display_parameters

subroutine log_parameters(m,f,o,i,c,a,s,fn)
implicit none
type(model2d),   intent(in) :: m
type(fd2d),      intent(in) :: f
type(inv),       intent(in) :: i
type(other),     intent(in) :: o
type(coord2d),   intent(in) :: c
type(aperture),  intent(in) :: a
type(source),    intent(in) :: s
type(fname),     intent(in) :: fn
open(99,file=fn%log_file)
call date_and_time(date,time_now)
write(99,*)"Date ",date
write(99,*)"Time ",time_now
write(99,*)
write(99,*)"nx=",m%nx, " nz=",m%nz, " dx=",m%dx
write(99,*)"nt=",f%nt," dt=",f%dt
if (s%isSource) then
   write(99,*)"Make source wavelet, frequency=",s%freq
else
   write(99,*)"Read in wavelt file, source file:",fn%source_file
endif
write(99,*)"Sort in:",o%sort_in
write(99,*)"Coordinats file type:",o%coordfile_type
write(99,*)"Coordinats file:",fn%coord_in_file
write(99,*)"ns=",c%npro," ngmax=",c%ngmax
write(99,*)
write(99,*)"x0=",a%x0," z0=",a%z0
write(99,*)"IsAper=",a%isAper
if (a%isAper) then
   write(99,*) "Aper_x=",a%aper_x
   write(99,*) "Aper_z=",a%aper_z
endif
write(99,*)
write(99,*)"Velocity file:",fn%vel_file
write(99,*)"Z0 file:",fn%z0_mig_file
write(99,*)
write(99,*)"VMIN=",m%vmin," ,VMAX=",m%vmax
write(99,*)
!if (o%wt_res.eq.1) write(*,*) "Residual Method: Hilbert Envelope "
!if (o%wt_res.eq.2) write(*,*) "Residual Method: Shihang's Method "
if (o%isIllum) then
   write(99,*)"Use the Illumination Compensation!"
   write(99,*)"Use the Prestack Illumination!"
else
   write(99,*)"No Use the Illumination Compensation!"
endif
write(99,*)
write(99,*)"Inversion type:",i%inv_type
if (i%inv_type.eq."CG") then
 write(99,*)"CG type=",i%cg_type
endif
write(99,*)"IsMig0:",i%isMig0
if (i%isMig0) then
 write(99,*)"Mig0_file:",fn%mig0_file
endif
write(99,*)"Iteration Number:",i%nit
write(99,*)
end subroutine log_parameters

end program fwi2d_ms 

