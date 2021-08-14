program a2d_rtm_h
use module_global, only : parfile
use module_sf90_mpi
use module_parser
use module_datatype
use module_acquisition2d
use module_read_para
use module_string
use module_source
use module_io
use module_aperture2d
use module_a2d
use module_utility
implicit none
type(fname)          :: fn
type(other)          :: o
type(aperture)       :: a
type(source)         :: s
type(fd2d)           :: f
type(model2d)        :: m
type(coord2d)        :: c_all,c
type(fdmod2d)        :: fdm
type(shot2d)         :: seis
type(image2d)        :: img
type(bc2d_general)   :: bc
type(wf2d)           :: wf
integer              :: ipro,npro_me,ih,ind
real                 :: t1, t2,f11,f22,f33,f44
logical              :: isError
integer, allocatable :: ipro_me(:)
real,    allocatable :: image(:,:,:),image_temp(:,:,:),image_illum(:,:),&
                        image_illum_temp(:,:,:),temp(:,:),img_h(:,:,:)
character(len=slen)   :: output

call cpu_time(t1)
call readCmd("par",parfile,"parfile")
call init_mpi
isError=.false.
! Read in parameters
call readparamfile(parfile,fn,o,a,s,f,m)
if (isError) goto 999
! Read acquisition geometry data
call read_coordfile_mpi(fn%coord_in_file,o%coordfile_type,c_all)
! Set up free surface
call allocate_and_initial(f%fs,m%nx)
call read_topofile_mpi(fn%topo_file, m%nx, f%npml, m%dx, a%z0, f%fs)
! Read velocity model
call allocate_and_initial(m%v,m%nz,m%nx)
call read_binfile_mpi(fn%vel_file, m%v, m%nz,m%nx)
if (f%fd_type.eq."SG") then                           ! modified by Bowen Guo
   call allocate_and_initial(m%den,m%nz,m%nx)         ! modified by Bowen Guo
   call read_binfile_mpi(fn%den_file,m%den,m%nz,m%nx) ! modified by Bowen Guo
endif
m%vmin=minval(m%v)
m%vmax=maxval(m%v)
if (rank.eq.0) then 
   write(*,*) "Suboffset - Reverse Time Migration!"
   call display_parameters(m,f,o,c_all,a,fn,s)
   call log_parameters(m,f,o,c_all,a,s,fn)
endif
! Setup Ricker source wavelet
call allocate_and_initial(s%sou,s%nw)
if (s%isSource) then
   call ricker(s)
   if (f%fd_type.eq."SG") call integrate(s%sou,s%nw)   ! modified by Bowen Guo
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
!write(*,*)"rank=",rank," ns_me=",ns_me," is_me=",is_me(1:ns_me)
!call flush(6)
!call allocate_and_initial(image_temp,m%nz,m%nx,o%nh)
allocate(image_temp(m%nz,m%nx,o%nh))
!call allocate_and_initial(img%img_h,m%nz,m%nx,o%nh)

if (rank.eq.0) then
   call allocate_and_initial(image,m%nz,m%nx,o%nh)
endif

if (o%isIllum) then
   !call allocate_and_initial(image_illum_temp,m%nz,m%nx,o%nh)
   allocate(image_illum_temp(m%nz,m%nx,o%nh))
endif


    if(rank.eq.0) write(*,*) "NH = ",o%nh
do ipro=1,npro_me
   f11 = o%fl
   f22 = f11+1.0;
   f33 = o%fh
   f44 = f33+2.5;
   if ((rank.eq.0).and.(o%is_csg_filt))   write(*,*) "!" ,f11,f22,f33,f44," Hz !"

   if (o%is_csg_filt)  call csg_filt(s%sou,s%nw,s%dt,f11,f22,f33,f44)
 
   call aperture_init(c,a,m,f,s,o,ipro,seis,img,fdm)
   call allocate_and_initial(img%img_h,m%nz,m%nx,o%nh)
   
   call filename(output, fn%csg_out, c%proid(ipro), ".bin")
   if (.not.f%isReadData) then
      call allocate_and_initial(seis%seis,seis%nt,seis%ng)
      call read_binfile(output,seis%seis,seis%nt,seis%ng)
   endif
   call date_and_time(date,time_now)

   if (o%is_csg_filt) call csg_filt(seis%seis,seis%nt,seis%dt,seis%ng,f11,f22,f33,f44)
   !write(*,*)time_now," RTM, ipro=", ipro," npro_me=",npro_me
   if (rank.eq.0) then
     write(*,*)" RTM, ipro=", ipro," npro_me=",npro_me
    ! write(*,*) size(img%img_h,3)
   endif
   if (f%fd_type.eq."2ND") then ! Center FD
      !call a2d_mig(seis,bc,wf,fdm,img,output)
      call a2d_mig(seis,bc,wf,fdm,img,output,o%nh)
   elseif (f%fd_type.eq."SG") then ! Staggerred FD
     ! call a2d_mig_SG(seis,bc,wf,fdm,img,output)  ! modified by Bowen Guo
      !call a2d_mig(seis,bc,wf,fdm,img,output)
      call a2d_mig(seis,bc,wf,fdm,img,output,o%nh)
   endif
   call output_illum(o%isIllum,o%isSaveIllum,fn%pre_illum_out,&
           c%proid(ipro),img%illum,img%ix1,img%ix2,img%nz,img%nx)
   call output_prestack_image(o%isPreStackImg,fn%pre_illum_img_out,&
           c%proid(ipro),img%img,img%ix1,img%ix2,img%nz,img%nx)

  ! call allocate_and_initial(temp,m%nz,m%nx,o%nh)
  ! call aperture_image_l2g(a%isAper,img%img,img%iz1,img%iz2,img%ix1,img%ix2,&
  !                         temp,m%nz,m%nx)
   
    image_temp=image_temp+img%img_h
 !  call deallocate_and_free(temp)
   if (o%isIllum) then
      do ih=1,o%nh
       !img%img_h(:,:,ih)=img%img_h(:,:,ih)/img%illum
       img%img_h(:,:,ih)=img%img_h(:,:,ih)/(img%illum+0.001*maxval(img%illum))
      enddo
  !    call allocate_and_initial(temp,m%nz,m%nx)
  !    call aperture_image_l2g(a%isAper,img%img,img%iz1,img%iz2,img%ix1,img%ix2,&
  !                            temp,m%nz,m%nx)
      image_illum_temp=image_illum_temp+img%img_h
 !     call deallocate_and_free(temp)
   endif
   call aperture_finalize(seis,img,fdm)
enddo

!if (rank.eq.0) then
!   call allocate_and_initial(image,m%nz,m%nx)
!endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)


 call sf90_sum(image_temp,image,m%nz,m%nx,o%nh)
if (o%isIllum) then
 call sf90_sum(image_illum_temp,image,m%nz,m%nx,o%nh)
endif
!call MPI_REDUCE(image_temp,image,m%nz*m%nx,MPI_REAL,MPI_SUM,&
!        0,MPI_COMM_WORLD,ierr)


if (rank.eq.0) then
   write(*,*) "Migration Finished, Output to Files..."
   call write_binfile(fn%img_file,image,m%nz,m%nx,o%nh)
   call deallocate_and_free(image)
endif

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
deallocate(image_temp)
if (o%isIllum) then
   if (rank.eq.0) then
      call allocate_and_initial(image_illum,m%nz,m%nx)
   endif
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(image_illum_temp,image_illum,m%nz*m%nx,&
            MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   if (rank.eq.0) then
      call write_binfile(fn%img_illum_file,image_illum,m%nz,m%nx)
      call deallocate_and_free(image_illum)
   endif
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   deallocate(image_illum_temp)
endif
call deallocate_type(f,m,s)
call deallocate_type(c)
deallocate(ipro_me)
999 continue
if (rank.eq.0) call log_end(99)
call stop_mpi
call cpu_time(t2)
if (rank == 0) write(*,*) 'Runtime = ', t2-t1
if (rank == 0) write(*,*) 'Runtime = ', (t2-t1)/nsize

contains

subroutine display_parameters(m,f,o,c,a,fn,s)
implicit none
type(model2d),       intent(in) :: m
type(fd2d),          intent(in) :: f
type(other),         intent(in) :: o
type(coord2d),       intent(in) :: c
type(aperture),      intent(in) :: a
type(fname),         intent(in) :: fn
type(source),        intent(in) :: s
write(*,*)"nx=",m%nx, " nz=",m%nz, " npml=",f%npml," dx=",m%dx
write(*,*)"nt=",f%nt," nt_out=",o%nt_out," dt=",f%dt
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
write(*,*)"Coordinats file type:",o%coordfile_type
write(*,*)"Coordinats file:",fn%coord_in_file
write(*,*)"Topography file:",fn%topo_file
write(*,*)
write(*,*)"VMIN=",m%vmin," ,VMAX=",m%vmax
write(*,*)"DX REQUIRED=VMIN/(15*FREQ)=",m%vmin/(15.0*s%freq)
if (m%dx>m%vmin/(15.0*s%freq)) then
   write(*,*) "WARRING!!! Stability condition MAY BE VIOLATED"
endif
write(*,*)"DT REQUIRED=0.6*DX/VMAX=",0.6*m%dx/m%vmax
if (f%dt>0.6*m%dx/m%vmax) then
   write(*,*) "WARRING!!! Stability condition MAY BE VIOLATED"
endif
write(*,*)
if (f%fd_type.eq."2ND") write(*,*)"Finite difference type: 2nd Order FD"
if (f%fd_type.eq."SG") write(*,*)"Finite difference type: Staggerred FD"
if (f%bc_type==1) write(*,*)"Boundary condition type: Absorbing"
if (f%bc_type==2) write(*,*)"Boundary condition type: Sponge"
write(*,*)"Finite Difference Accuracy=",f%fd_order
write(*,*)
write(*,*)"IsMakeSource=",s%isSource
if (s%isSource) write(*,*) "Frequency=",s%freq
if (.not.s%isSource) write(*,*) "Source file:",fn%source_file
write(*,*)
write(*,*)"IsFreeSurface=",f%isFS
if (.not.f%isFS) write(*,*) "IsDipole=",f%isDipole
write(*,*)"IsSaveBC=",f%isSaveBC
write(*,*)"IsSaveWF=",f%isSaveWF
write(*,*)"CSG_OUT File=",fn%csg_out
write(*,*)"Is Read Data during Calculation=",f%isReadData
write(*,*)"IsPreStackImg=",o%isPreStackImg
write(*,*)"IsIllum=",o%isIllum
write(*,*)"IMAGE_CONDITION=",f%ic
if (o%isIllum) write(*,*)"IsSaveIllum",o%isSaveIllum
end subroutine display_parameters

subroutine log_parameters(m,f,o,c,a,s,fn)
implicit none
type(model2d),   intent(in) :: m
type(fd2d),      intent(in) :: f         
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
write(99,*)"nx=",m%nx, " nz=",m%nz, " npml=",f%npml," dx=",m%dx
write(99,*)"nt=",f%nt," nt_out=",o%nt_out," dt=",f%dt
write(99,*)"ns=",c%npro," ngmax=",c%ngmax
write(99,*)
write(99,*)"x0=",a%x0," z0=",a%z0
write(*,*)"IsAper=",a%isAper
if (a%isAper) then
   write(99,*) "Aper_x=",a%aper_x
   write(99,*) "Aper_z=",a%aper_z
endif
write(99,*)
write(99,*)"Velocity file:",fn%vel_file
write(99,*)"Coordinats file type:",o%coordfile_type
write(99,*)"Coordinats file:",fn%coord_in_file
write(99,*)"Topography file:",fn%topo_file
write(99,*)
if (f%fd_type.eq."2ND") write(99,*)"Finite difference type: 2nd Order FD"
if (f%fd_type.eq."SG") write(99,*)"Finite difference type: Staggerred FD"
if (f%bc_type==1) write(99,*)"Boundary condition type: Absorbing"
if (f%bc_type==2) write(99,*)"Boundary condition type: Sponge"
write(99,*)"Finite Difference Accuracy=",f%fd_order
write(99,*)
write(99,*)"IsMakeSource=",s%isSource
if (s%isSource) write(99,*) "Frequency=",s%freq
if (.not.s%isSource) write(99,*) "Source file:",fn%source_file
write(99,*)
write(99,*)"IsFreeSurface=",f%isFS
if (.not.f%isFS) write(99,*) "IsDipole=",f%isDipole
write(99,*)"IsSaveBC=",f%isSaveBC
write(99,*)"IsSaveWF=",f%isSaveWF
write(99,*)"CSG_OUT File=",fn%csg_out
write(99,*)"Is Read Data during Calculation=",f%isReadData
write(99,*)"IsPreStackImg=",o%isPreStackImg
write(99,*)"IsIllum=",o%isIllum
if (o%isIllum) write(99,*)"IsSaveIllum",o%isSaveIllum
end subroutine log_parameters

end program a2d_rtm_h
