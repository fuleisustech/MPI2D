program a2d_modeling
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
use module_math        ! modified by Bowen Guo
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
type(bc2d_general)   :: bc
type(wf2d)           :: wf
integer              :: ipro,npro_me
real                 :: t1, t2
logical              :: isError
integer, allocatable :: ipro_me(:)
character(len=slen)   :: output
call cpu_time(t1)
call readCmd("par",parfile,"parfile")
call init_mpi
! Read in parameters
isError=.false.
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
m%vmin=minval(m%v);   m%vmax=maxval(m%v)
   !call allocate_and_initial(m%den,m%nz,m%nx)
   !call read_binfile_mpi(fn%den_file, m%den, m%nz,m%nx)
if (f%fd_type.eq."SG") then 
   call allocate_and_initial(m%den,m%nz,m%nx)
   call read_binfile_mpi(fn%den_file, m%den, m%nz,m%nx)
endif   
if (rank.eq.0) then 
   write(*,*) "Acoustic WE Modeling"         ! modified by Bowen Guo
   call display_parameters(m,f,o,c_all,a,fn,s)
   call log_parameters(m,f,o,c_all,a,s,fn)
endif
! Setup Ricker source wavelet
call allocate_and_initial(s%sou,s%nw)
if (s%isSource) then   
   call ricker(s)
   !call integrate(s%sou,s%nw) 
   s%sou = s%sou*1.0/(maxval(abs(s%sou)))
   if (f%fd_type.eq."SG") call integrate(s%sou,s%nw) ! modified by Bowen Guo
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
! Modeling kernal
do ipro=1,npro_me
   
   call aperture_init(c,a,m,f,s,o,ipro,seis,fdm)
   if (.not.f%isWriteData) then
      call allocate_and_initial(seis%seis,seis%nt,seis%ng)
   endif
   call date_and_time(date,time_now)
   if (rank.eq.0) then 
      write(*,*)time_now,", sid=",c%proid(ipro)," ipro=",ipro," npro_me=",npro_me
   endif
   call filename(output, fn%csg_out, c%proid(ipro), ".bin")
   if ((f%fd_type.eq."2ND").or.(f%fd_type.eq."Lax_wen")) then ! Center FD
      call a2d_mod(fdm,bc,wf,seis,output) 
   elseif (f%fd_type.eq."SG") then ! Staggerred FD
     ! call a2d_mod_SG(fdm,bc,wf,seis,output)  ! modified by Bowen Guo
      call a2d_mod(fdm,bc,wf,seis,output) 
   endif
   if (.not.f%isWriteData) then
      call write_binfile(output,seis%seis,seis%nt, seis%ng)
   endif
   call aperture_finalize(seis,fdm)
enddo
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
use module_datatype
implicit none
type(model2d),       intent(in) :: m
type(fd2d),          intent(in) :: f
type(other),         intent(in) :: o
type(coord2d),       intent(in) :: c
type(aperture),      intent(in) :: a
type(fname),         intent(in) :: fn
type(source),        intent(in) :: s
write(*,*)"nx=",m%nx, " nz=",m%nz, " npml=",f%npml," dx=",m%dx
write(*,*)"nt=",f%nt," dt=",f%dt," nt_out=",o%nt_out," dt_out=",o%dt_out
write(*,*)"ns=",c%npro," ngmax=",c%ngmax
write(*,*)
write(*,*)"x0=",a%x0," z0=",a%z0
write(*,*)"IsAper=",a%isAper
if (a%isAper) then
   write(*,*) "Aper_X=",a%aper_x
   write(*,*) "Aper_Z=",a%aper_z
endif
write(*,*)
write(*,*)"Velocity file:",fn%vel_file
write(*,*)"Coordinats file Type:",o%coordfile_type
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
if (f%fd_type.eq."Lax_wen") write(*,*)"Finite difference type: 2nd Order FD with Lax-Wendroff method"
if (f%fd_type.eq."SG") write(*,*)"Finite difference type: Staggerred FD"
if (f%bc_type==1) write(*,*)"Boundary condition type: Absorbing"
if (f%bc_type==2) write(*,*)"Boundary condition type: Sponge"
write(*,*)"Finite Difference Accuracy=",f%fd_order
write(*,*)
write(*,*)"IsMakeSource=",s%isSource
if (s%isSource) write(*,*) "Frequency=",s%freq
if (.not.s%isSource) write(99,*) "Source file:",fn%source_file
write(*,*)
write(*,*)"IsFreeSurface=",f%isFS
if (.not.f%isFS) write(*,*) "IsDipole=",f%isDipole
write(*,*)"CSG_OUT File=",fn%csg_out
write(*,*)"Is Write Data during Calculation=",f%isWriteData
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
write(99,*)"nt=",f%nt," dt=",f%dt," nt_out=",o%nt_out," dt_out=",o%dt_out
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
write(99,*)"Coordinats file Type:",o%coordfile_type
write(99,*)"Coordinats file:",fn%coord_in_file
write(99,*)"Topography file:",fn%topo_file
write(99,*)
if (f%fd_type.eq."2ND") write(99,*)"Finite difference type: 2nd Order FD"
if (f%fd_type.eq."Lax_wen") write(99,*)"Finite difference type: 2nd Order FD with Lax-Wendroff method"
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
write(99,*)"CSG_OUT File=",fn%csg_out
write(99,*)"Is Write Data during Calculation=",f%isWriteData
end subroutine log_parameters

end program a2d_modeling
