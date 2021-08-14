module module_acquisition2d
use module_global, only : I4, slen
use module_datatype
use module_sf90_mpi
use module_array
use module_utility
use module_string
implicit none
logical, private :: isAscii
integer, private :: proid0,proid1,proid,sid,gid,npro,ng,ngmax,&
  ipro,is,ig,itrace,ip,ip1,itmp1,ns_total,ng_total,rec0
real,    private :: xs, zs, xg, zg, t, zp, zp1, zp2, &
   theta_x1, theta_x2, xref_1, xref_2,tmp1,tmp2,tmp3,tmp4

interface read_coordfile_all
   module procedure read_gen_coordfile_all
   module procedure read_ms_coordfile_all
   module procedure read_pw1_coordfile_all
   module procedure read_pw2_coordfile_all
end interface read_coordfile_all

interface read_coordfile_all_mpi
   module procedure read_gen_coordfile_all_mpi
   module procedure read_ms_coordfile_all_mpi
   module procedure read_pw1_coordfile_all_mpi
   module procedure read_pw2_coordfile_all_mpi
end interface read_coordfile_all_mpi

interface write_coordfile_all
   module procedure write_gen_coordfile_all
   module procedure write_ms_coordfile_all
   module procedure write_pw1_coordfile_all
   module procedure write_pw2_coordfile_all
end interface write_coordfile_all

interface read_coordfile_single
   module procedure read_gen_coordfile_single
   module procedure read_ms_coordfile_single
   module procedure read_pw1_coordfile_single
   module procedure read_pw2_coordfile_single
end interface read_coordfile_single

interface write_coordfile_single
   module procedure write_gen_coordfile_single
   module procedure write_ms_coordfile_single
   module procedure write_pw1_coordfile_single
   module procedure write_pw2_coordfile_single
end interface write_coordfile_single

interface read_coordfile
   module procedure read_gen_coordfile
   module procedure read_eik_coordfile
   module procedure read_ms_coordfile
   module procedure read_pw1_coordfile
   module procedure read_pw2_coordfile
end interface read_coordfile

interface read_coordfile_mpi
   module procedure read_gen_coordfile_mpi
   module procedure read_eik_coordfile_mpi
   module procedure read_ms_coordfile_mpi
   module procedure read_pw1_coordfile_mpi
   module procedure read_pw2_coordfile_mpi
end interface read_coordfile_mpi

interface write_coordfile
   module procedure write_gen_coordfile
   module procedure write_ms_coordfile
   module procedure write_pw1_coordfile
   module procedure write_pw2_coordfile
end interface write_coordfile

interface read_topofile
   module procedure read_topofile
end interface read_topofile

interface read_topofile_mpi
   module procedure read_topofile_mpi
end interface read_topofile_mpi

interface read_pw1_par
   module procedure read_pw1_par
end interface read_pw1_par

interface read_pw1_par_mpi
   module procedure read_pw1_par_mpi
end interface read_pw1_par_mpi

interface read_pw2_par
   module procedure read_pw2_par
end interface read_pw2_par

interface read_pw2_par_mpi
   module procedure read_pw2_par_mpi
end interface read_pw2_par_mpi

interface assign_coord
   module procedure assign_gen_coord
   module procedure assign_ms_coord
   module procedure assign_pw1_coord
   module procedure assign_pw2_coord
end interface assign_coord

private get_isAscii,read_coordfile_all_core,write_coordfile_all_core

contains

!--------------------------------------------------------------------------
subroutine read_coordfile_all_core(coordfile,coordfile_type,&
   npro,proid,ngmax,ng,ntrace)
character(len=*),intent(in) :: coordfile,coordfile_type
integer, intent(out) :: npro,ngmax,ntrace
integer, allocatable, intent(inout) :: proid(:),ng(:)
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="read")
else
   open(10,file=coordfile,access="direct",recl=2*I4)
endif
! Determine the number of profiles
if (isAscii) then
   read(10,*,end=100) proid0, itmp1
else
   read(10,rec=1,err=100) proid0, itmp1
endif
npro = 1
do itrace=2,1000000000
   if (isAscii) then
      read(10,*,end=100) proid1, itmp1 
   else
      read(10,rec=itrace,err=100) proid1, itmp1
   endif
   if (proid1.ne.proid0) then
      proid0 = proid1;   npro = npro + 1
   endif
enddo
100 continue
rewind(10)
call allocate_and_initial(ng,proid,npro)
! Read all information 
do ipro=1,npro
   if (isAscii) then
      read(10,*,end=200) proid(ipro), ng(ipro)
   else
      read(10,rec=ipro,err=200) proid(ipro), ng(ipro)
   endif
enddo
200 continue
close(10)
ngmax=maxval(ng);   ntrace=sum(ng)
end subroutine read_coordfile_all_core

subroutine write_coordfile_all_core(coordfile,coordfile_type,npro,proid,ng)
character(len=*),intent(in) :: coordfile,coordfile_type
integer, intent(in) :: npro,proid(:),ng(:)
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form='formatted')
else
   open(10,file=coordfile,access="direct",recl=2*I4)
endif
do ipro=1,npro
   if (isAscii) then
      write(10,100)proid(ipro),ng(ipro)
   else
      write(10,rec=ipro)proid(ipro),ng(ipro)
   endif
enddo
100 format(I7," ",I7)
close(10)
end subroutine write_coordfile_all_core

subroutine read_gen_coordfile(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(coord2d),    intent(out)  :: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="read")
else
   open(10,file=coordfile,access="direct",recl=8*I4,action="read")
endif
! Determine the number of profiles
if (isAscii) then
   read(10,*,end=100) proid0, sid, gid, xs, zs, xg, zg, t
else
   read(10,rec=1,err=100) proid0, sid, gid, xs, zs, xg, zg, t
endif
npro = 1
do itrace=2,1000000000
   if (isAscii) then
      read(10,*,end=100) proid, sid, gid, xs, zs, xg, zg, t
   else
      read(10,rec=itrace,err=100) proid, sid, gid, xs, zs, xg, zg, t
   endif
   if (proid.ne.proid0) then
      proid0 = proid
      npro = npro + 1
   endif
enddo
100 continue
rewind(10)
c%npro = npro
call allocate_and_initial(c%ng,c%proid,npro)
! Determine the number of traces for each profile
if (isAscii) then
   read(10,*,end=100) proid0, sid, gid, xs, zs, xg, zg, t
else
   read(10,rec=1,err=100) proid0, sid, gid, xs, zs, xg, zg, t
endif
npro = 1;   ng = 1
do itrace=2,1000000000
   if (isAscii) then
      read(10,*,end=200) proid, sid, gid, xs, zs, xg, zg, t
   else
      read(10,rec=itrace,err=200) proid, sid, gid, xs, zs, xg, zg, t
   endif
   if (proid.ne.proid0) then
      proid0 = proid; c%ng(npro) = ng; npro = npro + 1; ng = 1
   else
      ng = ng + 1
   endif
enddo
200 continue
c%ng(npro) = ng
rewind(10)
c%ntrace=sum(c%ng(1:c%npro));   c%ngmax = maxval(c%ng(1:c%npro))
call allocate_and_initial(c%sid,c%gid,c%ngmax,c%npro)
call allocate_and_initial(c%xs,c%zs,c%xg,c%zg,c%t,c%ngmax,c%npro)
! Read all information 
itrace=0
do ipro=1,c%npro
   do ig=1,c%ng(ipro)
      if (isAscii) then
         read(10,*,end=300) c%proid(ipro), c%sid(ig,ipro), c%gid(ig,ipro), &
            c%xs(ig,ipro),c%zs(ig,ipro),c%xg(ig,ipro),c%zg(ig,ipro),c%t(ig,ipro)
      else
         itrace=itrace+1
         read(10,rec=itrace,err=300) c%proid(ipro), c%sid(ig,ipro), c%gid(ig,ipro), &
            c%xs(ig,ipro),c%zs(ig,ipro),c%xg(ig,ipro),c%zg(ig,ipro),c%t(ig,ipro)
      endif
   enddo
enddo
300 continue
close(10)
end subroutine read_gen_coordfile

subroutine read_gen_coordfile_mpi(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(coord2d),    intent(out)  :: c
if (rank == 0) call read_gen_coordfile(coordfile, coordfile_type, c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro,c%ngmax,c%ntrace)
if (rank > 0) then
   call allocate_and_initial(c%proid,c%ng,c%npro)
   call allocate_and_initial(c%sid,c%gid,c%ngmax,c%npro)
   call allocate_and_initial(c%xs,c%zs,c%xg,c%zg,c%t,c%ngmax,c%npro)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%ng,c%npro)
call sf90_bcast(c%sid,c%gid,c%ngmax,c%npro)
call sf90_bcast(c%xs,c%zs,c%xg,c%zg,c%t,c%ngmax,c%npro)
end subroutine read_gen_coordfile_mpi

subroutine write_gen_coordfile(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(coord2d),   intent(in):: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form='formatted')
else
   open(10,file=coordfile,access="direct",recl=8*I4)
endif
itrace=0
do ipro=1,c%npro
   do ig=1,c%ng(ipro)
      if (isAscii) then
         write(10,100)c%proid(ipro),c%sid(ig,ipro), c%gid(ig,ipro), c%xs(ig,ipro),&
            c%zs(ig,ipro), c%xg(ig,ipro), c%zg(ig,ipro), c%t(ig,ipro)
      else
         itrace=itrace+1
         write(10,rec=itrace)c%proid(ipro),c%sid(ig,ipro), c%gid(ig,ipro), &
            c%xs(ig,ipro),c%zs(ig,ipro), c%xg(ig,ipro), c%zg(ig,ipro), c%t(ig,ipro)
      endif
   enddo
enddo
100 format(I7," ",I7," ",I7,4(" ",F8.2)," ",F8.5)
close(10)
end subroutine write_gen_coordfile

subroutine read_gen_coordfile_single(coordfile,coordfile_type,c)
character(len=*), intent(in)    :: coordfile,coordfile_type
type(coord2d),    intent(inout) :: c
call allocate_and_initial(c%sid,c%gid,c%ng(1),1)
call allocate_and_initial(c%xs,c%zs,c%xg,c%zg,c%t,c%ng(1),1)
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(11,file=coordfile,form="formatted",action="read")
else
   open(11,file=coordfile,access="direct",recl=8*I4,action="read")
endif
do ig=1,c%ng(1)
   if (isAscii) then
      read(11,*)c%proid(1),c%sid(ig,1),c%gid(ig,1),c%xs(ig,1),&
         c%zs(ig,1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1)
   else
      read(11,rec=ig)c%proid,c%sid(ig,1),c%gid(ig,1),c%xs(ig,1),&
         c%zs(ig,1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1)
   endif
enddo
close(11)
end subroutine read_gen_coordfile_single

subroutine write_gen_coordfile_single(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(coord2d),   intent(in):: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form='formatted')
else
   open(10,file=coordfile,access="direct",recl=8*I4)
endif
do ig=1,c%ng(1)
   if (isAscii) then
      write(10,100)c%proid(1),c%sid(ig,1), c%gid(ig,1), c%xs(ig,1),&
         c%zs(ig,1), c%xg(ig,1), c%zg(ig,1), c%t(ig,1)
   else
      write(10,rec=ig)c%proid(1),c%sid(ig,1), c%gid(ig,1), c%xs(ig,1),&
         c%zs(ig,1), c%xg(ig,1), c%zg(ig,1), c%t(ig,1)
   endif
enddo
100 format(I7," ",I7," ",I7,4(" ",F8.2)," ",F8.5)
close(10)
end subroutine write_gen_coordfile_single

subroutine read_gen_coordfile_all(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(coord2d),    intent(out)  :: c
call read_coordfile_all_core(coordfile,coordfile_type,c%npro,&
   c%proid,c%ngmax,c%ng,c%ntrace)
end subroutine read_gen_coordfile_all

subroutine read_gen_coordfile_all_mpi(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(coord2d),    intent(out)  :: c
if (rank == 0)   call read_coordfile_all(coordfile, coordfile_type, c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro,c%ngmax,c%ntrace)
if (rank > 0)   call allocate_and_initial(c%proid,c%ng,c%npro)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%ng,c%npro)
end subroutine read_gen_coordfile_all_mpi

subroutine write_gen_coordfile_all(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(coord2d),   intent(in):: c
call write_coordfile_all_core(coordfile,coordfile_type,c%npro,c%proid,c%ng)
end subroutine write_gen_coordfile_all

subroutine assign_gen_coord(c_all,c,npro_me,ipro_me)
type(coord2d), intent(in)    :: c_all
type(coord2d), intent(inout) :: c
integer,       intent(in)    :: npro_me, ipro_me(:)
c%npro=npro_me
call allocate_and_initial(c%proid,c%ng,npro_me)
do ipro=1,npro_me
   c%proid(ipro)=c_all%proid(ipro_me(ipro))
   c%ng(ipro)=c_all%ng(ipro_me(ipro))
enddo
c%ngmax=maxval(c%ng);  c%ntrace=sum(c%ng)
call allocate_and_initial(c%sid,c%gid,c%ngmax,c%npro)
call allocate_and_initial(c%xs,c%zs,c%xg,c%zg,c%t,c%ngmax,c%npro)
do ipro=1,c%npro
   c%sid(1:c%ng(ipro),ipro)=c_all%sid(1:c%ng(ipro),ipro_me(ipro))
   c%gid(1:c%ng(ipro),ipro)=c_all%gid(1:c%ng(ipro),ipro_me(ipro))
   c%xs(1:c%ng(ipro),ipro)=c_all%xs(1:c%ng(ipro),ipro_me(ipro))
   c%zs(1:c%ng(ipro),ipro)=c_all%zs(1:c%ng(ipro),ipro_me(ipro))
   c%xg(1:c%ng(ipro),ipro)=c_all%xg(1:c%ng(ipro),ipro_me(ipro))
   c%zg(1:c%ng(ipro),ipro)=c_all%zg(1:c%ng(ipro),ipro_me(ipro))
   c%t(1:c%ng(ipro),ipro)=c_all%t(1:c%ng(ipro),ipro_me(ipro))
enddo
end subroutine assign_gen_coord

!--------------------------------------------------------------------------

subroutine read_ms_coordfile(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile, coordfile_type
type(ms_coord2d), intent(out)  :: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="read")
else
   open(10,file=coordfile,access="direct",recl=5*I4,action="read")
endif
! Read the number of supergather
rec0=0
if (isAscii) then
   read(10,*)c%npro
else
   read(10,rec=rec0+1)c%npro,tmp1,tmp2,tmp3,tmp4
endif
call allocate_and_initial(c%proid,c%ns,c%ng,c%npro)
rec0=rec0+1
do ipro=1,c%npro
   if (isAscii) then
      read(10,*)c%proid(ipro),c%ns(ipro),c%ng(ipro)
   else
      read(10,rec=ipro+rec0)c%proid(ipro),c%ns(ipro),c%ng(ipro),tmp1,tmp2
   endif
enddo
! Determine the maximum number of shots and geophones per supergather
c%nsmax=maxval(c%ns); c%ngmax=maxval(c%ng); c%ntrace=sum(c%ng)
call allocate_and_initial(c%sid,c%nsmax,c%npro)
call allocate_and_initial(c%xs,c%zs,c%d,c%p,c%nsmax,c%npro)
rec0=rec0+c%npro
ns_total=0
do ipro=1,c%npro
   do is=1,c%ns(ipro)
      if (isAscii) then
         read(10,*)c%sid(is,ipro),c%xs(is,ipro),c%zs(is,ipro),&
                   c%d(is,ipro),c%p(is,ipro)
      else
         ns_total=ns_total+1
         read(10,rec=rec0+ns_total)c%sid(is,ipro),c%xs(is,ipro),&
            c%zs(is,ipro),c%d(is,ipro),c%p(is,ipro)
      endif
   enddo
enddo
call allocate_and_initial(c%gid,c%ngmax,c%npro)
call allocate_and_initial(c%xg,c%zg,c%t,c%ngmax,c%npro)
rec0=rec0+ns_total
ng_total=0
do ipro=1,c%npro
   do ig=1,c%ng(ipro)
      if (isAscii) then
         read(10,*)c%gid(ig,ipro),c%xg(ig,ipro),c%zg(ig,ipro),c%t(ig,ipro)
      else
         ng_total=ng_total+1
         read(10,rec=rec0+ng_total)c%gid(ig,ipro),c%xg(ig,ipro),&
            c%zg(ig,ipro),c%t(ig,ipro),tmp1,tmp2
      endif
   enddo
enddo
close(10)
end subroutine read_ms_coordfile

subroutine read_ms_coordfile_mpi(coordfile,coordfile_type,c)
character(len=*), intent(in)   :: coordfile, coordfile_type
type(ms_coord2d), intent(out)  :: c
if (rank.eq.0)  call read_ms_coordfile(coordfile,coordfile_type,c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro,c%nsmax,c%ngmax,c%ntrace)
if (rank > 0) then
   call allocate_and_initial(c%proid,c%ns,c%ng,c%npro)
   call allocate_and_initial(c%sid,c%nsmax,c%npro)
   call allocate_and_initial(c%xs,c%zs,c%d,c%p,c%nsmax,c%npro)
   call allocate_and_initial(c%gid,c%ngmax,c%npro)
   call allocate_and_initial(c%xg,c%zg,c%t,c%ngmax,c%npro)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%ns,c%ng,c%npro)
call sf90_bcast(c%sid,c%nsmax,c%npro)
call sf90_bcast(c%xs,c%zs,c%d,c%p,c%nsmax,c%npro)
call sf90_bcast(c%gid,c%ngmax,c%npro)
call sf90_bcast(c%xg,c%zg,c%t,c%ngmax,c%npro)
end subroutine read_ms_coordfile_mpi

subroutine write_ms_coordfile(coordfile,coordfile_type,c)
character(len=*), intent(in) :: coordfile,coordfile_type
type(ms_coord2d), intent(in):: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="write")
else
   open(10,file=coordfile,access="direct",recl=5*I4,action="write")
endif
if (isAscii) then
   write(10,50)c%npro
else
   write(10,rec=1)c%npro,tmp1,tmp2,tmp3,tmp4
endif
50 format(I4)
rec0=rec0+1
do ipro=1,c%npro
   if (isAscii) then
      write(10,100)c%proid(ipro),c%ns(ipro),c%ng(ipro)
   else
      write(10,rec=1+ipro)c%proid(ipro),c%ns(ipro),c%ng(ipro),tmp1,tmp2
   endif
enddo
100 format(3(I7))
rec0=1+c%npro
ns_total=0
do ipro=1,c%npro
   do is=1,c%ns(ipro)
      if (isAscii) then
         write(10,200)c%sid(is,ipro),c%xs(is,ipro),c%zs(is,ipro),&
                      c%d(is,ipro),c%p(is,ipro)
      else
         ns_total=ns_total+1
         write(10,rec=rec0+ns_total)c%sid(is,ipro),c%xs(is,ipro),&
            c%zs(is,ipro),c%d(is,ipro),c%p(is,ipro)
      endif
   enddo
enddo
200 format(I7,2(" ",F10.2)," ",F8.4," ",F5.2)
rec0=rec0+ns_total
ng_total=0
do ipro=1,c%npro
   do ig=1,c%ng(ipro)
      if (isAscii) then
         write(10,300)c%gid(ig,ipro),c%xg(ig,ipro),c%zg(ig,ipro),c%t(ig,ipro)
      else
         ng_total=ng_total+1
         write(10,rec=rec0+ng_total)c%gid(ig,ipro),c%xg(ig,ipro),&
            c%zg(ig,ipro),c%t(ig,ipro),tmp1
      endif
   enddo
enddo
300 format(I7,2(" ",F10.2)," ",F8.4)
close(10)
end subroutine write_ms_coordfile

subroutine read_ms_coordfile_single(coordfile,coordfile_type,c)
character(len=*),  intent(in)    :: coordfile,coordfile_type
type(ms_coord2d),  intent(inout) :: c
call allocate_and_initial(c%ns,c%ng,1)
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(11,file=coordfile,form="formatted",action="read")
   read(11,*)c%ns(1),c%ng(1)
else
   open(11,file=coordfile,access="direct",recl=5*I4,action="read")
   read(11,rec=1)c%ns(1),c%ng(1),tmp1,tmp2,tmp3
endif
call allocate_and_initial(c%sid,c%ns(1),1)
call allocate_and_initial(c%gid,c%ng(1),1)
call allocate_and_initial(c%xs,c%zs,c%d,c%p,c%ns(1),1)
call allocate_and_initial(c%xg,c%zg,c%t,c%ng(1),1)
do is=1,c%ns(1)
   if (isAscii) then
      read(11,*)c%sid(is,1),c%xs(is,1),c%zs(is,1),c%d(is,1),c%p(is,1)
   else
      read(11,rec=is+1)c%sid(is,1),c%xs(is,1),c%zs(is,1),c%d(is,1),c%p(is,1)
   endif
enddo
do ig=1,c%ng(1)
   if (isAscii) then
      read(11,*)c%gid(ig,1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1)
   else
      read(11,rec=1+c%ns(1)+ig)c%gid(ig,1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1),tmp1
   endif
enddo
close(11)
end subroutine read_ms_coordfile_single

subroutine write_ms_coordfile_single(coordfile,coordfile_type,c)
character(len=*), intent(in) :: coordfile,coordfile_type
type(ms_coord2d), intent(in):: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="write")
   write(10,100)c%ns(1),c%ng(1)
else
   open(10,file=coordfile,access="direct",recl=5*I4,action="write")
   write(10,rec=1)c%ns(1),c%ng(1),tmp1,tmp2,tmp3
endif
100 format(2(I7))
ns_total=0
do is=1,c%ns(1)
   if (isAscii) then
      write(10,200)c%sid(is,1),c%xs(is,1),c%zs(is,1),c%d(is,1),c%p(is,1)
   else
      ns_total=ns_total+1
      write(10,rec=1+ns_total)c%sid(is,1),c%xs(is,1),c%zs(is,1),c%d(is,1),c%p(is,1)
   endif
enddo
200 format(I7,2(" ",F10.2)," ",F8.4," ",F5.2)
rec0=1+c%ns(1)   
ng_total=0
do ig=1,c%ng(1)
   if (isAscii) then
      write(10,300)c%gid(ig,1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1)
   else
      ng_total=ng_total+1
      write(10,rec=rec0+ng_total)c%gid(ig,1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1),tmp1
   endif
enddo
300 format(I7,2(" ",F10.2)," ",F8.4)
close(10)
end subroutine write_ms_coordfile_single

subroutine read_ms_coordfile_all(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(ms_coord2d),    intent(out)  :: c
call read_coordfile_all_core(coordfile,coordfile_type,c%npro,&
   c%proid,c%ngmax,c%ng,c%ntrace)
end subroutine read_ms_coordfile_all

subroutine read_ms_coordfile_all_mpi(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(ms_coord2d),    intent(out)  :: c
if (rank == 0)   call read_coordfile_all(coordfile, coordfile_type, c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro,c%ngmax,c%ntrace)
if (rank > 0)   call allocate_and_initial(c%proid,c%ng,c%npro)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%ng,c%npro)
end subroutine read_ms_coordfile_all_mpi

subroutine write_ms_coordfile_all(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(ms_coord2d),   intent(in):: c
call write_coordfile_all_core(coordfile,coordfile_type,c%npro,c%proid,c%ng)
end subroutine write_ms_coordfile_all

subroutine assign_ms_coord(c_all,c,npro_me,ipro_me)
type(ms_coord2d), intent(in)    :: c_all
type(ms_coord2d), intent(inout) :: c
integer,          intent(in)    :: npro_me, ipro_me(:)
c%npro=npro_me
call allocate_and_initial(c%proid,c%ns,c%ng,npro_me)
do ipro=1,npro_me
   c%proid(ipro)=c_all%proid(ipro_me(ipro))
   c%ns(ipro)=c_all%ns(ipro_me(ipro))
   c%ng(ipro)=c_all%ng(ipro_me(ipro))
enddo
c%nsmax=maxval(c%ns);   c%ngmax=maxval(c%ng);   c%ntrace=sum(c%ng)
call allocate_and_initial(c%sid,c%nsmax,c%npro)
call allocate_and_initial(c%xs,c%zs,c%d,c%p,c%nsmax,c%npro)
call allocate_and_initial(c%gid,c%ngmax,c%npro)
call allocate_and_initial(c%xg,c%zg,c%t,c%ngmax,c%npro)
do ipro=1,c%npro
   c%sid(1:c%ns(ipro),ipro) = c_all%sid(1:c%ns(ipro),ipro_me(ipro))
   c%xs(1:c%ns(ipro),ipro)  = c_all%xs(1:c%ns(ipro),ipro_me(ipro))
   c%zs(1:c%ns(ipro),ipro)  = c_all%zs(1:c%ns(ipro),ipro_me(ipro))
   c%d(1:c%ns(ipro),ipro)   = c_all%d(1:c%ns(ipro),ipro_me(ipro))
   c%p(1:c%ns(ipro),ipro)   = c_all%p(1:c%ns(ipro),ipro_me(ipro))
   c%gid(1:c%ng(ipro),ipro) = c_all%gid(1:c%ng(ipro),ipro_me(ipro))
   c%xg(1:c%ng(ipro),ipro)  = c_all%xg(1:c%ng(ipro),ipro_me(ipro))
   c%zg(1:c%ng(ipro),ipro)  = c_all%zg(1:c%ng(ipro),ipro_me(ipro))
   c%t(1:c%ng(ipro),ipro)   = c_all%t(1:c%ng(ipro),ipro_me(ipro))
enddo
end subroutine assign_ms_coord

!--------------------------------------------------------------------------

subroutine read_pw1_coordfile(coordfile,coordfile_type,c)
character(len=*),   intent(in)   :: coordfile,coordfile_type
type(pw1_coord2d),  intent(out)  :: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="read")
else
   open(10,file=coordfile,access="direct",recl=8*I4,action="read")
endif
! Determine the number of profiles
npro = 0
if (isAscii) then
   read(10,*,end=100) proid0, gid, zp1, theta_x1, xref_1, xg, zg, t
else
   read(10,rec=1,err=100) proid0, gid, zp1, theta_x1, xref_1, xg, zg, t
endif
npro = 1
do itrace=2,1000000000
   if (isAscii) then
      read(10,*,end=100) proid, gid, zp1, theta_x1, xref_1, xg, zg, t
   else
      read(10,rec=itrace,err=100) proid, gid, zp1, theta_x1, xref_1, xg, zg, t
   endif
   if (proid.ne.proid0) then
      proid0 = proid;   npro = npro + 1
   endif
enddo
100 continue
rewind(10)
c%npro = npro
call allocate_and_initial(c%ng,c%proid,npro)
! Determine the number of traces for each profile
if (isAscii) then
   read(10,*,end=200) proid0, gid, zp1, theta_x1, xref_1, xg, zg, t
else
   read(10,rec=1) proid0, gid, zp1, theta_x1, xref_1, xg, zg, t
endif
npro = 1;  ng = 1
do itrace=2,1000000000
   if (isAscii) then
      read(10,*,end=200) proid, gid, zp1, theta_x1, xref_1, xg, zg, t
   else
      read(10,rec=itrace) proid, gid, zp1, theta_x1, xref_1, xg, zg, t
   endif
   if (proid.ne.proid0) then
      proid0 = proid;    c%ng(npro) = ng
      npro = npro + 1;   ng = 1
   else
      ng = ng + 1
   endif
enddo
200 continue
c%ng(npro) = ng
rewind(10)
! Deternube the total number of traces for the whole dataset
c%ntrace=sum(c%ng(1:c%npro))
! Determine the maximum number of geophones per shot
c%ngmax = maxval(c%ng(1:c%npro))
call allocate_and_initial(c%gid,c%ngmax,npro)
call allocate_and_initial(c%zp1,c%theta_x1,c%xref_1,npro)
call allocate_and_initial(c%xg,c%zg,c%t,c%ngmax,npro)
! Read all information 
itrace=0
do ipro=1,c%npro
   do ig=1,c%ng(ipro)
      if (isAscii) then
         read(10,*,end=300)c%proid(ipro),c%gid(ig,ipro),c%zp1(ipro), &
            c%theta_x1(ipro), c%xref_1(ipro),c%xg(ig,ipro),c%zg(ig,ipro),c%t(ig,ipro)
      else
         itrace=itrace+1
         read(10,rec=itrace,err=300) c%proid(ipro), c%gid(ig,ipro), c%zp1(ipro), &
            c%theta_x1(ipro), c%xref_1(ipro),c%xg(ig,ipro),c%zg(ig,ipro),c%t(ig,ipro)
      endif
   enddo
enddo
300 continue
close(10)
end subroutine read_pw1_coordfile

subroutine read_pw1_coordfile_mpi(coordfile,coordfile_type,c)
character(len=*),   intent(in)   :: coordfile,coordfile_type
type(pw1_coord2d),  intent(out)  :: c
if (rank.eq.0)   call read_pw1_coordfile(coordfile,coordfile_type,c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro,c%ngmax,c%ntrace)
if (rank > 0) then
   call allocate_and_initial(c%proid,c%ng,c%npro)
   call allocate_and_initial(c%zp1,c%theta_x1,c%xref_1,c%npro)
   call allocate_and_initial(c%gid,c%ngmax,c%npro)
   call allocate_and_initial(c%xg,c%zg,c%t,c%ngmax,c%npro)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%ng,c%npro)
call sf90_bcast(c%zp1,c%theta_x1,c%xref_1,c%npro)
call sf90_bcast(c%gid,c%ngmax,c%npro)
call sf90_bcast(c%xg,c%zg,c%t,c%ngmax,c%npro)
end subroutine read_pw1_coordfile_mpi

subroutine write_pw1_coordfile(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(pw1_coord2d),   intent(in):: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="write")
else
   open(10,file=coordfile,access="direct",recl=8*I4,action="write")
endif
itrace=0
do ipro=1,c%npro
   do ig=1,c%ng(ipro)
      if (isAscii) then
         write(10,100)c%proid(ipro),c%gid(ig,ipro),c%zp1(ipro),c%theta_x1(ipro),&
            c%xref_1(ipro),c%xg(ig,ipro),c%zg(ig,ipro),c%t(ig,ipro)
      else
         itrace=itrace+1
         write(10,rec=itrace)c%proid(ipro),c%gid(ig,ipro),c%zp1(ipro),c%theta_x1(ipro),&
            c%xref_1(ipro),c%xg(ig,ipro),c%zg(ig,ipro),c%t(ig,ipro)
      endif
   enddo
enddo
100 format(2(I7),F8.2," ",F9.5,3(" ",F8.2)," ",F9.5)
close(10)
end subroutine write_pw1_coordfile

subroutine read_pw1_coordfile_single(coordfile,coordfile_type,c)
character(len=*),  intent(in)    :: coordfile,coordfile_type
type(pw1_coord2d), intent(inout) :: c
call allocate_and_initial(c%zp1,c%theta_x1,c%xref_1,1)
call allocate_and_initial(c%gid,c%ng(1),1)
call allocate_and_initial(c%xg,c%zg,c%t,c%ng(1),1)
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(11,file=coordfile,form="formatted",action="read")
else
   open(11,file=coordfile,access="direct",recl=8*I4,action="read")
endif
do ig=1,c%ng(1)
   if (isAscii) then
      read(11,*)c%proid(1),c%gid(ig,1),c%zp1(1),c%theta_x1(1),&
         c%xref_1(1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1)
   else
      read(11,rec=ig)c%proid(1),c%gid(ig,1),c%zp1(1),c%theta_x1(1),&
         c%xref_1(1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1)
   endif
enddo
close(11)
end subroutine read_pw1_coordfile_single

subroutine write_pw1_coordfile_single(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(pw1_coord2d),   intent(in):: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="write")
else
   open(10,file=coordfile,access="direct",recl=8*I4,action="write")
endif
do ig=1,c%ng(1)
   if (isAscii) then
      write(10,100)c%proid(1),c%gid(ig,1),c%zp1(1),c%theta_x1(1),&
         c%xref_1(1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1)
   else
      write(10,rec=ig)c%proid(1),c%gid(ig,1),c%zp1(1),c%theta_x1(1),&
         c%xref_1(1),c%xg(ig,1),c%zg(ig,1),c%t(ig,1)
   endif
enddo
100 format(2(I7),F8.2," ",F9.5,3(" ",F8.2)," ",F9.5)
close(10)
end subroutine write_pw1_coordfile_single

subroutine read_pw1_coordfile_all(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(pw1_coord2d),    intent(out)  :: c
call read_coordfile_all_core(coordfile,coordfile_type,c%npro,&
   c%proid,c%ngmax,c%ng,c%ntrace)
end subroutine read_pw1_coordfile_all

subroutine read_pw1_coordfile_all_mpi(coordfile, coordfile_type, c)
character(len=*),  intent(in)   :: coordfile,coordfile_type
type(pw1_coord2d), intent(out)  :: c
if (rank == 0)   call read_pw1_coordfile_all(coordfile, coordfile_type, c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro,c%ngmax,c%ntrace)
if (rank > 0)   call allocate_and_initial(c%proid,c%ng,c%npro)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%ng,c%npro)
end subroutine read_pw1_coordfile_all_mpi

subroutine write_pw1_coordfile_all(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(pw1_coord2d),   intent(in):: c
call write_coordfile_all_core(coordfile,coordfile_type,c%npro,c%proid,c%ng)
end subroutine write_pw1_coordfile_all

subroutine assign_pw1_coord(c_all,c,npro_me,ipro_me)
type(pw1_coord2d), intent(in)    :: c_all
type(pw1_coord2d), intent(inout) :: c
integer,           intent(in)    :: npro_me, ipro_me(:)
c%npro=npro_me
call allocate_and_initial(c%proid,c%ng,npro_me)
call allocate_and_initial(c%zp1,c%theta_x1,c%xref_1,npro_me)
do ipro=1,npro_me
   c%proid(ipro)    = c_all%proid(ipro_me(ipro))
   c%ng(ipro)       = c_all%ng(ipro_me(ipro))
   c%zp1(ipro)      = c_all%zp1(ipro_me(ipro))
   c%theta_x1(ipro) = c_all%theta_x1(ipro_me(ipro))
   c%xref_1(ipro)   = c_all%xref_1(ipro_me(ipro))
enddo
c%ngmax=maxval(c%ng);  c%ntrace=sum(c%ng)
call allocate_and_initial(c%gid,c%ngmax,c%npro)
call allocate_and_initial(c%xg,c%zg,c%t,c%ngmax,c%npro)
do ipro=1,c%npro
   c%gid(1:c%ng(ipro),ipro) = c_all%gid(1:c%ng(ipro),ipro_me(ipro))
   c%xg(1:c%ng(ipro),ipro)  = c_all%xg(1:c%ng(ipro),ipro_me(ipro))
   c%zg(1:c%ng(ipro),ipro)  = c_all%zg(1:c%ng(ipro),ipro_me(ipro))
   c%t(1:c%ng(ipro),ipro)   = c_all%t(1:c%ng(ipro),ipro_me(ipro))
enddo
end subroutine assign_pw1_coord

!--------------------------------------------------------------------------

subroutine read_pw2_coordfile(coordfile,coordfile_type,c)
character(len=*),   intent(in)   :: coordfile,coordfile_type
type(pw2_coord2d),  intent(out)  :: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="read")
else
   open(10,file=coordfile,access="direct",recl=9*I4,action="read")
endif
! Determine the number of profiles
npro = 0
if (isAscii) then
   read(10,*,end=100)proid0,gid,zp1,theta_x1,xref_1,zp2,theta_x2,xref_2,t
else
   read(10,rec=1,err=100)proid,gid,zp1,theta_x1,xref_1,zp2,theta_x2,xref_2,t
endif
npro = 1
do itrace=2,1000000000
   if (isAscii) then
      read(10,*,end=100)proid,gid,zp1,theta_x1,xref_1,zp2,theta_x2,xref_2,t
   else
      read(10,rec=itrace,err=100)proid,gid,zp1,theta_x1,xref_1,zp2,theta_x2,xref_2,t
   endif
   if (proid.ne.proid0) then
      proid0 = proid;   npro = npro + 1
   endif
enddo
100 continue
rewind(10)
c%npro = npro
call allocate_and_initial(c%ng,c%proid,npro)
! Determine the number of traces for each profile
if (isAscii) then
   read(10,*,end=200)proid0,gid,zp1,theta_x1,xref_1,zp2,theta_x2,xref_2,t
else
   read(10,rec=1,err=200)proid0,gid,zp1,theta_x1,xref_1,zp2,theta_x2,xref_2,t
endif
npro = 1;   ng = 1
do itrace=2,1000000000
   if (isAscii) then
      read(10,*,end=200)proid,gid,zp1,theta_x1,xref_1,zp2,theta_x2,xref_2,t
   else
      read(10,rec=itrace,err=200)proid,gid,zp1,theta_x1,xref_1,zp2,theta_x2,xref_2,t
   endif
   if (proid.ne.proid0) then
      proid0 = proid;   c%ng(npro) = ng;   npro = npro + 1;   ng = 1
   else
      ng = ng + 1
   endif
enddo
200 continue
c%ng(npro) = ng
rewind(10)
! Deternube the total number of traces for the whole dataset
c%ntrace=sum(c%ng(1:c%npro))
! Determine the maximum number of planes for all shots
c%ngmax = maxval(c%ng(1:c%npro))
call allocate_and_initial(c%gid,c%ngmax,npro)
call allocate_and_initial(c%zp1,c%theta_x1,c%xref_1,npro)
call allocate_and_initial(c%zp2,c%theta_x2,c%xref_2,c%t,c%ngmax,npro)
! Read all information 
itrace=0
do ipro=1,c%npro
   do ig=1,c%ng(ipro)
      if (isAscii) then
         read(10,*,end=300)c%proid(ipro),c%gid(ig,ipro),c%zp1(ipro),&
            c%theta_x1(ipro), c%xref_1(ipro),c%zp2(ig,ipro),c%theta_x2(ig,ipro),&
            c%xref_2(ig,ipro),c%t(ig,ipro)
      else
         itrace=itrace+1
         read(10,rec=itrace,err=300)c%proid(ipro),c%gid(ig,ipro),c%zp1(ipro), &
            c%theta_x1(ipro), c%xref_1(ipro),c%zp2(ig,ipro),c%theta_x2(ig,ipro),&
            c%xref_2(ig,ipro),c%t(ig,ipro)
      endif
   enddo
enddo
300 continue
close(10)
end subroutine read_pw2_coordfile

subroutine read_pw2_coordfile_mpi(coordfile,coordfile_type,c)
character(len=*),   intent(in)   :: coordfile,coordfile_type
type(pw2_coord2d),  intent(out)  :: c
if (rank.eq.0) call read_pw2_coordfile(coordfile,coordfile_type,c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro,c%ngmax,c%ntrace)
if (rank > 0) then
   call allocate_and_initial(c%proid,c%ng,c%npro)
   call allocate_and_initial(c%zp1,c%theta_x1,c%xref_1,c%npro)
   call allocate_and_initial(c%gid,c%ngmax,c%npro)
   call allocate_and_initial(c%zp2,c%theta_x2,c%xref_2,c%t,c%ngmax,c%npro)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%ng,c%npro)
call sf90_bcast(c%zp1,c%theta_x1,c%xref_1,c%npro)
call sf90_bcast(c%gid,c%ngmax,c%npro)
call sf90_bcast(c%zp2,c%theta_x2,c%xref_2,c%t,c%ngmax,c%npro)
end subroutine read_pw2_coordfile_mpi

subroutine write_pw2_coordfile(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(pw2_coord2d),   intent(in):: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form='formatted')
else
   open(10,file=coordfile,access="direct",recl=9*I4)
endif
itrace=0
do ipro=1,c%npro
   do ig=1,c%ng(ipro)
      if (isAscii) then
         write(10,100)c%proid(ipro),c%gid(ig,ipro),c%zp1(ipro),c%theta_x1(ipro),&
            c%xref_1(ipro),c%zp2(ig,ipro),c%theta_x2(ig,ipro),c%xref_2(ig,ipro),c%t(ig,ipro)
      else
         itrace=itrace+1
         write(10,rec=itrace)c%proid(ipro),c%gid(ig,ipro),c%zp1(ipro),c%theta_x1(ipro),&
            c%xref_1(ipro),c%zp2(ig,ipro),c%theta_x2(ig,ipro),c%xref_2(ig,ipro),c%t(ig,ipro)
      endif
   enddo
enddo
100 format(2(I7),F8.2," ",F9.5,2(" ",F8.2)," ",F8.5," ",F8.2," ",F8.5)
close(10)
end subroutine write_pw2_coordfile

subroutine read_pw2_coordfile_single(coordfile,coordfile_type,c)
character(len=*),  intent(in)    :: coordfile,coordfile_type
type(pw2_coord2d), intent(inout) :: c
call allocate_and_initial(c%zp1,c%theta_x1,c%xref_1,1)
call allocate_and_initial(c%gid,c%ng(1),1)
call allocate_and_initial(c%zp2,c%theta_x2,c%xref_2,c%t,c%ng(1),1)
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(11,file=coordfile,form="formatted",action="read")
else
   open(11,file=coordfile,access="direct",recl=9*I4,action="read")
endif
do ig=1,c%ng(1)
   if (isAscii) then
      read(11,*)c%proid(1),c%gid(ig,1),c%zp1(1),c%theta_x1(1),&
         c%xref_1(1),c%zp2(ig,1),c%theta_x2(ig,1),c%xref_2(ig,1),c%t(ig,1)
   else
      read(11,rec=ig)c%proid(1),c%gid(ig,1),c%zp1(1),c%theta_x1(1),&
         c%xref_1(1),c%zp2(ig,1),c%theta_x2(ig,1),c%xref_2(ig,1),c%t(ig,1)
   endif
enddo
close(11)
end subroutine read_pw2_coordfile_single

subroutine write_pw2_coordfile_single(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(pw2_coord2d),   intent(in):: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="write")
else
   open(10,file=coordfile,access="direct",recl=8*I4,action="write")
endif
do ig=1,c%ng(1)
   if (isAscii) then
      write(10,100)c%proid(1),c%gid(ig,1),c%zp1(1),c%theta_x1(1),&
         c%xref_1(1),c%zp2(ig,1),c%theta_x2(ig,1),c%xref_2(ig,1),c%t(ig,1)
   else
      write(10,rec=ig)c%proid(1),c%gid(ig,1),c%zp1(1),c%theta_x1(1),&
         c%xref_1(1),c%zp2(ig,1),c%theta_x2(ig,1),c%xref_2(ig,1),c%t(ig,1)
   endif
enddo
100 format(2(I7),F8.2," ",F9.5,2(" ",F8.2)," ",F8.5," ",F8.2," ",F8.5)
close(10)
end subroutine write_pw2_coordfile_single

subroutine read_pw2_coordfile_all(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(pw2_coord2d),    intent(out)  :: c
call read_coordfile_all_core(coordfile,coordfile_type,c%npro,&
   c%proid,c%ngmax,c%ng,c%ntrace)
end subroutine read_pw2_coordfile_all

subroutine read_pw2_coordfile_all_mpi(coordfile, coordfile_type, c)
character(len=*),  intent(in)   :: coordfile,coordfile_type
type(pw2_coord2d), intent(out)  :: c
if (rank == 0)   call read_pw2_coordfile_all(coordfile, coordfile_type, c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro,c%ngmax,c%ntrace)
if (rank > 0)   call allocate_and_initial(c%proid,c%ng,c%npro)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%ng,c%npro)
end subroutine read_pw2_coordfile_all_mpi

subroutine write_pw2_coordfile_all(coordfile,coordfile_type,c)
character(len=*),intent(in) :: coordfile,coordfile_type
type(pw2_coord2d),   intent(in):: c
call write_coordfile_all_core(coordfile,coordfile_type,c%npro,c%proid,c%ng)
end subroutine write_pw2_coordfile_all

subroutine assign_pw2_coord(c_all,c,npro_me,ipro_me)
type(pw2_coord2d), intent(in)    :: c_all
type(pw2_coord2d), intent(inout) :: c
integer,           intent(in)    :: npro_me, ipro_me(:)
c%npro=npro_me
call allocate_and_initial(c%proid,c%ng,npro_me)
call allocate_and_initial(c%zp1,c%theta_x1,c%xref_1,npro_me)
do ipro=1,npro_me
   c%proid(ipro)    = c_all%proid(ipro_me(ipro))
   c%ng(ipro)       = c_all%ng(ipro_me(ipro))
   c%zp1(ipro)      = c_all%zp1(ipro_me(ipro))
   c%theta_x1(ipro) = c_all%theta_x1(ipro_me(ipro))
   c%xref_1(ipro)   = c_all%xref_1(ipro_me(ipro))
enddo
c%ngmax=maxval(c%ng);  c%ntrace=sum(c%ng)
call allocate_and_initial(c%gid,c%ngmax,c%npro)
call allocate_and_initial(c%zp2,c%theta_x2,c%xref_2,c%t,c%ngmax,c%npro)
do ipro=1,c%npro
   c%gid(1:c%ng(ipro),ipro)      = c_all%gid(1:c%ng(ipro),ipro_me(ipro))
   c%zp2(1:c%ng(ipro),ipro)      = c_all%zp2(1:c%ng(ipro),ipro_me(ipro))
   c%theta_x2(1:c%ng(ipro),ipro) = c_all%theta_x2(1:c%ng(ipro),ipro_me(ipro))
   c%xref_2(1:c%ng(ipro),ipro)   = c_all%xref_2(1:c%ng(ipro),ipro_me(ipro))
   c%t(1:c%ng(ipro),ipro)        = c_all%t(1:c%ng(ipro),ipro_me(ipro))
enddo
end subroutine assign_pw2_coord

!--------------------------------------------------------------------------

subroutine read_eik_coordfile(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(eik_coord2d),    intent(out)  :: c
isAscii=get_isAscii(coordfile_type)
if (isAscii) then
   open(10,file=coordfile,form="formatted",action="read")
else
   open(10,file=coordfile,access="direct",recl=4*I4,action="read")
endif
! Determine the number of profiles
npro=0
do itrace=1,1000000000
   if (isAscii) then
      read(10,*,end=100) proid, sid, xs, zs
   else
      read(10,rec=itrace,err=100) proid, sid, xs, zs
   endif
   npro=npro+1
enddo
100 continue
rewind(10)
c%npro = npro
call allocate_and_initial(c%proid,c%sid,npro)
call allocate_and_initial(c%xs,c%zs,npro)
do itrace=1,c%npro
   if (isAscii) then
      read(10,*,end=200) c%proid(itrace),c%sid(itrace),c%xs(itrace),c%zs(itrace)
   else
      read(10,rec=itrace,err=200) c%proid(itrace),c%sid(itrace),c%xs(itrace),c%zs(itrace)
   endif
enddo
200 continue
close(10)
end subroutine read_eik_coordfile

subroutine read_eik_coordfile_mpi(coordfile, coordfile_type, c)
character(len=*), intent(in)   :: coordfile,coordfile_type
type(eik_coord2d),    intent(out)  :: c
if (rank == 0) call read_eik_coordfile(coordfile, coordfile_type, c)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%npro)
if (rank > 0) then
   call allocate_and_initial(c%proid,c%sid,c%npro)
   call allocate_and_initial(c%xs,c%zs,c%npro)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(c%proid,c%sid,c%npro)
call sf90_bcast(c%xs,c%zs,c%npro)
end subroutine read_eik_coordfile_mpi

!--------------------------------------------------------------------------------------

subroutine read_topofile(topofile, nx, npml, dx, z0, fs)
character(*), intent(in)  :: topofile
integer,      intent(in)  :: nx, npml
real,         intent(in)  :: dx, z0
integer,      intent(out) :: fs(:)
real                      :: topo(nx)
open(10,file=topofile,access="direct",status="old",recl=nx*I4,err=111,action="read")
read(10,rec=1) topo
close(10)
goto 222
111 continue
topo = z0
222 continue
! Set up free surface
fs = npml+1+int(topo/dx)
end subroutine read_topofile

subroutine read_topofile_mpi(topofile, nx, npml, dx, z0, fs)
character(*), intent(in)  :: topofile
integer,      intent(in)  :: nx, npml
real,         intent(in)  :: dx, z0
integer,      intent(out) :: fs(:)
if (rank == 0) then
   call read_topofile(topofile, nx, npml, dx, z0, fs)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(fs,nx)
end subroutine read_topofile_mpi

!--------------------------------------------------------------------------------------

logical function get_isAscii(coordfile_type)
character(len=*), intent(in) :: coordfile_type
if (coordfile_type.eq."ASCII") then
   get_isAscii=.true.
else
   get_isAscii=.false.
endif
end function get_isAscii

!--------------------------------------------------------------------------------------

subroutine read_pw1_par(pw1_par_file,pc)
character(len=*), intent(in)    :: pw1_par_file
type(pw1_coord2d),intent(inout) :: pc
open(10,file=pw1_par_file,form="formatted",action="read")
npro = 0
do ip=1,2000000
   read(10,*,end=100) proid, theta_x1, xref_1, zp1
   npro=npro+1
enddo
100 continue
rewind(10)
pc%npro=npro
call allocate_and_initial(pc%proid,npro)
call allocate_and_initial(pc%theta_x1,pc%xref_1,pc%zp1,npro)
do ip=1,pc%npro
   read(10,*,end=200) pc%proid(ip), pc%theta_x1(ip), pc%xref_1(ip), pc%zp1(ip)
enddo
200 continue
close(10)
end subroutine read_pw1_par

subroutine read_pw1_par_mpi(pw1_par_file,pc)
character(len=*), intent(in)    :: pw1_par_file
type(pw1_coord2d),intent(inout) :: pc
if (rank.eq.0)  call read_pw1_par(pw1_par_file,pc)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(pc%npro)
if (rank.ne.0) then
   call allocate_and_initial(pc%proid,pc%npro)
   call allocate_and_initial(pc%theta_x1,pc%xref_1,pc%zp1,pc%npro)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(pc%proid,pc%npro)
call sf90_bcast(pc%theta_x1,pc%xref_1,pc%zp1,pc%npro)
end subroutine read_pw1_par_mpi

subroutine read_pw2_par(pw2_par_file,pc)
character(len=*), intent(in)    :: pw2_par_file
type(pw2_coord2d),intent(inout) :: pc
open(10,file=pw2_par_file,form="formatted",action="read")
npro = 0
read(10,*,end=100) proid0, gid, theta_x1, xref_1, zp1, theta_x2, xref_2, zp2
npro = 1
do ip1=2,2000000
   read(10,*,end=100) proid, gid, theta_x1, xref_1, zp1, theta_x2, xref_2, zp2
   if (proid.ne.proid0) then
      npro=npro+1;   proid0=proid
   endif
enddo
100 continue
rewind(10)
pc%npro=npro
call allocate_and_initial(pc%proid,pc%ng,npro)
call allocate_and_initial(pc%theta_x1,pc%xref_1,pc%zp1,npro)
! Determine the number of pw2(ng) for each pw1
read(10,*,end=100) proid0, gid, theta_x1, xref_1, zp1, theta_x2, xref_2, zp2
npro = 1; ng = 1
do ip1=2,2000000
   read(10,*,end=200) proid, gid, theta_x1, xref_1, zp1, theta_x2, xref_2, zp2
   if (proid.ne.proid0) then
      proid0=proid
      pc%ng(npro) = ng
      npro = npro + 1
      ng = 1
   else
      ng = ng +1
   endif
enddo
200 continue
pc%ng(npro)=ng
rewind(10)
pc%ntrace=sum(pc%ng)
! Determine the maximum number of geophones per shot
pc%ngmax=maxval(pc%ng)
call allocate_and_initial(pc%theta_x2,pc%xref_2,pc%zp2,pc%ngmax,pc%npro)
do ip1=1,pc%npro
   do ig=1,pc%ng(ip1)
      read(10,*,end=300) pc%proid(ip1), pc%gid(ig,ip1), pc%theta_x1(ip1), &
         pc%xref_1(ip1), pc%zp1(ip1), pc%theta_x2(ig,ip1), pc%xref_2(ig,ip1), &
         pc%zp2(ig,ip1)
   enddo
enddo
300 continue
close(10)
end subroutine read_pw2_par

subroutine read_pw2_par_mpi(pw2_par_file,pc)
character(len=*), intent(in)    :: pw2_par_file
type(pw2_coord2d),intent(inout) :: pc
if (rank.eq.0)  call read_pw2_par(pw2_par_file,pc)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(pc%npro,pc%ngmax,pc%ntrace)
if (rank.ne.0) then
   call allocate_and_initial(pc%proid,pc%npro)
   call allocate_and_initial(pc%theta_x1,pc%xref_1,pc%zp1,pc%npro)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call sf90_bcast(pc%proid,pc%npro)
call sf90_bcast(pc%theta_x1,pc%xref_1,pc%zp1,pc%npro)
call sf90_bcast(pc%theta_x2,pc%xref_2,pc%zp2,pc%ngmax,pc%npro)
end subroutine read_pw2_par_mpi

subroutine coord_crg2pw1(c,pc)
type(coord2d),     intent(in)    :: c
type(pw1_coord2d), intent(inout) :: pc
pc%ngmax=c%npro
pc%ntrace=c%npro*pc%npro
call allocate_and_initial(pc%ng,pc%ngmax)
pc%ng=c%npro
call allocate_and_initial(pc%gid,pc%ngmax,pc%npro)
call allocate_and_initial(pc%xg,pc%zg,pc%t,pc%ngmax,pc%npro)
do ipro=1,pc%npro
   do ig=1,pc%ngmax
      pc%gid(ig,ipro)=c%gid(1,ig);pc%xg(ig,ipro)=c%xg(1,ig);pc%zg(ig,ipro)=c%zg(1,ig)
   enddo
enddo
end subroutine coord_crg2pw1

subroutine coord_csg2pw1(c,pc)
type(coord2d),     intent(in)    :: c
type(pw1_coord2d), intent(inout) :: pc
pc%ngmax=c%npro
pc%ntrace=c%npro*pc%npro
call allocate_and_initial(pc%ng,pc%ngmax)
pc%ng=c%npro
call allocate_and_initial(pc%gid,pc%ngmax,pc%npro)
call allocate_and_initial(pc%xg,pc%zg,pc%t,pc%ngmax,pc%npro)
do ipro=1,pc%npro 
   do ig=1,pc%ngmax ! == is=1,c%npro
      pc%gid(ig,ipro)=c%sid(1,ig);pc%xg(ig,ipro)=c%xs(1,ig);pc%zg(ig,ipro)=c%zs(1,ig)
   enddo
enddo
end subroutine coord_csg2pw1

end module module_acquisition2d
