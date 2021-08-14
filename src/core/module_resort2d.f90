module module_resort2d
use module_datatype
use module_array
use module_sf90_mpi
implicit none
real,    allocatable, private :: xs(:),zs(:),xg(:),zg(:),t(:)
integer, allocatable, private :: sid(:),gid(:),num_chan(:)

private resort_core, convert_1d_to_2d, convert_2d_to_1d

contains
!=======================================================================
subroutine csg2crg_resort(c1,c2,indx,sort_tol)
type(coord2d), intent(in)  :: c1
type(coord2d), intent(inout) :: c2
integer,       intent(out) :: indx(:)
real,          intent(in)  :: sort_tol
c2%ntrace=c1%ntrace
call allocate_and_initial(xs,zs,xg,zg,t,c1%ntrace)
call allocate_and_initial(sid,gid,num_chan,c1%ntrace)
call convert_2d_to_1d(c1)
call resort_core(c1%ntrace,xg,xs,c2%npro,num_chan,indx,sort_tol)
c2%ngmax=maxval(num_chan(1:c2%npro))
call allocate_and_initial(c2%proid,c2%ng,c2%npro)
call allocate_and_initial(c2%sid,c2%gid,c2%ngmax,c2%npro)
call allocate_and_initial(c2%xs,c2%zs,c2%xg,c2%zg,c2%t,c2%ngmax,c2%npro)
call convert_1d_to_2d(c2,indx)
call deallocate_and_free(xs,zs,xg,zg,t)   
call deallocate_and_free(sid,gid,num_chan)   
end subroutine csg2crg_resort
!=======================================================================

!=======================================================================
subroutine csg2cmp_resort(c1,c2,indx,sort_tol)
type(coord2d), intent(in)  :: c1
type(coord2d), intent(inout) :: c2
integer,       intent(out) :: indx(:)
real,          intent(in)  :: sort_tol
real,allocatable           :: mp(:)
c2%ntrace=c1%ntrace
call allocate_and_initial(xs,zs,xg,zg,t,mp,c1%ntrace)
call allocate_and_initial(sid,gid,num_chan,c1%ntrace)
call convert_2d_to_1d(c1)
mp=(xs+xg)/2.0
call resort_core(c1%ntrace,mp,xs,c2%npro,num_chan,indx,sort_tol)
c2%ngmax=maxval(num_chan(1:c2%npro))
call allocate_and_initial(c2%proid,c2%ng,c2%npro)
call allocate_and_initial(c2%sid,c2%gid,c2%ngmax,c2%npro)
call allocate_and_initial(c2%xs,c2%zs,c2%xg,c2%zg,c2%t,c2%ngmax,c2%npro)
call convert_1d_to_2d(c2,indx)
call deallocate_and_free(xs,zs,xg,zg,t,mp)   
call deallocate_and_free(sid,gid,num_chan)   
end subroutine csg2cmp_resort
!=======================================================================

!=======================================================================
subroutine cmp2csg_resort(c1,c2,indx,sort_tol)
type(coord2d), intent(in)  :: c1
type(coord2d), intent(inout) :: c2
integer,       intent(out) :: indx(:)
real,          intent(in)  :: sort_tol
real,allocatable           :: mp(:)
c2%ntrace=c1%ntrace
call allocate_and_initial(xs,zs,xg,zg,t,mp,c1%ntrace)
call allocate_and_initial(sid,gid,num_chan,c1%ntrace)
call convert_2d_to_1d(c1)
mp=(xs+xg)/2.0
call resort_core(c1%ntrace,xs,mp,c2%npro,num_chan,indx,sort_tol)
c2%ngmax=maxval(num_chan(1:c2%npro))
call allocate_and_initial(c2%proid,c2%ng,c2%npro)
call allocate_and_initial(c2%sid,c2%gid,c2%ngmax,c2%npro)
call allocate_and_initial(c2%xs,c2%zs,c2%xg,c2%zg,c2%t,c2%ngmax,c2%npro)
call convert_1d_to_2d(c2,indx)
call deallocate_and_free(xs,zs,xg,zg,t,mp)   
call deallocate_and_free(sid,gid,num_chan)   
end subroutine cmp2csg_resort
!=======================================================================

!=======================================================================
subroutine csg2cog_resort(c1,c2,indx,sort_tol)
type(coord2d), intent(in)  :: c1
type(coord2d), intent(inout) :: c2
integer,       intent(out) :: indx(:)
real,          intent(in)  :: sort_tol
real,allocatable           :: offset(:)
c2%ntrace=c1%ntrace
call allocate_and_initial(xs,zs,xg,zg,t,offset,c1%ntrace)
call allocate_and_initial(sid,gid,num_chan,c1%ntrace)
call convert_2d_to_1d(c1)
offset=(xs-xg)/2.0
call resort_core(c1%ntrace,offset,xs,c2%npro,num_chan,indx,sort_tol)
c2%ngmax=maxval(num_chan(1:c2%npro))
call allocate_and_initial(c2%proid,c2%ng,c2%npro)
call allocate_and_initial(c2%sid,c2%gid,c2%ngmax,c2%npro)
call allocate_and_initial(c2%xs,c2%zs,c2%xg,c2%zg,c2%t,c2%ngmax,c2%npro)
call convert_1d_to_2d(c2,indx)
call deallocate_and_free(xs,zs,xg,zg,t,offset)   
call deallocate_and_free(sid,gid,num_chan)   
end subroutine csg2cog_resort
!=======================================================================

!=======================================================================
subroutine cog2csg_resort(c1,c2,indx,sort_tol)
type(coord2d), intent(in)  :: c1
type(coord2d), intent(inout) :: c2
integer,       intent(out) :: indx(:)
real,          intent(in)  :: sort_tol
real,allocatable           :: offset(:)
c2%ntrace=c1%ntrace
call allocate_and_initial(xs,zs,xg,zg,t,offset,c1%ntrace)
call allocate_and_initial(sid,gid,num_chan,c1%ntrace)
call convert_2d_to_1d(c1)
offset=(xs-xg)/2.0
call resort_core(c1%ntrace,xs,offset,c2%npro,num_chan,indx,sort_tol)
c2%ngmax=maxval(num_chan(1:c2%npro))
call allocate_and_initial(c2%proid,c2%ng,c2%npro)
call allocate_and_initial(c2%sid,c2%gid,c2%ngmax,c2%npro)
call allocate_and_initial(c2%xs,c2%zs,c2%xg,c2%zg,c2%t,c2%ngmax,c2%npro)
call convert_1d_to_2d(c2,indx)
call deallocate_and_free(xs,zs,xg,zg,t,offset)   
call deallocate_and_free(sid,gid,num_chan)   
end subroutine cog2csg_resort
!=======================================================================

!=======================================================================
subroutine cog2cmp_resort(c1,c2,indx,sort_tol)
type(coord2d), intent(in)  :: c1
type(coord2d), intent(inout) :: c2
integer,       intent(out) :: indx(:)
real,          intent(in)  :: sort_tol
real,allocatable           :: offset(:),mp(:)
c2%ntrace=c1%ntrace
call allocate_and_initial(xs,zs,xg,zg,t,mp,offset,c1%ntrace)
call allocate_and_initial(sid,gid,num_chan,c1%ntrace)
call convert_2d_to_1d(c1)
offset=(xs-xg)/2.0
mp=(xs+xg)/2.0
call resort_core(c1%ntrace,mp,offset,c2%npro,num_chan,indx,sort_tol)
c2%ngmax=maxval(num_chan(1:c2%npro))
call allocate_and_initial(c2%proid,c2%ng,c2%npro)
call allocate_and_initial(c2%sid,c2%gid,c2%ngmax,c2%npro)
call allocate_and_initial(c2%xs,c2%zs,c2%xg,c2%zg,c2%t,c2%ngmax,c2%npro)
call convert_1d_to_2d(c2,indx)
call deallocate_and_free(xs,zs,xg,zg,t,mp,offset)   
call deallocate_and_free(sid,gid,num_chan)   
end subroutine cog2cmp_resort
!=======================================================================

!=======================================================================
subroutine cmp2cog_resort(c1,c2,indx,sort_tol)
type(coord2d), intent(in)  :: c1
type(coord2d), intent(inout) :: c2
integer,       intent(out) :: indx(:)
real,          intent(in)  :: sort_tol
real,allocatable           :: offset(:),mp(:)
c2%ntrace=c1%ntrace
call allocate_and_initial(xs,zs,xg,zg,t,mp,offset,c1%ntrace)
call allocate_and_initial(sid,gid,num_chan,c1%ntrace)
call convert_2d_to_1d(c1)
offset=(xs-xg)/2.0
mp=(xs+xg)/2.0
call resort_core(c1%ntrace,offset,mp,c2%npro,num_chan,indx,sort_tol)
c2%ngmax=maxval(num_chan(1:c2%npro))
call allocate_and_initial(c2%proid,c2%ng,c2%npro)
call allocate_and_initial(c2%sid,c2%gid,c2%ngmax,c2%npro)
call allocate_and_initial(c2%xs,c2%zs,c2%xg,c2%zg,c2%t,c2%ngmax,c2%npro)
call convert_1d_to_2d(c2,indx)
call deallocate_and_free(xs,zs,xg,zg,t,mp,offset)   
call deallocate_and_free(sid,gid,num_chan)   
end subroutine cmp2cog_resort

!=======================================================================
subroutine resort_core(ntrace,rec1,rec2,npro,num_chan,indx_out,sort_tol)
integer,   intent(in)  :: ntrace
real,      intent(in)  :: rec1(:), rec2(:), sort_tol
integer,   intent(out) :: npro, num_chan(:), indx_out(:)
integer,   allocatable :: indx(:),indx1(:),indx2(:),indx3(:)
real,      allocatable :: rec1_new(:),rec2_new(:)
integer                :: itrace, itrace1, itrace2, ii
call allocate_and_initial(indx,indx1,indx2,indx3,ntrace)
call allocate_and_initial(rec1_new,rec2_new,ntrace)
call indexx(ntrace,rec1,indx)
do itrace=1,ntrace
   rec1_new(itrace)=rec1(indx(itrace))
   rec2_new(itrace)=rec2(indx(itrace))   
enddo
npro=1
num_chan(npro)=1
itrace1=1
do itrace=2,ntrace
   if (abs(rec1_new(itrace)-rec1_new(itrace-1))>sort_tol) then
      itrace2=itrace-1
      indx1(1:num_chan(npro))=indx(itrace1:itrace2)
      call indexx(num_chan(npro),rec2_new(itrace1:itrace2),indx2(1:num_chan(npro)))
      do ii=1,num_chan(npro)
         indx3(ii)=indx1(indx2(ii))
      enddo
      do ii=itrace1,itrace2
         indx_out(ii)=indx3(ii-itrace1+1)
      enddo
      npro=npro+1
      num_chan(npro)=1
      itrace1=itrace
   else
      num_chan(npro)=num_chan(npro)+1
   endif
enddo  
if (num_chan(npro)>=1) then
   itrace2=ntrace
   indx1(1:num_chan(npro))=indx(itrace1:itrace2)
   call indexx(num_chan(npro),rec2_new(itrace1:itrace2),indx2(1:num_chan(npro)))
   do ii=1,num_chan(npro)
      indx3(ii)=indx1(indx2(ii))
   enddo
   do ii=itrace1,itrace2
      indx_out(ii)=indx3(ii-itrace1+1)
   enddo
endif
call deallocate_and_free(indx,indx1,indx2,indx3)
call deallocate_and_free(rec1_new,rec2_new)
end subroutine resort_core

subroutine convert_2d_to_1d(c1)
type(coord2d), intent(in)  :: c1
integer :: itrace, ipro, ig
itrace=0
do ipro=1,c1%npro
   do ig=1,c1%ng(ipro)
      itrace=itrace+1
      xs(itrace)=c1%xs(ig,ipro)
      zs(itrace)=c1%zs(ig,ipro)
      xg(itrace)=c1%xg(ig,ipro)
      zg(itrace)=c1%zg(ig,ipro)
      sid(itrace)=c1%sid(ig,ipro)
      gid(itrace)=c1%gid(ig,ipro)
      t(itrace)=c1%t(ig,ipro)
   enddo
enddo
end subroutine convert_2d_to_1d

subroutine convert_1d_to_2d(c2,indx)
integer,       intent(in)    :: indx(:)
type(coord2d), intent(inout) :: c2
integer                      :: itrace,ipro,ig
itrace=0
do ipro=1,c2%npro
   c2%proid(ipro)=ipro
   c2%ng(ipro)=num_chan(ipro)
   do ig=1,c2%ng(ipro)
      itrace=itrace+1
      c2%xs(ig,ipro)=xs(indx(itrace))
      c2%zs(ig,ipro)=zs(indx(itrace))
      c2%xg(ig,ipro)=xg(indx(itrace))
      c2%zg(ig,ipro)=zg(indx(itrace))
      c2%sid(ig,ipro)=sid(indx(itrace))
      c2%gid(ig,ipro)=gid(indx(itrace))
      c2%t(ig,ipro)=t(indx(itrace))
   enddo
enddo
end subroutine convert_1d_to_2d

end module module_resort2d
