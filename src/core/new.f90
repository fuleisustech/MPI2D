subroutine modeling(fdm,bc,wf,s,output)
type(fdmod2d),      intent(inout) :: fdm
type(shot2d),       intent(inout) :: s
type(bc2d_general), intent(inout) :: bc
type(wf2d),         intent(inout) :: wf
character(len=*),   intent(in)    :: output
integer                           :: fid,it_out,mod_type
call initial_variables(fdm,s)
if (allocated(bc%p_top).or.allocated(wf%pwf)) then
   if (allocated(bc%p_top)) then
      mod_type=2
   else
      mod_type=3
   endif
else
   mod_type=1
endif
if (isWriteData) then
   call MPI_FILE_OPEN(MPI_COMM_SELF, output, &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, &
         MPI_INFO_NULL, fid, ierr)
else
   s%seis = 0.0
endif
call setup_source_receiver(s)
call allocate_and_initial(p,p0,p1,nz_pml,nx_pml)
call wf_assign(p,pp,p0,pp0,p1,pp1)
if (bc_type==2) call allocate_and_initial(u1,nz_pml,nx_pml)
do it=1,nt
   if (bc_type==1) then
      call a2d_abc_kernal(pp,pp0,pp1,temp1,temp2,alpha,nz_pml,nx_pml,fd_order)
   else
      call a2d_sponge_kernal(pp,pp0,pp1,u1,c1t2,beta_dtdx,spgx,spgz,spcx1,spcx2,&
              spcz1,spcz2,npml,nz_pml,nx_pml,fd_order)
                                                                                        215,1         30%
endif
   call add_source(pp,beta_dt,fdm%s%sou(it),isFS,isDipole)
   if (isFS) then
      call add_FS(pp,npml,fs_thick)
   endif
   if (mod_type.eq.2) then
      call saveBC2d(pp,nz,nx,nzp,nxp,npml,bc%p_top(:,:,it),bc%p_bottom(:,:,it),&
               bc%p_left(:,:,it),bc%p_right(:,:,it),bc_len)
   elseif (mod_type.eq.3) then
      wf%pwf(:,:,it)=pp(npml:nzp+1,npml:nxp+1)
   endif
   !if (mod(it,1000)==1) call snapshot("p_",p,it,nz_pml,nx_pml)
   if (isWriteData) then
      if (dnt_mod.ne.1) then
         if (mod(it,dnt_mod).eq.1) then
            it_out=(it-1)/dnt_mod+1
g=1,s%ng
               disp=((ig-1)*s%nt+(it_out-1))*4
               call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
               call MPI_FILE_WRITE(fid, pp(igz(ig),igx(ig)),1, MPI_REAL, &
                         MPI_STATUS_IGNORE, ierr)
            enddo
         endif
      else
         do ig=1,s%ng
            disp=((ig-1)*s%nt+(it-1))*4
            call MPI_FILE_SEEK(fid,disp,MPI_SEEK_SET,ierr)
            call MPI_FILE_WRITE(fid, pp(igz(ig),igx(ig)),1, MPI_REAL, &
                      MPI_STATUS_IGNORE, ierr)
         enddo
      endif
   else
      call store_data(pp,s%ng,it,s%seis)
   endif
call wf_refresh(pp0,pp1,pp)
enddo  ! End of time-marching loop
if (mod_type.eq.2) then
   call copy_array(pp0,bc%p_nt_1,p1,bc%p_nt)
endif
if (isWriteData) then
   call MPI_FILE_CLOSE(fid,ierr)
endif
call finalize_variables()
call deallocate_and_free(p,p0,p1)
if (bc_type==2) call deallocate_and_free(u1)
end subroutine modeling

