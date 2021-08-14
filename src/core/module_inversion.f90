module module_inversion
implicit none
double precision, private :: a,b,beta

interface gradient
   module procedure gradient_2d
   module procedure gradient_LSMF_2d  ! Added by Bowen Guo
   module procedure gradient_3d
end interface gradient

interface step_length
   module procedure step_length_2d
   module procedure step_length_LSMF_2d !  Added by Bowen Guo
   module procedure step_length_3d
end interface step_length

interface trial_steplength
   module procedure trial_steplength2d
   module procedure trial_steplength3d
end interface trial_steplength

interface numerical_steplength
   module procedure numerical_steplength2d
   module procedure numerical_steplength3d
end interface numerical_steplength

contains

subroutine gradient_2d(inv_type,cg_type,isPreCondition,it,gk,gk1,prec_gk1,dk,dk1)
character(len=*), intent(in)  :: inv_type,cg_type
logical,          intent(in)  :: isPreCondition
integer,          intent(in)  :: it
real,             intent(in)  :: gk(:,:),gk1(:,:),prec_gk1(:,:),dk(:,:)
real,             intent(out) :: dk1(:,:)
if (inv_type.eq."CG") then
   if (it.eq.1) then
      if (isPreCondition) then
         dk1=prec_gk1
      else
         dk1=gk1
      endif
   else
      if (cg_type.eq."Fletcher-Reeves") then
         if (isPreCondition) then
            a=sum(gk1*prec_gk1)
         else
            a=sum(gk1**2.0)
         endif
      elseif (cg_type.eq."Polak-Rebiere") then
         if (isPreCondition) then
            a=sum((gk1-gk)*prec_gk1)
         else
            a=sum((gk1-gk)*gk1)
         endif
      endif
      b=sum(gk*gk)
      beta=a/b
      if (isPreCondition) then
         dk1=prec_gk1*beta*dk
         !dk1=prec_gk1+beta*dk
      else
         dk1=gk1+beta*dk
      endif
   endif
elseif (inv_type.eq."SD") then
   if (isPreCondition) then
      dk1=prec_gk1
   else
      dk1=gk1
   endif
endif
end subroutine gradient_2d
!========================================================================================
! Added by Bowen Guo
subroutine gradient_LSMF_2d(inv_type,cg_type,isPreCondition,it,gk_p,gk_s,gk1_p,gk1_s,&
                            prec_gk1_p,prec_gk1_s,dk_p,dk_s,dk1_p,dk1_s)
character(len=*), intent(in)  :: inv_type,cg_type
logical,          intent(in)  :: isPreCondition
integer,          intent(in)  :: it
real,             intent(in)  :: gk_p(:,:),gk_s(:,:),gk1_p(:,:),gk1_s(:,:),prec_gk1_p(:,:),prec_gk1_s(:,:),&
                                 dk_p(:,:),dk_s(:,:)
real,             intent(out) :: dk1_p(:,:),dk1_s(:,:)

if (inv_type.eq."CG") then 
   if (it.eq.1) then 
      if (isPreCondition) then 
         dk1_p=prec_gk1_p
         dk1_s=prec_gk1_s
      else
         dk1_p=gk1_p
         dk1_s=gk1_s
      endif
   else
      if (cg_type.eq."Fletcher-Reeves") then 
         if (isPreCondition) then
            a=sum(gk1_p*prec_gk1_p)+sum(gk1_s*prec_gk1_s)
         else
            a=sum(gk1_p**2.0)+sum(gk1_s**2.0)
         endif
      elseif (cg_type.eq."Polak-Rebiere") then 
         if (isPreCondition) then
            a=sum((gk1_p-gk_p)*prec_gk1_p)+sum((gk1_s-gk_s)*prec_gk1_s)
         else
            a=sum((gk1_p-gk_p)*gk1_p)+sum((gk1_s-gk_s)*gk1_s)
         endif
     endif
     b=sum(gk_p**2.0)+sum(gk_s**2.0)
     beta=a/b
     if (isPreCondition) then          !!! The preconditioning Matrix might be wrong
        dk1_p=prec_gk1_p*beta*dk_p
        dk1_s=prec_gk1_s*beta*dk_s
     else
        dk1_p=gk1_p+beta*dk_p
        dk1_s=gk1_s+beta*dk_s
     endif
   endif
elseif (inv_type.eq."SD") then
   if (isPreCondition) then
         dk1_p=prec_gk1_p
         dk1_s=prec_gk1_s
   else
         dk1_p=gk1_p
         dk1_s=gk1_s
   endif
endif
end subroutine gradient_LSMF_2d
!=================================================================================================

subroutine gradient_3d(inv_type,cg_type,isPreCondition,it,gk,gk1,prec_gk1,dk,dk1)
character(len=*), intent(in)  :: inv_type,cg_type
logical,          intent(in)  :: isPreCondition
integer,          intent(in)  :: it
real,             intent(in)  :: gk(:,:,:),gk1(:,:,:),prec_gk1(:,:,:),dk(:,:,:)
real,             intent(out) :: dk1(:,:,:)
if (inv_type.eq."CG") then
   if (it.eq.1) then
      if (isPreCondition) then
         dk1=prec_gk1
      else
         dk1=gk1
      endif
   else
      if (cg_type.eq."Fletcher-Reeves") then
         if (isPreCondition) then
            a=sum(gk1*prec_gk1)
         else
            a=sum(gk1**2.0)
         endif
      elseif (cg_type.eq."Polak-Rebiere") then
         if (isPreCondition) then
            a=sum((gk1-gk)*prec_gk1)
         else
            a=sum((gk1-gk)*gk1)
         endif
      endif
      b=sum(gk*gk)
      beta=a/b
      if (isPreCondition) then
         dk1=prec_gk1*beta*dk
      else
         dk1=gk1+beta*dk
      endif
   endif
elseif (inv_type.eq."SD") then
   if (isPreCondition) then
      dk1=prec_gk1
   else
      dk1=gk1
   endif
endif
end subroutine gradient_3d

!subroutine step_length_2d(inv_type,cg_type,isPreCondition,gk1,prec_gk1,dk1,gTLTLg,alpha)
subroutine step_length_2d(inv_type,isPreCondition,gk1,prec_gk1,dk1,gTLTLg,alpha)
!character(len=*), intent(in)  :: inv_type,cg_type
character(len=*), intent(in)  :: inv_type
logical,          intent(in)  :: isPreCondition
real,   intent(in) :: gk1(:,:),prec_gk1(:,:),dk1(:,:)
double precision, intent(in)  :: gTLTLg
double precision, intent(out):: alpha
if (inv_type.eq."CG") then
   a=sum(gk1*dk1)
  ! write(*,*) "The step length calculation has been initiated";call flush(6)
else
   if (isPreCondition) then
      a=sum(gk1*prec_gk1)
   else
      a=sum(gk1*gk1)
   endif
endif
!b=sum(Lg*Lg)
alpha=a/gTLTLg
!write(*,*)"alpha",alpha
end subroutine step_length_2d


!===========================================================================================
! Added by Bowen Guo
subroutine step_length_LSMF_2d(inv_type,isPreCondition,gk1_p,gk1_s,prec_gk1_p,prec_gk1_s,dk1_p,dk1_s,gTLTLg,alpha)
character(len=*), intent(in)  :: inv_type
logical,          intent(in)  :: isPreCondition
real,   intent(in)            :: gk1_p(:,:),gk1_s(:,:),prec_gk1_p(:,:),prec_gk1_s(:,:),&
                                 dk1_p(:,:),dk1_s(:,:)
double precision, intent(in)  :: gTLTLg
double precision, intent(out) :: alpha
if (inv_type.eq."CG") then 
   a=sum(gk1_p*dk1_p)+sum(gk1_s*dk1_s)
  ! write(*,*) "The step length calculation has been initiialed"
else
   if (isPreCondition) then
      a=sum(gk1_p*prec_gk1_p)+sum(gk1_s*prec_gk1_s)
   else
      a=sum(gk1_p*gk1_p)+sum(gk1_s*gk1_s)
   endif
endif

alpha=a/gTLTLg

end subroutine step_length_LSMF_2d




!===========================================================================================





!subroutine step_length_3d(inv_type,cg_type,isPreCondition,gk1,prec_gk1,dk1,gTLTLg,alpha)
subroutine step_length_3d(inv_type,isPreCondition,gk1,prec_gk1,dk1,gTLTLg,alpha)
!character(len=*), intent(in)  :: inv_type,cg_type
character(len=*), intent(in)  :: inv_type
logical,          intent(in)  :: isPreCondition
real,   intent(in) :: gk1(:,:,:),prec_gk1(:,:,:),dk1(:,:,:)
double precision, intent(in)  :: gTLTLg
double precision,   intent(out):: alpha
if (inv_type.eq."CG") then
   a=sum(gk1*dk1)
else
   if (isPreCondition) then
      a=sum(gk1*prec_gk1)
   else
      a=sum(gk1*gk1)
   endif
endif
!b=sum(Lg*Lg)
alpha=a/gTLTLg
end subroutine step_length_3d

!---------------------------------------------------------------------------------------
subroutine trial_steplength2d(m,dk1,stl_coe,alpha)
implicit none
real, intent(in)     :: m(:,:),dk1(:,:),stl_coe
double precision, intent(out)    :: alpha
real                 :: m_max,dk1_max
m_max=maxval(abs(m))
dk1_max=maxval(abs(dk1))
alpha=m_max*stl_coe/dk1_max
end subroutine trial_steplength2d

!---------------------------------------------------------------------------------------
subroutine trial_steplength3d(m,dk1,stl_coe,alpha)
implicit none
real, intent(in)     :: m(:,:,:),dk1(:,:,:),stl_coe
double precision, intent(out)    :: alpha
real                 :: m_max,dk1_max
m_max=maxval(abs(m))
dk1_max=maxval(abs(dk1))
alpha=m_max*stl_coe/dk1_max
!write(*,*)"m_max=",m_max," dk1_max=",dk1_max," alpha=",alpha
!call flush(6)
end subroutine trial_steplength3d

subroutine numerical_steplength2d(alpha,res0,res1,dk1,gama)
implicit none
double precision, intent(in)  :: alpha,res0,res1
real, intent(in)              :: dk1(:,:)
double precision, intent(out) :: gama
double precision              :: b
b=sum(dk1*dk1)
gama=-b/((res1-res0-b*alpha)/(alpha*alpha)*2.0)
end subroutine numerical_steplength2d

subroutine numerical_steplength3d(alpha,res0,res1,dk1,gama)
implicit none
double precision, intent(in)  :: alpha,res0,res1
real, intent(in)              :: dk1(:,:,:)
double precision, intent(out) :: gama
double precision              :: b
b=sum(dk1*dk1)
gama=-b/((res1-res0-b*alpha)/(alpha*alpha)*2.0)
end subroutine numerical_steplength3d

end module module_inversion
