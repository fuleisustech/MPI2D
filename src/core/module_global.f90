!
!  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
!  King Abdullah University of Science and Technology, All rights reserved.
!
!  Sponsors of CSIM are granted a non-exclusive, irrevocable royalty free
!  world-wide license to use this software and associated documentation files
!  (the "Software"), in their business, including the rights to modify and
!  distribute the Software to their affiliates, partners, clients or consultants
!  as necessary in the conduct of the sponsors business, all without accounting
!  to the King Abdullah University of Science and Technology, subject to the
!  following conditions:
!
!  The above copyright notice and this permission notice shall be included
!  in all copies or substantial portions of the Software.
!
!  Warranty Disclaimer:
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!  DEALINGS IN THE SOFTWARE.
!

!
!  file:     global.f90
!  @version: 1.0
!  @author:  Chaiwoot Boonyasiriwat
!  email:    chaiwoot@kaust.edu.sa
!  date:     December 2009
!  purpose:  module that defines global variables or constants
!

module module_global

implicit none

! Constants
real, parameter :: PI        = 3.14159265358979
real, parameter :: TWOPI     = 2.0*PI

! Timing variables
character(len=8)  :: date
character(len=10) :: time_now
 
! Memory variables
integer(kind=8), parameter :: MEM_UNIT=4*1024*1024

! Record length unit
integer, parameter :: I4 = 4

! Imaginary Unit Number         ! Bowen Guo
complex, parameter :: J = (0,1) ! Bowen Guo

! Length of string
integer, parameter :: slen=400

! Type of job
character(len=slen) :: jobtype, parfile

! FD coefficients
real, parameter :: C21 = -2
real, parameter :: C22 = 1.0
real, parameter :: C41 = -2.5
real, parameter :: C42 = 4.0/3.0
real, parameter :: C43 = -1.0/12.0
real, parameter :: C61 = -49.0/18.0
real, parameter :: C62 = 3.0/2.0
real, parameter :: C63 = -3.0/20.0
real, parameter :: C64 = 1.0/90.0
real, parameter :: C81 = -205.0/72.0
real, parameter :: C82 = 8.0/5.0
real, parameter :: C83 = -1.0/5.0
real, parameter :: C84 = 8.0/315.0
real, parameter :: C85 = -1.0/560.0

real, parameter :: S21 = 1.0
real, parameter :: S41 = 9.0/8.0
real, parameter :: S42 = -1.0/24.0
real, parameter :: S61 = 1.17187
real, parameter :: S62 = -6.51042E-2
real, parameter :: S63 = 4.68750E-3
real, parameter :: S81 = 1.19629
real, parameter :: S82 = -7.97526E-2
real, parameter :: S83 = 9.57031E-3
real, parameter :: S84 = -6.97545E-4

end module module_global

