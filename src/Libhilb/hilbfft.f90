subroutine hilbfft(a,nx,flag)
  implicit none
  integer :: nx,nf,flag
  real    :: a(1:nx)
  complex :: ca(1:nx)
  nf=(nx-1)/2+2
  ca(1:nx)=cmplx(a(1:nx),0.)
  call fft1d(ca,nx,-1)

  ca(nf:nx)=0.
  ca=2*ca
  call fft1d(ca,nx,1)
  if(flag)then
  a=abs(cmplx(a,imag(ca(1:nx))))
  else
  a(1:nx)=imag(ca(1:nx))
  endif
  return
end subroutine hilbfft
