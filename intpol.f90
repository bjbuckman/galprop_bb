!=============================================================================!
!=============================================================================!
!=============================================================================!
module spectra
  implicit none
  save
  integer, parameter :: n1=26,n2=100,n3=19,n4=22,n5=21,n6=4
  integer :: sym
  double precision, parameter :: m_p = 938.272d6        ! proton mass/eV 
  double precision, parameter :: E0=1.d10,d=1.5848,E_min=1d5,dn=0.1d0
  double precision xa(n1),ya(n2),fqgs(n1,n2,n6)
  double precision d_f,x0
  double precision, parameter :: Etrans=40.d9,width=0.01d0
end module spectra
!=============================================================================!
!=============================================================================!
!program photonspectra 
! changed to subroutine to avoid main clash AWS20130503
  subroutine  photonspectra
  use spectra, only : sym,m_p 
  implicit none
  integer i
  double precision E_p,E_g,T,s1,s2,s3,s4

  call init

! input: total energy E_p/eV of projectile nucleus, 
!        energy E_g/eV of photon
!        reac: 1 pp, 2 pHe, 3 Hep, 4 HeHe
!        sym = 1, Kamae symmetrised, sym=0 standard
! output: dsigma(E_p,E_g)/dE_g in mbarn/eV 
! parameters: transition energy Etrans between Kamae and QGSJET
!             width of transition


  T = 0.681d9
  E_p = m_p+T 
  write(*,*) 'T and E_p/GeV',T/1.d9,E_p/1.d9

  open(21,file='T0.681')                     ! plots in style of Chuck's
  sym = 0
  do i=1,10000,100
     E_g = 0.999**i * E_p
     call spectrum(E_p,E_g,1,s1)
     write(21,*) E_g/1.d6,s1
  end do
  write(21,*) 
  sym = 1
  do i=1,10000,100
     E_g = 0.999**i * E_p
     call spectrum(E_p,E_g,1,s1)
     write(21,*) E_g/1.d6,s1
  end do
  close(21)

! comparison of reac 1-4 at 100 GeV, setting s=0 or 1 at line 114: 
  E_p = 100.d9

  open(21,file='nuclei')                            
!  open(21,file='nuclei_QGS')                     
!  open(21,file='nuclei_Kamae')                   
  do i=1,10000,100
     E_g = 0.999**i * E_p
     call spectrum(E_p,E_g,1,s1)
     call spectrum(E_p,E_g,2,s2)
     call spectrum(E_p,E_g,3,s3)
     call spectrum(E_p,E_g,4,s4)
     write(21,*) real(E_g/1.d6),real(s1),real(s2),real(s3),real(s4)
  end do
  close(21)

end subroutine  photonspectra ! AWS20130503
!end program photonspectra      AWS20130503
!=============================================================================!
!=============================================================================!
subroutine spectrum(E_p,E_g,reac,spec)
  use spectra
  implicit none
  integer reac
  double precision :: E_p,E_g,E,spec,x,s,q,k
  double precision, parameter :: Esym=67.285d6
  double precision kamae_nd,spec_QGS

  if (E_g<Esym.and.sym==1) then
     E = Esym**2/E_g
  else
     E = E_g
  end if
  
  q = spec_QGS(E_p,E_g,reac) *1.d6/E_g

  select case (reac)
  case (1)                          ! pp
     x = log10(E_p)
     s = (tanh((x-x0)/(width*x0))+1.d0)/2.d0
     k = kamae_nd(E_p/1.d9,E/1.d9) *1.d6/E
  case (2)                          ! p+He
     x = log10(E_p)
     s = (tanh((x-x0)/(width*x0))+1.d0)/2.d0
     k = 3.03d0 * kamae_nd(E_p/1.d9,E/1.d9) *1.d6/E
  case (3)                          ! He+p
     x = log10(E_p/4.d0)
     s = (tanh((x-x0)/(width*x0))+1.d0)/2.d0
     k = 4.d0 * kamae_nd(E_p/4.d9,E/1.d9) *1.d6/E
  case (4)                          ! He+He
     x = log10(E_p/4.d0)
     s = (tanh((x-x0)/(width*x0))+1.d0)/2.d0
     k = 12.1d0 * kamae_nd(E_p/4.d9,E/1.d9) *1.d6/E 
  end select

!  s=1.d0
  spec = k*(1.d0-s) +q*s

end subroutine spectrum
!=============================================================================!
!=============================================================================!
double precision function spec_QGS(E_p,E_g,k)
  use spectra
  implicit none
  integer i,j,k
  double precision E_p,E_g
  double precision x,y,y1,y2,y3,y4,t,u,r

  if (E_p<1.5d10) then
     spec_QGS = 0.d0
     return
  end if

  x=log10(E_p)
  y=log10(E_g)

  i=int(log(E_p/E0)/log(d))
  j=int((y-d_f)/dn)

  y1 = fqgs(i,j,k)
  y2 = fqgs(i+1,j,k)
  y4 = fqgs(i,j+1,k)
  y3 = fqgs(i+1,j+1,k)

  t = (x-xa(i))/(xa(i+1)-xa(i))
  u = (y-ya(j))/(ya(j+1)-ya(j))
  r = (1.d0-t)*(1.d0-u)*y1 + t*(1.d0-u)*y2 + t*u*y3  + (1.d0-t)*u*y4 

  spec_QGS = 10.d0**r 

end function spec_QGS
!=============================================================================!
!=============================================================================!
