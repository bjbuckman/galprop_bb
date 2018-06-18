!=============================================================================!
!=========================changed to nm======================================!
subroutine init
  use spectra
  implicit none
  integer i,j 
  double precision x,y,f

  d_f = log10(E_min)-0.1d0         
  x0 = log10(Etrans)

  open(21,file='gam_ppR')
  do i=1,n1
     do j=1,n2
        read(21,*) x,y,f
        xa(i)=log10(x)
        ya(j)=log10(y)
        if (f>0.d0) then
           fqgs(i,j,1)=log10(f)
        else
           fqgs(i,j,1)=-20.d0
        end if
     end do
     read(21,*) 
  end do
  close(21)

  open(21,file='gam_pHe_nm')
  do i=1,n3
     do j=1,n2
        read(21,*) x,y,f
        xa(i)=log10(x)
        ya(j)=log10(y)
        if (f>0.d0) then
           fqgs(i,j,2)=log10(f)
        else
           fqgs(i,j,2)=-20.d0
        end if
     end do
     read(21,*) 
  end do
  close(21)


  open(21,file='gam_Hep_nm')
  do i=1,n4
     do j=1,n2
        read(21,*) x,y,f
        xa(i)=log10(x)
        ya(j)=log10(y)
        if (f>0.d0) then
           fqgs(i,j,3)=log10(f)
        else
           fqgs(i,j,3)=-20.d0
        end if
     end do
     read(21,*) 
  end do
  close(21)

  open(21,file='gam_HeHeR')
  do i=1,n5
     do j=1,n2
        read(21,*) x,y,f
        xa(i)=log10(x)
        ya(j)=log10(y)
        if (f>0.d0) then
           fqgs(i,j,4)=log10(f)
        else
           fqgs(i,j,4)=-20.d0
        end if
     end do
     read(21,*) 
  end do
  close(21)

end subroutine init
!=============================================================================!
!=============================================================================!
