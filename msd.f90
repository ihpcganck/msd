! 
!     Copyright (C) 2021  Chee Kwan Gan (ihpcganck@gmail.com)
!     Copyright (C) 2020  Chee Kwan Gan (ihpcganck@gmail.com)
! 
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
! 

module extra
  use commod
  implicit none

  character(len=*),parameter :: parfile='msd.par'

  integer,parameter :: absparu   =11
  integer,parameter :: ufile1    =12
  integer,parameter :: dfile1    =13
  integer,parameter :: zeropu    =14
  integer,parameter :: msdu      =15

  real(double),parameter :: msd_unit=1.0d-3*Angstrom**2.0d0

  integer,parameter :: nspaces=4

  type bstype
    integer :: ns,np,nb
    real(double),allocatable :: omega(:,:,:)
    real(double),allocatable :: dis(:)
    complex(double),allocatable :: eigenvec(:,:,:,:,:)
  end type bstype
contains

  subroutine msd(qind,natom,zerostruc,Tbegin,Tend,Tn,msdatom,nb,D1,evec1,w1,zeroordered,work,lwork,rwork,lrwork)
    integer :: nb,qind,natom,Tn
    type(supercell) :: zerostruc
    real(double) :: msdatom(natom,0:Tn)
    complex(double) :: evec1(nb,nb)
    integer :: lwork,lrwork
    complex(double) :: work(lwork)
    complex(double) :: D1(nb,nb)
    real(double) :: w1(nb)
    real(double) :: rwork(lrwork)
    real(double) :: zeroordered(nb)
    integer :: info
    integer :: i,j,k
    real(double) :: omega,Tbegin,Tend,Temper
    complex(double) :: e(3),dotpr
    real(double) :: r,fac,massj

    qind = qind

    evec1(1:nb,1:nb) = D1(1:nb,1:nb)
    call zheev('V','U',nb,evec1(1,1),nb,w1(1),work(1),lwork,rwork(1),info)
    if(info /= 0) then
      write(*,*) 'info = ',info
      stop 'err: zheev error.'
    endif

    do i = 1, nb
      zeroordered(i) = w1(i)
    enddo

    do i = 1, nb
      omega  = sqrt(abs(zeroordered(i)))*Rydoverh2invcm
      if(zeroordered(i) < zero) then
        zeroordered(i) = -omega
      else
        zeroordered(i) = omega
      endif
    enddo

     do i = 1, nb
       if(qind <= 10) then
         write(*,*) 'zeroordered(i) = ',zeroordered(i)
       endif

       do j = 1, natom
         e(1:3) = evec1(  (j-1)*3+1: 3*j,i)
         dotpr = zero
         do k = 1, 3
           dotpr = dotpr + conjg(e(k))*e(k)
         enddo
         do k = 0, Tn
           Temper = Tbegin + k*(Tend-Tbegin)/(Tn*1d0)
           massj = zerostruc%at(j)%mass*amu
           fac = real(dotpr)*hPlanck/(8d0*pi*pi*zeroordered(i)*invcm*massj)

           r = hPlanck*zeroordered(i)*invcm/(kBoltz*Temper)/two
           msdatom(j,k) = msdatom(j,k) + fac/tanh(r)/(msd_unit)
         enddo
       enddo
     enddo
  end subroutine msd

  subroutine ReadBSandEigenModes(ufile,bs)
    integer :: ufile
    type(bstype) :: bs
    integer :: natom,dis_ind,ndis,ns,np,nb,i,j
    character(len=300) :: tmpstr1

    read(ufile,*) tmpstr1,ns,np,nb
    bs%ns = ns
    bs%np = np
    bs%nb = nb
    write(*,*) 'bs%ns,bs%np,bs%nb=',ns,np,nb

    natom = nb/3
    if(nb /= natom*3) then
      write(*,'(A,I6,I6)') 'in ReadBSandEigenModes: not consistent: nb,natom =',nb,natom
      stop 'err: not divisible by 3.'
    endif
    write(*,*) 'natom = ',natom

    allocate(bs%omega(nb,np,ns))
    allocate(bs%eigenvec(3,natom,nb,np,ns))
    ndis = np*ns
    allocate(bs%dis(ndis))

    dis_ind = 0
    do i = 1, ns
      do j = 1, np
        dis_ind = dis_ind + 1
        read(ufile,*) bs%dis(dis_ind),bs%omega(1:nb,j,i)
      enddo
    enddo
    if(dis_ind /= ndis) then
      write(*,*) 'dis_ind,ndis = ',dis_ind,ndis
      stop 'err: ndis problem .'
    endif
    write(*,*) 'bs%omega(nb,np,ns)=',bs%omega(nb,np,ns)
  end subroutine ReadBSandEigenModes

  subroutine read_dyn(du,totalqpoints,natom,nb,dyn)
    integer :: du,totalqpoints,natom,nb
    complex(double) :: dyn(nb,nb,totalqpoints)
    integer :: qpt,tmpi,tmpj
    real(double) :: w1,w2,w3,w4,w5,w6
    complex(double) :: z1,z2,z3
    integer :: i,j
    character(len=DSL) :: tmpstr

    do qpt = 1, totalqpoints
      read(du,*) tmpstr
      read(du,*) tmpstr
      do i = 1, natom
        do j = 1, natom
          read(du,*) tmpi,tmpj
          if(tmpi /= i .and. tmpj /= j) then
            write(*,*) 'tmpi,i,tmpj,j=',tmpi,i,tmpj,j
            stop 'err: wrong index.'
          endif

          read(du,*) w1,w2,w3,w4,w5,w6
          z1 = dcmplx(w1,w2)
          z2 = dcmplx(w3,w4)
          z3 = dcmplx(w5,w6)
          dyn(  (i-1)*3+1, (j-1)*3+1:j*3,qpt) = (/z1,z2,z3/)
          read(du,*) w1,w2,w3,w4,w5,w6
          z1 = dcmplx(w1,w2)
          z2 = dcmplx(w3,w4)
          z3 = dcmplx(w5,w6)
          dyn(  (i-1)*3+2, (j-1)*3+1:j*3,qpt) = (/z1,z2,z3/)
          read(du,*) w1,w2,w3,w4,w5,w6
          z1 = dcmplx(w1,w2)
          z2 = dcmplx(w3,w4)
          z3 = dcmplx(w5,w6)
          dyn(  (i-1)*3+3, (j-1)*3+1:j*3,qpt) = (/z1,z2,z3/)
        enddo
      enddo
    enddo
  end subroutine read_dyn

end module extra

program sample
  use extra
  implicit none

  type(bstype) :: zerobs
  integer :: i,j,k,ns,np,nb,qind
  integer :: natom
  character(len=DSL) :: zero_epsfile
  character(len=DSL) :: zero_dfile
  character(len=DSL) :: zero_poscar
  character(len=DSL) :: dir,absparfile
  character(len=DSL) :: filestem
  character(len=DSL) :: inputformat

  integer :: lwork,lrwork
  complex(double),allocatable :: work(:)
  complex(double),allocatable :: evec1(:,:)

  integer :: totalqpoints
  complex(double),allocatable :: dyn1(:,:,:)
  integer :: m
  complex(double),allocatable :: D1(:,:)
  real(double),allocatable :: zeroordered(:)
  real(double),allocatable :: msdatom(:,:)
  real(double),allocatable :: w1(:),rwork(:)
  integer :: nq,Tn
  real(double) :: Tbegin,Tend
  type(dtype) :: variable,defaultv
  type(supercell) :: zerostruc
  real(double) :: massj,Temper
  character(len=DSL) :: outputfile

  write(*,*)
  write(*,*)
  call fargn(1)
  call fargv(1,dir)
  write(*,'(A)') 'dir is '//trim(dir)
  absparfile=trim(dir)//'/'//trim(parfile)
  write(*,'(A)') 'To open the file '//trim(absparfile)
  open(absparu,file=trim(absparfile),status='old',action='read')
  write(*,'(A)') trim(absparfile)//' is opened.'
  close(absparu)
  write(*,'(A)') trim(absparfile)//' is immediately closed.'

  call read_gendata(rtype,absparu,trim(absparfile),'Tbegin',variable,NONDEFAULTOPT,defaultv)
  Tbegin = variable%r
  call read_gendata(rtype,absparu,trim(absparfile),'Tend',variable,NONDEFAULTOPT,defaultv)
  Tend = variable%r
  call read_gendata(itype,absparu,trim(absparfile),'Tn',variable,NONDEFAULTOPT,defaultv)
  Tn = variable%i

  zero_epsfile=trim(dir)//'/1-bs.dat'

  write(*,*)
  open(ufile1,file=trim(zero_epsfile),status='old')
  write(*,*) 'reading zerobs:'
  call ReadBSandEigenModes(ufile1,zerobs)
  write(*,*)
  close(ufile1)

  ns = zerobs%ns
  np = zerobs%np
  nb = zerobs%nb

  natom = nb/3
  if(nb /= natom*3) then
    write(*,*) 'not consistent: nb,natom =',nb,natom
    stop 'err: not divisible by 3.'
  endif

  allocate(msdatom(natom,0:Tn))

  zero_dfile=trim(dir)//'/allq.dyn'
  write(*,'(A)') 'zero_dfile is '//trim(zero_dfile)

  totalqpoints=ns*np
  write(*,*) 'totalqpoints = ',totalqpoints
  allocate(dyn1(nb,nb,totalqpoints))

  open(dfile1,file=trim(zero_dfile),action='read',status='old')
  call read_dyn(dfile1,totalqpoints,natom,nb,dyn1(1,1,1))
  close(dfile1)

  allocate(D1(nb,nb))
  allocate(zeroordered(nb))

  allocate(w1(nb))
  lwork=2*nb
  lrwork=3*nb

  allocate(work(lwork))
  allocate(rwork(lrwork))
  allocate(evec1(nb,nb))

  nq = ns*np

  do i = 1, natom
    do j = 0, Tn
      msdatom(i,j) = zero
    enddo
  enddo
  qind = 0

  zero_poscar=trim(dir)//'/p.vasp'
  write(*,'(A)') 'zero_poscar is '//trim(zero_poscar)
  call read_struc(zeropu,zero_poscar,filestem,inputformat,zerostruc)

  do i = 1, natom
    massj = zerostruc%at(i)%mass
    write(*,*) 'i,mass = ',i,massj
  enddo

  do i = 1, ns

    do j = 1, np
      qind = qind + 1
      write(*,'(A,I8,A,I8)') 'qind = ',qind, '/',nq

      do k = 1, nb
        do m = 1, nb
          D1(m,k) = dyn1(m,k,qind)
        enddo
      enddo

      call msd(qind,natom,zerostruc,Tbegin,Tend,Tn,msdatom(1,0),nb,D1(1,1),evec1(1,1),w1(1),zeroordered(1),work(1),lwork,rwork(1),lrwork)
    enddo
  enddo

  do i = 1, natom
    do j = 0, Tn
      msdatom(i,j) = msdatom(i,j)/(nq*1.0d0)
    enddo
  enddo

  do i = 1, natom

    outputfile=trim(dir)//'/msd-atom-'//num2str(i,nspaces)//'-T-gives-msd.dat'

    open(unit=msdu,file=trim(outputfile),status='replace')
    do j = 0, Tn
      Temper = Tbegin + j*(Tend-Tbegin)/(Tn*1.0d0)
      write(msdu,*) Temper,msdatom(i,j)
    enddo
    close(msdu)
    write(*,'(A)') 'xmgrace '//trim(outputfile)
  enddo

end program sample
