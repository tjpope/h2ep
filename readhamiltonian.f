!*********************************************************************!
!*******  readhamiltonian  *******************************************!
!*********************************************************************!
        subroutine readhamiltonian
        use commandments
        use hamiltonian
        implicit none
        integer::nspin,columns,ncomp,i,j,k,l
        integer,allocatable,dimension(:)::z
        integer,allocatable,dimension(:,:)::q
        double precision::sr,si,hr,hi,kpoints(3)
        complex*16::im
        parameter(im=(0,1))
        character(200)::a,frmt
        n=0; columns=0; ef=0.0; nspin=0
10      read(*,'(a200)',end=20) a
        if(index(a,'name: nspin').gt.0) then
          read(*,*) a; read(*,*) nspin
        elseif(index(a,'name: FermiE').gt.0) then
          read(*,*) a; read(*,*) ef
        elseif (index(a,'name: iorb').gt.0) then
          read(*,*) a; read(*,'(a200)') a
          i = scan(a,':',BACK=.true.)
          read(a(i+1:),'(i20)') n
          read(*,'(a200)') a
          i = scan(a,':',BACK=.true.)
          read(a(i+1:),'(i10)') columns
          if (n.gt.100000) then
            write(*,'(a,/,i0)') 'ERROR:: n is too large:', n
            stop
          endif
          allocate(hamcomp(n),q(n,3),z(n),h(n,n),s(n,n))
          hamcomp=0; h=(0,0); s=(0,0)
          do i=1,n
            read(*,*) hamcomp(i),(q(i,j),j=1,3),z(i)
          enddo
          allocate(hatoms(hamcomp(n),2))
          j=1; k=1
          do i=1,n
            if(hamcomp(i).eq.hamcomp(i+1)) then
              j=j+1
              if(i+1.gt.n) then
                hatoms(k,1)=hamcomp(i); hatoms(k,2)=j; goto 10
              endif
            else
              hatoms(k,1)=hamcomp(i); hatoms(k,2)=j; j=1; k=k+1
            endif
          enddo
        elseif(index(a,'name: kpoints').gt.0) then
          read(*,*) a; read(*,*) a; read(*,*) a
          read(*,*) kpoints(1:3)
        elseif (index(a,'name: HS').gt.0) then
          read(*,*) a; read(*,'(a200)') a
          i = scan(a,':',BACK=.true.)
          read(a(i+1:),'(i10)') ncomp
          ncomp = int(ncomp/maxval(kpoints))
          read(*,'(a200)') a
          i = scan(a,':',BACK=.true.)
          read(a(i+1:),'(i10)') columns
          do l=1,ncomp
            read(*,*) k,i,j,sr,si,hr,hi
            h(i,j)=hr+im*hi; s(i,j)=sr+im*si
          enddo
        endif
        if(verbose) then
         open(60,file='read.ham')
         write(frmt,'("(",i0,"(f3.1,1x))")')n
         do i=1,n
          write(60,frmt) (real(h(i,j)),j=1,n)
         enddo
         close(60)
        endif
        goto 10
20      deallocate(hamcomp)
        return
        end subroutine readhamiltonian
!*********************************************************************!
