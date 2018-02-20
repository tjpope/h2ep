!*********************************************************************!
!       More Pictures From Transport                                  !
!*********************************************************************!
!       This program generates a bunch of extra figures to plot       !
!        from the transport. Respectively, they are,                  !
!        G(E_{F})/G_{0} -- Thermally averaged Transport ($2) tp.dat   !
!                           and iv.dat                                !
!        S(E_{F})       -- Thermopower ($3) tp.dat                    !
!        <x>            -- standard deviation of p(x) ($4) tp.dat     !
!        ke             -- thermal conductance ($5) tp.dat            !
!        ZT             -- Figure of Merit ($6) tp.dat                !
!        I/V            -- Current vs. Voltage ($3) iv.dat            !
!*********************************************************************!
!       The program reads the first and second column of a data file  !
!        as the energy and corresponding transmission probability     !
!        respectively. If there are blank lines, the calculation      !
!        will fail. For the best results, use my friend h2ep to       !
!        calculate the transmission probability.                      !
!*********************************************************************!
!       All the science is explained better in the liturature;        !
!                                                                     !
!         Algharagholy, Laith A., et al. "Tuning thermoelectric       !
!         properties of graphene/boron nitride heterostructures."     !
!         Nanotechnology 26.47 (2015): 475401.                        !
!                                                                     !
!       Just remember, the thermal stuff ignores phonons, which might !
!        be important. So be careful.                                 !
!        Also, the I/V calculations assume that the transport dosn't  !
!        change with bias voltage - which it sometimes does.          !
!*********************************************************************!
!       I'm easy to compile, just type (no need for any fancy flags): !
!                                                                     !
!         gfortran -o mpft mpft.f                                     !
!                                                                     !
!       Then to run me, put the transport data it t.dat and type:     !
!                                                                     !
!         ./mpft < t.dat                                              !
!                                                                     !
!*********************************************************************!
        module constants
        real::kT,kb,T,del,h,el,G0,zero,pi
        parameter(T=300,kb=8.6173324e-5,zero=1e-12)
        parameter(h=6.626068e-34,G0=7.7480917346e-5,el=1.60217646e-19)
        parameter(pi=3.14159265359)
        end module constants
!*********************************************************************!
        program mpft
        use constants
        implicit none
        integer::ne,ie,nfe,ife,ifs
        real::l0,l1,l2,l1t,l2t,ln,lnt,trblmkr,fermidist,dfe,de,Ef,
     +                    Efmin,G,S,Pe,ka,ZT,cnst,integrate,temperature
        real,allocatable,dimension(:)::tr,p,E,iv
        character(300)::hash
        logical::verbose
        call chat(verbose,temperature,del)
        if(temperature.lt.0.0) temperature=T
        kT=kb*temperature
        call genlog
        open(60,file='in.log')
        open(61,file='iv.dat')
        open(62,file='tp.dat')
10      read(60,*) hash
        if (index(hash,'#').ne.0) goto 10
        backspace(60)
        ne=0; ifs=0
20      read(60,*,end=30); ne=ne+1; goto 20
30      allocate(tr(ne),p(ne),E(ne),iv(ne))
        do ie=0,ne; backspace(60); enddo
        do ie=1,ne
         read(60,*) E(ie),tr(ie) 
         if(E(ie).gt.E(1)+del.and.ifs.eq.0) ifs=ie
        enddo
        close(60)
        if(verbose) write(61,'(3a20)') '#','V','T(E=V)','I'
        cnst=(4*el**2*pi)/h
        do ie=1,ne
         if (E(ie).lt.0) then
          iv(ie)=-integrate(tr,e,ne,E(ie),-E(ie))*cnst
         else
          iv(ie)=integrate(tr,e,ne,-E(ie),E(ie))*cnst
         endif
        enddo
        do ie=1,ne
         write(61,'(t5,3es20.10)') E(ie),tr(ie),iv(ie)
        enddo
        close(61)
        nfe=ne-(2*ifs)
        if(verbose) write(62,'(a1,6a20)') '#','Ef','l0','S','Pe','ka',
     +                                                             'ZT'
        do ife=1,nfe
         Ef=E(ifs+ife)
         do ie=1,ne; p(ie)=tr(ie)*fermidist(E(ie),Ef); enddo
         l0=ln(0,Ef,E,p,ne); l1=ln(1,Ef,E,p,ne); l2=ln(2,Ef,E,p,ne) 
         l1t=l1/l0; l2t=l2/l0; trblmkr=(l2t-(l1t**2))
         G=G0*l0; S=-l1t/T; ka=(G*trblmkr)/T; ZT=(l1t**2)/trblmkr
         write(62,'(t5,6es20.10)') Ef,l0,S,l1t,ka,ZT
        enddo
        close(62)
        end program mpft
!*********************************************************************!
        function ln(n,Ef,E,p,ne)
        use constants
        implicit none
        integer::n,ne,i
        real::ln,Ef,emin,emax,integrate
        real,dimension(ne)::p,E,f 
        emin=Ef-del; emax=Ef+del
        do i=1,ne; f(i)=(((E(i)-Ef)**n)*p(i)); enddo
        ln=integrate(f,E,ne,emin,emax)
        end function
!*********************************************************************!
        function integrate(f,x,n,xmin,xmax)
        integer::n,i
        real::xmin,xmax,integrate
        real,dimension(n)::x,f
        integrate=0
        do i=2,n
         if(x(i-1).ge.xmin.and.x(i).le.xmax) 
     +              integrate=integrate+(f(i-1)+f(i))*(x(i)-x(i-1))*0.5
        enddo
        end function        
!*********************************************************************!
        function fermidist(E,Ef)
        use constants
        implicit none
        real::E,Ef,fermidist,x
        x=exp((E-Ef)/kT); fermidist=x/(((x+1)**2)*kT)
        end function
!*********************************************************************!
        subroutine genlog
        implicit none
        integer::i
        character(200)::a
        open(60,file='in.log')
10      read(*,'(a200)',end=20) a
        write(60,'(a)') trim(adjustl(a))
        goto 10
20      close(60)
        return
        end subroutine genlog
!*********************************************************************!
        subroutine chat(verbose,temperature,del)
        implicit none
        integer::i
        real::temperature,del
        character(200),allocatable,dimension(:)::comlin        
        logical::verbose
        allocate(comlin(command_argument_count()))
        do i=1,command_argument_count()
         call getarg(i,comlin(i))
        enddo
        verbose=.false.; temperature=-1; del=0.25
        do i=1,command_argument_count()
         if(index(comlin(i),'-v').gt.0) then
          verbose=.true. 
         elseif(index(comlin(i),'-T').gt.0) then 
          read(comlin(i+1),'(f10.5)') temperature
         elseif(index(comlin(i),'-d').gt.0) then 
          read(comlin(i+1),'(f10.5)') del
         elseif(index(comlin(i),'--help').gt.0) then
          write(*,100); stop; 
         endif
        enddo
        return
100     format(t5,'mpft: instructions',/,/,t5,
     +   'mpft {command line input} < {transport file}',/,
     +   /,t5,'-T $f',t30,
     +   'Temperature used in calculation [default = 300]',
     +   /,t5,'-d $r',t30,
     +   'Width of integeral in calculation [default = 0.25]',
     +   /,t5,'-v',t30,
     +   'Let me tell you about the calculation [default .false.]',/)
        end subroutine chat
!*********************************************************************!
