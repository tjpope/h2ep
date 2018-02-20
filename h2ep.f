!*********************************************************************!
!*******  h2ep  ******************************************************!
!*********************************************************************!
!         Hamiltonian 2 Electronic Properties                         !
!*********************************************************************!
!         The tranmission probability is calculated using the non-    !
!         equilibrium Green’s function method. The hamiltonian from   !
!         SIESTA is needed in the GOLLUM format. Also, the size of    !
!         the leads (how many atoms). I assume the leads are          !
!         symmetric. If this isn't the case for you, go find another  !
!         code.                                                       !
!                                                                     !
!         I assume the hamiltonian is in this form:                   !
!                                                                     !
!              | XX  XX  XX  XX  XX  XX  XX |                         !
!              | XX  h0  h1  00  00  00  XX |                         !
!              | XX  h1  h0  gl  00  00  XX |                         !
!          H = | XX  00  gl  hs  gr  00  XX |                         !
!              | XX  00  00  gr  h0  h1  XX |                         !
!              | XX  00  00  00  h1  h0  XX |                         !
!              | XX  XX  XX  XX  XX  XX  XX |                         !
!                                                                     !
!         You need to tell me where the first h0 starts and ends.     !
!         I'll do the rest                                            !
!                                                                     !
!         For more info, check out the GOLLUM paper;                  !
!                                                                     !
!          Ferrer, Jaime, et al. "GOLLUM: a next-generation           !
!          simulation tool for electron, thermal and spin transport." !
!          New Journal of Physics 16.9 (2014): 093029.                !
!                                                                     !
!*********************************************************************!
!         This code is parallel and uses the OPENMP enviroment.       !
!         It also uses LAPACK and BLAS.                               !
!         It also uses stuff from fortran 2003 or later and WILL NOT  !
!         compile properly with previous versions.                    !
!*********************************************************************!
!         I only calculate the transmission probability. If you want  !
!         the thermopower or the IV characteristics, check out my     !
!         friend mpft.f - found in the src/ folder.                   !
!*********************************************************************!
!       I'm easy to compile, just let the makefile know where you     !
!       keep LAPACK and BLAS and type:                                !
!                                                                     !
!         make                                                        !
!                                                                     !
!       Then, for further instructions, type:                         !
!                                                                     !
!         ./h2ep --help                                               !
!                                                                     !
!*********************************************************************!
        program h2ep
        use my_types
        use commandments
        use hamiltonian
        use h0h1
        use omp_lib
        implicit none
        integer::num_threads,thread_num,per_thread,first,last,mem_store,mem_green,mem_invsn,ne_save
        double precision,allocatable,dimension(:,:)::transport
        call moses
        call readhamiltonian
        call cutuphamiltonian
        if(verbose) then
         mem_store=2*n*n*16+(34*l_lead%n0+2)*l_lead%n0*16+4*l_lead%n0
         mem_green=(18*l_lead%n0+14)*l_lead%n0*16+16*l_lead%n0*8+5*l_lead%n0
         mem_invsn=(5*n*n+16*l_lead%n0*l_lead%n0-4*n*l_lead%n0)*16
         ne_save=ne
        endif
!$omp parallel private(first,last)
!$omp single
        num_threads=omp_get_num_threads()
        if(ne.lt.num_threads) ne=num_threads
        per_thread=int(ne/num_threads)
        ne=num_threads*per_thread
        allocate(transport(ne,3))
        if(verbose) then
         write(*,'("# I used this much memory (...ish): ",i0," BYTES")') max(mem_store+mem_green,mem_store+mem_invsn)*num_threads
         write(*,'("# The number of threads was ",i0," and the number of points calculated was ",i0)') num_threads,ne
         if(ne_save.ne.ne)write(*,'(a)')"# I calculated more points than you asked for because the number of points was not divisible by the number of threads... you're welcome"
        endif
!$omp end single
        first=per_thread*(omp_get_thread_num())+1
        last=per_thread*(omp_get_thread_num()+1)
        call parallel_transport(l_lead,r_lead,ns,cl,scl,hs,ss,cr,scr,first,last,transport(first:last,:))
!$omp end parallel
        do epoint=1,ne; write(*,'(t5,2(es20.10,1x),f5.2)') transport(epoint,:); enddo
        contains
        subroutine parallel_transport(l_lead,r_lead,ns,cl,scl,hs,ss,cr,scr,first_point,last_point,t)
        use my_types
        use commandments
        use hamiltonian
        implicit none
        type(periodic) l_lead,r_lead
        type(periodic) l_leade, r_leade
        type(physics) l_leadp, r_leadp     
        integer::ns,first_point,last_point,point  
        double precision::e,e0,del,der
        double precision,dimension(last_point-first_point+1,3)::t
        complex*16::hs(ns,ns),ss(ns,ns),cl(l_lead%n0,ns),scl(l_lead%n0,ns),cr(r_lead%n0,ns),scr(r_lead%n0,ns)
        complex*16,allocatable,dimension(:,:)::green,gf
        logical::nochannels
        integer::epoint_thread
        allocate(l_leade%h0(l_lead%n0,l_lead%n0),l_leade%h1(l_lead%n0,l_lead%n0),r_leade%h0(r_lead%n0,r_lead%n0),r_leade%h1(r_lead%n0,r_lead%n0),l_leade%s0(l_lead%n0,l_lead%n0),l_leade%s1(l_lead%n0,l_lead%n0),r_leade%s0(r_lead%n0,r_lead%n0),r_leade%s1(r_lead%n0,r_lead%n0),l_leade%h0i(l_lead%n0,l_lead%n0),r_leade%h0i(r_lead%n0,r_lead%n0),l_leade%k(2*l_lead%n0,2*l_lead%n0),r_leade%k(2*r_lead%n0,2*r_lead%n0))
        der=dble((emax-emin)/(ne-1)); del=dble(emin-der)
        do epoint_thread=first_point,last_point
         point=epoint_thread-first_point+1
         nochannels=.false.
         e0=dble(der*epoint_thread+del); e=e0+ef
         call initialize_lead(l_lead,l_leade,e)
         call initialize_lead(r_lead,r_leade,e)
         allocate(l_leadp%gv(2*l_lead%n0),l_leadp%v_op(l_lead%n0,l_lead%n0),l_leadp%vectors(l_lead%n0,l_lead%n0*2),l_leadp%duals(l_lead%n0*2,l_lead%n0),l_leadp%oi(l_lead%n0),l_leadp%oo(l_lead%n0),r_leadp%gv(2*r_lead%n0),r_leadp%v_op(r_lead%n0,r_lead%n0),r_leadp%vectors(r_lead%n0,r_lead%n0*2),r_leadp%duals(r_lead%n0*2,r_lead%n0),r_leadp%oi(r_lead%n0),r_leadp%oo(r_lead%n0),green(l_lead%n0+r_lead%n0,l_lead%n0+r_lead%n0),gf(l_lead%n0+r_lead%n0,l_lead%n0+r_lead%n0))
         green=(0,0); l_leadp%n0=l_lead%n0; r_leadp%n0=r_lead%n0;
         call getgreens(e0,l_lead%n0,l_leade%h0,l_leade%h1,l_leade%s0,l_leade%s1,green(:l_lead%n0,:l_lead%n0),l_leadp%v_op,l_leadp%gv,l_leadp%nc,l_leadp%oo,l_leadp%oi,l_leadp%vectors,l_leadp%duals,.true.,verbose,nochannels)
         call getgreens(e0,r_lead%n0,r_leade%h0,r_leade%h1,r_leade%s0,r_leade%s1,green(l_lead%n0+1:,l_lead%n0+1:),r_leadp%v_op,r_leadp%gv,r_leadp%nc,r_leadp%oo,r_leadp%oi,r_leadp%vectors,r_leadp%duals,.false.,verbose,nochannels)
         t(point,2)=0.0; t(point,1)=e0; t(point,3)=dble(l_leadp%nc(1))
         if(.not.nochannels) then
          call inverteffham(l_lead%n0,r_lead%n0,green,ns,l_leade%h1,r_leade%h1,hs-e*ss,cl-e*scl,cr-e*scr,l_leade%k,r_leade%k,gf)
          call getscatteringmatrix(l_leadp,r_leadp,gf,t(point,2))
         endif
         deallocate(l_leadp%gv,r_leadp%gv,l_leadp%v_op,r_leadp%v_op,l_leadp%vectors,r_leadp%vectors,l_leadp%oi,l_leadp%oo,r_leadp%oi,r_leadp%oo,l_leadp%duals,r_leadp%duals,green,gf)
        enddo
        end subroutine parallel_transport
        end program h2ep
