!*********************************************************************!
!*******  moses  *****************************************************!
!*********************************************************************!
        subroutine moses
        use commandments
        implicit none
        integer::i,j
        character(200),allocatable,dimension(:)::comlin        
        allocate(comlin(command_argument_count())); do i=1,command_argument_count(); call getarg(i,comlin(i)); enddo
        ne=10; n0min=0; n0max=0; emin=-2.0; emax=2.0; verbose=.false.
        do i=1,command_argument_count()
         if (index(comlin(i),'-epoints').gt.0) then; read(comlin(i+1),'(i10)') ne
         elseif(index(comlin(i),'-emin').gt.0) then; read(comlin(i+1),'(f10.5)') emin
         elseif(index(comlin(i),'-emax').gt.0) then; read(comlin(i+1),'(f10.5)') emax
         elseif(index(comlin(i),'-h0h1').gt.0) then; j=scan(comlin(i+1),':',BACK=.true.); read(comlin(i+1)(1:j-1),'(i10)') n0min; read(comlin(i+1)(j+1:),'(i10)') n0max
         elseif(index(comlin(i),'-v').gt.0) then; verbose=.true.
         elseif(index(comlin(i),'--help').gt.0) then; write(*,100); stop; endif
        enddo
        if(n0min.eq.0.and.n0max.eq.0) stop'error in lead dimensions'
        if(n0min.ge.n0max) stop'error in lead dimensions'
        return
100     format(t5,'h2ep: instructions',/,/,t5,'h2ep {command line input} < {hamiltonian file} > {output file}',/,/,t5,'-epoints $i',t30,'number of energy-points used in calculation. [default = 10]',/,t5,'-emin $r / -emax $r',t30,'Minimum/Maximum values in the energy range [default = -2.0 / 2.0]',/,t5,'-h0h1 $i:$i',t30,'Initial and final atom number in h0 [no default]',/,t5,'-v',t30,'Let me tell you about the calculation [default .false.]',/)
        end subroutine moses
!*********************************************************************!
