!*********************************************************************!
!*******  desingularize  *********************************************!
!*********************************************************************!
        subroutine desingularize(n,mat)
        implicit none
        integer::n
        complex*16,dimension(n,n)::mat
        double precision,dimension(n,n)::random
        double precision::zero,delta
        parameter(zero=1d-9)
        delta=zero*maxval(abs(mat))
        call random_seed(); call random_number(random)
        mat=mat+delta*(-1.d0+2.d0*random)
        end subroutine desingularize
