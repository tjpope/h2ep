!*********************************************************************!
!*******  make_k  ****************************************************!
!*********************************************************************!
        subroutine make_k(n0,h0,h0i,h1,k)
        implicit none
        integer::n0
        complex*16,dimension(n0,n0)::h0,h1,h0i
        complex*16,dimension(2*n0,2*n0)::k
        k(1:n0,1:n0)=h0i
        k(n0+1:,n0+1:)=h0
        k(1:n0,n0+1:)=transpose(conjg(h1))
        k(n0+1:,1:n0)=h1
        end subroutine make_k
