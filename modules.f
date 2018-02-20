!*********************************************************************!
        module my_types
        type periodic
         integer::n0
         complex*16,allocatable,dimension(:,:)::h0,h1,s0,s1,h0i,s0i,k
        end type periodic
        type physics
         integer::n0
         integer,dimension(2)::nc
         integer,allocatable,dimension(:)::oi,oo
         complex*16,allocatable,dimension(:)::gv
         complex*16,allocatable,dimension(:,:)::v_op,vectors,duals,gr
        end type physics
        end module my_types
!*********************************************************************!
        module commandments
        integer::n0min,n0max,epoint,ne
        double precision::emin,emax
        logical::verbose
        end module 
!*********************************************************************!
        module hamiltonian
        integer::n
        integer,allocatable,dimension(:)::hamcomp
        integer,allocatable,dimension(:,:)::hatoms
        double precision::ef
        complex*16,allocatable,dimension(:,:)::h,s
        end module
!*********************************************************************!
        module h0h1
        use my_types, only: periodic
        type(periodic) l_lead,r_lead
        integer::ns
        complex*16,allocatable,dimension(:,:)::cl,scl,hs,ss,cr,scr
        end module
!*********************************************************************!
        module invert_m
        contains
        function invert(a,n) result(b)
        implicit none
        integer::n,info        
        integer,dimension(n)::ipiv
        complex*16,dimension(n,n)::a,b
        complex*16,dimension(n*n)::work_inv
        intent(in)a,n
        b=a
        call zgetrf(n,n,b,n,ipiv,info)
        call zgetri(n,b,n,ipiv,work_inv,n*n,info)
        end function invert
        end module invert_m
!*********************************************************************!
