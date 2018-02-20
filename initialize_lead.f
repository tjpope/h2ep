!*********************************************************************!
!*******  initialize_lead  *******************************************!
!*********************************************************************!
        subroutine initialize_lead(lead,leade,e)
        use my_types, only: periodic
        implicit none
        type(periodic) lead,leade
        double precision::e
        leade%h0=lead%h0-e*lead%s0; leade%s0=lead%s0
        leade%h1=lead%h1-e*lead%s1; leade%s1=lead%s1
        leade%h0i=lead%h0i-e*lead%s0i
        call desingularize(lead%n0,leade%h1)
        call make_k(lead%n0,leade%h0,leade%h0i,leade%h1,leade%k)
        end subroutine initialize_lead

