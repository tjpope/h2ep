!*********************************************************************!
!*******  cutuphamiltonian  ******************************************!
!*********************************************************************!
        subroutine cutuphamiltonian
        use commandments
        use hamiltonian
        use h0h1
        implicit none
        integer::i,j,n0minl,n0maxl,n0minr,n0maxr
!---------------------------------------------------------------------!
!            |     H0    H1    ..   |                                 !
! left  H =  |   H1^{d}  HI    ..   |                                 !
!            |     ..    ..    ..   |                                 !
!---------------------------------------------------------------------!
        n0minl=1; do i=1,n0min-1; n0minl=n0minl+hatoms(i,2); enddo
        n0maxl=0; do i=1,n0max; n0maxl=n0maxl+hatoms(i,2); enddo
        l_lead%n0=n0maxl+1-n0minl
        if(l_lead%n0.ge.int(n*0.5)) then; write(*,100) l_lead%n0, n; stop; endif
        allocate(l_lead%h0(l_lead%n0,l_lead%n0),l_lead%h1(l_lead%n0,l_lead%n0),l_lead%s0(l_lead%n0,l_lead%n0),l_lead%s1(l_lead%n0,l_lead%n0),l_lead%h0i(l_lead%n0,l_lead%n0),l_lead%s0i(l_lead%n0,l_lead%n0))
        l_lead%h0=h(n0minl:n0maxl,n0minl:n0maxl); l_lead%s0=s(n0minl:n0maxl,n0minl:n0maxl)
        l_lead%h1=h(n0minl:n0maxl,n0minl+l_lead%n0:n0maxl+l_lead%n0); l_lead%s1=s(n0minl:n0maxl,n0minl+l_lead%n0:n0maxl+l_lead%n0)
        l_lead%h0i=h(n0minl+l_lead%n0:n0maxl+l_lead%n0,n0minl+l_lead%n0:n0maxl+l_lead%n0); l_lead%s0i=s(n0minl+l_lead%n0:n0maxl+l_lead%n0,n0minl+l_lead%n0:n0maxl+l_lead%n0)
!---------------------------------------------------------------------!
!            |     ..    ..    ..   |                                 !
! right H =  |     ..    HI    H1   |                                 !
!            |     ..  H1^{d}  H0   |                                 !
!---------------------------------------------------------------------!
        n0maxr=1; do i=1,n0min-1; n0maxr=n0maxr+hatoms(i,2); enddo
        n0minr=0; do i=1,n0max; n0minr=n0minr+hatoms(i,2); enddo
        n0minr=n+1-n0minr; n0maxr=n+1-n0maxr; r_lead%n0=n0maxr+1-n0minr
        if(r_lead%n0.ge.int(n*0.5)) then; write(*,100) r_lead%n0, n; stop; endif
        allocate(r_lead%h0(r_lead%n0,r_lead%n0),r_lead%h1(r_lead%n0,r_lead%n0),r_lead%s0(r_lead%n0,r_lead%n0),r_lead%s1(r_lead%n0,r_lead%n0),r_lead%h0i(r_lead%n0,r_lead%n0),r_lead%s0i(r_lead%n0,r_lead%n0))
        r_lead%h0=h(n0minr:n0maxr,n0minr:n0maxr); r_lead%s0=s(n0minr:n0maxr,n0minr:n0maxr)
        r_lead%h1=h(n0minr:n0maxr,n0minr-r_lead%n0:n0maxr-r_lead%n0); r_lead%s1=s(n0minr:n0maxr,n0minr-r_lead%n0:n0maxr-r_lead%n0)
        r_lead%h0i=h(n0minr-r_lead%n0:n0maxr-r_lead%n0,n0minr-r_lead%n0:n0maxr-r_lead%n0); r_lead%s0i=s(n0minr-r_lead%n0:n0maxr-r_lead%n0,n0minr-r_lead%n0:n0maxr-r_lead%n0)
!---------------------------------------------------------------------!
!            |     ..    cl    ..   |                                 !
! scat  H =  |     ..    HS    ..   |                                 !
!            |     ..    cr    ..   |                                 !
!---------------------------------------------------------------------!
        ns=n-(l_lead%n0*2+r_lead%n0*2+n0minl*2)+2; i=n0maxl+1+l_lead%n0; j=n0maxl+ns+l_lead%n0
        allocate(hs(ns,ns),ss(ns,ns),cl(l_lead%n0,ns),scl(l_lead%n0,ns),cr(r_lead%n0,ns),scr(r_lead%n0,ns))
        hs=h(i:j,i:j); ss=s(i:j,i:j)
        cl=h(n0minl+l_lead%n0:n0maxl+l_lead%n0,i:j); scl=s(n0minl+l_lead%n0:n0maxl+l_lead%n0,i:j)
        cr=h(n0minr-r_lead%n0:n0maxr-r_lead%n0,i:j); scr=s(n0minr-r_lead%n0:n0maxr-r_lead%n0,i:j)
!---------------------------------------------------------------------!
        deallocate(h,s,hatoms)
        return
100        format("h0 and h1 dimensions are too large:",i0,/,"size of input hamiltonian:",i0)
        end subroutine cutuphamiltonian
!*********************************************************************!
