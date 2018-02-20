!*********************************************************************!
!*******  inverteffham  **********************************************!
!*********************************************************************!
        subroutine inverteffham(n0l,n0r,green,nscatt,h1l,h1r,hs,cl,cr,kl,kr,gf)
        use invert_m
        implicit none
!-------  input
        integer::n0reen,n0l,n0r,nl,nr,nscatt
        complex*16,dimension(n0l,n0l)::h1l
        complex*16,dimension(n0r,n0r)::h1r
        complex*16,dimension(n0l+n0l,n0l+n0l)::kl
        complex*16,dimension(n0r+n0r,n0r+n0r)::kr
        complex*16,dimension(nscatt,nscatt)::hs
        complex*16,dimension(n0l,nscatt)::cl
        complex*16,dimension(n0r,nscatt)::cr
!-------  build h_eff
        integer::neff
        complex*16,allocatable,dimension(:,:)::heff,coup
!-------  green's function
        complex*16,dimension(n0l+n0r,n0l+n0r)::green,gf
        nl=n0l+n0l; nr=n0r+n0r+nl; neff=nr+nscatt; n0reen=n0l+n0r
        allocate(heff(neff,neff),coup(n0reen,neff))
        heff=(0,0); coup=(0,0)
        heff(:nscatt,:nscatt)=hs
        heff(nscatt+1:nscatt+nl,nscatt+1:nscatt+nl)=kl
        heff(nscatt+nl+1:nscatt+nr,nscatt+nl+1:nscatt+nr)=kr
        heff(nscatt+1:nscatt+n0l,:nscatt)=cl
        heff(:nscatt,nscatt+1:nscatt+n0l)=transpose(conjg(cl))
        heff(nscatt+nl+1:nscatt+nl+n0r,:nscatt)=cr
        heff(:nscatt,nscatt+nl+1:nscatt+nl+n0r)=transpose(conjg(cr))
        coup(:n0l,nscatt+n0l+1:nscatt+nl)=h1l
        coup(n0l+1:,nscatt+nl+n0r+1:nscatt+nr)=h1r
        gf=invert(green-matmul(matmul(coup,invert(heff,neff)),transpose(conjg(coup))),n0reen)
        return
        end subroutine inverteffham
!*********************************************************************!
