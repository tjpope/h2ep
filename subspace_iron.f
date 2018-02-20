!*********************************************************************!
!*******  subspace_iron  *********************************************!
!*********************************************************************!
        subroutine subspace_iron(n,nc,chan,indx,vr,eig,ieig,k,h1)
        implicit none
!-------  in
        integer::n,n2,nc
        integer,dimension(n)::chan,indx
        complex*16,dimension(2*n)::eig,ieig,k
        complex*16,dimension(n,n)::h1
        complex*16,dimension(2*n,2*n)::vr
!-------  calc
        integer::nxs,nss,mss,i,j,is,js,jm,info
        integer,dimension(nc)::iss,jss
        double precision::zero
        double precision,allocatable,dimension(:)::rwrk,ws
        complex*16::im
        complex*16,dimension(n)::auxv
        complex*16,dimension(2*n,2*n)::vl
        complex*16,allocatable,dimension(:)::wrk
        complex*16,allocatable,dimension(:,:)::curr,vs
        parameter(im=(0,1),zero=1d-8); n2=2*n
        jss=chan(indx(1:nc)); nxs=nc; nss=1
        do while(nss.ne.0)
          mss=0; nss=0
          do i=1,nxs
            j=jss(i)
            if(abs(k(j)-k(jss(1))).lt.zero) then
              mss=mss+1; iss(mss)=j
            else
              nss=nss+1; jss(nss)=j
            endif
          enddo
          nxs=nss
          if(mss.gt.1)then
            allocate(curr(mss,mss),ws(mss),vs(mss,mss),wrk(2*mss),rwrk(2*mss))
            curr=0.d0
            do is=1,mss; do js=1,mss; do j=1,n; do jm=1,n
              curr(js,is)=curr(js,is)+im*conjg(vr(j,iss(js)))*(h1(j,jm)*eig(iss(1))-conjg(h1(jm,j))*ieig(iss(1)))*vr(jm,iss(is))
            enddo; enddo; enddo; enddo
            call zgeev('N','V',mss,curr,mss,ws,vl,n2,vs,mss,wrk,2*mss,rwrk,info)
            do is=1,mss
              auxv=0.d0
              do js=1,mss; do j=1,n2
                auxv(j)=auxv(j)+vr(j,jss(js))*vs(js,is)
              enddo; enddo
              vr(:,jss(is))=auxv
            enddo
            deallocate(rwrk,wrk,vs,ws,curr)
          endif
        enddo
        return
        end subroutine subspace_iron
!*********************************************************************!
