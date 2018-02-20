!*********************************************************************!
!*******  getgreens  *************************************************!
!*********************************************************************!
        subroutine getgreens(e0,n0,h0,h1,s0,s1,g,v_op,gv,nc,lindex,rindex,vectors,duals,left,verbose,nochannels)
        use invert_m
        implicit none
!-------  input 
        integer::n0,n2
        double precision::e0
        complex*16,dimension(n0,n0)::h0,h1,s0,s1
        logical::left,verbose,nochannels
!-------  adjust for overlap matirx
        integer::i,j
!-------  load a and b matrices
        complex*16,dimension(2*n0,2*n0)::a,b
!-------  solve eigenproblem
        integer::info
        double precision,dimension(16*n0)::rwork
        complex*16::vl
        complex*16,dimension(2*n0)::alpha,beta
        complex*16,dimension(4*n0)::work
        complex*16,dimension(2*n0,2*n0)::eigenvector
!-------  calculate eigenvalues
        double precision::pi
        complex*16::im        
        complex*16,dimension(2*n0)::eig,ieig,k
        parameter(im=(0,1),pi=3.14159265359)
!-------  calculate group velocity
        integer::l
        complex*16::norm
        complex*16,dimension(2*n0)::vg,gv
!-------  number of open channels
        integer::nchan
        integer,dimension(2)::nchanl,nchanr,nc
        double precision::zero
        parameter(zero=1d-7)
!-------  storing the eigenvectors
        integer,dimension(n0)::lchan,rchan,lindex,rindex
        complex*16,dimension(n0,n0)::ve
!-------  dual vectors
        complex*16,dimension(n0,2*n0)::vectors,chi
        complex*16,dimension(2*n0,n0)::duals
!-------  greens functions
        complex*16,dimension(n0,n0)::v_op,g
!---------------------------------------------------------------------!
!      | H_{0}-E  -H_{1}^{d} |      | H_{1}   0 |                     !
!  a = |                     |, b = |           |                     !
!      |    I          0     |      |   0     I |                     !
!---------------------------------------------------------------------!
        n2=n0*2; a=(0,0); b=(0,0); do i=n0+1,n2; b(i,i)=(1,0); enddo
        b(:n0,:n0)=h1; a(n0+1:,:n0)=b(n0+1:,n0+1:)
        a(:n0,:n0)=-h0; a(:n0,n0+1:)=-transpose(conjg(h1))
!---------------------------------------------------------------------!
! | E-H_{0}  -H_{1}^{d} |   | vec1 |         | H_{1}   0 |   | vec1 | !
! |                     | x |      | = eig x |           | x |      | !
! |    I          0     |   | vec2 |         |   0     I |   | vec2 | !
!---------------------------------------------------------------------!
        info=0; eigenvector=(0,0); alpha=(0,0); beta=(0,0)
        call zggev('N','V',n2,a,n2,b,n2,alpha,beta,vl,1,eigenvector,n2,work,2*n2,rwork,info)
!---------------------------------------------------------------------!
!   eig = alpha / beta                                                !
!---------------------------------------------------------------------!
        eig=alpha/beta; ieig=1.d0/eig; k=-1.d0*im*log(eig)/pi
!---------------------------------------------------------------------!
!   group velocity operator: Vg = i ( H_{1} * z - H_{1}^{d} / z )     !
!   group velocity: vg = <vector|Vg|vector>                           !
!   note: z=1 for an open channel                                     !
!---------------------------------------------------------------------!
        vg=(0,0);vectors=(0,0);chi=(0,0);ve=(0,0);v_op=(0,0);g=(0,0)
        do i=1,n2
         norm=0.d0
         do j=1,n0; do l=1,n0; norm=norm+conjg(eigenvector(j,i))*(s0(j,l)+s1(j,l)*eig(i)+conjg(s1(l,j))*ieig(i))*eigenvector(l,i); enddo; enddo
         eigenvector(:n0,i)=eigenvector(:n0,i)/sqrt(norm)
         do j=1,n0; do l=1,n0; vg(i)=vg(i)+im*dconjg(eigenvector(j,i))*(h1(j,l)*eig(i)-dconjg(h1(l,j))*ieig(i))*eigenvector(l,i); enddo; enddo
        enddo
!---------------------------------------------------------------------!
!  count the open channels                                            !
!                                                                     !
! open if |eig|=1 (left if vg<0, right if vg>0)                       !
! closed left if |eig|<1, closed right if |eig|>1                     !
!---------------------------------------------------------------------!
        nchanl=0;nchanr=0;nchan=0;lchan=0;rchan=0;lindex=0;rindex=0 
        do i=1,n2
         if(abs(abs(eig(i))-1.d0).le.zero) then
          nchan=nchan+1
          if(dble(vg(i)).lt.0.0d0) then
           nchanl(1)=nchanl(1)+1
           lindex(nchanl(1))=nchanl(1)+nchanl(2)
           lchan(nchanl(1)+nchanl(2))=i
          elseif(dble(vg(i)).gt.0.0d0) then
           nchanr(1)=nchanr(1)+1
           rindex(nchanr(1))=nchanr(1)+nchanr(2)
           rchan(nchanr(1)+nchanr(2))=i
          endif
         elseif(abs(eig(i)).gt.1.0d0) then
          nchanl(2)=nchanl(2)+1
          lchan(nchanl(2)+nchanl(1))=i
         elseif(abs(eig(i)).lt.1.0d0) then
          nchanr(2)=nchanr(2)+1
          rchan(nchanr(1)+nchanr(2))=i
         endif
        enddo
        nc(1)=nchanl(1); nc(2)=nchanr(1)
        if(minval(nc).lt.1) then
         nochannels=.true.; return
        endif
!---------------------------------------------------------------------!
!  collect the group velocities                                       !
!---------------------------------------------------------------------!
        gv=(0,0); gv(:n0)=vg(lchan); gv(n0+1:)=vg(rchan)
!---------------------------------------------------------------------!
!  Degenerate subspace ironing                                        !
!---------------------------------------------------------------------!
        call subspace_iron(n0,nc(1),lchan,lindex,eigenvector,eig,ieig,k,h1)
        call subspace_iron(n0,nc(2),rchan,rindex,eigenvector,eig,ieig,k,h1)
        do i=1,n2; eigenvector(:n0,i)=eigenvector(:n0,i)*abs(eigenvector(1,i))/eigenvector(1,i); enddo
!---------------------------------------------------------------------!
!  collect the vectors                                                !
!---------------------------------------------------------------------!
        vectors(:,:n0)=eigenvector(:n0,lchan); duals(:n0,:)=invert(vectors(:,:n0),n0)
        vectors(:,n0+1:)=eigenvector(:n0,rchan); duals(n0+1:,:)=invert(vectors(:,n0+1:),n0)
!---------------------------------------------------------------------!
!      xl+=<kL|eig|kL_bar>,  xr-=<kL|eig^-1|kL_bar>                   !
!    xl-=<kL|eig^-1|kL_bar>,   xr+=<kR|eig|kR_bar>                    ! 
!---------------------------------------------------------------------!
        do i=1,n0; ve(:,i)=vectors(:,i)*ieig(lchan(i)); enddo; chi(:,:n0)=matmul(ve,duals(:n0,:))
        do i=1,n0; ve(:,i)=vectors(:,n0+i)*ieig(rchan(i)); enddo; chi(:,n0+1:)=matmul(ve,duals(n0+1:,:))
!---------------------------------------------------------------------!
!  v_op=H1(<kL|eig^-1|kL_bar>-<kR|eig^-1|kR_bar>)                     !
!---------------------------------------------------------------------!
        v_op=matmul(transpose(conjg(h1)),(chi(:,n0+1:)-chi(:,:n0)))
!---------------------------------------------------------------------!
!  g^{-1} = h0 + h1 * <kL|eig^-1|kL_bar>                              !
!---------------------------------------------------------------------!
        g=h0+matmul(transpose(conjg(h1)),chi(:,:n0))
        return
100        format(t5,es20.5,5x,5i10)
        end subroutine getgreens
!*********************************************************************!
