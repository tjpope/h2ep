!*********************************************************************!
!*******  getscatteringmatrix  ***************************************!
!*********************************************************************!
        subroutine getscatteringmatrix(l_leadp,r_leadp,green,trans)
        use my_types, only: physics
        implicit none
        type(physics) l_leadp, r_leadp
        integer::i,j,is,info
        double precision,dimension(l_leadp%nc(1))::lgvl,lgvr
        double precision,dimension(r_leadp%nc(1))::rgvl,rgvr
        double precision::trans
        double precision,allocatable,dimension(:)::eig,rwork
        complex*16,dimension(max(l_leadp%n0,r_leadp%n0))::v_o,v_i
        complex*16,allocatable,dimension(:)::work
        complex*16,dimension(l_leadp%n0+r_leadp%n0,l_leadp%n0+r_leadp%n0)::middle,green,sgf
        complex*16,dimension(l_leadp%nc(1)+r_leadp%nc(1),l_leadp%nc(1)+r_leadp%nc(1))::s
        complex*16,dimension(l_leadp%nc(1),r_leadp%nc(1))::tl
!---------------------------------------------------------------------!
!       find the sqrt of the abs(group velcoities) for open channels  !
!---------------------------------------------------------------------!
        do i=1,l_leadp%nc(1); lgvl(i)=dsqrt(cdabs(l_leadp%gv(l_leadp%oo(i)))); enddo
        do i=1,l_leadp%nc(1); lgvr(i)=dsqrt(cdabs(l_leadp%gv(l_leadp%n0+l_leadp%oi(i)))); enddo
        do i=1,r_leadp%nc(1); rgvl(i)=dsqrt(cdabs(r_leadp%gv(r_leadp%oo(i)))); enddo
        do i=1,r_leadp%nc(1); rgvr(i)=dsqrt(cdabs(r_leadp%gv(r_leadp%n0+r_leadp%oi(i)))); enddo
        trans=0; s=(0,0)
!---------------------------------------------------------------------!
!       transmission from left lead                                   !
!                                                                     !
!       the reflection is found trivally, so here it's ignored        !
!       also, the transmision from the right is probably the same, so !
!       it's also ignored                                             !
!---------------------------------------------------------------------!
        tl=(0,0)
        sgf(:l_leadp%n0,:r_leadp%n0)=green(:l_leadp%n0,l_leadp%n0+1:)
        do i=1,l_leadp%nc(1); do j=1,r_leadp%nc(1)
         v_o(:l_leadp%n0)=l_leadp%duals(l_leadp%oo(i),:l_leadp%n0)*lgvl(i)
         v_i(:r_leadp%n0)=r_leadp%vectors(:r_leadp%n0,r_leadp%n0+r_leadp%oi(j))/rgvr(j)
         middle(:l_leadp%n0,:r_leadp%n0)=matmul(sgf(:l_leadp%n0,:r_leadp%n0),r_leadp%v_op(:r_leadp%n0,:r_leadp%n0))
         do is=1,r_leadp%n0; s(i,j)=s(i,j)-sum(v_o(:l_leadp%n0)*middle(:l_leadp%n0,is)*v_i(is)); enddo
        enddo; enddo
        middle(:l_leadp%nc(1),:l_leadp%nc(1))=matmul(s(:l_leadp%nc(1),:r_leadp%nc(1)),transpose(conjg(s(:l_leadp%nc(1),:r_leadp%nc(1)))))
        allocate(work(2*l_leadp%nc(1)),rwork(3*l_leadp%nc(1)),eig(l_leadp%nc(1)))
        call zheev('N','U',l_leadp%nc(1),middle(1:l_leadp%nc(1),1:l_leadp%nc(1)),l_leadp%nc(1),eig,work,2*l_leadp%nc(1),rwork,info)
        trans=sum(eig)
        deallocate(work,rwork,eig)
        end subroutine getscatteringmatrix
!*********************************************************************!
