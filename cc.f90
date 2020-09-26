subroutine ccsd
    use init
    use math
    implicit real*8(a-h,o-z)
    real*8,allocatable :: t1(:,:),t2(:,:,:,:)
    real*8,allocatable :: F_ccsd(:,:),W_ccsd(:,:,:,:)
    real*8,allocatable :: tau(:,:,:,:),tau_tilde(:,:,:,:)
    real*8,allocatable :: D1(:,:),D2(:,:,:,:)
    real*8 E_ccsd
    
    call init_postHF
    call ao2mo
    call build_spinorb_inte

    allocate(t1(nspinorb,nspinorb))
    allocate(t2(nspinorb,nspinorb,nspinorb,nspinorb))
    allocate(F_ccsd(nspinorb,nspinorb))
    allocate(W_ccsd(nspinorb,nspinorb,nspinorb,nspinorb))
    allocate(tau(nspinorb,nspinorb,nspinorb,nspinorb))
    allocate(tau_tilde(nspinorb,nspinorb,nspinorb,nspinorb))
    allocate(D1(nspinorb,nspinorb))
    allocate(D2(nspinorb,nspinorb,nspinorb,nspinorb))


    allocate(spinfock(nspinorb,nspinorb))
    spinfock=0
    do i=1,nspinorb
        do j=1,nspinorb
            do m=1,nele
                spinfock(i,j)=spinfock(i,j)+spinorb_inte(i,m,j,m)-spinorb_inte(i,m,m,j)
            end do
            if(mod(i,2)==mod(j,2))then
                spinfock(i,j)=spinfock(i,j)+T(index_spatial2spin(i),index_spatial2spin(j))+V(index_spatial2spin(i),index_spatial2spin(j))
            end if
        end do
    end do

    t1=0
    t2=0
    do i=1,nele
        do j=1,nele
            do k=nele+1,nspinorb
                do l=nele+1,nspinorb
                   t2(i,j,k,l)=(spinorb_inte(i,j,k,l)-spinorb_inte(i,j,l,k))/(E(index_spatial2spin(i))+E(index_spatial2spin(j))-E(index_spatial2spin(k))-E(index_spatial2spin(l)))
                end do
            end do
        end do
    end do

    E_ccsd=0
    do i=1,nele
        do j=1,nele
            do k=nele+1,nspinorb
                do l=nele+1,nspinorb
                    E_ccsd=E_ccsd+0.25*(spinorb_inte(i,j,k,l)-spinorb_inte(i,j,l,k))*t2(i,j,k,l)
                end do
            end do
        end do
    end do
    write(*,*) E_ccsd

    tau=0
    tau_tilde=0
    do i=1,nele
        do j=1,nele
            do ia=nele+1,nspinorb
                do ib=nele+1,nspinorb
                    tau(i,j,ia,ib)=t2(i,j,ia,ib)+t1(i,ia)*t1(j,ib)-t1(i,ib)*t1(j,ia)
                    tau_tilde(i,j,ia,ib)=t2(i,j,ia,ib)+0.5*(t1(i,ia)*t1(j,ib)-t1(i,ib)*t1(j,ia))
                end do
            end do
        end do
    end do

    F_ccsd=0
    do ia=nele+1,nspinorb
        do ie=nele+1,nspinorb
            F_ccsd(ia,ie)=(1-Kronecker_delta(ia,ie))*spinfock(ia,ie)
            do m=1,nele
                F_ccsd(ia,ie)=F_ccsd(ia,ie)-0.5*spinfock(m,ie)*t1(m,ia)
                do jf=nele+1,nspinorb
                    F_ccsd(ia,ie)=F_ccsd(ia,ie)+t1(m,jf)*(spinorb_inte(m,ia,jf,ie)-spinorb_inte(m,ia,ie,jf))
                    do n=1,nele            
                        F_ccsd(ia,ie)=F_ccsd(ia,ie)-0.5*tau_tilde(m,n,ia,jf)*(spinorb_inte(m,n,ie,jf)-spinorb_inte(m,n,jf,ie))
                    end do
                end do
            end do
        end do
    end do
    do m=1,nele
        do i=1,nele
            F_ccsd(m,i)=(1-Kronecker_delta(m,i))*spinfock(m,i)
            do ie=nele+1,nspinorb
                F_ccsd(m,i)=F_ccsd(m,i)+0.5*spinfock(i,ie)*t1(m,ie)
                do n=1,nele
                    F_ccsd(m,i)=F_ccsd(m,i)+t1(n,ie)*(spinorb_inte(m,n,i,ie)-spinorb_inte(m,n,ie,i))
                    do jf=nele+1,nspinorb            
                        F_ccsd(m,i)=F_ccsd(m,i)+0.5*tau_tilde(i,n,ie,jf)*(spinorb_inte(m,n,ie,jf)-spinorb_inte(m,n,jf,ie))
                    end do
                end do
            end do
        end do
    end do
    do m=1,nele
        do ie=nele+1,nspinorb
            F_ccsd(m,ie)=spinfock(m,ie)
            do n=1,nele
                do jf=nele+1,nspinorb
                    F_ccsd(m,ie)=F_ccsd(m,ie)+t1(n,jf)*(spinorb_inte(m,n,ie,jf)-spinorb_inte(m,n,jf,ie))
                end do
            end do
        end do
    end do


    W_ccsd=0
    do m=1,nele
        do n=1,nele
            do i=1,nele
                do j=1,nele
                    W_ccsd(m,n,i,j)=spinorb_inte(m,n,i,j)-spinorb_inte(m,n,j,i)
                    do ie=nele+1,nspinorb
                        W_ccsd(m,n,i,j)=W_ccsd(m,n,i,j)+t1(j,ie)*(spinorb_inte(m,n,i,ie)-spinorb_inte(m,n,ie,i))-t1(i,ie)*(spinorb_inte(m,n,j,ie)-spinorb_inte(m,n,ie,j))
                        do jf=nele+1,nspinorb
                            W_ccsd(m,n,i,j)=W_ccsd(m,n,i,j)+0.25*tau(i,j,ie,jf)*(spinorb_inte(m,n,ie,jf)-spinorb_inte(m,n,jf,ie))
                        end do
                    end do
                end do
            end do
        end do
    end do
    do ia=nele+1,nspinorb
        do ib=nele+1,nspinorb
            do ie=nele+1,nspinorb
                do jf=nele+1,nspinorb
                    W_ccsd(ia,ib,ie,jf)=spinorb_inte(ia,ib,ie,jf)-spinorb_inte(ia,ib,jf,ie)
                    do m=1,nele
                        W_ccsd(ia,ib,ie,jf)=W_ccsd(ia,ib,ie,jf)+t1(m,ib)*(spinorb_inte(ia,m,ie,jf)-spinorb_inte(ia,m,jf,ie))-t1(m,ia)*(spinorb_inte(ib,m,ie,jf)-spinorb_inte(ib,m,jf,ie))
                        do n=1,nele
                            W_ccsd(ia,ib,ie,jf)=W_ccsd(ia,ib,ie,jf)+0.25*tau(m,n,ia,ib)*(spinorb_inte(m,n,ie,jf)-spinorb_inte(m,n,jf,ie))
                        end do
                    end do
                end do
            end do
        end do
    end do
    do m=1,nele
        do ib=nele+1,nspinorb
            do ie=nele+1,nspinorb
                do j=i,nele
                    W_ccsd(m,ib,ie,j)=spinorb_inte(m,ib,ie,j)-spinorb_inte(m,ib,j,ie)
                    do jf=nele+1,nspinorb
                        W_ccsd(m,ib,ie,j)=W_ccsd(m,ib,ie,j)+t1(j,jf)*(spinorb_inte(m,ib,ie,jf)-spinorb_inte(m,ib,jf,ie))
                    end do
                    do n=1,nele
                        W_ccsd(m,ib,ie,j)=W_ccsd(m,ib,ie,j)-t1(n,ib)*(spinorb_inte(m,n,ie,j)-spinorb_inte(m,n,j,ie))
                        do jf=nele+1,nspinorb
                            W_ccsd(m,ib,ie,j)=W_ccsd(m,ib,ie,j)-(0.5*t2(j,n,jf,ib)+t1(j,jf)*t1(n,ib))*(spinorb_inte(m,n,ie,jf)-spinorb_inte(m,n,jf,ie))
                        end do
                    end do
                end do
            end do
        end do
    end do

    D1=0
    D2=0
    do i=1,nele
        do ia=nele+1,nspinorb
            D1(i,ia)=spinfock(i,i)-spinfock(ia,ia)
            do j=1,nele
                do ib=nele+1,nspinorb
                    D2(i,j,ia,ib)=spinfock(i,i)+spinfock(j,j)-spinfock(ia,ia)-spinfock(ib,ib)
                end do
            end do
        end do
    end do



end subroutine