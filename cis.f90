subroutine cis
    use init
    use math
    implicit real*8(a-h,o-z)


    call init_postHF
    call ao2mo
    call build_spinorb_inte

    n_single_exci=nele*(nspinorb-nele)
    allocate(cis_ham(n_single_exci,n_single_exci))
    allocate(cis_c(n_single_exci,n_single_exci))
    allocate(cis_e(n_single_exci))

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
    ! write(*,*) spinfock

    index1=0
    index2=0
    do i=1,nele
        do ia=nele+1,nspinorb
            index1=index1+1
            do j=1,nele
                do ib=nele+1,nspinorb
                    index2=index2+1
                    cis_ham(index1,index2)=spinfock(ia,ib)*Kronecker_delta(i,j)-spinfock(i,j)*Kronecker_delta(ia,ib)+spinorb_inte(ia,j,i,ib)-spinorb_inte(ia,j,ib,i)
                end do
            end do
            index2=0
        end do
    end do

    ! write(*,*) cis_ham

    call dia_symmat( n_single_exci,cis_ham,cis_e,cis_c)

    write(*,*) cis_e*27.21138602







end subroutine