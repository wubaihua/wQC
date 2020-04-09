subroutine cal_nucp(atom,natom,nucp)
    use def
    !use math
    use geo
    implicit none
    integer natom,i,j
    type(atomtype) :: atom(natom)
    real*8 nucp,d
    
    nucp=0
    do i=2,natom
        do j=1,i-1
            call geo_distance(atom,i,j,d)
            nucp=nucp+real(atom(i)%charge*atom(j)%charge)/d
        end do
    end do
    
end subroutine

    
subroutine RHF(nbas,nele,nucp,S,T,V,eri,D,E,C)
    use math
    implicit real*8(a-h,o-z)
    integer nbas,nele,i,j,m
    real*8 S(nbas,nbas),T(nbas,nbas),V(nbas,nbas),eri(nbas,nbas,nbas,nbas)
    real*8 H_core(nbas,nbas),E(nbas),C(nbas,nbas),D(nbas,nbas),Di(nbas,nbas)
    real*8 S_haf(nbas,nbas),Fock(nbas,nbas)
    real*8 E_ele,E_tot,nucp,E_toti,E_elei
    logical conv
    
    call mat_power(nbas,S,-0.5_8,S_haf)
    H_core=T+V
    
    Fock=matmul(matmul(transpose(S_haf),H_core),S_haf)

    write(*,*) "Fock="
    do i=1,nbas
        write(*,*) Fock(i,:)
    end do

    
    call dia_symmat(nbas,Fock,E,C)

    write(*,*) "C="
    do i=1,nbas
        write(*,*) C(i,:)
    end do

    write(*,*) "E="
    do i=1,nbas
        write(*,*) E(i)
    end do
    
    D=0
    do i=1,nbas
        do j=1,nbas
            do m=1,nele/2
                D(i,j)=D(i,j)+C(i,m)*C(j,m)
            end do
        end do
    end do


    write(*,*) "D="
    do i=1,nbas
        write(*,*) D(i,:)
    end do
    
    E_ele=0
    do i=1,nbas
        do j=1,nbas
            E_ele=E_ele+D(i,j)*(H_core(i,j)+Fock(i,j))
        end do
    end do
    
    E_tot=E_ele+nucp

    write(*,*) "E_ele=",E_ele
    write(*,*) "E_tot=",E_tot
    
    
    icyc=1
    conv=.false.
    
    do while(icyc<128)
        
        do i=1,nbas
            do j=1,nbas 
                Fock(i,j)=H_core(i,j)
                do k=1,nbas
                    do l=1,nbas
                     Fock(i,j)=Fock(i,j)+D(k,l)*(2.0*eri(i,j,k,l)-eri(i,k,j,l))
                    end do
                end do
            end do
        end do 
        
        call dia_symmat(nbas,Fock,E,C)
        
        C=matmul(S_haf,C)

        Di=0
        do i=1,nbas
            do j=1,nbas
                do m=1,nele/2
                    Di(i,j)=Di(i,j)+C(i,m)*C(j,m)
                end do
            end do
        end do
    
        E_elei=0
        do i=1,nbas
            do j=1,nbas
                E_elei=E_elei+Di(i,j)*(H_core(i,j)+Fock(i,j))
            end do
        end do
    
        E_toti=E_elei+nucp

        
        deltaE=E_toti-E_tot
        
        
        rmsd=0
        do i=1,nbas
            do j=1,nbas
                rmsd=rmsd+(Di(i,j)-D(i,j))**2
            end do
        end do
        rmsd=rmsd**0.5
        
        write(*,*) "step=",icyc
        write(*,*) "E_tot=",E_toti
        write(*,*) "delta_E=",deltaE
        write(*,*) "RMSD=",rmsd
        
        E_tot=E_toti
        E_ele=E_elei
        D=Di
        
       
        
        
        icyc=icyc+1
        
        
        if(deltaE<1.0E-6 .and. rmsd<1.0E-8)then
            conv=.true.
            exit
        end if

    end do
    
    if(.not. conv)then
        write(*,*) "ERROR!!!!!!!!!!!!!!"
        write(*,*) "SCF fails to convergeat at",128,"step."
        return
    end if

    write(*,*) "SCF convergence at",icyc,"step."
    write(*,*) "The RHF orbital eigenvalues are"

    do i=1,nbas
        if(i<=nele/2)then
            write(*,*) "occ.",E(i)
        else
            write(*,*) "virt.",E(i)
        end if
    end do






end subroutine