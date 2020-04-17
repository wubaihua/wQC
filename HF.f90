! This file is a part of wQC(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! Last update: 2020-4-17

!   HF: the part about Hartree-Fock method, including Restricted HF.


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

    
subroutine RHF(idout,nbas,nele,nucp,S,T,V,eri,D,E,C)
    use math
    implicit real*8(a-h,o-z)
    integer nbas,nele,i,j,m,idout
    real*8 S(nbas,nbas),T(nbas,nbas),V(nbas,nbas),eri(nbas,nbas,nbas,nbas)
    real*8 H_core(nbas,nbas),E(nbas),C(nbas,nbas),D(nbas,nbas),Di(nbas,nbas)
    real*8 S_haf(nbas,nbas),Fock(nbas,nbas),Fock_orth(nbas,nbas)
    real*8 E_ele,E_tot,nucp,E_toti,E_elei
    logical conv

    write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(idout,"(a)") "Restricted Hartree-Fock for closed-shell molecule:"
    
    call mat_power(nbas,S,-0.5_8,S_haf)
    H_core=T+V
    
    Fock=matmul(matmul(transpose(S_haf),H_core),S_haf)

    ! write(idout,*) "Fock="
    ! do i=1,nbas
    !     write(idout,*) Fock(i,:)
    ! end do

    
    call dia_symmat(nbas,Fock,E,C)

    ! write(idout,*) "C="
    ! do i=1,nbas
    !     write(idout,*) C(i,:)
    ! end do

    ! write(idout,*) "E="
    ! do i=1,nbas
    !     write(idout,*) E(i)
    ! end do
    
    D=0
    do i=1,nbas
        do j=1,nbas
            do m=1,nele/2
                D(i,j)=D(i,j)+C(i,m)*C(j,m)
            end do
        end do
    end do

    

    ! write(idout,*) "D="
    ! do i=1,nbas
    !     write(idout,*) D(i,:)
    ! end do
    
    E_ele=0
    do i=1,nbas
        do j=1,nbas
            E_ele=E_ele+D(j,i)*(H_core(i,j)+Fock(i,j))
        end do
    end do
    
    E_tot=E_ele+nucp

    ! write(idout,*) "E_ele=",E_ele
    ! write(idout,*) "E_tot=",E_tot
    
    
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

        ! write(idout,*) "Fock="
        ! do i=1,nbas
        !     write(idout,*) Fock(i,:)
        ! end do

        Fock_orth=matmul(matmul(transpose(S_haf),Fock),S_haf)

        ! write(*,*) "Fock="
        ! do i=1,nbas
        !     write(*,*) Fock(i,:)
        ! end do
        
        call dia_symmat(nbas,Fock_orth,E,C)
        
        C=matmul(S_haf,C)

        Di=0
        do i=1,nbas
            do j=1,nbas
                do m=1,nele/2
                    Di(i,j)=Di(i,j)+C(i,m)*C(j,m)
                end do
            end do
        end do

        ! do i=1,nbas
        !     do j=1,nbas
        !         if(Di(i,j)<1E-10)then
        !             Di(i,j)=0
        !         end if
        !     end do
        ! end do

        ! write(idout,*) "Di="
        ! do i=1,nbas
        !     write(idout,*) Di(i,:)
        ! end do

    
        E_elei=0
        do i=1,nbas
            do j=1,nbas
                E_elei=E_elei+Di(j,i)*(H_core(i,j)+Fock(i,j))
            end do
        end do

        ! E_elei=0
        ! do i=1,nele/2
            
        !     E_elei=E_elei+H_core(i,i)+E(i)
            
        ! end do
    
        E_toti=E_elei+nucp

        
        deltaE=E_toti-E_tot
        
        
        rmsd=0
        do i=1,nbas
            do j=1,nbas
                rmsd=rmsd+(Di(i,j)-D(i,j))**2
            end do
        end do
        rmsd=rmsd**0.5
        
        write(idout,"(a)") "//////////////////////////////////////////"
        write(idout,*) "SCF step=",icyc
        write(idout,*) "E_tot=",E_toti,"a.u."
        write(idout,*) "delta_E=",deltaE,"a.u."
        write(idout,*) "RMSD=",rmsd
        write(idout,"(a)") "//////////////////////////////////////////"

        E_tot=E_toti
        E_ele=E_elei
        D=Di
        
       
        
        
        icyc=icyc+1
        
        
        if(deltaE<1.0E-6 .and. rmsd<1.0E-8)then
            conv=.true.
            exit
        end if

    end do
    
    write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    if(.not. conv)then
        write(idout,*) "ERROR!!!!!!!!!!!!!!"
        write(idout,*) "SCF fails to convergeat at",128,"step."
        write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        return
    end if

    write(idout,*) "SCF convergence at",icyc,"step."
    write(idout,*) "The Electronic Energy is"
    write(idout,*) "E_ele=",E_ele,"a.u."
    write(idout,*) 
    write(idout,*) "The Total Energy is"
    write(idout,*) "E_tot=",E_tot,"a.u."
    write(idout,*) 
    write(idout,*) "The RHF orbital eigenvalues are (a.u.)"

    do i=1,nbas
        if(i<=nele/2)then
            write(idout,*) "occ.",E(i)
        else
            write(idout,*) "virt.",E(i)
        end if
    end do

    write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"




end subroutine







