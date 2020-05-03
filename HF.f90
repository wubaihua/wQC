! This file is a part of wQC(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! Last update: 2020-4-21

!   HF: the part about Hartree-Fock method, including Restricted HF,
! Restricted HF with DIIS speeding up and Unrestricted HF.


subroutine cal_nucp(idout,atom,natom,nucp)
    use def
    !use math
    use geo
    implicit none
    integer natom,i,j,idout
    type(atomtype) :: atom(natom)
    real*8 nucp,d
    
    write(idout,"(a)") "---------------------------------------------------"

    nucp=0
    do i=2,natom
        do j=1,i-1
            call geo_distance(atom,i,j,d)
            nucp=nucp+real(atom(i)%charge*atom(j)%charge)/d
        end do
    end do
    
    write(idout,*) "nuclear repulsion energy(a.u.)=",nucp

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

    write(idout,"(a)") "---------------------------------------------------"
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
        
        
        if(deltaE<1.0E-8 .and. rmsd<1.0E-8)then
            conv=.true.
            exit
        end if

    end do
    
    write(idout,"(a)") "---------------------------------------------------"

    if(.not. conv)then
        write(idout,*) "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(idout,*) "SCF fails to convergeat at",128,"steps."
        write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        call error_end(idout,"SCF ERROR")
    end if

    write(idout,*) "SCF convergence at",icyc-1,"steps."
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

    write(idout,"(a)") "---------------------------------------------------"




end subroutine





subroutine RHF_DIIS(idout,nbas,nele,nucp,S,T,V,eri,D,E,C,ndiis)
    use math
    implicit real*8(a-h,o-z)
    integer nbas,nele,i,j,m,idout,ndiis
    real*8 S(nbas,nbas),T(nbas,nbas),V(nbas,nbas),eri(nbas,nbas,nbas,nbas)
    real*8 H_core(nbas,nbas),E(nbas),C(nbas,nbas),D(nbas,nbas),Di(nbas,nbas)
    real*8 S_haf(nbas,nbas),Fock(nbas,nbas),Fock_orth(nbas,nbas)
    real*8 F_diis(ndiis,nbas,nbas),D_diis(ndiis,nbas,nbas),e_diis(ndiis,nbas,nbas)
    real*8 mat_diis(ndiis+1,ndiis+1),alpha_diis(ndiis+1),b_diis(ndiis+1)!,a(nbas,nbas)
    real*8 E_ele,E_tot,nucp,E_toti,E_elei
    logical conv

    write(idout,"(a)") "---------------------------------------------------"
    write(idout,"(a)") "Restricted Hartree-Fock for closed-shell molecule:"
    write(idout,*) "DIIS Optimizer used, size of DIIS space=",ndiis
    
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
        
        if(icyc<=ndiis)then
        
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

            F_diis(icyc,:,:)=Fock
            D_diis(icyc,:,:)=D
        else
            do i=1,ndiis
                e_diis(i,:,:)=matmul(matmul(F_diis(i,:,:),D_diis(i,:,:)),S)-matmul(S,matmul(D_diis(i,:,:),F_diis(i,:,:)))
            end do
            mat_diis=-1
            mat_diis(ndiis+1,ndiis+1)=0
            do i=1,ndiis
                do j=1,ndiis
                    call mat_dot_prod(nbas,e_diis(i,:,:),e_diis(j,:,:),mat_diis(i,j))
                end do
            end do
            b_diis=0
            b_diis(ndiis+1)=-1
            
            call solv_LES(ndiis+1,mat_diis,alpha_diis,b_diis)

            ! write(*,*) "sumalpha=:",sum(alpha_diis)
            ! a=0
            ! do i=1,ndiis
            !     a=a+alpha_diis(i)*e_diis(i,:,:)
            ! end do
            ! write(*,*) "a=",a
                
            Fock=0
            do i=1,ndiis
                Fock=Fock+alpha_diis(i)*F_diis(i,:,:)
            end do

            F_diis(1:ndiis-1,:,:)=F_diis(2:ndiis,:,:)
            F_diis(ndiis,:,:)=Fock

            D_diis(1:ndiis-1,:,:)=D_diis(2:ndiis,:,:)
            

        end if


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
        if(icyc>ndiis)D_diis(ndiis,:,:)=D
        
       
        
        
        icyc=icyc+1
        
        
        if(deltaE<1.0E-8 .and. rmsd<1.0E-8)then
            conv=.true.
            exit
        end if

    end do
    
    write(idout,"(a)") "---------------------------------------------------"

    if(.not. conv)then
        write(idout,*) "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(idout,*) "SCF fails to convergeat at",128,"steps."
        write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        call error_end(idout,"SCF ERROR")
    end if

    write(idout,*) "SCF convergence at",icyc-1,"steps."
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

    
    write(idout,"(a)") "---------------------------------------------------"




end subroutine



subroutine UHF(idout,nbas,nele_alpha,nele_beta,nucp,S,T,V,eri,D_alpha,D_beta,E_alpha,E_beta)
    use math
    implicit real*8(a-h,o-z)
    integer nbas,nele_alpha,nele_beta,i,j,m,idout,info
    real*8 S(nbas,nbas),T(nbas,nbas),V(nbas,nbas),eri(nbas,nbas,nbas,nbas)
    real*8 H_core(nbas,nbas),F_alpha(nbas,nbas),F_beta(nbas,nbas),E(nbas),C(nbas,nbas)
    real*8 E_alpha(nbas),E_beta(nbas),C_alpha(nbas,nbas),C_beta(nbas,nbas),D_alpha(nbas,nbas),D_beta(nbas,nbas),D_tot(nbas,nbas)
    real*8 S_haf(nbas,nbas),Fock(nbas,nbas),Fock_orth_alpha(nbas,nbas),Fock_orth_beta(nbas,nbas)
    real*8 Di_alpha(nbas,nbas),Di_beta(nbas,nbas),Di_tot(nbas,nbas)
    real*8 E_ele,E_tot,nucp,E_toti,E_elei,S_ab(nbas,nbas)
    logical conv


    write(idout,"(a)") "---------------------------------------------------"
    write(idout,"(a)") "Unrestricted Hartree-Fock for open-shell molecule:"
    !write(*,*) "nele_alpha1=",nele_alpha
    S_haf=0
    call mat_power(nbas,S,-0.5_8,S_haf)
    !write(*,*) "nbas=",nbas
    !write(*,*) S_haf
    H_core=T+V
    
    F_alpha=matmul(matmul(transpose(S_haf),H_core),S_haf)
    F_beta=matmul(matmul(transpose(S_haf),H_core),S_haf)
    ! write(*,*) "nele_alpha2=",nele_alpha
    ! write(*,*) "nele_beta2=",nele_beta
    ! write(*,*) "Fock="
    ! do i=1,nbas
    !     write(*,*) Fock(i,:)
    ! end do

    E_alpha(:)=0
    C_alpha(:,:)=0
    E_beta(:)=0
    C_beta(:,:)=0
    call dia_symmat(nbas,F_alpha,E_alpha,C_alpha)
    call dia_symmat(nbas,F_beta,E_beta,C_beta)
    ! C=Fock
    ! call dsyev('V','L',nbas,C,nbas,E,work,4*nbas,info)
    ! write(*,*) "nbas=",nbas
    ! write(*,*) E
    ! write(*,*) C
    ! E_alpha=E
    ! C_alpha=C
    ! E_beta=E_alpha
    ! C_beta=C_alpha
    ! write(*,*) E_alpha
    ! write(*,*) C_alpha
    ! write(*,*) E_beta
    ! write(*,*) C_beta
    ! write(*,*) "nele_alpha3=",nele_alpha
    ! write(*,*) "nele_beta3=",nele_beta
    ! write(idout,*) "C="
    ! do i=1,nbas
    !     write(idout,*) C(i,:)
    ! end do

    ! write(idout,*) "E="
    ! do i=1,nbas
    !     write(idout,*) E(i)
    ! end do
    ! write(*,*) "nbas=",nbas
    ! write(*,*) "nele_alpha=",nele_alpha
    D_alpha=0
    do i=1,nbas
        do j=1,nbas
            do m=1,nele_alpha
                D_alpha(i,j)=D_alpha(i,j)+C_alpha(i,m)*C_alpha(j,m)
            end do
        end do
    end do

    D_beta=0
    do i=1,nbas
        do j=1,nbas
            do m=1,nele_beta
                D_beta(i,j)=D_beta(i,j)+C_beta(i,m)*C_beta(j,m)
            end do
        end do
    end do

    D_tot=D_alpha+D_beta

    ! write(idout,*) "D="
    ! do i=1,nbas
    !     write(idout,*) D(i,:)
    ! end do
    
    E_ele=0
    do i=1,nbas
        do j=1,nbas
            E_ele=E_ele+0.5*(D_tot(i,j)*H_core(i,j)+F_alpha(i,j)*D_alpha(i,j)+F_beta(i,j)*D_beta(i,j))
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
                F_alpha(i,j)=H_core(i,j)
                F_beta(i,j)=H_core(i,j)
                do k=1,nbas
                    do l=1,nbas
                     F_alpha(i,j)=F_alpha(i,j)+D_tot(l,k)*eri(i,j,k,l)-D_alpha(l,k)*eri(i,l,k,j)
                     F_beta(i,j)=F_beta(i,j)+D_tot(l,k)*eri(i,j,k,l)-D_beta(l,k)*eri(i,l,k,j)
                    end do
                end do
            end do
        end do 

        ! write(idout,*) "Fock="
        ! do i=1,nbas
        !     write(idout,*) Fock(i,:)
        ! end do

        Fock_orth_alpha=matmul(matmul(transpose(S_haf),F_alpha),S_haf)
        Fock_orth_beta=matmul(matmul(transpose(S_haf),F_beta),S_haf)

        ! write(*,*) "Fock="
        ! do i=1,nbas
        !     write(*,*) Fock(i,:)
        ! end do
        
        call dia_symmat(nbas,Fock_orth_alpha,E_alpha,C_alpha)
        call dia_symmat(nbas,Fock_orth_beta,E_beta,C_beta)
        
        C_alpha=matmul(S_haf,C_alpha)
        C_beta=matmul(S_haf,C_beta)


        Di_alpha=0
        Di_beta=0
        do i=1,nbas
            do j=1,nbas
                do m=1,nele_alpha
                    Di_alpha(i,j)=Di_alpha(i,j)+C_alpha(i,m)*C_alpha(j,m)
                end do
                do m=1,nele_beta
                    Di_beta(i,j)=Di_beta(i,j)+C_beta(i,m)*C_beta(j,m)
                end do
            end do
        end do

        Di_tot=Di_alpha+Di_beta

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
                E_elei=E_elei+0.5*(D_tot(i,j)*H_core(i,j)+F_alpha(i,j)*D_alpha(i,j)+F_beta(i,j)*D_beta(i,j))
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
                rmsd=rmsd+(Di_tot(i,j)-D_tot(i,j))**2
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
        D_tot=Di_tot
        D_alpha=Di_alpha
        D_beta=Di_beta
        
       
        
        
        icyc=icyc+1
        
        
        if(deltaE<1.0E-8 .and. rmsd<1.0E-8)then
            conv=.true.
            exit
        end if

    end do
    
    write(idout,"(a)") "---------------------------------------------------"

    if(.not. conv)then
        write(idout,*) "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(idout,*) "SCF fails to convergeat at",128,"steps."
        write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        call error_end(idout,"SCF ERROR")
    end if

    write(idout,*) "SCF convergence at",icyc-1,"steps."
    write(idout,*) "The Electronic Energy is"
    write(idout,*) "E_ele=",E_ele,"a.u."
    write(idout,*) 
    write(idout,*) "The Total Energy is"
    write(idout,*) "E_tot=",E_tot,"a.u."
    write(idout,*) 
    write(idout,*) "The UHF alpha orbital eigenvalues are (a.u.)"

    do i=1,nbas
        if(i<=nele_alpha)then
            write(idout,*) "occ.",E_alpha(i)
        else
            write(idout,*) "virt.",E_alpha(i)
        end if
    end do

    write(idout,*) "The UHF beta orbital eigenvalues are (a.u.)"

    do i=1,nbas
        if(i<=nele_beta)then
            write(idout,*) "occ.",E_beta(i)
        else
            write(idout,*) "virt.",E_beta(i)
        end if
    end do

    write(idout,*)
    
    write(idout,"(a)") "Spin Calculation:"
    S_ab=matmul(matmul(transpose(C_alpha),S),C_beta)

    spin_contamination=nele_beta
    do i=1,nele_alpha
        do j=1,nele_beta
            spin_contamination=spin_contamination-(S_ab(i,j))**2
        end do
    end do
    write(idout,*) "spin contamination=",spin_contamination
    write(idout,*) "<S^2>=",spin_contamination+0.5*(nele_alpha-nele_beta)*(0.5*(nele_alpha-nele_beta)+1)



    write(idout,"(a)") "---------------------------------------------------"


end subroutine



