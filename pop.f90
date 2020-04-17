! wQC: a simple fortran Quantum Chemistry/Electronic Structure Program
!(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! Last update: 2020-4-17

! pop: the population analysis part.


subroutine pop_analy(idout,atom,nshl,nbas,natom,D_alpha,D_beta,S,shl_belong_to_atom,angl,MLK_charge,LDW_charge)
    use def
    use math
    implicit none
    integer,intent(in) :: nbas,natom,nshl,idout
    integer,intent(in) :: shl_belong_to_atom(nshl),angl(nshl)
    real*8,intent(in) :: S(nbas,nbas),D_alpha(nbas,nbas),D_beta(nbas,nbas)
    real*8,intent(inout) :: MLK_charge(natom),LDW_charge(natom)
    real*8 S_haf(nbas,nbas),DS(nbas,nbas),SDS(nbas,nbas),D(nbas,nbas)
    type(atomtype),intent(in) :: atom(natom)
    integer nbas_in_atom(natom),i,j,k,m,n,ibas,jbas

    write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(idout,"(a)") "Population Analysis"
    write(idout,*)
    !write(*,*) "test1"
    D=D_alpha+D_beta

    nbas_in_atom=0

    do i=1,natom
        do j=1,nshl
            if(shl_belong_to_atom(j)==i)then
                !write(*,*) "i,j,angl(j),shl_degen_sph(angl(j))=",i,j,angl(j),shl_degen_sph(angl(j))
                nbas_in_atom(i)=nbas_in_atom(i)+shl_degen_sph(angl(j))
            end if
        end do
    end do

    !write(*,*) nbas_in_atom

    MLK_charge=real(atom(:)%charge)
    LDW_charge=real(atom(:)%charge)
    
    call mat_power(nbas,S,0.5_8,S_haf)
    DS=matmul(D,S)
    SDS=matmul(matmul(S_haf,D),S_haf)
    !write(*,*) "nbas1=",nbas
    
    ! rho=0
    ! do i=1,nbas
    !     rho=rho+SDS(i,i)
    ! end do
    ! write(*,*) "rho=",rho

    ibas=1
    do i=1,natom
        ! write(*,*) "ibas=",ibas
        ! write(*,*) "ibas+nbas_in_atom(i)-1=",ibas+nbas_in_atom(i)-1
        do k=ibas,ibas+nbas_in_atom(i)-1
            !write(*,*) "i,k=",i,k
            MLK_charge(i)=MLK_charge(i)-DS(k,k)
            LDW_charge(i)=LDW_charge(i)-SDS(k,k)
            !write(*,*) "MLK(",i,")=",LDW_charge(i)
        end do
        ibas=ibas+nbas_in_atom(i)
    end do

    !   write(*,*) "test4"
    write(idout,"(a)") "The Atomic Charge Analysis:"
    write(idout,*) "index_of_atom  notation_of_atom  Mulliken_Charge  Löwdin_Charge"
    do i=1,natom
        write(idout,*) i,atom(i)%name,MLK_charge(i),LDW_charge(i)
    end do
    !write(*,*) "nbas2=",nbas
    ! ibas=1
    ! do i=1,natom
    !     do j=1,natom
    !         bond_order(i,:)=0
    !     end do
    ! end do
    ! ! write(*,*) bond_order
    ! do i=1,natom
    !     jbas=1
    !     do j=1,natom
    !         !if(j==i)cycle
    !         bond_order(i,j)=0
    !         do m=ibas,ibas+nbas_in_atom(i)-1
                
    !             do n=jbas,jbas+nbas_in_atom(j)-1
    !                 !write(*,*) "i,j,m,n=",i,j,m,n
                    
    !                 bond_order(i,j)=bond_order(i,j)+DS(m,n)*DS(n,m)
    !             end do
    !         end do
    !         jbas=jbas+nbas_in_atom(j)
    !     end do
    !     ibas=ibas+nbas_in_atom(i)
    ! end do
    ! ! write(*,*) "nbas3=",nbas
    ! do i=1,natom
    !     write(*,*) "i=",i
    !     bond_order(i,i)=0
    ! end do

    ! write(idout,"(a)") "The Mayer Bond Order Analysis:"
    ! write(idout,*) "  ",1,"--",natom
    ! do i=1,natom
    !     write(*,*) "i=",i
    !     write(idout,*) i,bond_order(:,i)
    ! end do
    
    ! write(*,*) "nbas4=",nbas
    
    
    
    
    write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    ! write(*,*) "nbas5=",nbas
    ! write(*,*) "natom5=",natom

end subroutine