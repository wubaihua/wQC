! This file is a part of wQC(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! Last update: 2020-5-3

! module init: initialization wQC. The arrays in wQC are defined here. 


module init
    use def

    character*200 filepath !input file load
    character*200 workpath !wQC load 
    character*200 string
    character*200 basispath !basis-set file load
    character*200 outpath !output file load
    character*20 basset !basis-set name
    !character, allocatable :: filename(:)
    integer chr !molecular charge
    integer spinmul !molecular spin multiplicity
    type(atomtype),allocatable :: atom(:) !atom information
    integer natom !number of atoms
    integer nele,& !total number of electrons
            nele_alpha,& !number of alpha electrons
            nele_beta !number of beta electrons

    integer nshl !number of shells
    integer nprim !number of primitive Gaussian functions
    integer nbas !number of basis functions
    integer,allocatable :: cntr_odr(:) !number of primitive gaussian function of each shell
    integer,allocatable :: angl(:) !angular momentum of shell
    integer,allocatable :: shl_belong_to_atom(:) !shell belongs to atom
    integer,allocatable :: sh_indx(:) !shell index in array
    integer,allocatable :: charge(:) !atom charge
    real*8,allocatable :: expnt(:),& !GTO exponent
                          coeff(:),& !GTO coefficient
                          geom(:,:)  !atom geom information
    real*8,allocatable :: S(:,:),& !overlap matrix
                          T(:,:),& !kinetic matrix 
                          V(:,:),& !potential matrix
                          eri(:,:,:,:) !electron repulsion integral (2-e integral)
    real*8,allocatable :: C(:,:),& !coefficient matrix of orbital
                          D(:,:),& !density matrix of orbital
                          E(:) !eigenvalue of orbital
    real*8,allocatable :: MLK_charge(:),& !Mulliken charge 
                          LDW_charge(:) !Lowdin charge                     
    real*8,allocatable :: D_alpha(:,:),& !density matrix of alpha electron orbital
                          D_beta(:,:),& !density matrix of beta electron orbital
                          E_alpha(:),& !eigenvalue of alpha electron orbital
                          E_beta(:) !eigenvalue of beta electron orbital

    integer norder !number of orders
    character*20,allocatable :: order(:) !all orders 
    
    real*8 nucp !nuclear repulsion energy




    contains

    subroutine init_atom()

        allocate(atom(natom)) 
        allocate(order(norder))

    end subroutine



    subroutine init_bas()

        allocate(cntr_odr(nshl))
        allocate(angl(nshl))
        allocate(shl_belong_to_atom(nshl))
        allocate(sh_indx(nshl))
        allocate(expnt(nprim))
        allocate(coeff(nprim))
        allocate(geom(3,natom))
        allocate(charge(natom))

        geom(1,:)=atom(:)%x
        geom(2,:)=atom(:)%y
        geom(3,:)=atom(:)%z
        charge=atom%charge


        

        
    end subroutine


    subroutine init_HF()

        allocate(S(nbas,nbas))
        allocate(T(nbas,nbas))
        allocate(V(nbas,nbas))
        allocate(eri(nbas,nbas,nbas,nbas))

        
        allocate(C(nbas,nbas))
        allocate(D(nbas,nbas))
        allocate(E(nbas))
        allocate(D_alpha(nbas,nbas))
        allocate(E_alpha(nbas))
        allocate(D_beta(nbas,nbas))
        allocate(E_beta(nbas))

        allocate(MLK_charge(natom))
        allocate(LDW_charge(natom))

    end subroutine

   


end module
