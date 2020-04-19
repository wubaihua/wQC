! wQC: a simple fortran Quantum Chemistry/Electronic Structure Program
!(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! This program is still unfinished. It may still have some bugs which
! can make wrong result. I will still develop and maintain this program. 

! Last update: 2020-4-17

! This file is the main program part of wQC.

program wQC
    use def
    use cint
    
    implicit real*8(a-h,o-z)
    character*200::filepath,string,basispath,outpath
    character, allocatable :: filename(:)
    integer:: chr,spinmul
    type(atomtype) atom(:)
    allocatable atom
    integer,allocatable :: cntr_odr(:),angl(:),shl_belong_to_atom(:),sh_indx(:),charge(:)
    real*8,allocatable :: expnt(:),coeff(:),geom(:,:)
    real*8,allocatable :: S(:,:),T(:,:),V(:,:),eri(:,:,:,:)
    real*8,allocatable :: C(:,:),D(:,:),E(:)
    real*8,allocatable :: MLK_charge(:),LDW_charge(:)
    real*8 nucp
    
    

    write(*,*)"Input the file path:"
    
    !read*, filepath
    filepath="methyl.inp     "
        
    
    open(10, file=trim(filepath), status='old')  

    open(15,file=filepath(1:len(trim(filepath))-4)//".out")

    call cpu_time(t1)

    call out_init(15,filepath)

    write(15,*) "load file successfully!"

    call get_natom(10,natom)
    rewind(10)
    allocate(atom(natom)) 

    call read_inp(10,natom,atom,chr,spinmul)
    nele=sum(atom%index)-chr
    write(15,*) "The number of atoms:" ,natom
    write(15,*) "The number of electrons:" ,nele
    write(15,*) "The charge of molecule:" ,chr
    if(spinmul==1)write(15,*)'Close Shell molecule'
    if(spinmul>1)write(15,*)'Open Shell molecule'
    nele_alpha=(nele-(spinmul-1))/2+(spinmul-1)
    nele_beta=(nele-(spinmul-1))/2
    write(15,*) "The number of Alpha electrons:",nele_alpha
    write(15,*) "The number of Beta electrons:",nele_beta
    close(10)
    
    basispath="basis/"//"sto-3g"//".gbs"
    open(20,file=basispath,status="old")
    call get_bas_para(20,nshl,nprim,nbas,atom,natom)
    write(15,*) 'nshl=',nshl
    write(15,*) 'nprim=',nprim
    
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
    
    call read_bas(20,nshl,nprim,nbas,atom,natom,cntr_odr,angl,shl_belong_to_atom,sh_indx,expnt,coeff)
    
    !write(15,*) "nbas=",nbas
    ! nbas2=nbas
    ! nele2=nele
    allocate(S(nbas,nbas))
    allocate(T(nbas,nbas))
    allocate(V(nbas,nbas))
    allocate(eri(nbas,nbas,nbas,nbas))
    
    call cal_eint(nbas,natom,nprim,nshl,cntr_odr,charge,angl,shl_belong_to_atom,sh_indx,expnt,coeff,geom,S,T,V,eri)
    
    ! write(*,*) "nbas1=",nbas
    ! write(*,*) "nele1=",nele
    ! write(*,*) "nshl1=",nshl
    ! write(*,*) "natom1=",natom

	!read(*,*) istep
	!if(istep==1)then
 !      call HF(natom,atom)
 !   else
 !       write(*,*) 'wrong'
 !   end if
    
    call cal_nucp(atom,natom,nucp)
    write(15,*) "nucp=",nucp
    allocate(C(nbas,nbas))
    allocate(D(nbas,nbas))
    allocate(E(nbas))
    write(*,*) "nele_alpha=",nele_alpha
    !call RHF(15,nbas,nele,nucp,S,T,V,eri,D,E,C)
    !call RHF_DIIS(15,nbas,nele,nucp,S,T,V,eri,D,E,C,6)
    call UHF(15,nbas,nele_alpha,nele_beta,nucp,S,T,V,eri,D_alpha,D_beta,E_alpha,E_beta)
    
    allocate(MLK_charge(natom))
    allocate(LDW_charge(natom))
    ! write(*,*) "nbas=",nbas
    ! write(*,*) "nele=",nele

    !call pop_analy(15,atom,nshl,nbas,natom,D,D,S,shl_belong_to_atom,angl,MLK_charge,LDW_charge)
    !call pop_analy(15,atom,nshl,nbas,natom,D_alpha,D_beta,S,shl_belong_to_atom,angl,MLK_charge,LDW_charge)
    ! nbas=nbas2
    ! nele=nele2
    ! write(*,*) "nbas=",nbas
    ! write(*,*) "nele=",nele
    ! write(*,*) "nshl=",nshl
    ! write(*,*) "natom=",natom

    !call MP2(15,nbas,nele,E,C,eri,E_mp2)

    call cpu_time(t2)
    write(15,*) "Job Time:",t2-t1,"Seconds"
    write(15,"(a)") "End wQC"
    
    
    

end program



