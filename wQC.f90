! wQC: a simple fortran Quantum Chemistry/Electronic Structure Program
!(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! This program is still unfinished. It may still have some bugs which
! can make wrong result. I will still develop and maintain this program. 



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
    real*8 nucp
    
    write(*,*)"Input the file path:"
    
    read*, filepath
    !filepath="H2.txt     "
        
    
    open(10, file=trim(filepath), status='old')  

    open(15,file=filepath(1:len(trim(filepath))-4)//".out")

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
    close(10)
    
    basispath="basis/"//"def2svp"//".gbs"
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
    
    write(15,*) "nbas=",nbas
    
    allocate(S(nbas,nbas))
    allocate(T(nbas,nbas))
    allocate(V(nbas,nbas))
    allocate(eri(nbas,nbas,nbas,nbas))
    
    call cal_eint(nbas,natom,nprim,nshl,cntr_odr,charge,angl,shl_belong_to_atom,sh_indx,expnt,coeff,geom,S,T,V,eri)
    

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
    call RHF(15,nbas,nele,nucp,S,T,V,eri,D,E,C)
    
    
    
    
    

end program



