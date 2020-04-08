program MTChem
    use def
    use cint
    
    implicit real*8(a-h,o-z)
    character*200::filepath,string
    integer:: chr,spinmul
    integer,external :: GetFileN
    type(atomtype) atom(:)
    allocatable atom
    integer,allocatable :: cntr_odr(:),angl(:),shl_belong_to_atom(:),sh_indx(:),charge(:)
    real*8,allocatable :: expnt(:),coeff(:),geom(:,:)
    real*8,allocatable :: S(:,:),T(:,:),V(:,:),eri(:,:,:,:)
    real*8,allocatable :: C(:,:),D(:,:),E(:)
    real*8 nucp
    
    write(*,*)"Input the file path:"
    
!11  read*, filepath
    filepath="H2O.txt"
        
    
    open(10, file=trim(filepath), status='old')    
    write(*,*) "load file successfully!"
    call get_natom(10,natom)
    rewind(10)
    allocate(atom(natom)) 

    call read_inp(10,natom,atom,chr,spinmul)
    nele=sum(atom%index)-chr
    write(*,*) "The number of atoms:" ,natom
    write(*,*) "The number of electrons:" ,nele
    write(*,*) "The charge of molecule:" ,chr
    if(spinmul==1)write(*,*)'Close Shell molecule'
    if(spinmul>1)write(*,*)'Open Shell molecule'
    close(10)
    
    open(20,file="sto-3g.gbs",status="old")
    call get_bas_para(20,nshl,nprim,nbas,atom,natom)
    write(*,*) 'nshl=',nshl
    write(*,*) 'nprim=',nprim
    
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
    
    write(*,*) "nbas=",nbas
    
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
    write(*,*) "nucp=",nucp
    allocate(C(nbas,nbas))
    allocate(D(nbas,nbas))
    allocate(E(nbas))
    call RHF(nbas,nele,nucp,S,T,V,eri,D,E,C)
    
    
    
    
    
pause
end program



