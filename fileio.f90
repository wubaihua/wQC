! This file is a part of wQC(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! filoio: the file input/output part of wQC. 



subroutine get_ninp(idinp,idout,natom,norder)
    implicit real*8(a-h,o-z)
    integer idinp,natom,norder,ierror,ifound
    character*200 c200

    norder=0
    
    do while(.true.)
        read(idinp,*,iostat=ierror) c200
        if(index(c200,"#")/=0)norder=norder+1
        if(index(c200,"*")/=0)exit
    end do

   
    if(ierrir/=0)then
        write(idout,*) "Error: Molecule infomation not found!" 
        return
    end if
    

    natom=0
    read(idinp,*)
    do while(.true.)
        read(idinp,*,iostat=ierror) c200
       
        if(c200(1:1)=="*")exit
        natom=natom+1
    end do
    
    !write(*,*) "natom=",natom

end subroutine


subroutine read_inp(idinp,natom,atom,chr,spinmul,norder,order)
    use def
    implicit real*8(a-h,o-z)
    type(atomtype) atom(natom)
    integer:: chr,spinmul
    character*20 order(norder)
    character*200 c200

    i=1
    do while(.true.)
        read(idinp,*,iostat=ierror) c200
        if(index(c200,"#")/=0) then
            backspace(idinp)
            read(idinp,*) order(i)
            order(i)(1:19)=order(i)(2:20)
            i=i+1
        end if
        if(index(c200,"*")/=0)exit
    end do

    
    read(idinp,*) chr, spinmul
    do i=1,natom
        read(idinp,*) atom(i)%name,atom(i)%x,atom(i)%y,atom(i)%z
        do j=1,109
            if(name2ind(j)==atom(i)%name)then
                atom(i)%index=j
                atom(i)%charge=j
                exit
            end if
        end do
    end do
    


end subroutine
    





subroutine get_bas_para(idbas,nshl,nprim,nbas,atom,natom)
    use def
    implicit real*8(a-h,o-z)
    integer :: idbas,nshl,nbas,iatom,nprim,ifound,natom
    type(atomtype) :: atom(natom)
    character*20 :: c20
    character*2 :: shell
    character*200 :: c200
    nshl=0
    nbas=0
    nprim=0
    !write(*,*) "natom=",natom
  
    
    do i=1,natom
        rewind(idbas)
        do while(index(c200,trim(adjustl(atom(i)%name))//"     0")/=1)
            read(idbas,"(a)",iostat=ierror) c200
        end do
        
        do while(.true.)
            !read(idbas,*) c200
            !write(*,*) "c200=",c200
            read(idbas,*) shell,np,i2 
            !write(*,*) "shell=",shell
            ! write(*,*) "np=",np
            ! write(*,*) "i2=",i2
            !nshl=nshl+1
            if(shell=='SP')then
                nshl=nshl+2
                nprim=nprim+2*np
            else
                nshl=nshl+1
                nprim=nprim+np
            end if
            do j=1,np
                read(idbas,*)
            end do
            read(idbas,*) c200
            if(c200(1:4)=="****")then
                exit
            else
                backspace(idbas)
            end if
        end do
    end do

end subroutine


subroutine read_bas(idbas,idout,nshl,nprim,nbas,atom,natom,cntr_odr,angl,shl_belong_to_atom,sh_indx,expnt,coeff)
    use def
    implicit real*8(a-h,o-z)
    integer :: idbas,nshl,nbas,iatom,nprim,ifound,natom,idout
    integer cntr_odr(nshl),angl(nshl),shl_belong_to_atom(nshl),sh_indx(nshl)
    real*8 expnt(nprim),coeff(nprim)
    type(atomtype) :: atom(natom)
    character*20 :: c20
    character*2 :: shell
    character*200 :: c200
    
    !write(*,*) "natom=",natom
    
    ishl=1
    iprim=1
    
    do i=1,natom
        rewind(idbas)
        do while(index(c200,trim(adjustl(atom(i)%name))//"     0")/=1)
            read(idbas,"(a)",iostat=ierror) c200
        end do
        do while(.true.)
            read(idbas,*) shell,np,i2 
            select case(shell)
            case("S ")
                angl(ishl)=0
                cntr_odr(ishl)=np
                shl_belong_to_atom(ishl)=i
                ishl=ishl+1
                do j=1,np
                    read(idbas,*) expnt(iprim),coeff(iprim)
                    iprim=iprim+1
                end do
            case("P ")
                angl(ishl)=1
                cntr_odr(ishl)=np
                shl_belong_to_atom(ishl)=i
                ishl=ishl+1
                do j=1,np
                    read(idbas,*) expnt(iprim),coeff(iprim)
                    iprim=iprim+1
                end do
            case("D ")
                angl(ishl)=2
                cntr_odr(ishl)=np
                shl_belong_to_atom(ishl)=i
                ishl=ishl+1
                do j=1,np
                    read(idbas,*) expnt(iprim),coeff(iprim)
                    iprim=iprim+1
                end do
            case("F ")
                angl(ishl)=3
                cntr_odr(ishl)=np
                shl_belong_to_atom(ishl)=i
                ishl=ishl+1
                do j=1,np
                    read(idbas,*) expnt(iprim),coeff(iprim)
                    iprim=iprim+1
                end do
            case("SP")
                angl(ishl)=0
                cntr_odr(ishl)=np
                shl_belong_to_atom(ishl)=i
                ishl=ishl+1
                do j=1,np
                    read(idbas,*) expnt(iprim),coeff(iprim),x2
                    iprim=iprim+1
                end do
                angl(ishl)=1
                cntr_odr(ishl)=np
                shl_belong_to_atom(ishl)=i
                ishl=ishl+1
                do j=1,np
                    backspace(idbas)
                end do
                do j=1,np
                    read(idbas,*) expnt(iprim),x2,coeff(iprim)
                    iprim=iprim+1
                end do  
            end select
            !cntr_odr(ishl)=np
            !shl_belong_to_atom(ishl)=i
            !ishl=ishl+1    
            !do j=1,np
            !    read(idbas,*) expnt(iprim),coeff(iprim)
            !    iprim=iprim+1
            !end do
            read(idbas,*) c200
            if(c200(1:4)=="****")then
                exit
            else
                backspace(idbas)
            end if
        end do
    end do
    
    sh_indx(1)=1
    do i=2,nshl
        !write(*,*) shl_degen_sph(angl(i))
        sh_indx(i)=sh_indx(i-1)+shl_degen_sph(angl(i-1))
    end do
    
    nbas=0
    do i=1,nshl
        nbas=nbas+shl_degen_sph(angl(i))
    end do

    write(idout,*) 'The number of shells=',nshl
    write(idout,*) 'The number of primitive Gaussian functions=',nprim
    write(idout,*) 'The number of basis functions=',nbas
    
    !write(*,*) "cntr_odr=",cntr_odr
    !write(*,*) "angl=",angl
    !write(*,*) "shl_belong_to_atom",shl_belong_to_atom
    !write(*,*) "expnt=",expnt
    !write(*,*) "coeff=",coeff
    !write(*,*) "sh_indx=",sh_indx

end subroutine





!This function is taken from source code of Multiwfn by Dr.Tian Lu.    
subroutine loclabel(fileid,label,ifound,irewind,maxline)
integer fileid,ierror
integer,optional :: ifound,irewind,maxline
character*200 c200
CHARACTER(len=*) ::  label
if ((.not.present(irewind)).or.(present(irewind).and.irewind==1)) rewind(fileid)
if (.not.present(maxline)) then
	do while(.true.)
		read(fileid,"(a)",iostat=ierror) c200
		if (index(c200,label)/=0) then
			backspace(fileid)
			if (present(ifound)) ifound=1 !Found result
			return
		end if
		if (ierror/=0) exit
	end do
else
	do iline=1,maxline
		read(fileid,"(a)",iostat=ierror) c200
		if (index(c200,label)/=0) then
			backspace(fileid)
			if (present(ifound)) ifound=1 !Found result
			return
		end if
		if (ierror/=0) exit
	end do
end if
if (present(ifound)) ifound=0
end subroutine







subroutine out_init(idout,filepath)
    implicit none
    integer idout
    character*200 filepath
    character nowdate*20,nowtime*20
    
    write(idout,"(a)") "##########################################################################"
    write(idout,"(a)") "#  wQC: a simple fortran Quantum Chemistry/Electronic Structure Program  #"
    write(idout,"(a)") "#  Author:   Baihua Wu   (wubaihua@pku.edu.cn)                           #"    
    write(idout,"(a)") "##########################################################################"
    write(idout,*)
    write(idout,*) "wQC starts at"
    call date_and_time(nowdate,nowtime)
    write(idout,"('Date: ',a,'-',a,'-',a,'   Time: ',a,':',a,':',a)") nowdate(1:4),nowdate(5:6),nowdate(7:8),nowtime(1:2),nowtime(3:4),nowtime(5:6)
    write(idout,*)
    write(idout,"(a)") "wQC input file:"
    write(idout,"(a)") trim(filepath)



end subroutine






subroutine error_end(idout,massage)
    implicit none
    integer idout
    character(len=*) massage
    character nowdate*20,nowtime*20

    write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(idout,"(a)") massage
    write(idout,"(a)") "wQC ERROR termination at"
    call date_and_time(nowdate,nowtime)
    write(idout,"('Date: ',a,'-',a,'-',a,'   Time: ',a,':',a,':',a)") nowdate(1:4),nowdate(5:6),nowdate(7:8),nowtime(1:2),nowtime(3:4),nowtime(5:6)
    close(idout)
    stop
end subroutine


subroutine normal_end(idout)
    implicit none
    integer idout
    character nowdate*20,nowtime*20

    write(idout,"(a)") "----------------------------------------"
    write(idout,"(a)") "wQC NORMAL termination at"
    call date_and_time(nowdate,nowtime)
    write(idout,"('Date: ',a,'-',a,'-',a,'   Time: ',a,':',a,':',a)") nowdate(1:4),nowdate(5:6),nowdate(7:8),nowtime(1:2),nowtime(3:4),nowtime(5:6)
    close(idout)
    stop
end subroutine
