



subroutine get_basset(idout,norder,order,basset)
    implicit none
    integer i,idout,norder
    character*20,intent(in) :: order(norder)
    character*20,intent(inout) :: basset
    logical ifbas

    write(idout,"(a)") "---------------------------------------------------"
    do i=1,norder
        inquire(file='basis/'//trim(adjustl(order(i)))//'.gbs',exist=ifbas)
        !write(*,*) 'basis/'//trim(adjustl(order(i)))//'.gbs',ifbas
        if(ifbas)then
            basset=order(i)
            write(idout,*) "basis-set:",basset
            return
        end if
    end do
    if(.not.ifbas)then
        write(idout,"(a)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(idout,"(a)") "WORRING: basis-set input not found!"
        write(idout,"(a)") "Be sure the basis-set file(.gbs) exists in /basis/"
        write(idout,"(a)") "Using the default basis-set: sto-3g"
        basset="sto-3g"
    end if

end subroutine