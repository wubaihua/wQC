module geo
    use def
    use constant
    use math
contains  

subroutine molgeo(natom,atom,distance,angle,diangle)
    use def
    implicit real*8(a-h,o-z)
    type(atomtype),intent(in) :: atom(:)
    real*8 distance(natom,natom),angle(natom,natom,natom),diangle(natom,natom,natom,natom)
    
    
    do i=1,natom
        do j=1,natom
            call geo_distance(atom,i,j,distance(i,j))
            do m=1,natom
                call geo_angle(atom,i,j,k,angle(i,j,m))
                do n=1,natom
                    call geo_diangle(atom,i,j,m,n,diangle(i,j,m,n))
                end do
            end do    
        end do
    end do
    
end subroutine
    
    
    
subroutine geo_distance(atom,i,j,d)
    implicit none
    integer i,j
    real*8 d
    type(atomtype),intent(in) :: atom(:)
    
    d=sqrt((atom(i)%x-atom(j)%x)**2+(atom(i)%y-atom(j)%y)**2+(atom(i)%z-atom(j)%z)**2)
    
    
end subroutine



subroutine geo_angle(atom,i,j,k,an)
    implicit none
    integer i,j,k
    real*8 an,L1,L2,prod
    type(atomtype),intent(in) :: atom(:)
    
    prod=(atom(i)%x-atom(j)%x)*(atom(k)%x-atom(j)%x)+(atom(i)%y-atom(j)%y)*(atom(k)%y-atom(j)%y)+(atom(i)%z-atom(j)%z)*(atom(k)%z-atom(j)%z)
    L1=sqrt((atom(i)%x-atom(j)%x)**2+(atom(i)%y-atom(j)%y)**2+(atom(i)%z-atom(j)%z)**2)
    L2=sqrt((atom(k)%x-atom(j)%x)**2+(atom(k)%y-atom(j)%y)**2+(atom(k)%z-atom(j)%z)**2)

    an=acos(prod/(L1*L2))*180/pi
end subroutine

subroutine geo_diangle(atom,i,j,m,n,dian)
    implicit none
    integer i,j,m,n
    real*8 dian
    type(atomtype),intent(in) :: atom(:)
    real*8 a(3),b(3),c(3),d(3)
    
    a(1)=atom(i)%x
    a(2)=atom(i)%y
    a(3)=atom(i)%z
    b(1)=atom(j)%x
    b(2)=atom(j)%y
    b(3)=atom(j)%z
    c(1)=atom(m)%x
    c(2)=atom(m)%y
    c(3)=atom(m)%z
    d(1)=atom(n)%x
    d(2)=atom(n)%y
    d(3)=atom(n)%z
    
    call math_diangle(a,b,c,d,dian)    


end subroutine

end module



