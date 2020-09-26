module math
    use constant
    use LAPACK95
contains

!Calculate the distance between two points r1 and r2.
subroutine math_distance(r1,r2,d)
    implicit none
    real*8,intent(in) :: r1(:),r2(:)
    real*8 d
    integer nsize,i
    
    nsize=size(r1)
    d=0
    do i=1,nsize
        d=d+(r1(i)-r2(i))**2
    end do
    
    d=sqrt(d)
end subroutine
!calcualte the norm or \vec x, i.e. |\vec x|.
subroutine math_vecnorm(x,norm)
    implicit none
    real*8,intent(in) :: x(:)
    real*8 norm
    integer nsize,i
    norm=0
    do i=1,nsize
        norm=norm+x(i)**2
    end do
    norm =sqrt(norm)


end subroutine

!calculate \vec x \cdot \vec y 
subroutine math_dotprod(x,y,p)
    implicit none
    real*8,intent(in):: x(:),y(:)
    real*8 p
    integer nsize,i
    
    p=0
    nsize=size(x)
    do i=1,nsize
        p=p+x(i)*y(i)
    end do
    
end subroutine

!calculate \vec x \times \vec y,only in dim=3
subroutine math_timeprod(x,y,p)
    implicit none
    real*8 x(3),y(3)
    real*8 p(3)
    integer nsize,i
    
    p(1)=x(2)*y(3)-x(3)*y(2)
    p(2)=-x(1)*y(3)+x(3)*y(1)
    p(3)=x(1)*y(2)-x(2)*y(1)
    
end subroutine



!Calculate the angle of x-y-z.
subroutine math_angle(x,y,z,theta)
    implicit none
    real*8,intent(in) :: x(:),y(:),z(:)
    real*8 theta,L1,L2,prod
    integer nsize,i
    
    nsize=size(x)
    prod=0
    
    do i=1,nsize
        prod=prod+(x(i)-y(i))*(z(i)-y(i))
    end do
    call math_distance(x,y,L1)
    call math_distance(z,y,L2)
    
    theta=acos(prod/(L1*L2))*180/pi   
end subroutine

!Calculate the dihedral angle between a-b-c-d,only in dim=3.
subroutine math_diangle(a,b,c,d,dia)    
    implicit none
    real*8 a(3),b(3),c(3),d(3),dia
    real*8 ab(3),bc(3),cd(3),g(3),h(3)
    real*8 p,norm1,norm2
    ab=a-b
    bc=b-c
    cd=c-d
    
    call math_timeprod(bc,ab,g)
    call math_timeprod(bc,cd,h)
    call math_dotprod(g,h,p)
    call math_vecnorm(g,norm1)
    call math_vecnorm(h,norm2)
    
    dia=acos(p/(norm1*norm2))*180/pi
end subroutine

!Generate two random number x1,x2 satisfy Gaussian distribution N(miu,sigma).
subroutine box_muller(x1,x2,sigma,miu)
    implicit none
    real*8, parameter :: pi = 3.14159265358979323846
    real*8, parameter :: hbar    = 1.0
    real*8 x1,x2,u1,u2,sigma,miu
    call RANDOM_NUMBER(u1)
    call RANDOM_NUMBER(u2)
    x1=sqrt(-2*log(u1))*cos(2*pi*u2)
    x2=sqrt(-2*log(u1))*sin(2*pi*u2)
    !x1=(x1-miu)/sigma
    !x2=(x2-miu)/sigma
    x1=x1*sigma+miu
    x2=x2*sigma+miu
end subroutine
   
function Kronecker_delta(i,j)
    integer Kronecker_delta,i,j
    if(i==j)then
        Kronecker_delta=1
    else
        Kronecker_delta=0
    end if
    
end function



function heaviside(x)
    implicit none
    real*8 x
    integer heaviside
    if(x<0)then
        heaviside=0
    else
        heaviside=1
    end if
    

end function




! subroutine heaviside(x,h)
!     implicit none
!     real*8  h,x
!     if(x<0)then
!         h=0
!     else
!         h=1
!     end if
    

! end subroutine



! subroutine window1_pccp(gamma,n1,n2,w1,w2)
!     implicit none
    
!     real*8 n1,n2,w1,w2,gamma,h1,h2
    
!     call heaviside(gamma-abs(n1-1),h1)
!     w1=1/(2*gamma)*h1
!     call heaviside(gamma-abs(n1),h2)
!     w2=1/(2*gamma)*h2


! end subroutine

!Diagonal symmetric matrix
subroutine dia_symmat(n,A,E,C)
    implicit none
    integer n,info
    real*8 A(n,n),E(n),C(n,n),work(3*n+1)
    
    C=A
    call dsyev('V','L',n,C,n,E,work,3*n+1,info)

end subroutine

!calculate matrix A to the power of r : A^r 
subroutine mat_power(n,A,r,AR)
    implicit none
    integer n,i
    real*8 A(n,n),r,AR(n,n),U(n,n),E(n),Er(n,n)
    
    call dia_symmat(n,A,E,U)
    Er=0
    do i=1,n
        Er(i,i)=E(i)**r
    end do
    
    AR=matmul(matmul(U,Er),transpose(U))
    
end subroutine

!calculate e to the power of symmetric matrix A: e^A
subroutine exp_mat(n,A,eA)
    implicit none
    integer n,i
    real*8 A(n,n),E(n),U(n,n),eE(n,n),eA(n,n)

    call dia_symmat(n,A,E,U)
    eE=0
    do i=1,n
        eE(i,i)=exp(E(i))
    end do

    eA=matmul(matmul(U,eE),transpose(U))

end subroutine

!calcualte the dot product of two matrix A,B : x=A \cdot B
subroutine mat_dot_prod(n,A,B,x)
    implicit none
    integer i,n,j
    real*8,intent(in) :: A(n,n),B(n,n)
    real*8,intent(inout) :: x
    
    x=0
    do i=1,n
        do j=1,n
            x=x+A(i,j)*B(i,j)
        end do
    end do

end subroutine
function mat_dotprod(m,n,A,B)
    implicit none
    integer i,m,n,j
    real*8,intent(in) :: A(m,n),B(m,n)
    real*8 x,mat_dotprod
    
    x=0
    do i=1,m
        do j=1,n
            x=x+A(i,j)*B(i,j)
        end do
    end do

    mat_dotprod=x

end function


!save the linear equation set :Ax=b
subroutine solv_LES(n,A,x,b)
    implicit none
    integer n,ipiv(n),info
    real*8,intent(in) :: A(n,n),b(n)
    real*8 A1(n,n)
    real*8,intent(inout) :: x(n)

    A1=A
    x=b

    call dgesv( n, 1, A1, n, ipiv, x, n, info )

end subroutine











    
end module