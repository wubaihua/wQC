! This file is a part of wQC(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! Last update: 2020-4-18

! PT: the part about many-Body perturbation theory(MBPT) method, including 
! Moller-Plesset 2-order perturbation(MP2).



subroutine MP2(idout,nbas,nele,E,C,eri,E_mp2)
    implicit none
    integer idout,nbas,nele
    integer sigma,mu,lambda,nu
    integer p,q,r,s
    integer i,j,ia,ib
    real(kind=8) :: E(nbas),C(nbas,nbas),eri(nbas,nbas,nbas,nbas),E_mp2
    real(kind=8) :: erimo(nbas,nbas,nbas,nbas)
    real(kind=8) sum1,sum2,sum3

    
    write(idout,"(a)") "---------------------------------------------------"
    write(idout,"(a)") "Moller-Plesset 2-order perturbation(MP2) Method"
    

    sum1=0.0
    sum2=0.0
    sum3=0.0
    erimo=0.0

   

    do p=1,nbas
        do q=1,nbas
            do r=1,nbas
                do s=1,nbas
                    do mu=1,nbas
                        do nu=1,nbas
                            do lambda=1,nbas
                                do sigma=1,nbas
                                    sum1=sum1+C(sigma,s)*eri(mu,nu,lambda,sigma)
                                end do
                                sum2=sum2+C(lambda,r)*sum1
                                sum1=0
                            end do
                            sum3=sum3+C(nu,q)*sum2
                            sum2=0
                        end do
                        erimo(p,q,r,s)=erimo(p,q,r,s)+C(mu,p)*sum3
                        sum3=0
                    end do
                end do
            end do
        end do
    end do

    
    E_mp2=0
    do i=1,nele/2
        do j=1,nele/2
            do ia=nele/2+1,nbas
                do ib=nele/2+1,nbas
                    E_mp2=E_mp2+erimo(i,ia,j,ib)*(2.0*erimo(i,ia,j,ib)-erimo(i,ib,j,ia))/(E(i)+E(j)-E(ia)-E(ib))
                end do
            end do
        end do
    end do

    

    write(idout,*) "MP2 correction energy (a.u.):"
    write(idout,*) "delta E(MP2)=",E_mp2

    write(idout,"(a)") "---------------------------------------------------"
   

end subroutine
