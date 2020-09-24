subroutine ao2mo
    use init
    implicit real*8(a-h,o-z)

    write(idout,"(a)") "---------------------------------------------------"
    write(idout,"(a)") "Transform the 2-e integrals to MO basis:"
    

    sum1=0.0
    sum2=0.0
    sum3=0.0
    erimo=0.0

   

    do p=1,nbas
        do q=1,nbas
            do r=1,nbas
                do s1=1,nbas
                    do mu=1,nbas
                        do nu=1,nbas
                            do lambda=1,nbas
                                do sigma=1,nbas
                                    sum1=sum1+C(sigma,s1)*eri(mu,nu,lambda,sigma)
                                end do
                                sum2=sum2+C(lambda,r)*sum1
                                sum1=0
                            end do
                            sum3=sum3+C(nu,q)*sum2
                            sum2=0
                        end do
                        erimo(p,q,r,s1)=erimo(p,q,r,s1)+C(mu,p)*sum3
                        sum3=0
                    end do
                end do
            end do
        end do
    end do


end subroutine


subroutine build_spinorb_inte
    use init
    implicit real*8(a-h,o-z)
    integer,external :: index_spatial2spin

    do i=1,nspinorb
        do j=1,nspinorb
            do k=1,nspinorb
                do l=1,nspinorb

                    if(mod(i,2)/=mod(k,2))then
                        spinorb_inte(i,j,k,l)=0
                    else if(mod(j,2)/=mod(l,2))then
                        spinorb_inte(i,j,k,l)=0
                    else
                        spinorb_inte(i,j,k,l)=erimo(index_spatial2spin(i),index_spatial2spin(k),index_spatial2spin(j),index_spatial2spin(l))
                    end if

                end do
            end do
        end do
    end do

    ! write(*,*) spinorb_inte

end subroutine


function index_spatial2spin(i)
    integer i,j,index_spatial2spin

    if(mod(i,2)==0)then
        index_spatial2spin=i/2
    else
        index_spatial2spin=(i+1)/2
    end if
end function
