! This file is a part of wQC.

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! module cint: using Libcint(J.Compet.Chem.2015,36,1664–1671) to calculate 
! single and double electronic integral and build the overlap, kinetic and 
! potential matrix.

! The code in this file refers to St Maxwell's website
! (https://gensoukyo.me/use-libcint/). Author thanks for his contribution.
! Author also wants to thanks Dr. Qiming Sun, the author of Libcint, for 
! some advices about using Libcint.


module cint
    implicit none
    integer, dimension(:,:), allocatable :: atm
    integer, dimension(:,:), allocatable :: bas
    real*8, dimension(:), allocatable :: env
  
    integer, parameter :: CHARGE_OF  = 1
    integer, parameter :: PTR_COORD  = 2
    integer, parameter :: NUC_MOD_OF = 3
    integer, parameter :: PTR_ZETA   = 4
    integer, parameter :: ATM_SLOTS  = 6
    integer, parameter :: ATOM_OF    = 1
    integer, parameter :: ANG_OF     = 2
    integer, parameter :: NPRIM_OF   = 3
    integer, parameter :: NCTR_OF    = 4
    integer, parameter :: KAPPA_OF   = 5
    integer, parameter :: PTR_EXP    = 6
    integer, parameter :: PTR_COEFF  = 7
    integer, parameter :: BAS_SLOTS  = 8
    integer, parameter :: PTR_ENV_START = 20  
contains

    subroutine set_libcint_input(natm,nprm,nshl,cntr_odr,charge,angl,shl_belong_to_atom,sh_indx,expnt,coeff,geom)
        integer :: off, prim_off
        integer :: iatom, ishl, iprim, ioff
        integer :: natm,nprm,nshl
        integer :: charge(natm),cntr_odr(nshl),angl(nshl),shl_belong_to_atom(nshl),sh_indx(nshl)
        real*8 :: expnt(nprm),coeff(nprm),geom(3,natm)
        real*8, external :: CINTgto_norm
        
        allocate(atm(ATM_SLOTS,natm))
        allocate(bas(BAS_SLOTS,nshl))
        allocate(env(PTR_ENV_START+3*natm+2*nprm))
  
        off = PTR_ENV_START
        do iatom = 1, natm
            atm(CHARGE_OF,iatom) = charge(iatom)
            atm(PTR_COORD,iatom) = off
            atm(NUC_MOD_OF,iatom) = 1
            env(off+1:) = geom(:,iatom)
            off = off + 3
        end do
  
        prim_off = 1
        do ishl = 1, nshl
            bas(ATOM_OF,ishl) = shl_belong_to_atom(ishl) - 1
            bas(ANG_OF,ishl) = angl(ishl)
            bas(NPRIM_OF,ishl) = cntr_odr(ishl)
            bas(NCTR_OF,ishl) = 1
  
            bas(PTR_EXP,ishl) = off
            ioff = 0
            do iprim = prim_off, prim_off + cntr_odr(ishl)-1
                env(off+1+ioff) = expnt(iprim)
                ioff = ioff + 1
            end do
  
            off = off + cntr_odr(ishl)
  
            bas(PTR_COEFF,ishl) = off
            ioff = 0
            do iprim = prim_off, prim_off + cntr_odr(ishl)-1
                env(off+1+ioff) = coeff(iprim) * CINTgto_norm(angl(ishl), expnt(iprim))
                ioff = ioff + 1
            end do
  
            off = off + cntr_odr(ishl)
  
            prim_off = prim_off + cntr_odr(ishl)
        end do
  
    end subroutine set_libcint_input
    
    
    subroutine normalize(natm,nprm,nshl,cntr_odr,charge,angl,shl_belong_to_atom,sh_indx,expnt,coeff,geom)
        integer :: ishl, iprim
        integer :: di, dj
        real*8, dimension(:,:), allocatable :: buf1e
        integer, dimension(2) :: shls
        integer :: natm,nprm,nshl
        integer :: charge(natm),cntr_odr(nshl),angl(nshl),shl_belong_to_atom(nshl),sh_indx(nshl)
        real*8 :: expnt(nprm),coeff(nprm),geom(3,natm)
        integer, external :: CINTcgto_spheric
  
        do ishl = 1, nshl
            shls(1) = ishl - 1
            shls(2) = ishl - 1
            di = CINTcgto_spheric(ishl-1, bas)
            dj = CINTcgto_spheric(ishl-1, bas)
            allocate(buf1e(di,dj))
            call cint1e_ovlp_sph(buf1e,shls,atm,natm,bas,nshl,env)
            do iprim = 1, cntr_odr(ishl)
                env(bas(PTR_COEFF,ishl)+iprim) = env(bas(PTR_COEFF,ishl)+iprim) / sqrt(buf1e(1,1))
            end do
            deallocate(buf1e)
        end do
  
    end subroutine normalize
    
    
    subroutine cal_eint(nbas,natm,nprm,nshl,cntr_odr,charge,angl,shl_belong_to_atom,sh_indx,expnt,coeff,geom,S,T,V,eri)
        implicit none
        integer :: idbas,nshl,nbas,iatom,nprm,ifound,natm
        integer cntr_odr(nshl),angl(nshl),shl_belong_to_atom(nshl),sh_indx(nshl),charge(natm)
        real*8 expnt(nprm),coeff(nprm),geom(3,natm)
        integer :: i, j, k, l, di, dj, dk, dl, x, y, z, w
        real*8  S(nbas,nbas)
        real*8  T(nbas,nbas)
        real*8  V(nbas,nbas)
        real*8  eri(nbas,nbas,nbas,nbas)
        real*8, dimension(:,:), allocatable :: buf1e
        real*8, dimension(:,:,:,:), allocatable :: buf2e
        integer, dimension(4) :: shls
        integer, external :: CINTcgto_spheric
        
        call set_libcint_input(natm,nprm,nshl,cntr_odr,charge,angl,shl_belong_to_atom,sh_indx,expnt,coeff,geom)
    
        call normalize(natm,nprm,nshl,cntr_odr,charge,angl,shl_belong_to_atom,sh_indx,expnt,coeff,geom)
    
    
        do i = 1, nshl
            shls(1) = i - 1
            di = CINTcgto_spheric(i-1, bas)
            do j = 1, nshl
                shls(2) = j - 1
                dj = CINTcgto_spheric(j-1, bas)
                allocate(buf1e(di,dj))
                x = sh_indx(i); y = sh_indx(j)
                call cint1e_ovlp_sph(buf1e,shls,atm,natm,bas,nshl,env)
                S(x:,y:) = buf1e(:,:)
                call cint1e_kin_sph(buf1e,shls,atm,natm,bas,nshl,env)
                T(x:,y:) = buf1e(:,:)
                call cint1e_nuc_sph(buf1e,shls,atm,natm,bas,nshl,env)
                V(x:,y:) = buf1e(:,:)
                deallocate(buf1e)
            end do
        end do
        
        
        do i = 1, nshl
            shls(1) = i - 1
            di = CINTcgto_spheric(i-1, bas)
            do j = 1, nshl
                shls(2) = j - 1
                dj = CINTcgto_spheric(j-1, bas)
                do k = 1, nshl
                    shls(3) = k - 1
                    dk = CINTcgto_spheric(k-1, bas)
                    do l = 1, nshl
                        shls(4) = l - 1
                        dl = CINTcgto_spheric(l-1, bas)
                        allocate(buf2e(di,dj,dk,dl))
                        x = sh_indx(i); y = sh_indx(j)
                        z = sh_indx(k); w = sh_indx(l)
                        call cint2e_sph(buf2e,shls,atm,natm,bas,nshl,env,0_8)
                        eri(x:,y:,z:,w:) = buf2e(:,:,:,:)
                        deallocate(buf2e)
                    end do
                end do
            end do
        end do
    
 
    
    end subroutine
    
    
    
    
end module
