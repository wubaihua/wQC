module def
    type atomtype !自定义的记录原子信息的类型
    character*2 name !原子名
    integer index !原子序号。当使用了赝势时原子核电荷将小于原子序号
    real*8 x,y,z,charge !原子坐标(Bohr)及原子核电荷
    end type

    type primtype !自定义的记录GTF信息的专用类型
    integer center,functype !GTF所属原子以及GTF的类型
    real*8,allocatable:: e(:)!指数

    real*8,allocatable::  C(:,:)
    end type

character*2 :: name2ind(0:109)=(/ "Bq","H ","He", &   !原子名和原子序号的转换表。0号是虚原子
"Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
"Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
"Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
"Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
"Rf","Db","Sg","Bh","Hs","Mt" /) !104~109


character*2 :: shellindex(1:6)=(/'S ','P ','SP','D','F','G'/)

!下面三行是GTF类型与GTF的x,y,z上的指数的转换表。例如XY型GTF的x,y,z的指数分别是1,1,0，和下面三行的第8列对应，所以这样的GTF的functype是8。
integer :: type2ix(35)=(/ 0,1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1, 0,0,0,0,0,1,1,1,1,2,2,2,3,3,4 /)
integer :: type2iy(35)=(/ 0,0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1, 0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
integer :: type2iz(35)=(/ 0,0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1, 4,3,2,1,0,3,2,1,0,2,1,0,1,0,0 /)
integer :: nmo=0,nprims=0,ncenter=0 !轨道数、GTF数、原子数
integer wfntype !0/1/2代表波函数是R/U/ROHF波函数，3/4分别代表闭壳层和开壳层后HF波函数
real*8 :: nelec=0,naelec=0,nbelec=0 !总电子数、alpha和beta电子数
type(atomtype),allocatable :: a(:) !记录原子信息的数组
type(primtype),allocatable :: b(:) !记录GTF信息的数组
real*8,allocatable :: MOocc(:),MOene(:) !轨道占据数和轨道能量
integer,allocatable :: MOtype(:) !记录轨道类型。0/1/2分别代表无自旋/alpha/beta型轨道
!real*8,allocatable :: CO(:,:) !轨道展开系数矩阵。CO(i,j)代表第j个GTF在第i号轨道中的展开系数，系数中已经把收缩系数、归一化系数全都包含进去了。
!index of shell : 0--S,  1--P,  2--D,  3--F     
integer :: shl_degen_sph(0:3)=(/1,3,5,7/) !Degeneracy of shell in Spherical harmonic Gaussian function(5D,7F) 
integer :: shl_degen_cart(0:3)=(/1,3,6,10/) !Degeneracy of shell in Cartesian Gaussian function(6D,10F) 

contains


 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end module
