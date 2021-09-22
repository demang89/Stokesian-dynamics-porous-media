module HI
implicit none
real(8), parameter :: pi=3.1415926d0
integer :: dm, natom     !dimension, total no. of particles in simulation box
integer :: n_x, n_y, n_z !No. of nanopost particles in each dimension
integer, allocatable :: ngrid(:) !No. of grids for resistance matrix calculation in each dimension
integer, allocatable :: rn(:) !No. of real space vectors for resistance matrix calculation 
integer, allocatable :: kn(:) !No. of reciprocal space vectors for resistance matrix calculation
real(8), allocatable :: cell_size(:) ! Unit cell size in each dimension
real(8), allocatable :: box_size(:) !Simulation box size in each dimension
real(8), allocatable :: ibox_size(:) ! Inverse of simulation box size
real(8), allocatable :: gridsize(:) !Gridsize for resistance matrix calculation
real(8), allocatable :: igridsize(:) !Inverse of gridsize
real(8) :: alpha, rc, kc !Ewald sum parameters
real(8) :: zeta, s !Confinement parameter, s-s distance between nanoposts
real(8) :: rnp, rp !nanopost and particle radii
real(8) :: gx1,gx2,gx3,fx0,fx1,fx2,fx3,fx4,fx5,fx6,fx7,fx8,fx9,fx10,fx11 !Coefficients of R11 matrix
real(8) :: gy2,gy3,fy0,fy1,fy2,fy3,fy4,fy5,fy6,fy7,fy8,fy9,fy10,fy11 !Coefficients of R12 matrix
real(8) :: c11_0, c11_1, c11_2, c11_3, c11_4, c11_5, c11_6 !Fitting parameters for inverse M11 matrix
real(8) :: d11_0, d11_1, d11_2, d11_3, d11_4, d11_5, d11_6 !Fitting parameters for inverse M11 matrix
real(8) :: c12_0, c12_1, c12_2, c12_3, c12_4, c12_5, c12_6 !Fitting parameters for inverse M12 matrix
real(8) :: d12_0, d12_1, d12_2, d12_3, d12_4, d12_5, d12_6 !Fitting parameters for inverse M12 matrix
real(8), allocatable :: rpos(:,:) !Nanopost particle position array
real(8) :: vmax, theta !fluid velocity magnitude and direction angle in degree 
integer(8) :: nmax, t_grid !No. of integration timesteps
real(8) :: dt !integration timestep value
logical :: longrange, shortrange !logical parameters for inclusion of long-range and short-range HI
contains
!########################################################################

!Subroutine to read model data and nanopost particle positions
subroutine read_data()
implicit none
integer :: i
open(unit=14,file='post_pos.txt',action='read')
read(14,'(2f20.8)') rp, rnp
read(14,'(2I10)') dm, natom

allocate(rpos(dm,natom-1))
allocate(cell_size(dm))
allocate(box_size(dm))
allocate(ngrid(dm))
allocate(gridsize(dm))
allocate(rn(dm))
allocate(kn(dm))
allocate(ibox_size(dm))
allocate(igridsize(dm))

read(14,*) cell_size(:)
read(14,*) box_size(:)
read(14,*) ngrid(:)
read(14,*) gridsize(:)
read(14,*) alpha, rc, kc
read(14,*) rn(:)
read(14,*) kn(:)

do i = 1, natom-1
        read(14,*) rpos(:,i)
enddo
ibox_size = 1.d0 / box_size
igridsize = 1.d0 / gridsize
return
end subroutine read_data
!##########################################################################

!Subroutine to read fitting parameters of inverse of pair-wise mobility matrix
subroutine read_m2p_coeff()
implicit none
open(unit=14,file='r11_coeff.txt',action='read')
read(14,*) c11_0
read(14,*) c11_1
read(14,*) c11_2
read(14,*) c11_3
read(14,*) c11_4
read(14,*) c11_5
read(14,*) c11_6

read(14,*) d11_0
read(14,*) d11_1
read(14,*) d11_2
read(14,*) d11_3
read(14,*) d11_4
read(14,*) d11_5
read(14,*) d11_6

close(14)

open(unit=14,file='r12_coeff.txt',action='read')
read(14,*) c12_0
read(14,*) c12_1
read(14,*) c12_2
read(14,*) c12_3
read(14,*) c12_4
read(14,*) c12_5
read(14,*) c12_6

read(14,*) d12_0
read(14,*) d12_1
read(14,*) d12_2
read(14,*) d12_3
read(14,*) d12_4
read(14,*) d12_5
read(14,*) d12_6

close(14)

return
end subroutine read_m2p_coeff
!#####################################################################

!Subroutine for pair-wise Rotne_prager matrix
subroutine Rotneprager(mij,rij,a1,a2)
implicit none
real(8),intent(in) :: rij(dm), a1, a2
real(8),intent(out) :: mij(dm,dm)
real(8) :: a3, dr1, dr2, idr1, idr3, udr(dm)
real(8) :: rmat(dm,dm), umat(dm,dm)
integer :: i
               
umat=0.d0
forall(i=1:dm) umat(i,i)=1.d0

a3=a1*(0.5d0*(a1**2+a2**2))

dr2=sum(rij*rij)

if(dr2/=0.d0) then
        dr1=dsqrt(dr2)
        idr1=1.d0/dr1;  idr3=idr1**3
        udr=rij*idr1
        rmat=spread(udr,dim=2,ncopies=dm)*spread(udr,dim=1,ncopies=dm)
        mij=0.75d0*a1*idr1*(umat+rmat) +0.5d0*a3*idr3*(umat-3.d0*rmat)
else
        mij=umat
endif
end subroutine Rotneprager
!#############################################################################

!Subroutine for 3D pair-wise far-field mobility matrix using Ewald sum
subroutine RPEwald3d(mij,rij,a1,a2)
implicit none
real(8),intent(in) :: rij(3),a1,a2
real(8),intent(out) :: mij(3,3)
real(8) :: a3,e1, e2, e3, e4, e5, e7, v1, v2, volume, rc2, kc2,isqrpi
real(8) :: const, Kpre1, Kpre2, Kpre3, dr1, dr2, dr4, rmat(3,3), umat(3,3)
real(8) :: udr(3),idr1,idr2,idr3,var1,var2,var3,var4, var5, var6, rnij(3),kij(3)
integer :: i,j,k

rc2=rc*rc; kc2=kc*kc

isqrpi=pi**(-0.5d0)

umat=reshape((/1,0,0,0,1,0,0,0,1/),shape(umat))
Kpre1 = 2.d0 * pi / box_size(1) 
Kpre2 = 2.d0 * pi / box_size(2)
Kpre3 = 2.d0 * pi / box_size(3)

a3=a1*(0.5d0*(a1**2+a2**2))
e1=alpha; e2=e1*e1; e3=e2*e1; e4=e3*e1; e5=e4*e1; e7=e5*e2
const = 1.d0 - 6.d0*isqrpi*e1*a1 + 40.d0/3.d0*isqrpi*e3*a1**3
mij=0.d0
volume = box_size(1) * box_size(2) * box_size(3)
do i = -rn(1), rn(1)
        do j = -rn(2), rn(2)
                do k = -rn(3), rn(3)
                        rnij(1)=dble(i)*box_size(1)+rij(1)
                        rnij(2)=dble(j)*box_size(2)+rij(2)
                        rnij(3)=dble(k)*box_size(3)+rij(3)
                        dr2=sum(rnij*rnij)
                        if(dr2/=0.0 .and. dr2<=rc2) then
                                dr1=dsqrt(dr2); dr4=dr2*dr2
                                idr1=1.0/dr1; idr2=idr1**2; idr3=idr1*idr2
                                udr=rnij*idr1
                                rmat=spread(udr,dim=2,ncopies=3)*spread(udr,dim=1,ncopies=3)
                                var1=erfc(e1*dr1)
                                var2=exp(-e2*dr2)
                                var3=0.75d0*a1*idr1+0.5d0*a3*idr3
                                var4=(4.d0*e7*a3*dr4 + 3.d0*e3*a1*dr2-20.d0*e5*a3*dr2-4.5d0*e1*a1 + 14.d0*e3*a3 + e1*a3*idr2)*isqrpi
                                var5 = 0.75d0*a1*idr1-1.5d0*a3*idr3
                                var6=(-4.d0*e7*a3*dr4-3.d0*e3*a1*dr2+16.d0*e5*a3*dr2+1.5d0*e1*a1-2.d0*e3*a3 -3.d0*e1*a3*idr2)*isqrpi
                                v1=var3*var1 + var4*var2
                                v2=var5*var1 + var6*var2
                                mij = mij + (v1*umat + v2*rmat)
                        endif
                enddo
        enddo
enddo

do i = -kn(1), kn(1)
        do j = -kn(2), kn(2)
                do k = -kn(3), kn(3)
                        kij(1)=dble(i)*Kpre1
                        kij(2)=dble(j)*Kpre2
                        kij(3)=dble(k)*Kpre3
                        dr2=sum(kij*kij)
                        if(dr2/=0.d0 .and. dr2<=kc2) then
                                dr1=dsqrt(dr2)
                                idr1=1.0/dr1
                                udr=kij*idr1
                                rmat=spread(udr,dim=2,ncopies=3)*spread(udr,dim=1,ncopies=3)
                                var1=0.25d0*dr2/e2
                                var2=exp(-var1)
                                var3 = a1-a3*dr2/3.d0
                                var4 = 1.d0 + var1 + 2.d0*var1*var1
                                var5=6.d0*pi/dr2/volume
                                var6=dcos(dot_product(kij,rij))
                                mij = mij +((umat-rmat)*var2*var3*var4*var5*var6)
                        endif
                enddo
        enddo
enddo
if(sum(rij*rij)==0.d0) then
        mij = mij + const*umat
endif
return
end subroutine RPEwald3d
!################################################################################

!Subroutine for 2D pair-wise far-field mobility matrix using Ewald sum
subroutine RPEwald2d(mij,rij,a1,a2)
implicit none
real(8),intent(in) :: rij(2),a1,a2
real(8),intent(out) :: mij(2,2)
real(8) :: a3,e1, e2, e3, v1, v2, rc2, kc2,isqrpi, area
real(8) :: const, Kpre1,Kpre2, dr1, dr2, dr4, rmat(2,2), umat(2,2),kij(2)
real(8) :: udr(2),idr1,idr2,idr3,var1,var2,var3,var4, var5, var6, rnij(2)
integer :: i,j,k

isqrpi = pi**(-0.5d0)
rc2=rc*rc; kc2=kc*kc

umat=reshape((/1,0,0,1/),shape(umat))

Kpre1 = 2.d0 * pi / box_size(1)
Kpre2 = 2.d0 * pi / box_size(2)

a3=a1*(0.5d0*(a1**2+a2**2))

e1=alpha; e2=e1*e1; e3=e2*e1

const = 1.0d0 - 1.5d0 *isqrpi * a1 * e1 - 2.d0/3.d0 * isqrpi *a1**3 *e3

mij=0.d0!-3.d0*a1/e1/isqrpi/L**2*umat

area = box_size(1) * box_size(2)
do i = -rn(1), rn(1)
        do j = -rn(2), rn(2)
                rnij(1)=dble(i)*box_size(1)+rij(1)
                rnij(2)=dble(j)*box_size(2)+rij(2)
                dr2=sum(rnij*rnij)
                if(dr2/=0.d0 .and. dr2<=rc2) then
                        dr1=dsqrt(dr2)
                        idr1=1.0/dr1; idr2=idr1**2; idr3=idr1*idr2
                        udr=rnij*idr1
                        rmat=spread(udr,dim=2,ncopies=2)*spread(udr,dim=1,ncopies=2)
                        var1=erfc(e1*dr1)
                        var2=exp(-e2*dr2)
                        var3=0.75d0*a1*idr1+0.5d0*a3*idr3
                        var4= e1*a3*idr2 *isqrpi
                        var5 = 0.75d0*a1*idr1-1.5d0*a3*idr3
                        var6=(1.5d0*e1*a1 -2.d0*e3*a3 -3.d0*e1*a3*idr2)*isqrpi
                        v1=var3*var1 + var4*var2
                        v2=var5*var1 + var6*var2
                        mij = mij + (v1*umat + v2*rmat)
                endif
        enddo
enddo

do i = -kn(1), kn(1)
        do j = -kn(2), kn(2)
                kij(1)=dble(i)*Kpre1
                kij(2)=dble(j)*Kpre2
                dr2=sum(kij*kij)
                if(dr2/=0.d0 .and. dr2<=kc2) then
                        dr1=dsqrt(dr2)
                        idr1=1.0/dr1
                        udr=kij*idr1
                        rmat=spread(udr,dim=2,ncopies=2)*spread(udr,dim=1,ncopies=2)
                        var1=0.5d0*dr1/e1
                        var2=erfc(var1)
                        var3=exp(-var1**2)
                        var4=1.5d0*pi*idr1/area
                        var5=dcos(dot_product(rij,kij))
                        v1= 2.0d0*a1*var2
                        v2=(a1-2.d0/3.d0*a3*dr2)*var2 +a1*dr1/e1*isqrpi*var3
                        mij = mij + var4*var5*(v1*umat - v2*rmat)
                endif
        enddo
enddo
if(sum(rij*rij)==0.d0) then
        mij = mij + const*umat
endif
return
end subroutine RPEwald2d
!#########################################################################

!Subroutine for matrix intversion using cholesky decomposition technique
subroutine matrix_inv(mpbc,rpbc)
real(8), intent(in)  :: mpbc(dm*natom,dm*natom)
real(8), intent(out) :: rpbc(dm*natom,dm*natom)
integer :: i, j
rpbc=0.d0

do i=1,dm*natom
        do j=1,i
                if(j==i) then
                        rpbc(i,j)=dsqrt(mpbc(i,j)-sum(rpbc(j,1:j-1)**2))
                else
                        rpbc(i,j)=(mpbc(i,j)-sum(rpbc(j,1:j-1)*rpbc(i,1:j-1)))/rpbc(j,j)
                endif
        enddo
enddo

do i=1,dm*natom
        rpbc(i,i)=1.d0/rpbc(i,i)
        do j=i+1,dm*natom
                rpbc(j,i)=-sum(rpbc(j,i:j-1)*rpbc(i:j-1,i))/rpbc(j,j)
        enddo
enddo
rpbc=matmul(transpose(rpbc),rpbc)
return
end subroutine matrix_inv
!#############################################################################

!Subroutine for initial random seed
SUBROUTINE init_random_seed()
implicit none
integer :: i, n, clock
integer, allocatable:: seed(:)
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE init_random_seed
!###############################################################################

!Subroutine to compute coefficients of pair-wise resistance matrix
subroutine lubcoefficient()
implicit none
real(8) :: la
real(8) :: la2, la3, la4, la5, la6, la7, la8, la9, la10
la = rnp / rp
la2 = la*la; la3 = la2*la; la4 = la3*la; la5 = la4*la
la6 = la5*la; la7 = la6*la; la8 = la7*la; la9 = la8*la; la10 = la9*la
gx1 = 2.d0 * la2 / (1.d0+la)**3
gx2 = la * (1.d0+7.d0*la+la2) / 5.d0 / (1.d0+la)**3
gx3 = (1.d0 + 18.d0*la - 29.0d0*la2 + 18.d0*la3 + la4) / (1.d0+la)**3 / 42.d0
fx0 = 1.d0; fx1 = 3.d0*la; fx2 = 9.d0*la
fx3 = -4.d0*la + 27.d0*la2 - 4.d0*la3
fx4 = -24.d0*la + 81.d0*la2 + 36.d0*la3
fx5 = 72.d0*la2 + 243.d0*la3 + 72.d0*la4
fx6 = 16.d0*la + 108.d0*la2 + 281.d0*la3 + 648.d0*la4 + 144.d0*la5
fx7 = 288.d0*la2 + 1620.d0*la3 + 1515.d0*la4 + 1620.d0*la5 + 288.d0*la6
fx8 = 576.d0*la2 + 4848.d0*la3 + 5409.d0*la4 + 4524.d0*la5 + 3888.d0*la6 + 576.d0*la7
fx9 = 1152.d0*la2 + 9072.d0*la3 + 14752.d0*la4 + 26163.d0*la5 + 14752.d0*la6 + 9072.d0*la7 + 1152.d0*la8
fx10 = 2304*la2 + 20736*la3 + 42804*la4 + 115849*la5 + 76176*la6 + 39264*la7 + 20736*la8 + 2304*la9
fx11 = 4608*la2 + 46656*la3 + 108912*la4 + 269100*la5 + 319899*la6 + 269100*la7 + 108912*la8 + 46656*la9 + 4608*la10


gy2 = 4.d0/15.d0 * la * (2.d0+la+2.d0*la2) / (1.d0+la)**3
gy3 = 2.d0/375.d0 * (16.d0-45.d0*la + 58.d0*la2 - 45.d0*la3 + 16.d0*la4) / (1.d0+la)**3
fy0 = 1.d0; fy1 = 1.5d0*la; fy2 = 9.d0/4.d0*la
fy3 = 2.d0*la + 27.d0/8.d0*la2 + 2.d0*la3
fy4 = 6.d0*la + 81.d0/16.d0*la2 + 18.d0*la3
fy5 = 63.d0/2.d0*la2 + 243.d0/32.d0*la3 + 63.d0/2.d0*la4
fy6 = 4.d0*la + 54.d0*la2 + 1241.d0/64.d0*la3 + 81.d0*la4 + 72.d0*la5
fy7 = 144.d0*la2 + 1053.d0/8.d0*la3 + 19083.d0/128.d0*la4 + 1053.d0/8.d0*la5 + 144.d0*la6
fy8 = 279.d0*la2 + 4261.d0/8.d0*la3 + 126369.d0/256.d0*la4 - 117.d0/8.d0*la5 + 648.d0*la6 + 288.d0*la7
fy9 = 576.d0*la2 + 1134.d0*la3 + 60443.d0/32.d0*la4 + 766179.d0/512.d0*la5 + 60443.d0/32.d0*la6 + 1134.d0*la7 + 576*la8
fy10= 1152.0*la2 + 7857.0/4.0*la3 + 98487/16*la4 + 10548393/1024*la5 + 67617/8*la6 - 351/2*la7 + 3888*la8 + 1152*la9
fy11= 2304.d0*la2 + 7128.d0*la3 + 22071.d0/2.d0*la4 + 2744505.d0/128.d0*la5 &
        & + 95203835.d0/2048*la6 + 2744505.d0/128.d0*la7 + 22071.d0/2.d0*la8 + 7129.d0*la9 + 2304.d0*la10

fx1 = fx1/((1.d0+la)*2.d0)**1 - gx1 - 2.d0*gx2 - 4.d0*gx3
fx2 = fx2/((1.d0+la)*2.d0)**2 - gx1 - 2.d0*gx2/2.d0 - 4.d0*gx3/2.d0/2.d0
fx3 = fx3/((1.d0+la)*2.d0)**3 - gx1 - 2.d0*gx2/3.d0 + 4.d0*gx3/3.d0
fx4 = fx4/((1.d0+la)*2.d0)**4 - gx1 - 2.d0*gx2/4.d0 + 4.d0*gx3/4.d0/2.d0
fx5 = fx5/((1.d0+la)*2.d0)**5 - gx1 - 2.d0*gx2/5.d0 + 4.d0*gx3/5.d0/3.d0
fx6 = fx6/((1.d0+la)*2.d0)**6 - gx1 - 2.d0*gx2/6.d0 + 4.d0*gx3/6.d0/4.d0
fx7 = fx7/((1.d0+la)*2.d0)**7 - gx1 - 2.d0*gx2/7.d0 + 4.d0*gx3/7.d0/5.d0
fx8 = fx8/((1.d0+la)*2.d0)**8 - gx1 - 2.d0*gx2/8.d0 + 4.d0*gx3/8.d0/6.d0
fx9 = fx9/((1.d0+la)*2.d0)**9 - gx1 - 2.d0*gx2/9.d0 + 4.d0*gx3/9.d0/7.d0
fx10 = fx10/((1.d0+la)*2.d0)**10 - gx1 - 2.d0*gx2/10.d0 + 4.d0*gx3/10.d0/8.d0
fx11 = fx11/((1.d0+la)*2.d0)**11 - gx1 - 2.d0*gx2/11.d0 + 4.d0*gx3/11.d0/9.d0

fy1 = fy1/((1.d0+la)*2.d0)**1 - 2.d0*gy2/1.d0 - 4.d0*gy3
fy2 = fy2/((1.d0+la)*2.d0)**2 - 2.d0*gy2/2.d0 - 4.d0*gy3/2.d0/2.d0
fy3 = fy3/((1.d0+la)*2.d0)**3 - 2.d0*gy2/3.d0 + 4.d0*gy3/3.d0
fy4 = fy4/((1.d0+la)*2.d0)**4 - 2.d0*gy2/4.d0 + 4.d0*gy3/4.d0/2.d0
fy5 = fy5/((1.d0+la)*2.d0)**5 - 2.d0*gy2/5.d0 + 4.d0*gy3/5.d0/3.d0
fy6 = fy6/((1.d0+la)*2.d0)**6 - 2.d0*gy2/6.d0 + 4.d0*gy3/6.d0/4.d0
fy7 = fy7/((1.d0+la)*2.d0)**7 - 2.d0*gy2/7.d0 + 4.d0*gy3/7.d0/5.d0
fy8 = fy8/((1.d0+la)*2.d0)**8 - 2.d0*gy2/8.d0 + 4.d0*gy3/8.d0/6.d0
fy9 = fy9/((1.d0+la)*2.d0)**9 - 2.d0*gy2/9.d0 + 4.d0*gy3/9.d0/7.d0
fy10 = fy10/((1.d0+la)*2.d0)**10  - 2.d0*gy2/10.d0 + 4.d0*gy3/10.d0/8.d0
fy11 = fy11/((1.d0+la)*2.d0)**11  - 2.d0*gy2/11.d0 + 4.d0*gy3/11.d0/9.d0
return
end subroutine lubcoefficient
!########################################################################
!Subroutine to generate gaussian random number
real(8) function randomnormal()
integer :: iset
real(8) :: u1, u2, rsq, first,second, frac
save iset
save second
data iset/0/
if(iset.eq.0) then
        rsq=0.0d0
        do while((rsq >= 1.0d0) .or. (rsq == 0.0d0))
                call random_number(u1)
                call random_number(u2)
                u1=2.0d0*u1-1.0d0
                u2=2.0d0*u2-1.0d0
                rsq=u1*u1+u2*u2
        enddo
        frac = dsqrt(-2.0d0 * dlog(rsq)/rsq)
        second = frac * u1
        first = frac * u2
        iset = 1
else
        first=second
        iset=0
endif
randomnormal=first
return
end function randomnormal
!############################################################################

!Subroutine to compute R11 component (R^{2P}_{11}-(M^{RP}_{11})^{-1}) of lubrication correction matrix
subroutine lub11correction(rij,rlub,a1,a2)
implicit none
real(8), intent(in) :: rij(dm),a1,a2
real(8), intent(out) :: rlub(dm,dm)
real(8) :: csum, dsum, dr, idr, rmat(dm,dm), umat(dm,dm), udr(dm)
real(8) :: var1,var2,var3, rtot, rmin, tsq, tsh
real(8) :: idr1, idr2, idr3, idr4, idr5, idr6, p1, p2
real(8) :: ids1, ids2, ids4, ids6, ids8, ids10,q1
integer :: i

rtot=a1+a2
rmin=rtot+0.0001d0

dr=dsqrt(sum(rij*rij))
udr=rij/dr
rmat=spread(udr,dim=2,ncopies=dm)*spread(udr,dim=1,ncopies=dm)

umat=0.d0
forall(i=1:dm) umat(i,i)=1.d0

if(dr<=rmin) then
dr=rmin
endif

idr1=1.0d0/dr; idr2=idr1**2; idr3=idr1*idr2; idr4=idr1*idr3;
idr5=idr1*idr4; idr6=idr1*idr5
csum = c11_0 + c11_1*idr1 + c11_2*idr2 + c11_3*idr3 + c11_4*idr4 + c11_5*idr5 + c11_6*idr6
dsum = d11_0 + d11_1*idr1 + d11_2*idr2 + d11_3*idr3 + d11_4*idr4 + d11_5*idr5 + d11_6*idr6

ids1=rtot/dr; ids2=ids1**2; ids4=ids2*ids2; ids6=ids2*ids4
ids8=ids2*ids6; ids10=ids2*ids8

var1=1.0d0-ids2; var2=dlog(var1)

p1=gx1/var1 - gx2*var2 - gx3*var1*var2 + fx0 - gx1
tsq = (p1 + fx2*ids2 + fx4*ids4 + fx6*ids6 + fx8*ids8 + fx10*ids10)

p2 = -gy2*var2 - gy3*var1*var2  + fy0
tsh = (p2 + fy2*ids2 + fy4*ids4 + fy6*ids6 + fy8*ids8 + fy10*ids10)

rlub = (a1*(tsq-tsh)-dsum)*rmat + (a1*tsh-csum)*umat
return
end subroutine lub11correction
!##########################################################################

!Subroutine to compute R12 component (R^{2P}_{12}-(M^{RP}_{12})^{-1}) of lubrication correction matrix
subroutine lub12correction(rij,r12,a1,a2)
implicit none
real(8), intent(in) :: rij(dm),a1,a2
real(8), intent(out) :: r12(dm,dm)
real(8) :: csum, dsum, dr, idr, rmat(dm,dm), umat(dm,dm), udr(dm)
real(8) :: var1,var2,var3, rtot, rmin, tsq, tsh
real(8) :: idr1, idr2, idr3, idr4, idr5,idr6, p1, p2, p3, p4
real(8) :: ids1, ids2, ids3, ids5, ids7, ids9, ids11
integer :: i


rtot = a1+a2
rmin = rtot+0.0001d0

dr=dsqrt(sum(rij*rij))
udr=rij/dr
rmat=spread(udr,dim=2,ncopies=dm)*spread(udr,dim=1,ncopies=dm)


if(dr<=rmin) then
dr=rmin
endif

idr1=1.0d0/dr;  idr2=idr1**2; idr3=idr1*idr2; idr4=idr1*idr3;
idr5=idr1*idr4; idr6=idr5*idr1

csum = c12_0 + c12_1*idr1 + c12_2*idr2 + c12_3*idr3 + c12_4*idr4 + c12_5*idr5 + c12_6*idr6
dsum = d12_0 + d12_1*idr1 + d12_2*idr2 + d12_3*idr3 + d12_4*idr4 + d12_5*idr5 + d12_6*idr6

umat=0.d0
forall(i=1:dm) umat(i,i)=1.d0

ids1=rtot/dr; ids2=ids1**2; ids3=ids1*ids2; ids5=ids2*ids3
ids7=ids2*ids5; ids9=ids2*ids7; ids11=ids2*ids9

var1=1.0d0-ids2; var3=dlog((1.0d0+ids1)/(1.0d0-ids1))

p1 = ids1*gx1/var1 + gx2*var3 + gx3*var1*var3 + 2.0*gx3*ids1
tsq = -(p1 + fx1*ids1 + fx3*ids3 + fx5*ids5 + fx7*ids7 + fx9*ids9 + fx11*ids11)

p2 = gy2*var3 + gy3*var1*var3 + 2.0*gy3*ids1
tsh = -(p2 + fy1*ids1 + fy3*ids3 + fy5*ids5 + fy7*ids7 + fy9*ids9 + fy11*ids11)

r12 = (a1*(tsq-tsh)-dsum)*rmat+ (a1*tsh-csum)*umat
return
end subroutine lub12correction
!##################################################################

!Subroutine to compute lower triangular and inverse matrix
subroutine choleskydecomposition(rlub,lt,mij)
implicit none
real(8), intent(in) :: rlub(dm,dm)
real(8), intent(out) :: lt(dm,dm),mij(dm,dm)
integer :: i, j
lt=0.d0
do i=1,dm
        do j=1,i
                if(j==i) then
                        lt(i,j)=dsqrt(rlub(i,j)-sum(lt(j,1:j-1)**2))
                else
                        lt(i,j)=(rlub(i,j)-sum(lt(j,1:j-1)*lt(i,1:j-1)))/lt(j,j)
                endif
        enddo
enddo
mij=lt
do i=1,dm
        mij(i,i)=1.0/mij(i,i)
        do j=i+1,dm
                mij(j,i)=-sum(mij(j,i:j-1)*mij(i:j-1,i))/mij(j,j)
        enddo
enddo
mij=matmul(transpose(mij),mij)
return
end subroutine choleskydecomposition
!##################################################################

!Subroutine for trilinear interpolation
subroutine trilinearinterpolation(r_grid,r_frac,resm,r11)
implicit none
integer,intent(in) :: r_grid(dm)
real(8), intent(in) :: resm(0:ngrid(1), 0:ngrid(2), 0:ngrid(3), 3, 3)
real(8), intent(in) :: r_frac(dm)
real(8), intent(out) :: r11(3,3)
real(8) :: c000(3,3), c100(3,3), c001(3,3), c101(3,3), c010(3,3), c110(3,3),c011(3,3), c111(3,3)
real(8) :: c00(3,3), c10(3,3), c01(3,3), c11(3,3), c1(3,3), c0(3,3)
c000 = resm(r_grid(1)  , r_grid(2)  , r_grid(3)  ,:,:)
c100 = resm(r_grid(1)+1, r_grid(2)  , r_grid(3)  ,:,:)
c010 = resm(r_grid(1)  , r_grid(2)+1, r_grid(3)  ,:,:)
c001 = resm(r_grid(1)  , r_grid(2)  , r_grid(3)+1,:,:)
c110 = resm(r_grid(1)+1, r_grid(2)+1, r_grid(3)  ,:,:)
c101 = resm(r_grid(1)+1, r_grid(2)  , r_grid(3)+1,:,:)
c011 = resm(r_grid(1)  , r_grid(2)+1, r_grid(3)+1,:,:)
c111 = resm(r_grid(1)+1, r_grid(2)+1, r_grid(3)+1,:,:)
c00 = (1.d0-r_frac(1))*c000 + r_frac(1)*c100
c01 = (1.d0-r_frac(1))*c001 + r_frac(1)*c101
c10 = (1.d0-r_frac(1))*c010 + r_frac(1)*c110
c11 = (1.d0-r_frac(1))*c011 + r_frac(1)*c111
c0  = (1.d0-r_frac(2))*c00  + r_frac(2)*c10
c1  = (1.d0-r_frac(2))*c01  + r_frac(2)*c11
r11 = (1.d0-r_frac(3))*c0   + r_frac(3)*c1
return
end subroutine trilinearinterpolation
!####################################################################

!Subroutine for bilinear interpolation
subroutine bilinearinterpolation(r_grid,r_frac,resm,r11)
implicit none
integer,intent(in) ::  r_grid(dm)
real(8), intent(in) :: resm(0:ngrid(1), 0:ngrid(2),2,2)
real(8), intent(in) :: r_frac(dm)
real(8), intent(out) :: r11(2,2)
real(8) :: c00(2,2), c10(2,2), c01(2,2), c11(2,2), c1(2,2), c0(2,2)
c00 = resm(r_grid(1)  , r_grid(2)  ,:, :)
c10 = resm(r_grid(1)+1, r_grid(2)  ,:, :)
c01 = resm(r_grid(1)  , r_grid(2)+1,:, :)
c11 = resm(r_grid(1)+1, r_grid(2)+1,:, :)
c0  = (1.d0-r_frac(1))*c00 + r_frac(1)*c10
c1  = (1.d0-r_frac(1))*c01 + r_frac(1)*c11
r11 = (1.d0-r_frac(2))*c0  + r_frac(2)*c1
return
end subroutine bilinearinterpolation
!######################################################################1

!Subroutine to compute overlapping of diffusing particle with nanopost particles
subroutine periodic_HSInteraction(ri,rtot2,dhs2)
implicit none
real(8),intent(in) :: ri(dm), rtot2
real(8),intent(out) :: dhs2
real(8) :: r_temp(dm), rij(dm), dr2
dhs2 = 1.d0
r_temp = modulo(ri, cell_size)
rij = r_temp - 0.5d0*cell_size
dr2 = dot_product(rij,rij)
if(dr2 .lt. rtot2) then
        dhs2 = -1.d0
endif
return
end subroutine periodic_HSInteraction
!##########################################################

!Apply periodic boundary conditions
subroutine pbc(rij)
implicit none
real(8), intent(inout) :: rij(dm)
rij = rij - box_size * dnint(rij * ibox_size)
end subroutine pbc
!#########################################################

!Compute grid location for far-field resistance matrix
subroutine grid_location(ri,r_grid,r_frac)
implicit none
real(8), intent(in) :: ri(dm)
integer, intent(out) :: r_grid(dm)
real(8), intent(out) :: r_frac(dm)
real(8) :: rij(dm)
rij = modulo(ri, cell_size)
r_grid = idint(rij * igridsize)
r_frac = (rij - dble(r_grid) * gridsize) * igridsize
return
end subroutine grid_location
!################################################################
end module HI
