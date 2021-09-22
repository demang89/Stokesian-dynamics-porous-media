program main
use HI
implicit none
integer :: i, j, k, p, m, n
real(8) :: irp, irnp
real(8), allocatable :: mij(:,:), rsum(:,:), rself(:,:), rij(:), mpbc(:,:), rpbc(:,:)
character*20 :: filename

!Read nanopost particle positions and other data
call read_data()

irp = 1.d0/rp; irnp=1.d0/rnp !Inverse of particle and nanopost radii

allocate(rij(dm)) !Allocate size for separation vector
allocate(mij(dm,dm)) !Allocate size for pair-wise far-field mobility matrix
allocate(rsum(dm,dm)) !Allocate size for Rsum matrix
allocate(rself(dm,dm)) !Allocate size for Rself matrix
allocate(mpbc(dm*natom,dm*natom)) !Allocate size for grand farfield mobility matrix
allocate(rpbc(dm*natom,dm*natom)) !Allocate size for grand resistance matrix

!Initialize grand matrix to zero
mpbc = 0.d0
rpbc = 0.d0

!Open file to write gridwise Rsum and Rself matrix
open(unit=28, file='rself.txt', action='write')
open(unit=56, file='rsum.txt', action='write')

!If model dimension in 3D
if(dm==3) then
        rij=0.0d0
        !Calculate particle self mobility matrix and feed into grand mobility matrix
        call RPEwald3d(mij,rij,rp,rp)
        mpbc(1:dm,1:dm)=mij*irp
        !Calculate nanopost particles self mobility matrix and feed into grand mobility matrix
        call RPEwald3d(mij,rij,rnp,rnp)
        do i=2,natom
                mpbc(dm*(i-1)+1:dm*i, dm*(i-1)+1:dm*i) = mij*irnp
        enddo
        !Compute pairwise mobility matrix for nanopost particles and feed into grand mobility matrix 
        do i=2,natom
                do j=i+1,natom
                        rij=rpos(:,i-1)-rpos(:,j-1)
                        rij = rij - box_size * dnint(rij*ibox_size)
                        call RPEwald3d(mij,rij,rnp,rnp)
                        mpbc(dm*(i-1)+1:dm*i,dm*(j-1)+1:dm*j) = mij*irnp
                        mpbc(dm*(j-1)+1:dm*j,dm*(i-1)+1:dm*i) = mij*irnp
                enddo
        enddo
        !Compute Rsum and Rself matrix on grid points
        !It will generate nan values at grid points overlapping with nanopost
        do i = 0, ngrid(1)
                do j = 0, ngrid(2)
                        do k = 0, ngrid(3)
                                !Compute pairwise mobility matrix for particle
                                !and nanopost particle and feed into grand mobility matrix
                                do p = 2, natom
                                        rij(1) = dble(i) * gridsize(1) - rpos(1,p-1)
                                        rij(2) = dble(j) * gridsize(2) - rpos(2,p-1)
                                        rij(3) = dble(k) * gridsize(3) - rpos(3,p-1)
                                        rij = rij - box_size * dnint(rij*ibox_size)
                                        call RPEwald3d(mij,rij,rp,rnp)
                                        mpbc(1:dm,dm*(p-1)+1:dm*p) = mij*irp
                                        mpbc(dm*(p-1)+1:dm*p,1:dm) = mij*irp
                                enddo
                                !Compute inverse of grand mobility matrix
                                call matrix_inv(mpbc,rpbc)
                                !Compue Rsum and Rself matrix
                                do m=1,dm
                                        do n=1,dm
                                                rself(m,n)=rpbc(m,n)
                                                rsum(m,n)=sum(rpbc(m,n:dm*natom:dm))
                                        enddo
                                enddo
                                !Write Rsum and Rself matrix to files
                                write(28,'(9f16.5)')rself(:,:)
                                write(56,'(9f16.5)')rsum(:,:)
                        enddo
                enddo
        enddo
!If model dimension is 2D
elseif(dm==2) then
        rij=0.0d0
        !Calculate particle self mobility matrix and feed into grand mobility matrix
        call RPEwald2d(mij,rij,rp,rp)
        mpbc(1:dm,1:dm)=mij*irp
        !Calculate nanopost particles self mobility matrix and feed into grand mobility matrix
        call RPEwald2d(mij,rij,rnp,rnp)
        do i=2,natom
                mpbc(dm*(i-1)+1:dm*i, dm*(i-1)+1:dm*i) = mij*irnp
        enddo
        !Compute pairwise mobility matrix for nanopost particles and feed into grand mobility matrix 
        do i=2,natom
                do j=i+1,natom
                        rij=rpos(:,i-1)-rpos(:,j-1)
                        rij = rij - box_size * dnint(rij*ibox_size)
                        call RPEwald2d(mij,rij,rnp,rnp)
                        mpbc(dm*(i-1)+1:dm*i,dm*(j-1)+1:dm*j) = mij*irnp
                        mpbc(dm*(j-1)+1:dm*j,dm*(i-1)+1:dm*i) = mij*irnp
                enddo
        enddo
        !Compute Rsum and Rself matrix on grid points
        do i = 0, ngrid(1)
                do j = 0, ngrid(2)
                        !Compute pairwise mobility matrix for particle
                        !and nanopost particle and feed into grand mobility matrix
                        do p = 2, natom
                                rij(1) = dble(i) * gridsize(1) - rpos(1,p-1)
                                rij(2) = dble(j) * gridsize(2) - rpos(2,p-1)
                                rij = rij - box_size * dnint(rij(1)*ibox_size)
                                call RPEwald2d(mij,rij,rp,rnp)
                                mpbc(1:dm,dm*(p-1)+1:dm*p) = mij*irp
                                mpbc(dm*(p-1)+1:dm*p,1:dm) = mij*irp
                        enddo
                        !Compute inverse of grand mobility matrix
                        call matrix_inv(mpbc,rpbc)
                        !Compue Rsum and Rself matrix
                        do m=1,dm
                                do n=1,dm
                                        rself(m,n)=rpbc(m,n)
                                        rsum(m,n)=sum(rpbc(m,n:dm*natom:dm))
                                enddo
                        enddo
                        !Write Rsum and Rself matrix to files
                        write(28,'(4f16.5)')rself(:,:)
                        write(56,'(4f16.5)')rsum(:,:)
                enddo
        enddo

endif
close(28)
close(56)
end program main

