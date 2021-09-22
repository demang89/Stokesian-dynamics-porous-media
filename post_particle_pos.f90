program main
use HI
implicit none
real(8) :: xx, yy, zz 
integer :: i, j, k

open(unit=14, file='post_particle.txt', action='read')
read(14,*)
read(14,*) n_x, n_y, n_z
read(14,*)
read(14,*) zeta, rp, rnp

open(unit=28, file='post_pos.txt', action='write')

! If packing geometry is 3D square
if(n_z .ne. 1) then
        dm = 3; natom = n_x*n_y*n_z + 1
        allocate(cell_size(dm))
        allocate(box_size(dm))
        allocate(ngrid(dm))
        allocate(gridsize(dm))
        allocate(rn(dm))
        allocate(kn(dm))

        write(28,'(2f20.8)') rp, rnp

        write(28,'(2I10)') dm, natom

        s = (2.d0 * rp / zeta + 2.d0 * rnp)
        cell_size(1) = s
        cell_size(2) = s
        cell_size(3) = 2.d0 * rnp

        write(28,'(3f20.8)') cell_size(:)

        box_size(1) = cell_size(1) * dble(n_x)
        box_size(2) = cell_size(2) * dble(n_y)
        box_size(3) = cell_size(3) * dble(n_z)

        write(28,'(3f20.8)') box_size(:)

        ngrid = idnint(cell_size) * 10

        write(28,'(3I10)') ngrid(:)

        gridsize = cell_size / ngrid

        write(28,'(3f20.8)') gridsize(:)

        alpha = 2.d0 / (box_size(1) * box_size(2) * box_size(3))**(1.d0/3.d0)
        rc = 14.d0 / alpha
        kc = 2.d0 * alpha**2 * rc

        write(28,'(3f20.8)') alpha, rc, kc

        rn=idnint(rc/box_size); kn=idnint(kc*box_size/2.d0/pi)

        write(28,'(3I10)') rn(:)
        write(28,'(3I10)') kn(:)

        do i = 1, n_x
                do j = 1, n_y
                        do k = 1, n_z
                                xx = dble(i-0.5d0) * cell_size(1)
                                yy = dble(j-0.5d0) * cell_size(2)
                                zz = dble(k-0.5d0) * cell_size(3)
                                write(28,'(3f20.8)') xx, yy, zz
                        enddo
                enddo
        enddo


! If packing geometry is 2D square
else

        dm = 2; natom = n_x*n_y*n_z + 1
        allocate(cell_size(dm))
        allocate(box_size(dm))
        allocate(ngrid(dm))
        allocate(gridsize(dm))
        allocate(rn(dm))
        allocate(kn(dm))

        write(28,'(2f20.8)') rp, rnp

        write(28,'(2I10)') dm, natom

        s = (2.d0 * rp / zeta + 2.d0 * rnp)
        cell_size(1) = s
        cell_size(2) = s

        write(28,'(2f20.8)') cell_size(:)

        box_size(1) = cell_size(1) * dble(n_x)
        box_size(2) = cell_size(2) * dble(n_y)

        write(28,'(2f20.8)') box_size(:)

        ngrid = idnint(cell_size) * 10

        write(28,'(2I10)') ngrid(:)

        gridsize = cell_size / ngrid

        write(28,'(2f20.8)') gridsize(:)

        alpha = 2.0 / (box_size(1) * box_size(2))**0.5d0
        rc = 14.d0 / alpha
        kc = 2.d0 * alpha**2 * rc

        write(28,'(3f20.8)') alpha, rc, kc

        rn=idnint(rc/box_size); kn=idnint(kc*box_size/2.d0/pi)

        write(28,'(2I10)') rn(:)
        write(28,'(2I10)') kn(:)

        do i = 1, n_x
                do j = 1, n_y
                        xx = dble(i-0.5d0) * cell_size(1)
                        yy = dble(j-0.5d0) * cell_size(2)
                        write(28,'(3f20.8)') xx, yy
                enddo
        enddo
endif

close(14)
close(28)
end program main

