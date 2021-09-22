program main
use HI
implicit none
integer :: i, j, k, seq, reason
integer(8) :: p
real(8) :: const, mc, rtot, rtot2, dr2
real(8) :: x1, x2, x3, x4, dhs2
integer :: r_grid(3)
real(8) :: r_frac(3)
real(8) :: ract(3), rtemp(3), rij(3)
real(8) :: var(3), Fr(3), Ft(3), vinf(3), Ff(3), Fc(3)
real(8) :: rself(3,3), rsum(3,3), r11(3,3), r12(3,3), umat(3,3)
real(8) :: lt(3,3), lt_temp(3,3), mij(3,3), mij_temp(3,3)
real(8) :: test(9)
real(8), allocatable :: rout(:,:,:,:,:), rsout(:,:,:,:,:)
character*20 :: filename

read(*,*) seq !Trajectory sequence 

!Read data for running particle simulations
open(unit=14, file='traj.txt', action='read')
read(14,*)
read(14,*) vmax, theta !Fluid velocity magnitude and flow direction angle
read(14,*)
read(14,*) nmax, t_grid, dt !No. of timesteps, traj record timesteps, and timestep value
read(14,*)
read(14,*) longrange !Whether far-field HI to be included
read(14,*)
read(14,*) shortrange !Whether short-range HI to be included
close(14)

call read_data() !Read nanopost particle position and other data
call read_m2p_coeff() !Read fitting parameters of inverse of mobility matrix
call lubcoefficient() !Compute coefficients of pair-wise resistance matrix
call init_random_seed() !Generate initial random seed

!Identity matrix
umat=0.d0
forall(i=1:3) umat(i,i)=rp

rtot=rp+rnp; rtot2=rtot*rtot
const=dsqrt(2.0d0/dt); mc=200.0d0
vinf=(/dble(vmax)*dcos(Pi/180.d0*dble(theta)),dble(vmax)*dsin(Pi/180.d0*dble(theta)),0.d0/)

allocate(rsout(0:ngrid(1), 0:ngrid(2), 0:ngrid(3), dm, dm)) !Allocate array for Rsum matrix
allocate( rout(0:ngrid(1), 0:ngrid(2), 0:ngrid(3), dm, dm)) !Allocate array for Rself matrix

!Read pre-computed Rsum and Rself matrix
open(unit=14, file='rself.txt', action='read')
open(unit=28, file='rsum.txt', action='read')
do i=0,ngrid(1)
        do j=0,ngrid(2)
                do k=0,ngrid(3)
                        read(14,*,iostat=reason) test
                        rout(i,j,k,1:dm,1:dm)=reshape(test,(/dm,dm/))
                        if(reason > 0) then
                                print *, 'something is wrong in tensor.txt'
                                stop
                        elseif(reason < 0) then
                                print *, 'not enough data in tensor.txt'
                                stop
                        endif
                        read(28,*,iostat=reason) test
                        rsout(i,j,k,1:dm,1:dm)=reshape(test,(/dm,dm/))
                        if(reason > 0) then
                                print *, 'something is wrong in sumtensor.txt'
                                stop
                        elseif(reason < 0) then
                                print *, 'not enough data in sumtensor.txt'
                                stop
                        endif
                enddo
        enddo
enddo
close(14)
close(28)

!Trajectory filename
write(filename,'("file",i0,".txt")') seq
open(unit=14,file=filename,action="write")

!Initial random position in unit cell
dhs2= -1.d0
do while(dhs2 .eq. -1.d0)
        call random_number(x1)
        call random_number(x2)
        call random_number(x3)
        ract = (/dble(x1)*cell_size(1), dble(x2)*cell_size(2), dble(x3)*cell_size(3)/)
        call periodic_HSInteraction(ract, rtot2, dhs2)
end do
write(14,1000) 0, ract(1), ract(2), ract(3)

!If there is no external flow
if(vmax .eq. 0.d0) then
!If both long and short range HI are included
        if(longrange .and. shortrange) then
                do p=1,nmax
                        call grid_location(ract, r_grid, r_frac)
                        call trilinearinterpolation(r_grid, r_frac, rout, rself)
                        do j=1,natom-1
                                 rij = ract - rpos(:,j); call pbc(rij); dr2 = sum(rij*rij)
                                 if(dr2<50.0d0)then
                                         call lub11correction(rij, r11, rp, rnp)
                                         rself = rself + r11
                                 endif
                        enddo
                        call choleskydecomposition(rself,lt,mij)
                        dhs2=-1.d0
                        do while(dhs2 .eq. -1.d0)
                                var(1)=const*randomnormal()
                                var(2)=const*randomnormal()
                                var(3)=const*randomnormal()
                                Fr=matmul(lt,var)
                                rtemp=ract + (matmul(mij,Fr)*dt/mc)
                                call grid_location(rtemp,r_grid,r_frac)
                                call trilinearinterpolation(r_grid,r_frac,rout,rself)
                                do j=1,natom-1
                                        rij=rtemp-rpos(:,j); call pbc(rij); dr2=sum(rij*rij)
                                        if(dr2<50.0d0)then
                                                call lub11correction(rij,r11,rnp,rp)
                                                rself = rself + r11
                                        endif
                                enddo
                                call choleskydecomposition(rself,lt_temp,mij_temp)
                                mij_temp = 0.5d0*mc*(mij_temp-mij)
                                rtemp = ract + matmul(mij,Fr)*dt + matmul(mij_temp,Fr)*dt
                                call periodic_HSInteraction(rtemp,rtot2,dhs2)
                        enddo
                        ract=rtemp
                        if(modulo(p,t_grid)==0) then
                                write(14,1000) p, ract(1), ract(2), ract(3)
                        endif
                enddo

        !If only long range HI is included
        elseif(longrange .and. (.not.(shortrange))) then
                do p=1,nmax
                        call grid_location(ract, r_grid, r_frac)
                        call trilinearinterpolation(r_grid, r_frac, rout, rself)
                        call choleskydecomposition(rself,lt,mij)
                        dhs2=-1.d0
                        do while(dhs2 .eq. -1.d0)
                                var(1)=const*randomnormal()
                                var(2)=const*randomnormal()
                                var(3)=const*randomnormal()
                                Fr=matmul(lt,var)
                                rtemp=ract + (matmul(mij,Fr)*dt/mc)
                                call grid_location(rtemp,r_grid,r_frac)
                                call trilinearinterpolation(r_grid,r_frac,rout,rself)
                                call choleskydecomposition(rself,lt_temp,mij_temp)
                                mij_temp = 0.5d0*mc*(mij_temp-mij)
                                rtemp = ract + matmul(mij,Fr)*dt + matmul(mij_temp,Fr)*dt
                                call periodic_HSInteraction(rtemp,rtot2,dhs2)
                        enddo
                        ract=rtemp
                        if(modulo(p,t_grid)==0) then
                                write(14,1000) p, ract(1), ract(2), ract(3)
                        endif
                enddo

        !If only short range HI is included
        elseif((.not.(longrange)) .and. shortrange) then
                do p=1,nmax
                        rself = umat; rsum = umat
                        do j=1,natom-1
                                 rij = ract - rpos(:,j); call pbc(rij); dr2 = sum(rij*rij)
                                 if(dr2<50.0d0)then
                                         call lub11correction(rij, r11, rp, rnp)
                                         rself = rself + r11
                                 endif
                        enddo
                        call choleskydecomposition(rself,lt,mij)
                        dhs2=-1.d0
                        do while(dhs2 .eq. -1.d0)
                                var(1)=const*randomnormal()
                                var(2)=const*randomnormal()
                                var(3)=const*randomnormal()
                                Fr=matmul(lt,var)
                                rtemp=ract + (matmul(mij,Fr)*dt/mc)
                                do j=1,natom-1
                                        rij=rtemp-rpos(:,j); call pbc(rij); dr2=sum(rij*rij)
                                        if(dr2<50.0d0)then
                                                call lub11correction(rij,r11,rnp,rp)
                                                rself = rself + r11
                                        endif   
                                enddo
                                call choleskydecomposition(rself,lt_temp,mij_temp)
                                mij_temp = 0.5d0*mc*(mij_temp-mij)
                                rtemp = ract + matmul(mij,Fr)*dt + matmul(mij_temp,Fr)*dt
                                call periodic_HSInteraction(rtemp,rtot2,dhs2)
                        enddo
                        ract=rtemp
                        if(modulo(p,t_grid)==0) then
                                write(14,1000) p, ract(1), ract(2), ract(3)
                        endif
                enddo

        !If both long and short range HI are not included
        else
                do p=1,nmax
                        dhs2=-1.d0
                        do while(dhs2 .eq. -1.d0)
                                var(1)=const*randomnormal()
                                var(2)=const*randomnormal()
                                var(3)=const*randomnormal()
                                rtemp = ract + var * dt
                                call periodic_HSInteraction(rtemp,rtot2,dhs2)
                        enddo
                        ract=rtemp
                        if(modulo(p,t_grid)==0) then
                                write(14,1000) p, ract(1), ract(2), ract(3)
                        endif
                enddo
        endif 
!If there is external flow
else
!If both long and short range HI are included
        if(longrange .and. shortrange) then
                do p=1,nmax
                        call grid_location(ract, r_grid, r_frac)
                        call trilinearinterpolation(r_grid, r_frac, rout, rself)
                        call trilinearinterpolation(r_grid, r_frac, rsout, rsum)
                        do j=1,natom-1
                                 rij = ract - rpos(:,j); call pbc(rij); dr2 = sum(rij*rij)
                                 if(dr2<50.0d0)then
                                         call lub11correction(rij, r11, rp, rnp)
                                         call lub12correction(rij, r12, rp, rnp)
                                         rsum = rsum + r12 + r11
                                         rself = rself + r11
                                 endif
                        enddo
                        Ff=matmul(rsum,vinf)
                        call choleskydecomposition(rself,lt,mij)
                        dhs2=-1.d0
                        do while(dhs2 .eq. -1.d0)
                                var(1)=const*randomnormal()
                                var(2)=const*randomnormal()
                                var(3)=const*randomnormal()
                                Fr=matmul(lt,var)
                                Ft=Fr+Ff
                                rtemp=ract + (matmul(mij,Ft)*dt/mc)
                                call grid_location(rtemp,r_grid,r_frac)
                                call trilinearinterpolation(r_grid,r_frac,rout,rself)
                                do j=1,natom-1
                                        rij=rtemp-rpos(:,j); call pbc(rij); dr2=sum(rij*rij)
                                        if(dr2<50.0d0)then
                                                call lub11correction(rij,r11,rnp,rp)
                                                rself = rself + r11
                                        endif
                                enddo
                                call choleskydecomposition(rself,lt_temp,mij_temp)
                                mij_temp = 0.5d0*mc*(mij_temp-mij)
                                rtemp = ract + matmul(mij,Ft)*dt + matmul(mij_temp,Ft)*dt
                                call periodic_HSInteraction(rtemp,rtot2,dhs2)
                        enddo
                        ract=rtemp
                        if(modulo(p,t_grid)==0) then
                                write(14,1000) p, ract(1), ract(2), ract(3)
                        endif
                enddo

        !If only long range HI is included
        elseif(longrange .and. (.not.(shortrange))) then
                do p=1,nmax
                        call grid_location(ract, r_grid, r_frac)
                        call trilinearinterpolation(r_grid, r_frac, rout, rself)
                        call trilinearinterpolation(r_grid, r_frac, rsout, rsum)
                        Ff=matmul(rsum,vinf)
                        call choleskydecomposition(rself,lt,mij)
                        dhs2=-1.d0
                        do while(dhs2 .eq. -1.d0)
                                var(1)=const*randomnormal()
                                var(2)=const*randomnormal()
                                var(3)=const*randomnormal()
                                Fr=matmul(lt,var)
                                Ft=Fr+Ff
                                rtemp=ract + (matmul(mij,Ft)*dt/mc)
                                call grid_location(rtemp,r_grid,r_frac)
                                call trilinearinterpolation(r_grid,r_frac,rout,rself)
                                call choleskydecomposition(rself,lt_temp,mij_temp)
                                mij_temp = 0.5d0*mc*(mij_temp-mij)
                                rtemp = ract + matmul(mij,Ft)*dt + matmul(mij_temp,Ft)*dt
                                call periodic_HSInteraction(rtemp,rtot2,dhs2)
                        enddo
                        ract=rtemp
                        if(modulo(p,t_grid)==0) then
                                write(14,1000) p, ract(1), ract(2), ract(3)
                        endif
                enddo

        !If only short range HI is included
        elseif((.not.(longrange)) .and. shortrange) then
                do p=1,nmax
                        rself = umat; rsum = umat
                        do j=1,natom-1
                                 rij = ract - rpos(:,j); call pbc(rij); dr2 = sum(rij*rij)
                                 if(dr2<50.0d0)then
                                         call lub11correction(rij, r11, rp, rnp)
                                         call lub12correction(rij, r12, rp, rnp)
                                         rsum = rsum + r12 + r11
                                         rself = rself + r11
                                 endif
                        enddo
                        Ff=matmul(rsum,vinf)
                        call choleskydecomposition(rself,lt,mij)
                        dhs2=-1.d0
                        do while(dhs2 .eq. -1.d0)
                                var(1)=const*randomnormal()
                                var(2)=const*randomnormal()
                                var(3)=const*randomnormal()
                                Fr=matmul(lt,var)
                                Ft=Fr+Ff
                                rtemp=ract + (matmul(mij,Ft)*dt/mc)
                                do j=1,natom-1
                                        rij=rtemp-rpos(:,j); call pbc(rij); dr2=sum(rij*rij)
                                        if(dr2<50.0d0)then
                                                call lub11correction(rij,r11,rnp,rp)
                                                rself = rself + r11
                                        endif   
                                enddo
                                call choleskydecomposition(rself,lt_temp,mij_temp)
                                mij_temp = 0.5d0*mc*(mij_temp-mij)
                                rtemp = ract + matmul(mij,Ft)*dt + matmul(mij_temp,Ft)*dt
                                call periodic_HSInteraction(rtemp,rtot2,dhs2)
                        enddo
                        ract=rtemp
                        if(modulo(p,t_grid)==0) then
                                write(14,1000) p, ract(1), ract(2), ract(3)
                        endif
                enddo

        !If both long and short range HI are not included
        else
                do p=1,nmax
                        dhs2=-1.d0
                        do while(dhs2 .eq. -1.d0)
                                var(1)=const*randomnormal()
                                var(2)=const*randomnormal()
                                var(3)=const*randomnormal()
                                rtemp = ract + (var + vinf) * dt
                                call periodic_HSInteraction(rtemp,rtot2,dhs2)
                        enddo
                        ract=rtemp
                        if(modulo(p,t_grid)==0) then
                                write(14,1000) p, ract(1), ract(2), ract(3)
                        endif
                enddo
        endif 
endif
1000 format(i11,x,f20.8,x,f20.8,x,f20.8)
close(14)
deallocate(rout, rsout, rpos, cell_size, box_size, ngrid, gridsize, rn, kn)
deallocate(ibox_size, igridsize)
end program main
!############################################################
