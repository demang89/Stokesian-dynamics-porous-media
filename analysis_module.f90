module analysis_parameter
real(8), parameter :: pi = 3.1415926d0
integer :: dm, ntraj, nlines !system dimension, no. of trajectories, no. of lines in a trajectory
integer :: tsteps !Number of time-steps to calculate MSDs
integer :: del !number of lines for time-interval \Delta{t}=1.0
integer :: size_1 !Number of lines in a trajectory at time-interval \Delta{t}=1.0
integer :: size_2 ! Number of time-steps to calculate trajectory metrics
real(8) :: dt  !time-interval at which trajectory was recorded
real(8) :: theta !flow direction angle
integer :: flow
real(8), allocatable :: fx(:,:), fy(:,:) ! arrays to store trajectory in x and y directions
contains
!####################
!Subroutine to read input file
subroutine read_input()
implicit none
open(unit=100,file="analysis.txt",action="read")
read(100,*)
read(100,*) dm, ntraj, nlines
read(100,*)
read(100,*) tsteps, size_2
read(100,*)
read(100,*) del
read(100,*)
read(100,*) dt
read(100,*)
read(100,*) flow
read(100,*)
read(100,*) theta
close(100)
return
end subroutine read_input
!##########################################################
!Subroutine to read particle trajectories
subroutine read_traj()
implicit none
integer :: i, j
real(8) :: x1, x2, x3, x4
real(8) :: p1, p2
real(8) :: temp(dm)
character*20 :: filename

allocate(fx(ntraj,nlines))
allocate(fy(ntraj,nlines))

p1 = dcos(pi*theta/180.d0)
p2 = dsin(pi*theta/180.d0)

do i=1,ntraj
        write(filename,'("file",i0,".txt")') (i)
        open(unit=100,file=filename,action="read")
        do j=1,nlines
                read(100,*) x1, temp(:)
                fx(i,j) =  temp(1) * p1 + temp(2) * p2
                fy(i,j) =  temp(2) * p1 - temp(1) * p2
        enddo
        close(100)
enddo
return
end subroutine read_traj
!############################
!Subroutine to compute average particle velocity
subroutine average_velocity(velx,velxstd,vely,velystd)
implicit none
real(8), intent(out) :: velx, vely, velxstd, velystd
integer :: i
real(8) :: vx(ntraj), vy(ntraj)

do i = 1, ntraj
        vx(i) = (fx(i,nlines) - fx(i,1)) / dt / dble(nlines-1)
        vy(i) = (fy(i,nlines) - fy(i,1)) / dt / dble(nlines-1)
enddo

velx = sum(vx) / dble(ntraj)
vely = sum(vy) / dble(ntraj)
velxstd=dsqrt((sum(vx*vx)-dble(ntraj)*velx**2)/dble(ntraj-1)/dble(ntraj))
velystd=dsqrt((sum(vy*vy)-dble(ntraj)*vely**2)/dble(ntraj-1)/dble(ntraj))
return
end subroutine average_velocity
!###############################
!Subroutine to compute in-plane time-dependent diffusivity under no flow conditions
subroutine MSD_calc_noflow()
implicit none
integer :: i, j, k
real(8) :: x1, x2
real(8) :: msd(tsteps,ntraj),msdavg(tsteps),msdstd(tsteps)

do i=1,ntraj
        do j=1,tsteps
               k=j
               x1=sum((fx(i,k+1:nlines)-fx(i,1:nlines-k))**2)
               x2=sum((fy(i,k+1:nlines)-fy(i,1:nlines-k))**2)
               msd(j,i)=(x1+x2)/size(fx(1,k+1:nlines))/dt/dble(k)/4.d0
        enddo
enddo

msdavg=sum(msd,dim=2)/dble(ntraj)
msdstd=sqrt((sum(msd*msd,dim=2)-dble(ntraj)*msdavg**2)/dble(ntraj-1)/dble(ntraj))

open(unit=50,file='diff_avg.txt',action="write")
do i=1,tsteps
        write(50,'(3f20.8)') dt*i, msdavg(i),msdstd(i)
enddo
close(50)
return
end subroutine MSD_calc_noflow
!#########################################
!Subroutine to compute time-dependent dispersion coefficients under flow conditions
subroutine MSD_calc_flow(velx,velxstd,vely,velystd)
implicit none
real(8), intent(in) :: velx, velxstd, vely, velystd
integer :: i, j, k
real(8) :: x1, x2
real(8) :: msdx(tsteps,ntraj),msdxavg(tsteps),msdxstd(tsteps)
real(8) :: msdy(tsteps,ntraj),msdyavg(tsteps),msdystd(tsteps)

do i=1,ntraj
        do j=1,tsteps
               k=j
               x1=sum((fx(i,k+1:nlines)-fx(i,1:nlines-k)-velx*dt*dble(k))**2)
               x2=sum((fy(i,k+1:nlines)-fy(i,1:nlines-k)-vely*dt*dble(k))**2)
               msdx(j,i)=x1/size(fx(i,k+1:nlines))/dt/dble(k)/2.d0
               msdy(j,i)=x2/size(fx(i,k+1:nlines))/dt/dble(k)/2.d0

        enddo
enddo
msdxavg=sum(msdx,dim=2)/dble(ntraj)
msdxstd=sqrt((sum(msdx*msdx,dim=2)-dble(ntraj)*msdxavg**2)/dble(ntraj-1)/dble(ntraj))
msdyavg=sum(msdy,dim=2)/dble(ntraj)
msdystd=sqrt((sum(msdy*msdy,dim=2)-dble(ntraj)*msdyavg**2)/dble(ntraj-1)/dble(ntraj))

open(unit=50,file='disp_avg.txt',action="write")
write(50,'(4f20.8)') velx, velxstd, vely, velystd
do i=1,tsteps
        write(50,'(5f20.8)') dt*i, msdxavg(i),msdxstd(i),msdyavg(i),msdystd(i)
enddo
close(50)
return
end subroutine MSD_calc_flow
!######################################################################
!Subroutine to compute cosine correlation angle
subroutine correlation_func()
implicit none
integer :: i, j, k
real(8) :: x
real(8) :: vx(size_1 - 1), vy(size_1 - 1), vmag(size_1 - 1)
real(8) :: corr(size_2,ntraj), corr_avg(size_2), corr_std(size_2) 
integer :: length(size_1)

length = (/(i,i=1,nlines,del)/)

do i=1,ntraj
        do j = 1, size_1-1
                vx(j) = fx(i,length(j+1)) - fx(i,length(j))
                vy(j) = fy(i,length(j+1)) - fy(i,length(j)) 
                vmag(j) = dsqrt(vx(j)*vx(j) + vy(j)*vy(j))
        enddo
        do j = 1, size_2
                x=0
                do k= 1, size_1-1-j
                        x = x + (vx(k) * vx(k+j) + vy(k) * vy(k+j)) / (vmag(k) * vmag(k+j))
                enddo
                corr(j,i)=x/dble(size_1-j-1)
        enddo
enddo

corr_avg=sum(corr,dim=2)/dble(ntraj)
corr_std=sqrt((sum(corr*corr,dim=2)-dble(ntraj)*corr_avg**2)/dble(ntraj-1)/dble(ntraj))

open(unit=50,file='correlation_avg.txt',action="write")
do i=1,size_2
        write(50,'(3f20.8)') dble(i), corr_avg(i),corr_std(i)
enddo
close(50)
end subroutine correlation_func
!##################################################################
!Subroutine to compute tortuosity
subroutine tort_calc()
implicit none
integer :: i, j, k
real(8) :: x
real(8) :: disp(size_1-1)
real(8) :: tort(size_2,ntraj), tort_avg(size_2), tort_std(size_2)
integer :: length(size_1)

length = (/(i,i=1,nlines,del)/)

do i = 1, ntraj
        do j = 1, size_1-1
                disp(j) = dsqrt((fx(i,length(j+1))-fx(i,length(j)))**2+(fy(i,length(j+1))-fy(i,length(j)))**2)
        enddo
        do j = 1, size_2
                x = 0.0
                do k= 1, size_1-1-j
                        x = x + sum(disp(k:k+j)) / dsqrt((fx(i,length(k+j))-fx(i,length(k)))**2+(fy(i,length(k+j))-fy(i,length(k)))**2)
                end do
                tort(j,i) = x / dble(size_1-1-j)
        enddo
enddo

tort_avg=sum(tort,dim=2)/dble(ntraj)
tort_std=sqrt((sum(tort*tort,dim=2)-dble(ntraj)*tort_avg**2)/dble(ntraj-1)/dble(ntraj))

open(unit=50,file='tort_avg.txt',action="write")
do i=1,size_2
        write(50,'(3f20.8)') dble(i), tort_avg(i), tort_std(i)
enddo
close(50)
return
end subroutine tort_calc
!###############################
end module analysis_parameter
