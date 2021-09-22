program main
use analysis_parameter
implicit none
integer :: i
real(8) :: velx, velxstd, vely, velystd

!Read input file 'analysis.txt'
call read_input()

!Trajectory point index to be used to trajectory metric parameters
size_1 = size((/(i,i=1,nlines,del)/))

call read_traj() !Read trajectories

!Compute average particle velocity and dispersion coefficients in case of flow
if(flow) then
    call average_velocity(velx,velxstd,vely,velystd)
    call MSD_calc_flow(velx,velxstd,vely,velystd)
else
!Compute diffusivity in case of no flow
    call MSD_calc_noflow()
endif

!Compute cosine of correlation angle and tortuosity
call correlation_func()
call tort_calc()

end program main
