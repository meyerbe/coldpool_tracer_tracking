! 2018 Sep 08 O Henneb
!
! special vrsion for Bettina 
! Dec 2018 
!
! read rain tracking output txt for COGs (loop trough every line)
! read 2D velocity field (time-stepwise)
! read rain track cells for precipitation boundaries (timestepwise)
! SUBROUTINE neigh idetifies precipitation boundaries
! SUBROUTINE set_tracer sets tracer at the boundaries
! SUBROUTINE ....
! SUBROUTINE ... 
! 
! 2018 Oct:
! change initial tracer placement from precip boundary to circles
! calculate effective radius r = sqrt(A)/pi 
! Problem: tracers which are inside the precip area never chatch up
! alternatives: 
!    set tracer at outremost circle
!    distribute tracer equally around outline
!      * sort boundary points by angle 
!      * calculate distance from one gp to another and sum it up ->
!        circulferance    
!      * divide by circum by no of desired tracer to get distance 
!
! 2018 Nov:
! change order of CALL routine: 
! output before update of tracer -> output of windvector match the position of
! tracer at output now, velocities before subtimestepping
!
! 2018 Dec:
! include option to track tracer base don radial velocity instead of full vector
! output of radial and tangential velocity also included
!
! 2018 Dec:
! allow precipitation tracking starting with smaller precipitation area, but set
! tracer not untill a larger size is reached which is given in jobfile  
! number of precipitation cells is given as input now. But can be done nicer
!
! OUTPUT: 
! traced(CP-ID, int tracer ID, property)
! traced(:,:, 1)   x pos of tracer 
! traced(:,:, 2)   y pos of tracer
! traced(:,:, 3)   x gp 
! traced(:,:, 4)   y gp
! traced(:,:, 5)   distance to COG
! traced(:,:, 6)   timestep 
! traced(:,:, 7)   CP age defined by onset of precip ??? check this 
! traced(:,:, 8)   angle to COG
! traced(:,:, 9)   CP ID
! traced(:,:,10)   gloabl tracer number 
! traced(:,:,11)   active 0-1 
! traced(:,:,12)   start time of tracer 
! traced(:,:,13)   u 
! traced(:,:,14)   v
! traced(:,:,15)   distance from tracer to COG in x direction 
! traced(:,:,16)   merger 
! traced(:,:,17)   distance from tracer to COG in x direction
! traced(:,:,18)   precip is ongoing
! traced(:,:,19)   vrad
! traced(:,:,20)   vtan
PROGRAM cp_tracking

USE netcdf
USE cp_parameters !, ONLY: dsize_x, dsize_y, max_tracer_CP, dt ! &
!                          resolution, dt ! resolution in m! dt in sec

IMPLICIT NONE
REAL, ALLOCATABLE    :: vel(:,:,:)               ! velocity field
REAL, ALLOCATABLE    :: COMx(:), COMy(:)         ! store COG
REAL, ALLOCATABLE    :: rmax(:)                  ! area of precip
!INTEGER, ALLOCATABLE :: IDstart(:)               ! first timestep of precip 
INTEGER, ALLOCATABLE :: already_tracked(:)       ! memory of cell counter
REAL, ALLOCATABLE    :: traced(:,:,:)            ! tracked information, CP ID x internal tarcer ID x properties
REAL, ALLOCATABLE    :: traced_dummy(:,:)        ! dummy to sort tracer by angle
INTEGER,ALLOCATABLE  :: cpio(:,:)                ! to identify merger (dim1) and start time (dim2) splitting events have start time 0 to avoid them
INTEGER              :: ID                       ! local ID from rain track
INTEGER              :: i                        ! running index
INTEGER              :: ierr                     ! error index for reading files
REAL                 :: xCOG, yCOG               ! buffer variable for read COGs
INTEGER              :: onset                    ! first cell reaches treshold
INTEGER              :: max_no_of_cells          ! maximum number of CPs (can be retrieved from the number of rain cells
!INTEGER              :: max_tracer_CP            ! max no of tracers per CP
INTEGER              :: max_tracers              ! max no of commulative tracers 
INTEGER,ALLOCATABLE  :: tracpo(:,:)              ! keeps track of first two indices in traced (CP ID, internal tracer ID) for every tracer 
INTEGER              :: srv_header_input(8)
INTEGER              :: timestep
INTEGER              :: count_tracer             ! counts the internal tracer in an individual CP
INTEGER              :: tracking_end
!INTEGER              :: dsize_y,dsize_x
!INTEGER              :: n_sub                    ! number of sub-timesteps
!INITIALIZE some values
namelist /INPUTgeneral/  dt, res 
namelist /INPUTtracer/ max_tracer_CP, max_age, rad, n_sub
namelist /INPUTIO/ odir
open(100,file='job/namelist.dat')
read(100,nml=INPUTgeneral)
read(100,nml=INPUTIO)
read(100,nml=INPUTtracer)

!get number of tracked precip events
!CALL getarg(max_no_of_cells)
OPEN(1111,FILE='na.txt')
READ(1111,*) max_no_of_cells
max_tracers = max_no_of_cells*max_tracer_CP
write(*,*) 'max no of cps:', max_no_of_cells
 CALL get_dim(tracking_end,dsize_y,dsize_x,trim(odir) //'/input/uv_alltimes.nc') 
write(*,*) tracking_end,dsize_x,dsize_y
! allocate fields
ALLOCATE(cpio(max_no_of_cells,3))
ALLOCATE(traced(max_no_of_cells,max_tracer_CP,20))
ALLOCATE(traced_dummy(max_tracer_CP,20))
!ALLOCATE(IDstart(max_no_of_cells))
ALLOCATE(COMx(max_no_of_cells))
ALLOCATE(COMy(max_no_of_cells))
ALLOCATE(rmax(max_no_of_cells))
ALLOCATE(vel(dsize_x,dsize_y,2))
ALLOCATE(already_tracked(max_no_of_cells))
ALLOCATE(tracpo(2,max_tracers))

!read data for CP setting and start of tracer
 i =1
 OPEN(1,FILE=trim(odir) // '/input/mergingCPs.txt',FORM='formatted',ACTION='read',IOSTAT=ierr) 
 write(*,*) trim(odir) // '/input/mergingCPs.txt'
 IF ( ierr == 0) then
   DO
     READ(1,*,END=400) cpio(i,1), cpio(i,3), COMx(i), COMy(i), rmax(i)
     i = i +1
   END DO
 ELSE 
     write(*,*) 'Beim Oeffnen der Datei ist ein Fehler Nr.', ierr,' aufgetreten'
 END IF
400 CONTINUE

!OPEN(4,FILE=trim(odir) // '/input/input_u.srv',FORM='unformatted', ACTION='read')
!OPEN(5,FILE=trim(odir) // '/input/input_v.srv',FORM='unformatted',ACTION='read')
OPEN(40,FILE=trim(odir) // '/output/coldpool_tracer_out.txt',FORM='formatted', ACTION='write')

!INITIALIZE some values
already_tracked(:) = 0
traced(:,:,:) = 0.
count_tracer = 1 ! counts individual pixels !OCH was ist mit pixeln gemeint? der
!IDstart(:) = 0

! read when first precip is tracked
 onset = minval(cpio(:,3),1)


!n_sub = 1
write(*,*) "setting # sub-timesteps to ", n_sub


 write(*,*) "start main loop"
 DO timestep = onset,tracking_end
   write(*,*) 'timestep',timestep, 'onset', onset
! read velocity files
   CALL read_nc_2D (vel(:,:,1),'u',trim(odir) // '/input/uv_alltimes.nc',timestep)
   CALL read_nc_2D (vel(:,:,2),'v',trim(odir) // '/input/uv_alltimes.nc',timestep)

   !READ (4,END=200) srv_header_input
   !READ (4) vel(:,:,1)
   !READ (5,END=200) srv_header_input
   !READ (5) vel(:,:,2)

   ! read velocity files until start of tracking is reached, if so, also read tracking
   !   max_no_of_cells              ! max number of CPs
   !   traced                       ! tracked information, (CP ID x internal tarcer ID x properties) (no time dimension)
   !   tracpo                       ! keeps track of first two indices in traced (CP ID, internal tracer ID) for every tracer

   ! only called when timestep==cpio(:,3); if so, initialize new tracers
   CALL initCircle(max_no_of_cells, timestep, traced, COMx,COMy,&
                  rmax, cpio, max_tracers, tracpo, count_tracer, dsize_x, dsize_y)
   ! interpolate velocity field at t=timestep onto tracer position in traced (position at t=timestep-1 !??)
   !   and save onto traced(:,:,13:14)
   CALL velocity_interpol(vel(:,:,1), vel(:,:,2), timestep, traced, count_tracer, &
                  max_no_of_cells, tracpo, max_tracers, dsize_x, dsize_y, n_sub)
   ! compute polar coordinates (distance and angle) of tracers based on traced (at t=timestep-1!?) (>> traced(:,:,i),i=5,8,15,17)
   CALL geometry(traced,COMx,COMy,already_tracked,max_no_of_cells,dsize_x,dsize_y)
   ! compute radial and tangential velocity based on updated velocity (v(t=timestep,x(t=timestep-1))
   CALL radvel(traced,already_tracked,max_no_of_cells)
   ! write output: x(timestep-1), v(timestep,x(timestep-1))
   CALL write_output(traced,max_tracers,count_tracer,timestep,tracpo,&
                  max_no_of_cells,COMx,COMy)
   ! sub-time stepping: interpolate velocity onto tracer position v(timestep,x(timestep-1+n/5*dt_sub)) (for first subtime-step same as in velocity_interpol)
   !   and update velocity v(timestep, x(timestep-1+4/5*dt_sub) >> (traced(:,:,13:14))
   !   and update tracer position to where it will be advected by v(timestep, x(timestep)) (>> traced(:,:,1:4))
   CALL update_tracer(vel(:,:,1),vel(:,:,2),timestep,traced, count_tracer, &
                  max_no_of_cells,tracpo,max_tracers,dsize_x,dsize_y,n_sub)
 END DO
 WRITE(*,*) "finished main loop" 
! 200 CONTINUE 

CONTAINS
! -----------------------------------------------------------------------
!subroutines for reading and writing netcdf
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
! input 2D nc data
! -----------------------------------------------------------------------
  SUBROUTINE read_nc_2D (poutput,varname, filename,ctime)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  real, dimension(:,:), intent(inout) :: poutput
  integer :: ncId, rhVarId, nx, ny, nz, nt, ctime
  integer, dimension(nf90_max_var_dims) :: dimIDs
!  real, allocatable, dimension(:,:,:) ::  zvar
    write(*,*) varname, ctime
    CALL check(nf90_open(filename, nf90_NoWrite, ncid))
    CALL check(nf90_inq_varid(ncid,varname, rhVarId))
    CALL check(nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nt))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nx))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
    CALL check(nf90_get_var(ncid, rhVarId, poutput, start =(/1,1,ctime/),&
                                                    count= (/ny, nx, 1/)))
    CALL check(nf90_close(ncid))
  END SUBROUTINE read_nc_2D


! -----------------------------------------------------------------------
! input 2D nc data
! -----------------------------------------------------------------------
  SUBROUTINE get_dim(nt,ny,nx, filename)

  IMPLICIT NONE
  character(*), intent(in) :: filename
  integer, intent(out)     :: nt, nx, ny
  integer                  :: ncId, rhVarId
  integer, dimension(nf90_max_var_dims) :: dimIDs
  write(*,*) "get dimensions of input file: ", filename
    CALL check(nf90_open(filename, nf90_NoWrite, ncid))
    CALL check(nf90_inq_varid(ncid,"u", rhVarId))
    CALL check(nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nt))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nx))
    CALL check(nf90_close(ncid))
  END SUBROUTINE  get_dim


!---------------------------------------------------
! CHECK NETCD DATA
!--------------------------------------------------
  SUBROUTINE check(status)

    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  END SUBROUTINE check

!--------------------------------------------------------------------------------------
! Calculate  circle around COG dependent on size for initial tracer placement
! replaces routine neighbours
!--------------------------------------------------------------------------------------
!SUBROUTINE initCircle(max_no_of_cells,ts,IDstart,traced, COMx,COMy, &
!                      rmax,cpio,max_tracers,tracpo,count_tracer)
SUBROUTINE initCircle(max_no_of_cells,ts,traced, COMx,COMy, &
                      rmax,cpio,max_tracers,tracpo,count_tracer,dsize_x, dsize_y)
   USE cp_parameters, ONLY : max_tracer_CP
  
  INTEGER, INTENT(IN)       :: max_no_of_cells, max_tracers, ts
  REAL, INTENT(IN)          :: rmax(max_no_of_cells),&
                               COMx(max_no_of_cells), COMy(max_no_of_cells)
  INTEGER, INTENT(IN)       :: cpio(max_no_of_cells,3)
  INTEGER, INTENT(INOUT)    :: count_tracer
  INTEGER, INTENT(INOUT)    :: tracpo(2,max_tracers)
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: pi, inc !, reff 
  INTEGER                   :: i, j
  INTEGER, INTENT(IN)       :: dsize_x, dsize_y
  
  pi = 2.*asin(1.)
  inc = 2*pi/max_tracer_CP
  do i = 1,max_no_of_cells,1 ! loop trough precip cells
   if (cpio(i,3) == ts) then ! set new circle when precip begins
     !reff = sqrt(area(i))/pi 
      do j = 1,max_tracer_CP,1 ! loop trough tracer number
        tracpo(:,count_tracer)=(/i,j/)
        count_tracer           = count_tracer + 1
        traced(i,j, 1) = MOD(COMx(i) + rmax(i)*cos(j*inc)-1.,float(dsize_x))+1.
        traced(i,j, 2) = MOD(COMy(i) + rmax(i)*sin(j*inc)-1.,float(dsize_y))+1.
        traced(i,j, 6) = ts ! current time
        traced(i,j, 7) = 0   ! age 
        traced(i,j, 16) = 1
        traced(i,j,12) = cpio(i,3)
        !traced(i,j, 8) = j*inc 
        traced(i,j, 9) = i 
        traced(i,j,10) = count_tracer 
        traced(i,j,11) = 1   ! active tracer 
        !traced(i,j,13) = 0   ! u 
        !traced(i,j,14) = 0   ! v
      end do
   end if
  end do
END SUBROUTINE 


!--------------------------------------------------------------------------------------
! INTERPOLATE velocities 
!-----------------------------------------------------------------------------------
SUBROUTINE velocity_interpol(velx,vely,timestep,traced, count_tracer,max_no_of_cells,tracpo,max_tracers,dsize_x,dsize_y,n_sub)
USE cp_parameters, ONLY :  res, dt, max_tracer_CP

  INTEGER, INTENT(IN)       :: timestep,max_no_of_cells, max_tracers
  REAL, INTENT(IN)          :: velx(dsize_x,dsize_y), &
                               vely(dsize_x,dsize_y)
  INTEGER, INTENT(IN)       :: count_tracer !,max_tracer_CP
  !INTEGER, INTENT(IN)       :: cp_field(dsize_x,dsize_y)
  INTEGER                   :: tracer_ts, it
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: ix_new, iy_new, vx_intp, vy_intp
  REAL                      :: ix, iy, ixt, iyt    !t are half level (interface)
!  INTEGER                   :: ix_round, iy_round, &
  INTEGER                   :: ix_new_round, iy_new_round, start_time
  INTEGER                   :: ix_l, iy_l, ix_r, iy_r
  REAL                      :: wgt_x, wgt_y, wgt_xt, wgt_yt
  INTEGER, INTENT(IN)       :: tracpo(2,max_tracers)
  INTEGER                   :: sub_dt   ! subtimstepping
  INTEGER, INTENT(IN)       :: dsize_x, dsize_y
  !INTEGER                   :: n_sub
  INTEGER, INTENT(IN)       :: n_sub
  !n_sub = 5
  sub_dt = dt/n_sub ! update evry min
  write(*,*) 'velocity interpolation: dt ', dt, 'sub_dt ', sub_dt
  it = 1 ! counter trough traced
  DO WHILE (it .LT. count_tracer) ! count tracer are all tracers set until here

    ! determining the first time step of the event
     start_time=INT(traced(tracpo(1,it),tracpo(2,it),6))
     ! determining how many timesteps have passed since then
     tracer_ts =timestep-start_time+1

     ! getting the previous positions
     ix=traced(tracpo(1,it),tracpo(2,it),1)
     iy=traced(tracpo(1,it),tracpo(2,it),2) ! position in gridpoints

     !for u-Wind defined on xt (half level in x direction)
     ixt = ix-0.5

     wgt_xt  = MOD(ixt,1.)
     wgt_y  = MOD(iy,1.)

     ix_l = MOD(INT(ixt-wgt_xt)-1+dsize_x,dsize_x)+1
     iy_l = MOD(INT(iy-wgt_y)-1+dsize_y,dsize_y)+1

     ix_r = MOD(ix_l +dsize_x,dsize_x)+1
     iy_r = MOD(iy_l +dsize_y,dsize_y)+1

     vx_intp = velx(iy_l,ix_l)*(1-wgt_xt)*(1-wgt_y) &
             + velx(iy_l,ix_r)*(1-wgt_xt)*(  wgt_y) &
             + velx(iy_r,ix_l)*(  wgt_xt)*(1-wgt_y) &
             + velx(iy_r,ix_r)*(  wgt_xt)*(  wgt_y)

     !for v-Wind defined on xt (half level in x direction)
     iyt = iy-0.5

     wgt_x  = MOD(ix,1.)
     wgt_yt  = MOD(iyt,1.)

     ix_l = MOD(INT(ix-wgt_x)-1+dsize_x,dsize_x)+1
     iy_l = MOD(INT(iyt-wgt_yt)-1+dsize_y,dsize_y)+1

     ix_r = MOD(ix_l +dsize_x,dsize_x)+1
     iy_r = MOD(iy_l +dsize_y,dsize_y)+1

     vy_intp = vely(iy_l,ix_l)*(1-wgt_x)*(1-wgt_yt) &
             + vely(iy_l,ix_r)*(1-wgt_x)*(  wgt_yt) &
             + vely(iy_r,ix_l)*(  wgt_x)*(1-wgt_yt) &
             + vely(iy_r,ix_r)*(  wgt_x)*(  wgt_yt)

     !save the velocities of the tracer
     traced(tracpo(1,it),tracpo(2,it),13) = vx_intp
     traced(tracpo(1,it),tracpo(2,it),14) = vy_intp
     it = it+1
   END DO
  RETURN
END SUBROUTINE velocity_interpol 

!--------------------------------------------------------------------------------------
! UPDATE TRACER along the horizontal wind field
!-----------------------------------------------------------------------------------
SUBROUTINE update_tracer(velx,vely,timestep,traced, count_tracer,max_no_of_cells,tracpo,max_tracers,dsize_x,dsize_y,n_sub)
!OCH TO DO: dt and resolution should be parameter read by USE from module
USE cp_parameters, ONLY : res, dt, max_tracer_CP

  INTEGER, INTENT(IN)       :: timestep,max_no_of_cells, max_tracers
  REAL, INTENT(IN)          :: velx(dsize_x,dsize_y), &
                               vely(dsize_x,dsize_y)
  INTEGER, INTENT(IN)       :: count_tracer !,max_tracer_CP
  !INTEGER, INTENT(IN)       :: cp_field(dsize_x,dsize_y)
  INTEGER                   :: tracer_ts, it
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: ix_new, iy_new, vx_intp, vy_intp
  REAL                      :: ix, iy, ixt, iyt    !t are half level (interface)
!  INTEGER                   :: ix_round, iy_round, &
  INTEGER                   :: ix_new_round, iy_new_round, start_time
  INTEGER                   :: ix_l, iy_l, ix_r, iy_r
  REAL                      :: wgt_x, wgt_y, wgt_xt, wgt_yt
  INTEGER, INTENT(IN)       :: tracpo(2,max_tracers)
  INTEGER                   :: sub_dt   ! subtimstepping
  INTEGER, INTENT(IN)       :: dsize_x, dsize_y
  INTEGER, INTENT(IN)       :: n_sub
  !INTEGER                   :: n_sub
  !n_sub = 5
  sub_dt =dt/n_sub ! update evry min
  write(*,*) 'update tracer: sub_dt ', sub_dt, ', res (dx) ', res
  it = 1 ! counter trough traced
  DO WHILE (it .LT. count_tracer) ! count tracer are all tracers set until here
    ! determining the first time step of the event
     start_time=INT(traced(tracpo(1,it),tracpo(2,it),6))
     ! determining how many timesteps have passed since then
     tracer_ts =timestep-start_time+1

     ! getting the previous positions
     ix=traced(tracpo(1,it),tracpo(2,it),1)
     iy=traced(tracpo(1,it),tracpo(2,it),2) ! position in gridpoints

     !for u-Wind defined on xt (half level in x direction)
     do i = 1,n_sub
       ixt = ix-0.5

       wgt_xt  = MOD(ixt,1.)
       wgt_y  = MOD(iy,1.)

       ix_l = MOD(INT(ixt-wgt_xt)-1+dsize_x,dsize_x)+1
       iy_l = MOD(INT(iy-wgt_y)-1+dsize_y,dsize_y)+1 

       ix_r = MOD(ix_l +dsize_x,dsize_x)+1   
       iy_r = MOD(iy_l +dsize_y,dsize_y)+1 

       vx_intp = velx(iy_l,ix_l)*(1-wgt_xt)*(1-wgt_y) &
               + velx(iy_l,ix_r)*(1-wgt_xt)*(  wgt_y) &
               + velx(iy_r,ix_l)*(  wgt_xt)*(1-wgt_y) &
               + velx(iy_r,ix_r)*(  wgt_xt)*(  wgt_y)

       !for v-Wind defined on xt (half level in x direction)
       iyt = iy-0.5

       wgt_x  = MOD(ix,1.)
       wgt_yt  = MOD(iyt,1.)

       ix_l = MOD(INT(ix-wgt_x)-1+dsize_x,dsize_x)+1
       iy_l = MOD(INT(iyt-wgt_yt)-1+dsize_y,dsize_y)+1    

       ix_r = MOD(ix_l +dsize_x,dsize_x)+1   
       iy_r = MOD(iy_l +dsize_y,dsize_y)+1

       vy_intp = vely(iy_l,ix_l)*(1-wgt_x)*(1-wgt_yt) &
               + vely(iy_l,ix_r)*(1-wgt_x)*(  wgt_yt) &
               + vely(iy_r,ix_l)*(  wgt_x)*(1-wgt_yt) &
               + vely(iy_r,ix_r)*(  wgt_x)*(  wgt_yt)

       !save the velocities of the tracer
       traced(tracpo(1,it),tracpo(2,it),13) = vx_intp
       traced(tracpo(1,it),tracpo(2,it),14) = vy_intp
       ! update to new location
       ix_new = MOD((ix + sub_dt*vx_intp/res)-1.+FLOAT(dsize_x),FLOAT(dsize_x))+1.
       iy_new = MOD((iy + sub_dt*vy_intp/res)-1.+FLOAT(dsize_y),FLOAT(dsize_y))+1.
       ! uebergabe
       ix = ix_new
       iy = iy_new
     end do
     !save the velocities of the tracer
     traced(tracpo(1,it),tracpo(2,it),13) = vx_intp
     traced(tracpo(1,it),tracpo(2,it),14) = vy_intp
!     ix_new = MOD((ix + dt*vx_intp/res)-1.+FLOAT(dsize_x),FLOAT(dsize_x))+1.
!     iy_new = MOD((iy + dt*vy_intp/res)-1.+FLOAT(dsize_y),FLOAT(dsize_y))+1.
     ix_new_round = nint(ix_new) ! maybe mod required
     iy_new_round = nint(iy_new)

     ! save new values for next loop
     traced(tracpo(1,it),tracpo(2,it),1) = ix_new
     traced(tracpo(1,it),tracpo(2,it),2) = iy_new
     traced(tracpo(1,it),tracpo(2,it),3) = float(ix_new_round) 
     traced(tracpo(1,it),tracpo(2,it),4) = float(iy_new_round)

     traced(tracpo(1,it),tracpo(2,it),6) = start_time !timestep
     traced(tracpo(1,it),tracpo(2,it),7) = tracer_ts  !age
     traced(tracpo(1,it),tracpo(2,it),10) = it
     it = it + 1
   ENDDO
  RETURN
END SUBROUTINE update_tracer

! --------------------------------------------------------------------------------------
! CALCULATE ANGLE AND RADIUS FROM COG OF EACH TRACER 
!-----------------------------------------------------------------------------------
SUBROUTINE geometry(traced,COMx,COMy,already_tracked,max_no_of_cells,dsize_x,dsize_y)
USE cp_parameters, ONLY : max_tracer_CP

  INTEGER, INTENT(IN)       :: max_no_of_cells, dsize_x, dsize_y
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: DELTAx(max_no_of_cells,max_tracer_CP), &
                               DELTAy(max_no_of_cells,max_tracer_CP), pi
  INTEGER                   :: dx(max_no_of_cells), dy(max_no_of_cells) ! shift everything to have CP in thecenter
  INTEGER                   :: i,j
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)
  REAL, INTENT(IN)          :: COMx(max_no_of_cells),COMy(max_no_of_cells)
   pi = 2.*asin(1.)
  ! move COG into middle of center first 
  dx = INT(dsize_x)/2 - INT(COMx)
  dy = INT(dsize_y)/2 - INT(COMy)
  
  DO i = 1,max_no_of_cells ! loop trough every cp
!   IF (already_tracked(i) .gt. 0) THEN ! do only sth if there are already tracerfor the CP
    do j = 1,max_tracer_CP,1
    !DO j=1,already_tracked(i)
      ! shift to center
      traced(i,j,1) = mod(traced(i,j,1)+dx(i)-1+dsize_x,float(dsize_x))+1
      traced(i,j,2) = mod(traced(i,j,2)+dy(i)-1+dsize_y,float(dsize_y))+1

      DELTAx(i,j) = traced(i,j,1) -dsize_x/2. !(COMx(i)+dx(i))-1.
      DELTAy(i,j) = traced(i,j,2) -dsize_y/2. !(COMy(i)+dy(i))-1 !dsize_y/2.
      traced(i,j,15) = DELTAx(i,j)
      traced(i,j,17) = DELTAy(i,j)

      IF      (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .gt. 0) THEN
       traced(i,j,8) = atan(DELTAy(i,j)/DELTAx(i,j))

      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .gt. 0) THEN  ! 2nd
       traced(i,j,8) = atan(abs(DELTAx(i,j))/DELTAy(i,j)) + 1./2.*pi

      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .lt. 0) THEN  ! 3nd
       traced(i,j,8) = atan(DELTAy(i,j)/DELTAx(i,j)) + pi

      ELSE IF (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .lt. 0) THEN  ! 4nd
       traced(i,j,8) = atan(abs(DELTAx(i,j))/abs(DELTAy(i,j))) +3./2.* pi

      ELSE IF (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .eq. 0)THEN
       traced(i,j,8) = 0.

      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .eq. 0) THEN
       traced(i,j,8) =  pi

      ELSE IF (DELTAx(i,j) .eq. 0 .and. DELTAy(i,j) .gt. 0) THEN
       traced(i,j,8) = 1./2. *pi
      ELSE IF (DELTAx(i,j) .eq. 0 .and. DELTAy(i,j) .lt. 0) THEN
       traced(i,j,8) =  3./2.* pi
      END IF
      traced(i,j,5) = sqrt(DELTAy(i,j)**2.+DELTAx(i,j)**2.)
      ! shift back 
      traced(i,j,1) = mod(traced(i,j,1)-dx(i)-1+float(dsize_x),float(dsize_x))+1
      traced(i,j,2) = mod(traced(i,j,2)-dy(i)-1+float(dsize_y),float(dsize_y))+1
     END DO
!    END IF
   END DO
END SUBROUTINE geometry

! --------------------------------------------------------------------------------------
! Calculate radial velocities 
!-----------------------------------------------------------------------------------
SUBROUTINE radvel(traced, already_tracked,max_no_of_cells)
USE cp_parameters, ONLY : max_tracer_CP
  INTEGER, INTENT(IN)       :: max_no_of_cells
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,20)
  REAL                      :: pi
  INTEGER                   :: i,j
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)

  pi = 2.*asin(1.)
  DO i = 1,max_no_of_cells ! loop trough every cp
!   IF (already_tracked(i) .gt. 0) THEN ! do only sth if there are already
!   tracerfor the CP
    do j = 1,max_tracer_CP,1
      traced(i,j,19)   = traced(i,j,13)*cos(traced(i,j,8)) + traced(i,j,14) * sin(traced(i,j,8))
      traced(i,j,20)   = traced(i,j,14)*cos(traced(i,j,8)) - traced(i,j,13) * sin(traced(i,j,8))
    END DO
!    END IF
  END DO
 
END SUBROUTINE radvel

! --------------------------------------------------------------------------------------
! SORT traced 
!-----------------------------------------------------------------------------------
SUBROUTINE sort(traced_dummy,jc)

  INTEGER, INTENT(IN)      :: jc
  REAL, INTENT(INOUT)      :: traced_dummy(jc,20)
  INTEGER                  :: j,k
  REAL                     :: a(20)
 ! sort by angle for every cpID:
  ! initialize
    do j=1,jc !   count_tracer!  
      a=traced_dummy(j,:)  !save value at j 
      do k=j-1,1,-1 ! go backwarts from value below current
        ! if  val before current value is smaller than actual value
        if (traced_dummy(k,8)<= a(8)) goto 10
        ! set this value at the current position 
        traced_dummy(k+1,:)=traced_dummy(k,:)

        ! DELTA X abspeichern
      end do
      k=0
      10 traced_dummy(k+1,:)=a
    end do
! write output sorted by angle for every cold pool
 j=1
  ! updating previous tracers
  DO WHILE (j .LT. jc)

    WRITE(50,200) INT(traced_dummy(j,6)), INT(traced_dummy(j,9)), & !timestep,CP ID
                  traced_dummy(j,1), traced_dummy(j,2),traced_dummy(j,3), &!position
                  traced_dummy(j,4), traced_dummy(j,5), traced_dummy(j,10)!, &
    j = j+1
  end do
 200 FORMAT   (2(2X,I4) ,   3(2X,F11.5), 2X,F5.3, 2(2X,F11.3))

END SUBROUTINE sort

SUBROUTINE write_output(traced,max_tracers,count_tracer,timestep,tracpo,&
                       max_no_of_cells,COMx,COMy)
  WRITE(*,*) "write output"
USE  cp_parameters, ONLY :max_tracer_CP, max_age
  INTEGER, INTENT(IN)       :: count_tracer,max_tracers, timestep, &
                               max_no_of_cells !, max_tracer_CP 
  INTEGER                   :: it
  INTEGER,INTENT(IN)        :: tracpo(2,max_tracers)
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells, max_tracer_CP,20)
  REAL, INTENT(IN)          :: COMx(max_no_of_cells), COMy(max_no_of_cells)
  it=1

  150 FORMAT    (2(4X,I4),    & !timestep, age
                2X,I6, 4X,I4, & !tracer and CP ID
                2(2X,F10.5),  & !pos 1-2
                2(4X,I4),     & !rounded
                2(2X,F7.3),   & !distance and angle
                2(2X,F8.2),   & !velocity               
                2(2X,F8.2),   & !radial and tangential velocity 
                2(2X,F8.2),   & !x dist
                2(2X,F8.2))      !COG
                !1(2X,I1),     & ! merger
                !1(2X,I4))       ! precip ID
  ! updating previous tracers
  DO WHILE (it .LT. count_tracer)
    !IF (traced(tracpo(1,it),tracpo(2,it),11)  .eq. 1.) THEN  !trace only if tracer is active
    !      IF(INT(traced(tracpo(1,it),tracpo(2,it),7)) .le. max_age)then   ! output only up to  3hours
            WRITE(40,150) INT(timestep), INT(traced(tracpo(1,it),tracpo(2,it),7)),& !timestep, age
                        it, INT(traced(tracpo(1,it),tracpo(2,it),9)),& !tracer and CP ID
                        traced(tracpo(1,it),tracpo(2,it),1), traced(tracpo(1,it),tracpo(2,it),2),& !pos 1-2
                        INT(traced(tracpo(1,it),tracpo(2,it),3)), INT(traced(tracpo(1,it),tracpo(2,it),4)),& !rounded
                        traced(tracpo(1,it),tracpo(2,it),5), traced(tracpo(1,it),tracpo(2,it),8),& !distance and angle
                        traced(tracpo(1,it),tracpo(2,it),13), traced(tracpo(1,it),tracpo(2,it),14),& ! u, v Wind component
                        traced(tracpo(1,it),tracpo(2,it),19),& ! radial vel
                        traced(tracpo(1,it),tracpo(2,it),20),& ! tangential vel
                        traced(tracpo(1,it),tracpo(2,it),15),traced(tracpo(1,it),tracpo(2,it),17),& ! x and y distance  
                        COMx(tracpo(1,it)), COMy(tracpo(1,it))                                    ! center
                        !INT(traced(tracpo(1,it),tracpo(2,it),16)), INT(traced(tracpo(1,it),tracpo(2,it),9))   ! merger, prec ID
    !END IF

    !      END IF

    it = it+1
  END DO

END SUBROUTINE write_output

END PROGRAM cp_tracking
 
