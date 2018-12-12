MODULE cp_parameters
implicit none

INTEGER               :: dsize_x, dsize_y  ! domain size
INTEGER               :: max_age           ! max # of timesteps for tracking one CP
INTEGER               :: max_tracer_CP     ! # of tracer per CP
REAL                  :: res               ! resolution in m
REAL                  :: dt                ! time step in sec
CHARACTER(LEN=100)    :: odir
logical               :: rad               ! If True tracer follow radial vel
END MODULE cp_parameters
