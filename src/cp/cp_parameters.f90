MODULE cp_parameters
implicit none

INTEGER               :: dsize_x, dsize_y  ! domain size
INTEGER               :: max_age           ! max # of timesteps for tracking one CP
INTEGER               :: max_tracer_CP     ! # of tracer per CP
REAL                  :: res               ! resolution in m
REAL                  :: dt                ! time step in sec
CHARACTER(LEN=200)    :: odir
logical               :: rad               ! If True tracer follow radial vel
INTEGER               :: n_sub             ! # subtimesteps
END MODULE cp_parameters
