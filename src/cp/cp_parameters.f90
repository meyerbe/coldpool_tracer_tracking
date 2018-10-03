MODULE cp_parameters
implicit none

INTEGER               :: dsize_x
INTEGER               :: dsize_y 
INTEGER               :: max_tracer_CP
REAL                  :: res  ! resolution in m
REAL                  :: dt ! time step in sec
CHARACTER(LEN=100)    :: odir
END MODULE cp_parameters
