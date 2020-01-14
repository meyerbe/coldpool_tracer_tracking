Run tracer in 'job' repository:

>> use sh-script to call run_raintrack_bet.job


(1) run_raintrack_bet.job
> run tracer algorithm for one simulation
- adapt script: 
	JOBNAME: path to fields-data
        CIRCLES: coordinates and radius
- if necessary, compile using ./compile_bet.job (uses cp_parameters.mod)

(2) script_analyse_CP.sh & run_raintrack_loop.job
> run tracer algorithm for a set of simulations defined in script_analyse_CP.sh
- adapt sh-script: 
	- define path that is used as JOBNAME for job-script
- adapt job-script: 
	- JOBNAME should be an input
	- change CIRCLES (coordinates, radius) 




PLOTTING: 
- plot_tracer_analysis_all.py





actual FORTRAN files: 
- in src
- compile, using compile.sh
