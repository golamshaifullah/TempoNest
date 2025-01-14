! Include file for example MultiNest program 'Gaussian Rings' (see arXiv:0704.3704 & arXiv:0809.3437)

module params
implicit none

! Toy Model Parameters

	!dimensionality
      	integer sdim
      	parameter(sdim=2)
      
      	!no. of modes to generate
      	integer sModes 
	parameter(sModes=2)
      
      	!width of the Gaussian profile of each ring
	double precision sw(sModes)
	data sw /0.1d0,0.1d0/
      
      	!width of the rings
      	double precision sr(sModes)
	data sr /2.d0,2.d0/
      
      	!Center of the rings. 
      	!Centers in the 1st dimensions are set here while in 
      	!all the other dimensions, they are set to 0 in main.f90
      	double precision sc(sModes,sdim)
	data sc(1,1),sc(2,1) /-3.5d0,3.5d0/
      
      	!priors on the parameters
      	!uniform priors (-6,6) are used for all dimensions & are set in main.f90
      	double precision spriorran(sdim,2)
      


! Parameters for MultiNest
	
      	!whether to do use Nested Importance Sampling
	logical nest_IS
 	parameter(nest_IS=.true.)
	
      	!whether to do multimodal sampling
	logical nest_mmodal 
 	parameter(nest_mmodal=.false.)
	
      	!sample with constant efficiency
	logical nest_ceff
 	parameter(nest_ceff=.true.)
	
      	!max no. of live points
      	integer nest_nlive
	parameter(nest_nlive=10000)
      
      	!tot no. of parameters, should be sdim in most cases but if you need to
      	!store some additional parameters with the actual parameters then
      	!you need to pass them through the likelihood routine
	integer nest_nPar 
	parameter(nest_nPar=sdim)
      
      	!seed for MultiNest, -ve means take it from sys clock
	integer nest_rseed 
	parameter(nest_rseed=-1)
      
      	!evidence tolerance factor
      	double precision nest_tol 
      	parameter(nest_tol=0.5)
      
      	!enlargement factor reduction parameter
      	double precision nest_efr
      	parameter(nest_efr=0.05d0)
      
      	!root for saving posterior files
      	character*100 nest_root
	parameter(nest_root='chains/1-')
	
	!after how many iterations feedback is required & the output files should be updated
	!note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	integer nest_updInt
	parameter(nest_updInt=1000)
	
	!null evidence (set it to very high negative no. if null evidence is unknown)
	double precision nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
      	!max modes expected, for memory allocation
      	integer nest_maxModes 
      	parameter(nest_maxModes=10)
      
      	!no. of parameters to cluster (for mode detection)
      	integer nest_nClsPar
      	parameter(nest_nClsPar=2)
      
      	!whether to resume from a previous run
      	logical nest_resume
      	parameter(nest_resume=.true.)
      
      	!whether to write output files
      	logical nest_outfile
      	parameter(nest_outfile=.true.)
      
      	!initialize MPI routines?, relevant only if compiling with MPI
	!set it to F if you want your main program to handle MPI initialization
      	logical nest_initMPI
      	parameter(nest_initMPI=.true.)
      
      	!points with loglike < nest_logZero will be ignored by MultiNest
      	double precision nest_logZero
      	parameter(nest_logZero=-huge(1d0))
      
      	!max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
	!has done max no. of iterations or convergence criterion (defined through nest_tol) has been satisfied
      	integer nest_maxIter
      	parameter(nest_maxIter=0)
	
	!parameters to wrap around (0 is F & non-zero T)
	integer nest_pWrap(sdim)
	
      	!feedback on the sampling progress?
      	logical nest_fb 
      	parameter(nest_fb=.true.)
!=======================================================================


end module params
