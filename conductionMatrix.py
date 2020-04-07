def conductionMatrix(lambdaa,dt,dx,Nx):
# This program creates a matrix to be used in calculating conduction of heat between 
# different	walls, ceilings, floors, and windows of a given geometry. This is hardcoded to
# map out heat conduction flow throughout these different features in the geometry of the 
# Y2E2 building only. 






################################# Imports and Setup ############################################





	# Import numpy package
	import numpy as np

	# Pre allocate A matrix
	##### This may or may not change for res house -- DEPENDENT 
	A = np.zeros([Nx,Nx])





################################## Fill out Matrix ############################################




	# Fill out first and last row of A matrix to bind it
	##### This may or may not change for res house -- DEPENDENT 
	A[0,0]     = 1-lambdaa
	A[0,1]     = lambdaa
	A[Nx-1,Nx-2] = lambdaa
	A[Nx-1,Nx-1]   = 1-lambdaa

	# Fill out all other rows of A
	##### Nx may change based on personal preference for new problem res house. -- DEPENDENT
	for i in range(Nx-2):
	    A[i+1,i+1]   = 1-2*lambdaa
	    A[i+1,i] = lambdaa
	    A[i+1,i+2] = lambdaa

	# Return the matrix
	##### May change based on upstream data potentially changing. DEPENDENT
	return A

