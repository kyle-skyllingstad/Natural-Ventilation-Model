def initCondition(dt,hconv,hconv_side,Ainf,Fqi,Ttmi,Texp,Tout,qi):
# This function creates a thermal mass temperature to be used in calculation of Y2E2 air
# temperatures. This thermal mass temperature forms the initial condition upon which
# secondLoop.py will iterate through time when called by the main host function 
# newBoxModel.py.





################################# Imports and Setup ##########################################




	# Import numpy package
	import numpy as np


	# Import conductionMatrix function
	from conductionMatrix import conductionMatrix


	# Data (new values based on Fluent output 11/18/2015)
	##### Some of these will change significantly for a residential house.
	va        = 2880                   # Volume of air inside the building [m^3]
	rhoa      = 1.225                  # Density of air [kg/m^3]
	ca        = 1005                   # Specific heat of air [J/KgK]
	A1        = 1.61                   # Area of windows on floor 1 [m^2]
	A2        = 1.755                  # Area of windows on floor 2 [m^2]
	A3        = 1.755                  # Area of windows on floor 3 [m^2]
	h1        = 14.41-3.07			   # Height correction 1
	h2        = 14.41-7.79 			   # Height correction 2
	h3        = 14.41-12.34            # Stack height for windows[m]
	hinf      = 6.68                   # Infiltration equivalent stack height
	t_concr   = 0.10                   # Thickness of thermal mass floors[m]
	t_gyps    = 0.025                  # Thickness of thermal mass ceiling[m]
	Afloor1   = 184.65                 # Area of floor1 + atrium floor1 [m^2]
	Afloor2   = 186.41                 # Area of floor2 + floor office2 [m^2]
	Afloor3   = 180.86                 # Area of floor3 [m^2]
	Aceil1    = 168.91                 # Area of floor1 + atrium floor1 [m^2]
	Aceil2    = 186.41                 # Area of floor2 + floor office2 [m^2]
	Aceil3    = 180.86                 # Area of floor3 [m^2]
	Aoffice   = 555.25+455.21+489.52   # Area of side walls office
	Aout      = 21.83+41.18+49.32      # Area of side walls outside
	rho_concr = 2300                   # Density of Thermal Mass concrete [kg/m^3]
	rho_gyps  = 2320                   # Density of Thermal Mass gypsum[kg/m^3]
	k_concr   = 0.8                    # Conductivity of Thermal Mass concrete [W/mK]
	k_gyps    = 0.17                   # Conductivity of Thermal Mass gypsum [W/mK]
	cp_concr  = 750                    # Specific heat of Thermal Mass concrete [J/KgK]
	cp_gyps   = 1138                   # Specific heat of Thermal Mass gypsum [J/KgK]
	Tbase     = Texp[0]                # Basement temperature
	Toffice   = Texp[0]                # Office temperature
	Troof     = Tout[0]                # Roof temperature
	E         = 0                      # Thermal emissivity of the envelope components
	sigma     = 5.67e-8                # Stefan-Boldtzmann constant [W/m^2K^4]
	CC        = 0                      # cloudiness for clear sky


	# Conduction floors/ceilings
	##### These will change for residential house.
	Nwalls = 4   # floor1
	             # ceiling1-floor2
	             # ceiling2-floor3
	             # side walls office
	       

	# Properties of the surfaces.
	##### These will change for res house.
	keff   = (t_concr+t_gyps)/(t_concr/k_concr+t_gyps/k_gyps)		# k effective
	cpeff  = (t_concr*cp_concr+t_gyps*cp_gyps)/(t_concr+t_gyps)		# cp eff
	rhoeff = (t_concr*rho_concr+t_gyps*rho_gyps)/(t_concr+t_gyps)	# rho eff


	# Properties for four floors/ceilings, put in matrix for the 4 surfaces (nwalls). 
	##### These will change for res house bc Nwalls no longer equal to 4.
	rho = [rho_concr, rho_concr, rho_concr, rho_gyps]
	cp  = [cp_concr, cp_concr, cp_concr, cp_gyps]
	k   = [k_concr, k_concr, k_concr, k_gyps]
	ttm = [t_concr, t_concr, t_concr, t_gyps]


	# Number of gridpoints in the wall/surface.
	# 10 elements, t is the thickness. 
	##### These may or may not change for res house, based on preferences in the test.
	dx = np.zeros([4])
	dx[0] = ttm[0]/10 
	dx[1] = ttm[1]/10 
	dx[2] = ttm[2]/10 
	dx[3] = ttm[3]/10 


	# Pre allocate x-vector to count along wall/surface.
	##### This will change for res house -- DEPENDENT
	x = np.zeros([11,Nwalls]) # Hardcoded the 11 from the dimension of x in Matlab.


	# Build series of x-vectors for each wall/surface.
	##### This will change for res house -- DEPENDENT
	for i in range(Nwalls):
		x[:,i] = np.arange(0,ttm[i] + 0.0001,dx[i])
	    # Old code for clarification and visualization
	    # Old was    x[:,i] = 0:dx[i]:ttm[i]
	    
	# Background code/text.
	# Number of points in the wall/surface. 
	# number of grid points is given by the quantity Nx - 1 elements.             
	# Nx = np.size(x[:],1) for some reason, this did not give the correct value


	##### May or may not change for res house, thinking this is also probably dependent
	##### but different for res house -- MANUAL CHANGE REQUIRED. 
	Nx = 11; # This is hardcoded to match size of x.


	# Pre allocate alpha and lambda (named lambdaa b/c lambda is a function in Python)
	# These are vars seen in the eqns in the paper, they are constants.
	##### This will change for res house -- DEPENDENT
	alpha = np.zeros([len(rho)])
	lambdaa = alpha


	# Compute matrix for conduction
	# Get conduction constants preallocated right before.
	##### This will change for res house. DEPENDENT
	for i2 in range(len(rho)):
		alpha[i2] = k[i2]/(rho[i2]*cp[i2])
		lambdaa[i2] = alpha[i2]*dt/((dx[i2])**2) # The coeffiicient in front of the 
		# next-point temp value in the finite differencing scheme. 


	# Pre allocate A matrix for conduction.
	##### This will change for res house -- DEPENDENT
	A = np.zeros([Nx,Nx,Nwalls]) # Hardcoded Nx value for simplicity and increased speed


	# Carry out filling in A with lambdas as in the eqn in the paper.
	##### This will change for res house -- DEPENDENT
	for i in range(Nwalls):
	    A[:,:,i] = conductionMatrix(lambdaa[i], dt, dx[i], Nx)

	
	# Boundary conditions - the vector that the A matrix is multiplied by to give the 
	# temperature calculations downstream.
	##### This will change for res house. No ceiling 2, etc. 
	A[0,0,0]   = 1-2*lambdaa[0] # basement
	A[0,0,1]   = 1-2*lambdaa[1] # ceiling1
	A[0,0,2]   = 1-2*lambdaa[2] # ceiling2
	A[Nx-1,Nx-1,3] = 1-2*lambdaa[3] # side walls office


	# Internal internal heat flux calculation - formula from paper
	##### This will change for res house bc no floor2, floor3, etc. 
	##### this should actually be in internal heat flux but fits here for convenience.
	Qi = Fqi*(Afloor1+Afloor2+Afloor3)*qi


	# Old code and methodology for reference
	# Solar radiation heat flux  
	# qsol = 1230*exp(-1./(3.8*sin(alt+1.6)))





############################### Initial thermal mass loop ####################################




	# First loop for thermal mass initial temperature 
	# Initialize Night Flush Conditions   
	##### These will change for res house.
	Ta  = Texp[0] # Initial Temperature of Air
	T0  = [Tbase, Ttmi, Ttmi, Ttmi]
	TNx = [Ttmi, Ttmi, Ttmi, Toffice]

	# Pre allocate thermal mass temp.
	##### This will change for res house -- DEPENDENT 
	Ttm = np.zeros([Nx,Nwalls]) #hardcoded Nx value here.

	# Calculate thermal mass temperatures.
	##### This will change for res house. Probably dependent
	for i in range(Nwalls):
	    Ttm[:,i] = T0[i] + (TNx[i]-T0[i])*x[:,i]/ttm[i] # Initial Temperature of Thermal Mass
	


	# Calculate initial heat fluxes.
	# Run a pre-simulation to fix and optimize initial conditions.
	# The wall temperature distribution is unstable, so we use short time pre simulation 
	# to settle it down. 
	##### These will likely change for res house. B/c using pressure-driven flow instead of 
	##### buoyancy-driven.
	for t in range(len(qi)-1):

	    # Initial heat flux due to infiltration
	    ##### Assume infiltration is zero for the first modification.
	    To = Tout[t]
	    if (Ta > To):    
	        vinf = Ainf*(2*9.81*hinf*(1-To/Ta))**0.5   
	    else: 
	        vinf = 0.0
	    
	    # Calculate initial heat flux from infiltration -- eqn from paper
	    ##### This will change for res house -- DEPENDENT
	    Qinf = rhoa*ca*vinf*(Ta-To)  
	    

	    # Initial convective heat fluxes
	    ##### These will change for res house.
	    qconv_1f  = hconv*(Ttm[Nx-1,0] - Ta) # first floor
	    qconv_1c  = hconv*(Ttm[0, 1] - Ta) # first ceiling
	    qconv_2f  = hconv*(Ttm[Nx-1,1] - Ta) # second floor
	    qconv_2c  = hconv*(Ttm[0, 2] - Ta) # second ceiling        
	    qconv_3f  = hconv*(Ttm[Nx-1,2] - Ta) # third floor
	    qconv_off = hconv_side*(Ttm[0, 3] - Ta) # side walls office to inside

	    
	    # Multiply little q's by areas to get total intitial convective heat flux.
	    ##### These will change for res house.
	    Qconvi = Afloor1*qconv_1f  + Aceil1*qconv_1c + \
	             Afloor2*qconv_2f  + Aceil2*qconv_2c + \
	             Afloor3*qconv_3f  + Aoffice*qconv_off 
	     

		         ## Walls boundary conditions
	    # x = 0


	    # Initialize and fill out b matrix. This specifies the walls' boundary conditions.
	    # The variable b here contains the first and last element of each wall for the 
	    # boundary conditions. There are 4 walls. Hence b is 4 walls x 2 endpoints for each.
	    ##### This will change for res house.
	    b = np.zeros([Nx,len(lambdaa)])

	    b[0,0]  =  lambdaa[0]*Tbase
	    b[0,1]  =  lambdaa[1]*Ta # Old included: + qrad_1c)
	    b[0,2]  =  lambdaa[2]*Ta # Old included: + qrad_2c)
	    b[0,3]  = -lambdaa[3]*dx[3]/k[3]*(qconv_off)# old included: + qrad_off)

	    # Equality x = tm
	    b[Nx-1,0] = -lambdaa[0]*dx[0]/k[0]*(qconv_1f) # Old included: + qrad_1f)
	    b[Nx-1,1] = -lambdaa[1]*dx[1]/k[1]*(qconv_2f) # Old included: + qrad_2f)
	    b[Nx-1,2] = -lambdaa[2]*dx[2]/k[2]*(qconv_3f) # Old included: + qrad_3f)
	    b[Nx-1,3] =  lambdaa[3]*Toffice



	    # Thermal mass temperature calculation
	    # Temperature inside each of the four walls, at the 10 grid areas for each.  
	    # Initial conditions for all 10 areas, for all 4 walls.
	    ##### This will change for res house -- DEPENDENT
	    for i in range(Nwalls):
	    	
	    	# Old code for comparison and validation
	    	# for isecond in range(Nx):
	    		# Ttm[isecond,i] = A[isecond,isecond,i]*Ttm[isecond,i] + b[isecond,i]

	    	##### This will change for res house. However, b/c b and other things changed. 
	    	##### DEPENDENT
	    	Ttm[:,i] = np.dot(A[:,:,i],Ttm[:,i]) + b[:,i]

	    	# Old code for clarification and visualization
	        # Old was   Ttm[:,i] = A[:,:,i]*Ttm[:,i] + b[:,i]
	    

	    ## Air temperature 
	    # Initial case
	    ##### This will change for res house. DEPENDENT
	    Ta = Ta + (Qconvi + Qi[t])*dt/(va*rhoa*ca)
	    

	# Output variable thermal mass temperature
 	##### This will change for res house. DEPENDENT
	return Ttm

