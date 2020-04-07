def secondLoop(dt,Ttm_init,hconv,hconv_side,Cd01,Cd02,Cd03,Ainf,Fqi,Texp,Tout,qi):	
# This function loops through theoretical calculations of air temperatures at
# Stanford's Y2E2 building until they converge on their "guessed" values. This
# function is useful in providing these predicted values to be synthesized and
# post processed in the main host function, newBoxModel.py.





################################## Imports and Setup #########################################




	#import numpy package.
	import numpy as np


	# import conductionMatrix fcn.
	from conductionMatrix import conductionMatrix




	##### This section will change for a residential house.
	##### For a residential house, we will need to change these values. Set the volume
	##### of the house to a set reasonable value 100 m2 x 2.7 m so 270 m3 for now.
	##### It has two windows, and their openings are 1m x 0.5m for now as is standard today.
	##### Here we have wind-driven flow, not buoyancy-driven as in the Y2E2 building. 
	##### all the h's go away -- instead need to specify pressure coefficients. 



	# Data (new values based on Fluent output 11/18/2015)
	##### these will change for res house.
	va        = 2880                   # Volume of air inside the building [m^3]
	rhoa      = 1.225                  # Density of air [kg/m^3]
	ca        = 1005                   # Specific heat of air [J/KgK]
	A1        = 1.61                   # Area of windows on floor 1 [m^2]
	A2        = 1.755                  # Area of windows on floor 2 [m^2]
	A3        = 1.755                  # Area of windows on floor 3 [m^2]
	h1        = 14.41-3.07 
	h2        = 14.41-7.79 
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
	Tbase     = Texp[1]                # Basement temperature
	Toffice   = Texp[1]                # Office temperature
	Troof     = Tout[1]                # Roof temperature
	E         = 0                      # Thermal emissivity of the envelope components
	sigma     = 5.67e-8                # Stefan-Boldtzmann constant [W/m^2K^4]
	CC        = 0                      # cloudiness for clear sky


	## Conduction
	##### These will change for res house. Fewer ceilings/floors.
	Nwalls = 4  # floor1
	             # ceiling1-floor2
	             # ceiling2-floor3
	             # side walls office
	             

	# Properties (effective): k effective, cp eff, rho eff
	##### These will change for res house. May be dependent or not depending on formulae used.
	keff   = (t_concr+t_gyps)/(t_concr/k_concr+t_gyps/k_gyps)
	cpeff  = (t_concr*cp_concr+t_gyps*cp_gyps)/(t_concr+t_gyps)
	rhoeff = (t_concr*rho_concr+t_gyps*rho_gyps)/(t_concr+t_gyps)

	# Properties put into vectors for the 4 walls/surfaces
	##### These will change for res house.
	rho = [rho_concr, rho_concr, rho_concr, rho_gyps]
	cp  = [cp_concr, cp_concr, cp_concr, cp_gyps]
	k   = [k_concr, k_concr, k_concr, k_gyps]
	ttm = [t_concr, t_concr, t_concr, t_gyps]


	# Number of gridpoints in the wall
	# The variable t is thickness, 10 elements in this thickness. 
	##### These may or may not change based on preference.
	dx = np.zeros([len(ttm)])

	dx[0] = ttm[0]/10 
	dx[1] = ttm[1]/10
	dx[2] = ttm[2]/10
	dx[3] = ttm[3]/10


	# Pre allocate x counter vector.
	x = np.zeros([11,Nwalls])   # Hardcoded the 11 from the dimension of x in Matlab.

	# Fill out x vector for each wall/surface. 
	for i in range(Nwalls):
		x[:,i] = np.arange(0,ttm[i] + 0.0001,dx[i]) ##### Spatial vectors for each wall type
	    
	    # Old code for comparison and visualization purposes
	    # Old was    x[:,i] = 0:dx[i]:ttm[i]
	
	# Number of grid points is mapped by Nx - 1 elements.
	# Nx = np.size(x[:],1) did not give values to the desired accuracy.
	##### This may or may not change for res house depending on preferences -- DEPENDENT
	Nx = 11; # Hardcoded to match size of x.


	# Pre allocate constants lambda and alpha as in the equations in literature. 
	# This variable is named lambdaa because lambda is a python function.
	##### This will change for res house -- DEPENDENT
	alpha = np.zeros([Nwalls])

	# Pre-allocate lambda to be the same size as alpha. It is NOT equal. 
	lambdaa = alpha 


	# Compute matrix for conduction -- subject to future changes.
	# Get conduction constants preallocated right before.
	##### This will change for res house -- probably dependent
	for i2 in range(len(rho)):
		alpha[i2] = k[i2]/(rho[i2]*cp[i2]) # thermal diffusivity
		lambdaa[i2] = alpha[i2]*dt/((dx[i2])**2) # the coeffiicient in front of the next-point 
		# temperature value in the finite differencing scheme. 

	# Pre allocate A matrix
	##### this will change for res house -- DEPENDENT
	A = np.zeros([Nx,Nx,Nwalls])	


	# Carry out filling in conduction matrix A with lambdas as in the equation in the paper.
	##### This will change for res house -- DEPENDENT
	for i in range(Nwalls):

		# Call conductionMatrix function to create conduction matrix A
	    A[:,:,i] = conductionMatrix(lambdaa[i], dt, dx[i], Nx)
	

	# Boundary conditions - like the vector that the A matrix is multiplied by to give the
	# temperatures 
	##### This will change for res house. 
	A[0,0,0]   = 1-2*lambdaa[0] # basement
	A[0,0,1]   = 1-2*lambdaa[1] # ceiling1
	A[0,0,2]   = 1-2*lambdaa[2] # ceiling2
	A[Nx-1,Nx-1,3] = 1-2*lambdaa[3] # side walls office


	## Internal heat flux calculations begin.
	##### This will change for res house. 
	Qi = Fqi*(Afloor1+Afloor2+Afloor3)*qi

	## Initial condition specification
	##### This will change for res house. Probably dependent 
	Ta         = np.zeros(np.size(qi))
	Ta[0]      = Texp[0]
	Ttm        = np.zeros([Nx,Nwalls,len(qi)])
	Ttm[:,:,0] = Ttm_init
	b          = np.zeros([Nx,Nwalls])





################################# Calculation of Heat Fluxes #################################




	# Calculate heat fluxes.
	##### These will change for res house. That is becuase we will be using pressure-driven 
	##### flow instead of infiltration.
	for t in range(len(qi)-1):

	    # Heat flux out due to natural ventilation 
	    # Night Flush Enable = 23 C, Disable = 20 C
	    ##### This will change for res house, see below
	    To = Tout[t]
	    if (Ta[t] > To):  

	    	# Natural ventilation heat flux
	    	##### This would change if switched to wind-driven flow and not buoyancy-driven. 
	    	##### Look at chapter 4 in the natural ventilation book. Link in the CEE 251A 
	    	##### class for the book. Purely wind driven - page 103. CD coefficient, not Cz. 
	    	##### The variable vnv is the flow rate that changes.
	        vnv = Cd01*A1*(2*9.81*h1*(1-To/Ta[t]))**0.5 + \
	              Cd02*A2*(2*9.81*h2*(1-To/Ta[t]))**0.5 + \
	              Cd03*A3*(2*9.81*h3*(1-To/Ta[t]))**0.5     
	    else: 
	        vnv = 0.0
	    
	    # Calculate natural ventilation heat transfer -- areas times little q's.
	    ##### This will change for res house. DEPENDENT
	    Qnv = rhoa*ca*vnv*(Ta[t]-To)  
	    


	    ## Heat flux due to infiltration
	    ##### This will change for res house. 
	    if (Ta[t] > To):   
	    	##### Assume infiltration is zero for now. set this to 0.  
	        vinf = Ainf*(2*9.81*hinf*(1-To/Ta[t]))**0.5   
	    else: 
	        vinf = 0.0
	    
	    # Calculate heat flux due to infiltration.
	    ##### This will change for res house -- DEPENDENT
	    Qinf = rhoa*ca*vinf*(Ta[t]-To)  



	    ## Convective heat fluxes
	    ##### This will change for res house -- 1 floor, one ceiling, and some walls.
	    qconv_1f  = hconv*(Ttm[Nx-1,0,t] - Ta[t]) # first floor
	    qconv_1c  = hconv*(Ttm[0, 1,t] - Ta[t]) # first ceiling
	    qconv_2f  = hconv*(Ttm[Nx-1,1,t] - Ta[t]) # second floor
	    qconv_2c  = hconv*(Ttm[0, 2,t] - Ta[t]) # second ceiling        
	    qconv_3f  = hconv*(Ttm[Nx-1,2,t] - Ta[t]) # third floor
	    qconv_off = hconv_side*(Ttm[0, 3,t] - Ta[t]) # side walls office to inside


	    # Convective heat flux -- area times little q's.
	    ##### This will change for res house 
	    Qconvi = Afloor1*qconv_1f  + Aceil1*qconv_1c + \
	             Afloor2*qconv_2f  + Aceil2*qconv_2c + \
	             Afloor3*qconv_3f  + Aoffice*qconv_off



	    ## Walls boundary conditions
	    ##### This will change for res house. 
	    # x = 0
	    b[0,0]  =  lambdaa[0]*Tbase
	    b[0,1]  =  lambdaa[1]*Ta[t]#  + qrad_1c)
	    b[0,2]  =  lambdaa[2]*Ta[t]#  + qrad_2c)
	    b[0,3]  = -lambdaa[3]*dx[3]/k[3]*(qconv_off)# + qrad_off)

	    # x = tm
	    b[Nx-1,0] = -lambdaa[0]*dx[0]/k[0]*(qconv_1f)# + qrad_1f)
	    b[Nx-1,1] = -lambdaa[1]*dx[1]/k[1]*(qconv_2f)# + qrad_2f)
	    b[Nx-1,2] = -lambdaa[2]*dx[2]/k[2]*(qconv_3f)# + qrad_3f)
	    b[Nx-1,3] =  lambdaa[3]*Toffice


	    # Pre-allocate average thermal mass temperature.
	    ##### This will change for res house -- DEPENDENT
	    Ttm_avg = np.zeros([Nwalls,len(qi)-1])


	    # Thermal mass temperature calculation. 
	    ##### This will change for res house -- DEPENDENT
	    for i in range(Nwalls):
	    	# Old code for comparison and visualization purposes
	    	# for isecond in range(Nx):

	    		# Old code for comparison and visualization purposes
	    		# Ttm[isecond,i,t+1] = A[isecond,isecond,i]*Ttm[isecond,i,t] + b[isecond,i]

	    	##### This will change for res house. However, b/c b and other things changed. 
	    	##### DEPENDENT
	        Ttm[:,i,t+1] = np.dot(A[:,:,i],Ttm[:,i,t]) + b[:,i];
	        Ttm_avg[i,t] = np.average(Ttm[:,i,t])
	    

	    ## Air temperature
	    ##### This will change for res house. DEPENDENT
	    Ta[t+1] = Ta[t] + (Qconvi - Qnv + Qi[t])*dt/(va*rhoa*ca)


	##### This will change for res house. DEPENDENT
	return Ta, Ttm_avg

