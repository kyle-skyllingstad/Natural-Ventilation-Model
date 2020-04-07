def internalHeatFlux(day,hours,deg):
# This function calculates heat flux coefficients and it is called in the theoretical 
# calculation of Y2E2 building temperatures. 





################################# Imports and Setup ##########################################




	# Import numpy package
	import numpy as np

	# days  = vector of numbers between 1 and 7 relative to the days of the week
	# hours = vector of numbers between 1 and 24 relative to the hours of the day
	# deg   = degree of the polynomial interpolation
	
	# Old code for comparison     
	# qi_coeff = matrix [gradeXdays]

	## INTERNAL HEAT FLUX ESTIMATION (Data from table 9 Hult)


	# USAGE FRACTION


	# Internl heat flux weekday lights and occupants. This corresponds to W/ft2 units.
	# Each entry orresponds to an hour of the day, hence 24 entries. 
	# Each entry is essentially a weighted value.
	# These entries/values are different on the weekends compared to the week days. 
	##### These will be different for a residential house compared to Y2E2.
	LeO_week = [0.19, 0.13, 0.08, 0.08, 0.08, 0.08, 0.19, 0.30,\
				0.55, 0.70, 0.85, 0.89, 0.93, 0.97, 1.00, 1.00,\
				0.90, 0.80, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25]
	   
	# Same methodology for the equipment inside.
	##### also different for residential home.    
	E_week = [0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50,\
	          0.62, 0.75, 0.88, 1.00, 1.00, 1.00, 1.00, 1.00,\
	          1.00, 0.83, 0.66, 0.50, 0.50, 0.50, 0.50, 0.50]


	# Saturday here -- lights and occupants.	          
	# There is a difference on saturday -- a weekend here.     
	LeO_sat = [0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.12,\
	           0.16, 0.21, 0.25, 0.30, 0.30, 0.30, 0.30, 0.30,\
	           0.3, 0.275, 0.25, 0.225, 0.2, 0.16, 0.12, 0.08]
	  

	# Same thing for the equipment inside on saturday.      
	##### This will change for residential home.  
	E_sat = 0.5*np.ones(np.size(E_week))                            
	   

	# Sunday
	# Lights and occupants.
	##### This will change for res home.
	LeO_sun = [0.10, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.12,\
	           0.16, 0.21, 0.25, 0.30, 0.30, 0.30, 0.30, 0.30,\
	           0.3, 0.275, 0.25, 0.225, 0.2, 0.18, 0.16, 0.13]
	      
	# Equipment 
	##### This will change for res house -- DEPENDENT       
	E_sun = E_sat                                            
	      

	# USAGE PEAK [W/m^2] -- convert from feet to meters for design specifications
	##### This will need to be modified for a residential house. 
	L_peak = 1/(0.3048**2) 
	O_peak = 0.5/(0.3048**2) 
	E_peak = 0.7/(0.3048**2) 





################################# Heat Flux Calculation ########################################




	# Start calculating heat coefficient and heat flux -- for internal heat calculation
	# Create heat and heat coefficient counter vectors
	##### This may or may not change for res house -- up to us -- DEPENDENT
	qi = np.zeros((len(hours),1))
	qi_coeff = np.zeros((deg+1,1))


	# HEAT FLUX [W/m^2]
	##### These will change for res house -- DEPENDENT
	for j in range(1):

	# Old code for comparison
	# old was: for j in range(len(day)): changed code from this to in range(1) 
	# This originally tried to adjust for different days.
	    for i in range(len(hours)):

	    	# Set hours to be looked at
	    	inputhours = int(hours[i])

	    	# Reset for next day
	        if (hours[i] == 0):
	            hours[i] = 24

	        # Old code for comparison and validation
	       	# print(hours[i])


	       	# Calculate internal heat fluxes based on usage ratios and peak values
	       	##### These will change for a residential house. DEPENDENT
	        if (day < 6):
	            qi[i,j] = (L_peak+O_peak)*LeO_week[inputhours-1] + (E_peak)*E_week[inputhours-1] # weekdays
	        elif (day == 6):
	            qi[i,j] = (L_peak+O_peak)*LeO_sat[inputhours-1]  + (E_peak)*E_sat[inputhours-1]  # Saturday
	        else:
	            qi[i,j] = (L_peak+O_peak)*LeO_sun[inputhours-1]  + (E_peak)*E_sun[inputhours-1]  # Sunday
	        
	    
	    # Interpolation coefficients. 
	    # These are in a curvefit that is fit to the heat flux data points in the second for loop.
	    ##### these will change for residential houses. DEPENDENT
	    qi_coeff[:,j] = np.polyfit(hours, qi[:,j], deg)
	
	   

	# Output variable.
	##### this changes for residential house. DEPENDENT
	return qi_coeff

	