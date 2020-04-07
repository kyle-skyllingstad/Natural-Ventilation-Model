def experimentalData(num,dt):
# This function is called by the main host function newBoxModel.py to read in experimental
# air temperature data from Stanford's Y2E2 building to be compared to the theoretical
# Y2E2 air temperature predictions.





################################# Imports and Setup ##########################################




	# import Python add-on packages
	import numpy as np
	from scipy import interpolate
	import scipy.io
	import csv


	# import internalHeatFlux Subfunction
	from internalHeatFlux import internalHeatFlux

	# Given the date in number
	# Gives the experimental data needed by the box model



	# Open 2013 experimental data -- use csv reader module
	##### This may or may not change for res house depending on data sampling.
	with open('oneyear_data_K.csv', 'rU') as file1:
		reader1 = csv.reader(file1)
		oneyeardatalist = list(reader1)
	
	# Put data into usable list
	oneyeardatalist = zip(*oneyeardatalist)


	# Oneyear_data_K = scipy.io.loadmat('oneyear_data_K.mat') does not open data in the 
	# desired format.
	# The old function, for comparison purposes, was:
	# load ../../skyllink/Documents/oneyear_data_K
	# sorted(mat_contents.keys())
	# print(oneyear_data_K)

	# Extract and format time values.
	time_year  = [float(i) for i in oneyeardatalist[3]]

	# Get one year's exp data for each of the three floors in Y2E2. NOT 3 different years.
	##### Keep floor 1 for oneyear data for floor of Catherine's house, but no other floors.
	Texp_year1 = [float(i) for i in oneyeardatalist[16]]
	Texp_year2 = [float(i) for i in oneyeardatalist[17]]
	Texp_year3 = [float(i) for i in oneyeardatalist[18]]


	# Average three floors' experimental data; it is a fair estimation.
	##### Keep Texp_year as if all three floors were averaged
	Texp_year  = map(sum, zip(Texp_year1,Texp_year2,Texp_year3))
	Texp_year = [x/3 for x in Texp_year]


	# Get outside T for each of the three floors.
	##### Keep Tout 2013 of first floor and discard the rest as for Texp.  
	Tout_year1 = [float(i) for i in oneyeardatalist[19]]
	Tout_year2 = [float(i) for i in oneyeardatalist[20]]
	Tout_year3 = [float(i) for i in oneyeardatalist[21]]


	# Average three floors out; it is a fair estimation as for Texp.
	##### Keep Tout_year as if all three floors were averaged, just like exp T was.
	Tout_year  = map(sum, zip(Tout_year1,Tout_year2,Tout_year3))
	Tout_year = [x/3 for x in Tout_year]


	# Old code for testing and verification purposes
	# Old was: time_year  = oneyear_data_K.Time
	# Old was: Texp_year1 = oneyear_data_K.D1
	# Old was: Texp_year2 = oneyear_data_K.D2
	# Old was: Texp_year3 = oneyear_data_K.D3
	# Old was: Texp_year  = (oneyear_data_K.D1 + oneyear_data_K.D2 + oneyear_data_K.D3)/3
	# Old was: Tout_year  = (oneyear_data_K.O1 + oneyear_data_K.O2 + oneyear_data_K.O3)/3






############################### Data Extraction Section ######################################





	# Extract data by month and day of month
	##### Again, this may or may not change for res house b/c data choice is our preference.
	##### Month and day and year could change or not, up to testers.
	if (num == 909): # 9 September 

		# !!!!! YOU WILL NEED TO CHANGE THIS PATH TO RUN THIS CODE !!!!!
	    # !!!!! Make this the path where you save the two necessary csv files to run, !!!!!
	    # !!!!! todo_03 and oneyear_data_K !!!!!
	    path = '../../skyllink/Documents/Job Sample Code/Gorle'
	    date = '09/09'

	    # Old code for testing
	    # angles
	    # angles = xlsread(fullfile(path,'myAngles.xlsx'));

	    # 4pm-7pm (All windows closed)
	    Texp_pre  = Texp_year[6039:6043]
	    Tout_pre  = Texp_year[6039:6043]
	    hours_pre = time_year[6039:6043] 

	    # Old code for testing  
	    # alt_pre   = angles[17:20,2]*pi/180

	    # 7pm-11pm (All windows opened)
	    Texp  = Texp_year[6042:6047]
	    Texp1 = Texp_year1[6042:6047]
	    Texp2 = Texp_year2[6042:6047]
	    Texp3 = Texp_year3[6042:6047]    
	    Tout  = Tout_year[6042:6047]
	    hours = time_year[6042:6047]   
	    day   = 1 # Monday

	    # Old code for testing
	    # alt   = angles[20:24,2]

	elif (num == 430):
	    
	    # 30 April 
	    # !!!!! YOU WILL NEED TO CHANGE THIS PATH TO RUN THIS CODE !!!!!
	    # !!!!! Make this the path where you save the two necessary csv files to run, !!!!!
	    # !!!!! todo_03 and oneyear_data_K !!!!!
	    path = '../../skyllink/Documents/Job Sample Code/Gorle' 
	    date = '04/30'

	    # Old code for testing 
	    # angles
	    # angles = xlsread(fullfile(path,'myAngles.xlsx'));

	    # 5pm-8pm (All windows closed)
	    Texp_pre  = Texp_year[2872:2876]
	    Tout_pre  = Texp_year[2872:2876]
	    hours_pre = time_year[2872:2876]  

	    # Old code for testing
	    # alt_pre   = angles[17:20,2]*pi/180

	    # 8pm-112m (All windows opened)
	    Texp  = Texp_year[2875:2880]
	    Texp1 = Texp_year1[2875:2880]
	    Texp2 = Texp_year2[2875:2880]
	    Texp3 = Texp_year3[2875:2880]
	    Tout  = Tout_year[2875:2880]
	    hours = time_year[2875:2880]
	    day   = 2 # Tuesday   
	    
	elif (num == 503):
	    
	    # 3 May 
	    # For old code for testing: !!!!! YOU WILL NEED TO CHANGE THIS PATH LIKE THE OTHERS
	    path = '../postProcessing/2013/5-03'
	    date = '05/03'
	    # Old code for testing
	    # angles
	    # angles = xlsread(fullfile(path,'myAngles.xlsx')); 

	    # 5pm-8pm (All windows closed)
	    Texp_pre  = Texp_year[2944:2948]
	    Tout_pre  = Texp_year[2944:2948]
	    hours_pre = time_year[2944:2948]  

	    # Old code for testing
	    # alt_pre   = angles[17:20,2]*pi/180

	    # 8pm-112m (All windows opened)
	    Texp  = Texp_year[2947:2952]
	    Texp1 = Texp_year1[2947:2952]
	    Texp2 = Texp_year2[2947:2952]
	    Texp3 = Texp_year3[2947:2952]
	    Tout  = Tout_year[2947:2952]
	    hours = time_year[2947:2952]
	    day   = 5 # Friday
	        
	elif (num == 703):
	    
	    # 3 July
	    # For old code for testing: !!!!! YOU WILL NEED TO CHANGE THIS PATH LIKE THE OTHERS
	    path = '../postProcessing/2013/7-03'
	    date = '07/03'

	    # Old code for testing
	    # angles
	    # angles = xlsread(fullfile(path,'myAngles.xlsx'))

	    # 5pm-8pm (All windows closed)
	    Texp_pre  = Texp_year[4408:4412]
	    Tout_pre  = Texp_year[4408:4412]
	    hours_pre = time_year[4408:4412]  

	    # Old code for testing 
	    # alt_pre   = angles[17:20,2]*pi/180

	    # 8pm-12am (All windows opened)
	    Texp  = Texp_year[4411:4416]
	    Texp1 = Texp_year1[4411:4416]
	    Texp2 = Texp_year2[4411:4416]
	    Texp3 = Texp_year3[4411:4416]
	    Tout  = Tout_year[4411:4416]
	    hours = time_year[4411:4416]   
	    day   = 3 # Wednesday 
	




############################## Process and Formulate Data ####################################




	## Nightflush
	# Time step and simulation time (~4 hours of natural ventilation)


	# Make hour span vector for temps
	##### May or may not change for res house -- up to us -- DEPENDENT
	Tspan_h = np.arange(0,len(hours)-1 + 0.0001,dt/3600.0)

	# Old code for testing
	# Old was: Tspan_h = 0:dt/3600:len(hours)-1
	

	# Make hour span counter vector
	##### May or may not change for res house -- up to us -- DEPENDENT
	hspan   = np.linspace(19, 23, len(Tspan_h))


	# Experimental data 1 night (~4 hours of natural ventilation)
	##### May or may not change for res house -- up to us -- DEPENDENT
	Texpin = np.linspace(0,len(hours)-1,len(hours))
	Texpinout = np.linspace(0,len(hours)-1,len(hours))


	# Make Texp input vectors.
	##### This will change for res house   b/c keeping only f1rst floor -- DEPENDENT
	Texpin1 = Texpin
	Texpin2 = Texpin
	Texpin3 = Texpin


	# Splinefit Texp to data to get curves. One curvefit for each floor of Y2E2.
	##### This will change to one floor for Catherine's house test case because
	##### only 1 floor will be kept -- DEPENDENT
	# Prep spline fit interpolation fit
	tck = interpolate.splrep(Texpin, Texp)
	tck1 = interpolate.splrep(Texpin1, Texp1)
	tck2 = interpolate.splrep(Texpin2, Texp2)
	tck3 = interpolate.splrep(Texpin3, Texp3)

	# Fit the spline
	Texp = interpolate.splev(Tspan_h, tck)
	Texp1 = interpolate.splev(Tspan_h, tck1)
	Texp2 = interpolate.splev(Tspan_h, tck2)
	Texp3 = interpolate.splev(Tspan_h, tck3)


	# Old code for testing
	# Old was: Texp = interp1(Texpin, Texp, Tspan_h,'spline') #!!!!!!
	# Old was: Texp1 = interp1(Texpin1, Texp1, Tspan_h,'spline') #!!!!!!
	# Old was: Texp2 = interp1(Texpin2, Texp2, Tspan_h,'spline') #!!!!!!
	# Old was: Texp3 = interp1(Texpin3, Texp3, Tspan_h,'spline') #!!!!!!



	# Outside polynomial temperature. Fit Tout curve based on Tout data.
	##### May or may not change based on upstream data fit -- DEPENDENT
	Tout_coeff = np.polyfit(Texpin, np.transpose(Tout), 3) 
	Tout         = np.polyval(Tout_coeff, Tspan_h) 


	# Old code for testing
	# Solar radiation angles
	# alt = interp1(0:len(hours)-1, alt, Tspan_h,'spline') #!!!!!!


	# Internal heat flux
	# Create counter vector for I think a time interval.
	##### May or may not change for res house -- up to us
	int1 = np.linspace(19,23,5)

	# Call internalHeatFlux over the specified time interval to get heat transfer 
	##### Coeff -- I think depends on geometry, so will change for Catherine's house.
	##### This will have to be modified to fit the case where it is a residential home
	##### -- DEPENDENT
	qi_coeff = internalHeatFlux(day, int1, 3)
	qi         = np.polyval(qi_coeff[:,0], hspan) 



	# Initial condition
	# Create initial cond temp counter vector based on hours.
	##### May or may not change for res house -- up to us -- DEPENDENT
	Tspan_h_pre = np.arange(0,len(hours_pre)-1 + 0.0001,dt/3600.0)

	# Old code for testing
	# Old was: Tspan_h_pre = 0:dt/3600:len(hours_pre)-1 


	# Create init cond counter vector for time hours. 
	##### May or may not change for res house -- up to us -- DEPENDENT
	hspan_pre   = np.linspace(hours_pre[0],hours_pre[len(hours_pre)-1],len(Tspan_h_pre)) 


	# Experimental data for 1 night 
	# Splinefit init cond exp T curvefit to init cond exp T data. 
	##### May or may not change for res house -- up to us -- DEPENDENT
	Texpprein = np.linspace(0,len(hours_pre)-1,len(hours_pre))
	tckpre = interpolate.splrep(Texpprein, Texp_pre)
	Texp_pre = interpolate.splev(Tspan_h_pre, tckpre)
	
	# Old code for testing
	# Old was: Texp_pre = interp1(0:len(hours_pre)-1, Texp_pre, Tspan_h_pre,'spline')


	# Outside polynomial temperature
	# Fit init cond outside temp curvefit to outside T data.
	##### May or may not change for res house -- up to us -- DEPENDENT
	Tout_coeff = np.polyfit(Texpprein, np.transpose(Tout_pre), 2) 
	Tout_pre     = np.polyval(Tout_coeff, Tspan_h_pre) #!!!!!!!

	# Old code for testing
	# Solar radiation angles
	# alt_pre = interp1(0:len(hours_pre)-1, alt_pre, Tspan_h_pre,'spline');

	# Internal heat flux
	##### This will change for residential home - new geometry. -- DEPENDENT
	qi_coeff = internalHeatFlux(day, np.transpose(hours_pre), 3) 
	qi_pre     = np.polyval(qi_coeff[:,0], hspan_pre) #!!!!!!!


	##### qi will likely have to change, but keep as is for now. Pre method might 
	##### be chaning as well.

	# These will change for res house. -- DEPENDENT
	return Texp, Tout, qi, Texp_pre, Tout_pre, qi_pre, Tspan_h, path, date

	# Old variable outputs for testing and verification purposes
	#return Texp
	#return Tout
	#return qi
	#return Texp_pre
	#return Tout_pre
	#return qi_pre
	#return Tspan_h
	#return path
	#return date

