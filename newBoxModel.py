# This is the main host program of the Y2E2 ventilation model bundle.

# The goal of this bundle is to predict air temperatures inside Stanford's Y2E2 
# building using thermal modeling, and then compare these predicted theoretical
# air temperatures to experimental air temperatures recorded from sensors
# mounted just on the outside of the building.

# This program, newBoxModel.py, is the main host function that runs as long as 
# all other functions are included in the python bundle. A few data collections
# are necessary for the computations as well. Sections of this program and the
# others are labeled whether or not they would require modification in the case
# where this bundle were to be adapted to theoretically predict air temperatures
# inside and around residential homes instead of Stanford's Y2E2 building.





############################## Experimental Data and Setup ###################################




# Import python add-ons
import numpy as np
import math
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt

# Import other python functions
from experimentalData import experimentalData 
from initCondition import initCondition
from secondLoop import secondLoop



# Prcoess xperimental data
##### This will change as a result of residential house.

# The marker num is the date (4/30)
##### This could change for residential house
num = 430

# Time interval (s)
##### This could chanve for residential house
dt  = 15

##### This will change for res house. Depends on upstream information though -- DEPENDENT
# Call experimentalData function to carry this out.
[Texp, Tout, qi, Texp_pre, Tout_pre, qi_pre, Tspan_h, path, date] = experimentalData(num, dt)



## Theoretical Thermal Box model Input Parameters
# Uncertain parameters
    # Heat transfer coefficient h [W/m2K]
    # Discharge coefficient Cd,0 [-]
    # Infiltration heat flux [W]
    # Internal loads [-]
    # Thermal mass initial temperature [K]


##### Set initial values for these parameters in the heat equations. Set these values to the 
##### values prescribed below.





##################### Uncertainty quantification Analysis Section ############################




##### This may or may not change for res house.
with open('todo_03.csv', 'rU') as file2:
	reader2 = csv.reader(file2)
	todo_03list = list(reader2)

todo_03list = zip(*todo_03list)

# Hardcode array of 6 columns
todo_03columnlen = np.size(todo_03list)/6
todo = np.zeros((todo_03columnlen,6))

todo[:,0]  = [float(i) for i in todo_03list[0]]
todo[:,1]  = [float(i) for i in todo_03list[1]]
todo[:,2]  = [float(i) for i in todo_03list[2]]
todo[:,3]  = [float(i) for i in todo_03list[3]]
todo[:,4]  = [float(i) for i in todo_03list[4]]
todo[:,5]  = [float(i) for i in todo_03list[5]]



# Load and carry out todo stuff. It is not happening in this current model.
##### This may or may not change for res house -- DEPENDENT.
# todo = load(fullfile(path,'todo_03.dat'))


# Use set coefficients and terms, instead of UQ Analysis for this case
##### This will change for res house.
dt   = 15
Nuq  = 5
Nrun = todo_03columnlen
Nqoi = 17

uq_hconv_fc = 2.5     #todo(:,2)
uq_Cd0      = 0.5     #todo(:,3)
uq_Fqi      = 1       #todo(:,4)
uq_Ttmi     = Texp[0] #todo(:,5)
uq_hconv_s  = 3     #todo(:,6)





############################ Temperature Calculation Section #################################




# Loop through UQ params and implement them.
##### Right now, we are using set valyes and not really looking at the UQ analysis
##### quite yet. Loop therefore has length of 1. 
##### This will change for res house -- DEPENDENT

# Begin the loop which is one for now.
for iuq in range(1):#:Nrun
    ## Uncertain heat transfer model parameters 
    # Print(['UQ case ',np.num2str(iuq)])


    ##### This will change for res house -- DEPENDENT
    hconv_fc = uq_hconv_fc # was uq_hconv_fc[iuq]
    hconv_s  = uq_hconv_s # was uq_hconv_s[iuq]
    Cd01     = uq_Cd0 # was .....[iuq]
    Cd02     = uq_Cd0 # was .....[iuq]
    Cd03     = uq_Cd0 # was .....[iuq]
    Ainf     = 0 
    Fqi      = uq_Fqi # was .....[iuq]
    Ttmi     = uq_Ttmi # was .....[iuq]


    # First loop for the initial temperature of the thermal mass to set up time iteration.
    ##### This will change for res house -- DEPENDENT
    # Call initCondition function to carry this out.
    Ttm_init = initCondition(dt, hconv_fc, hconv_s, Ainf, Fqi, Ttmi, Texp_pre, Tout_pre, qi_pre);
    

    ## Second loop to iterate across time to calculate theoretical Y2E2 air temperatures.
    ##### This will change for res house -- DEPENDENT
    # Call secondLoop function to carry this out.
    [Ta, Ttm_avg] = secondLoop(dt, Ttm_init, hconv_fc, hconv_s, Cd01, Cd02, Cd03, Ainf, Fqi, Texp, Tout, qi);
    

# Calculate errors
##### This will change for res house -- DEPENDENT 
maxerror = max(abs(Texp-Ta))
maxerror = round(maxerror,4) # Cleaner value to be displayed.





#################################### Figure Section ##########################################




# Plot Temp comparison.
##### This will change for res house -- DEPENDENT 
fig1 = plt.figure("Perturb Plot")
plt.plot(Tspan_h[0::57], Ta[0::57], 'ko')
plt.plot(Tspan_h[0::57], Texp[0::57], 'ro')

# Old code for modification purposes
# plt.legend(('nominal','perturbed','measured'))
# text(1, 293.2, sprintf('max error = %6.4f [K]', max(abs(Texp-Ta))), 'fontsize', 18)
# set(gca,'xtick',[0 max(Tspan_h)],'FontSize',25) !!!!!!!!!!!!
# set(gca,'ytick',[289 298],'FontSize',25) !!!!!!!!!!!!!

# Format currently-used figure
plt.text(0.4, 293.2, 'max error = ')
plt.text(1, 293.2, str(maxerror))
plt.xlabel('Time [h]')
plt.ylabel('Air Temp. [K]')

# Old title code.
# title('Perturb $$h_{side}$$ by $$20\%$$', 'interpreter', 'latex')

# Formatting plot
plt.title('Perturb $h_{side}$ by 20%')
plt.axis([0, 4, 293, 296])
plt.legend(('Air Temp.','Experimental Temp.'))
plt.show() # Shows all plots in this script
plt.close('all') # Close all plots open in this script that are not showed 

# Old preference to be turned on/off.
# plt.set(gca,'fontsize', 18)

