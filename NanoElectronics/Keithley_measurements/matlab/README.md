# Keithley_OC_IV

## Matlab codes
______________________

### Matlab is used to plot Conductance vs. distance histograms (stretching)

#### Prepacking:
Open folder Python_Filip in matlab folder explorer.
Load the script **"MATFILE_1_Creator.m"** in editor.
Change the MeasurementPath to the folder where the "Opening curves" folder of the measurement is contained.
Run the script
This will produce the packed Matlab files with all the data. They are contained in the same folder as the opening curves after running the script.
(This procedure can be also used when evaluating the files with IV measurements, MAtlab file for IV data will also be created)

#### Speed calibration (extracting the coeficient):
---------------------------------------------
Load the produced Opening.m into Matlab.
Open the **"MATFILE_2Dhist_Opening_V3_Simmons_Part1.m"** in editor. This script contains multiple segments.
Set the correct parameters in the first segment (this needs some playing around sometimes to achieve wanted results)
Run the segments until(including) "%% Estimate the slope of x vs TS".
This should produce the plot of curve data converted to distance vs time data. 
Additionally the line is fitted through this data, and coefficient is extracted from this line.

#### Plotting the conductance vs distance:
-------------------------------------
Take this coefficient and insert it into **"MATFILE_2Dhist_Opening_V3_SimmonsTEST.m"** ("dx_o").
Modify other settings in the first segment if needed.
Run the script.
This will produce the plots.
To change the colormaping, "Colormaps_2DHisto_Salen.m" can be used.
While keeping previous plot open, run the segments from the colormaps script, and they will change colors accordingly.
If you are happy with the results this can be enough, but generally it is not.

#### Data correction (!!!):
-----------------------
The Keithley 6430 produces regime switching artefacts.
Before running **"Colormaps_2DHisto_Salen"** (if it was already run before, data should be reimported), scripts **"timestamp_interpolator.m"** or **"selective_TS_interpolator.m"** can be run. (demands some tinkering with how much interpolation is desired)
After this script "Colormaps_2DHisto_Salen" can be run.
Save the matlab workspace as a m-file.


#### Additional info:
------------------
Functions that are used in creating of the initial mat files from zip files are found in:
**"MCBJ_FUNCTIONS_MATFILE.m"**
Functions definitions for extracting the speed coeficient and possible fitting to the SLM(obsolete) can be found in:
simmonscurrent.m, SimmonsModel.m and SLM.m



## Python Codes:
_____________

### Opening histograms:
-------------------

#### Option 1:
Open the **"Opening_ploter_Single.py"** and change the "empty_all = loadmat('/home/filipk/Desktop/Linux pc/Novi folder artur/Shifted data stamps/matlab_mn_shifted-80.mat')"
to the file where workspace opening data from Matlab was saved.
Running this script will plot opening histograms for one molecule. It will be normalized to the maximum of the bin intensity.

#### Option 2:
Open the **"Opening_ploter2_counts_single.py"**, this is same as Option 1, except that normalization is done on the selected cound number. To change maximum number of points per bin(maximum of the colorbar) modify the next parts:
ax.pcolormesh(X,Y,np.transpose(h2),norm=mcolors.PowerNorm(1),cmap=colormap,vmin=75,vmax=1000) ..... change vmax to some other number
also change labels of the colorbar by modyfying the ticks capture in:
cbaraa=[r'$\mathbf{0}$',r'$\mathbf{200}$',r'$\mathbf{400}$',r'$\mathbf{600}$',r'$\mathbf{800}$',r'$\mathbf{1000}$']
in other words, for vmax =2000 do this 
cbaraa=[r'$\mathbf{0}$',r'$\mathbf{400}$',r'$\mathbf{800}$',r'$\mathbf{1200}$',r'$\mathbf{1600}$',r'$\mathbf{2000}$']

#### Option 3:
-- > **"Opening_ploter2_counts_gnod.py"**
Modify loadmat location for you data.
This code needs to be modified for single compound by commenting some code!!!
This is the most recent version of the plotting code.
This routine uses two different count limits for two different region( before 1G0 and after)
Currently it is 2000 for before and 1000 for tunneling curves.
Redefine tick fonts and other proporties for you specific case.
Also modify savefig location.

### IV evaluation
--------------
Open **"Final_hist_sum_clean_v3.py"**. This program is used for plotting the E0, Gamma and Sym parameters against the transmission.
paths_* contains location of the c++ fitted files (more precisely file GOF800)
Under that....
Variable paths_all contains files paths that are going to be plotted.
names_all contains names of the compounds (this will influence some parts in function that plots the data, so things should be modified acordingly)
This file dumps in the end pickled 1d E0 bin edges and counts per bin. This file is then imported in "A_lorentzian_hist.py".
**"A_lorentzian_hist.py"** plots 1d histograms in a stacked manner and fits some function to it.

**"IVPloter_holdon_compiled_GOF.py"** plots GOF condition selected IV curves as a accumulation or all of them if the filtering is overriden.


## MISC:
____

Clasifier --- produces QLD reduced data from the fitted curves

