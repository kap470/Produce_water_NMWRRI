Authors: Kevin Perez, Babak Bahaddin
New Mexico Water Resources Research Institute 
Produce water model version 1.0. 
Written  in Fortran 90.
Compiled with gfortran.

***************************************
RUN THE MODEL
***************************************
Run from any terminal the file ./Program.exe

***************************************
INPUT FILES
***************************************
1. Injection_wells.csv
2. Oil_wells.csv
3. Input_file.txt

***************************************
VALUES FOR Input_file.txt
***************************************
1. ow
Number of oil wells in the file    
2. iw
Number of Injection wells in the file   
3. tini  
Initial time of the simulation
4. tend  
End time of the simulation
5. dt    
size of each time step
6. Ncols
Number of columns of the spatial domain 
7. Nrows
Number of rows of the spatil domain 
8. Xcor
X Lower left coordinate of the spatial domian  
9. Ycor  
Y Lower left coordinate of the spatial domian
10. delta
Size of each pixel in the spatial domain
11. nan
No-data value   
12. cost  
Unitary value of transport produce water
13. OWR   
Rate at which the outflow of Oil wells stock becomes produce water
14. OWS0  
Initial value of oil wells stock (same value for all oil wells)
15. INJS0 
Initial value of injection wells stock (same value for all injection wells)

***************************************
OUTPUT FILES
***************************************
1. OW_OUTF.csv
File with outflow values for each oil well and for each time step 
2. OW_STOCK.csv
File with stock values for each oil well and for each time step 
3. OW_PW.csv
File with produce water values for each oil well and for each time step 
4. IJW_CAP.csv
File with stock values for each injection well and for each time step 
5. IJW_INF.csv
File with inflow values for each injection well and for each time step 
6. INJ_RC.csv
File with remaining capacity values for each injection well and for each time step 
7. Matrix_OilWells_InjWells_debug.csv
This matrix is usfeul to debug the code 
8. IJW_MASK.asc
Raster map with locations of injection wells (View in any GIS software)
9. OW_Mask.asc
Raster map with locations of oil wells (View in any GIS software)
10. VIsualization_Results.xlsm
Excel macro to view the results
