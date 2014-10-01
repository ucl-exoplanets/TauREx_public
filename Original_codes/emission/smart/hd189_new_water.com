0                               Pressure not variable
0                               Temperature not variable
3						format list directed
../emis/atm/HD189_burrows3.atm
11						# records to skip at TOF
1,2						columns of p and t
1.						pressure scaling factor
1.89E+03    					surface temperature
1                                               # of gas absorbers
1                                               AFGL gas code for gas 1
0                                               Gas not variable
1                                               number of coeff. types
3                                               type of abs coef (line)
../trans/bt2_1-39_1000.txt
3                                               format  (list directed)
../emis/atm/HD189_base.atm
11                                              # records to skip at TOF
1                                               vert coordinate (pres)
1,4                                             columns of p and rmix
1                                               type - volume mix ratio
1.0,1.E-07                                          p, rmix scaling factors
18.0                                            molecular weight
0                                               # of aerosol modes
0						include lambert surface
0                                          surface reflectance non-variable
3                                               format (list directed)
albedo_0.001.dat                                
5                                               # records to skip at TOF
1
2,3                                             columns of wl and albedo
1                                               wavenumber grid (wl)
1.0                                             wl scaling factors
0., 0.                                             albedo scaling factors
0.033                                           distance from sun (au)
20.8                                          	surface gravity
9.2E+04                                          radius of planet
2.3                                             mean mol. weight of atm
1                                               # of Rayleigh scatterers
5                                               Rayleigh Scatter index
1.                                              volume mixing ratio air
8                                               number of d_ostreams
3		                        	sources - thermal
3
3
../emis/Kurucz1cm-1_susim_atlas2.dat
12
1
2
1,2
1
60., 0.
0.01
1
1
0
1
2				units (w/m**2/micron)
300.0,2500.			min, max wavenumber
2				type of grid (slit)
2				triangular
1.				slit hwhm
0.25				output sampling
0.25				tau error criterion
0.15				pi0 error criterion
0.15				g error criteron
0.02				albedo error criterion
1				output format (ascii)
spec/em189_new_water
2



                                                                                                                                    






