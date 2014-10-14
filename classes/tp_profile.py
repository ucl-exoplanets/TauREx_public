################################################
#class profile
#
# Reads in all relevant data and performs pre-processing 
# such as sorting, grid-interpolations etc. 
#
# Input: -parameter object
#        -data object
#
#
# Output: - T-P profiles
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013 
#
################################################

#loading libraries     
from base import base   
import numpy as np
import pylab as pl



class tp_profile(base):


    def __init__(self,params,data):
        
        self.params       = params
        self.transmission = self.params.fit_transmission
        self.emission     = self.params.fit_emission
        
        #constants
        self.boltzmann = 1.3806488e-23# m2 kg s-2 K-1
        
        #derived values
        self.scaleheight = self.get_scaleheight(params.planet_temp, params.planet_grav, params.planet_mu)
     
        
        
        if params.tp_var_atm:
            self.nlayers = int(params.tp_atm_levels)
            self.ngas    = int(data.ngas)
            self.pta     = self.setup_pta_grid()
            self.P       = self.pta[:,0]
            self.T       = self.pta[:,1]
            self.Z       = self.pta[:,2]
            self.X       = np.zeros((self.ngas,self.nlayers))
            self.X      += 1e-5  #setting up initial mixing ratios
            self.rho     = self.get_rho()
            
        else:           
            self.nlayers = data.nlayers
            self.pta     = data.pta
            self.P       = self.pta[:,0]
            self.T       = data.pta[:,1]
            self.Z       = self.pta[:,2]
            self.X       = data.X     
            self.ngas    = int(data.ngas)  
            self.rho     = self.get_rho(T=self.T,P=self.P)
            
        self.num_T_params = 3 #number of free temperature parameters
        
#         pl.figure(104)
#         pl.plot(self.T,np.log(self.P))
#         pl.show()
#         exit()
        
        
        if self.params.fit_emission or self.params.fit_transmission:
            self.setup_prior_bounds()
            
        if self.params.fit_emission:
            self.PARAMS,self.TPindex, self.TPcount = self.setup_parameter_grid(emission=True)    
#             self.PARAMS[4] = 1200.0
            
        if self.params.fit_transmission:
            self.PARAMS,self.TPindex, self.TPcount = self.setup_parameter_grid(transmission=True)

#         self.PARAMS[0] = 1e-5
#         self.PARAMS[1] = 1e-5
#         self.PARAMS[2] = 1400
#         self.PARAMS[3] = 1400
#         self.PARAMS[4] = 1200


#         self.bounds = [(0.0, 0.01), (0.0, 0.01), (1000.0, 1800.0), (1000.0, 1800.0), (1000.0, 1800.0), (50000.0, 500000.0), (50.0, 150.0)]
#         self.Tpriors = [(1000.0, 1800.0), (1000.0, 1800.0), (1000.0, 1800.0)]
#         self.Ppriors = [(50000.0, 500000.0), (50.0, 150.0)]
#         self.Xpriors = [(0.0, 0.01), (0.0, 0.01)]
#         print self.bounds
#         
#          
#         self.PARAMS,self.TPindex, self.TPcount = self.setup_parameter_grid(emission=True)  
#         print self.PARAMS
#         PARAMS2 = self.PARAMS
#         PARAMS2[2] = 1400
#         PARAMS2[3] = 1400
#         PARAMS2[4] = 1200
# #         
# #         print PARAMS2
# #         
#         self.T,self.P,self.X = self.TP_profile(PARAMS=PARAMS2)
#         self.rho = self.get_rho(T=self.T,P=self.P)
#             
#         pl.figure(104)
#         pl.plot(T,np.log(P))
#         
#         rho = self.get_rho(T=T, P=P)
#         pl.figure(105)
#         pl.plot(rho)
#         pl.show()
#         exit()
#         pl.figure(200)
#         pl.plot(T,P)
#         pl.yscale('log')
#         pl.gca().invert_yaxis()
#         pl.show()
#         exit()
            
                

    #class methods 
            
    def get_scaleheight(self,T_aver,surf_g,mmw):
        
        return (self.boltzmann*T_aver)/(mmw*surf_g)
        
        
    # @profile #line-by-line profiling decorator
    def setup_pta_grid(self):
        #calculate pressure, temperature, altitude grid
        
        MAX_P    = self.params.tp_max_pres
        N_SCALE  = self.params.tp_num_scale
        N_LAYERS = self.nlayers
        
        max_z    = N_SCALE * self.scaleheight
        
#         dz       = max_z / N_LAYERS
        
        #generatinng altitude-pressure array
        PTA_arr = np.zeros((N_LAYERS,3))
        PTA_arr[:,2] = np.linspace(0,max_z,num=N_LAYERS)
        PTA_arr[:,0] = MAX_P * np.exp(-PTA_arr[:,2]/self.scaleheight)
        PTA_arr[:,1] = self.params.planet_temp

        
        return PTA_arr
        
    def setup_prior_bounds(self):
        #partially to be moved to parameter file i guess
        
        self.Xpriors = [self.params.fit_X_low, self.params.fit_X_up] #lower and upper bounds for column density
        self.Tpriors = [self.params.planet_temp - self.params.fit_T_low, self.params.planet_temp + self.params.fit_T_up] #lower and upper bounds for temperature
        
        #BE REALLY CAREFUL WITH THIS PARAMETER. THIS NEEDS SANITY CHECKING 
        self.Ppriors = [[5e4,5e5],[50.0,150.0]] #lower and upper bounds for individual TP transistion (e.g. tropopause height, mesospehere height) in Pa
        
        #setting up bounds for downhill algorithm
        #this may be merged into setup_parameter_grid() later but unsure of how complex this is going to be right now
        bounds = []
        for i in xrange(self.ngas):
            bounds.append((self.Xpriors[0],self.Xpriors[1]))
        if self.transmission:
            bounds.append((self.Tpriors[0],self.Tpriors[1]))
        if self.emission:
            for i in xrange(self.num_T_params):
                bounds.append((self.Tpriors[0],self.Tpriors[1]))
            for i in xrange(self.num_T_params-1):
                bounds.append((self.Ppriors[i][0],self.Ppriors[i][1]))
            
        self.bounds = bounds
        
        

    def setup_parameter_grid(self, transmission=False, emission=False):
        #I AM SAMPSAN I AM SAMPSAN
        
        #PARAMS [abundances (ngas), temperatures (num_T_params), pressures (num_T_params-1)]
        #COUNT [no. abundances, no. temperatures, no. pressures]
        #INDEX index to PARAMS for slicing 
        
        COUNT  = []
        PARAMS = []
        INDEX  = []
        
        #setting up mixing ratios for individual gases
        Xmean = np.mean(self.Xpriors)
        cgas = 0
        for i in xrange(self.ngas):
            PARAMS.append(Xmean)
            cgas += 1
        
        COUNT.append(cgas)
        
        #setting up temperature parameters
        T_mean = np.mean(self.Tpriors)
        num_T_params = self.num_T_params
        
        ctemp = 0; cpres = 0
        if transmission:
            ctemp +=1
            PARAMS.append(self.params.planet_temp)
        if emission:
            for i in xrange(num_T_params):
                PARAMS.append(T_mean)
                ctemp += 1
            for i in xrange(num_T_params-1):
                PARAMS.append(np.mean(self.Ppriors[i]))
                cpres += 1
                
        COUNT.append(ctemp)
        COUNT.append(cpres)

        cumidx = COUNT[0]
        INDEX.append(cumidx)
        for i in xrange(1,len(COUNT)):
            cumidx += COUNT[i]
            INDEX.append(cumidx)
    
#         print PARAMS[:INDEX[0]], PARAMS[INDEX[0]:INDEX[1]], PARAMS[INDEX[1]:]
        
        return PARAMS, INDEX, COUNT
                
        
    
    def TP_profile(self,PARAMS, T=None, P=None):
    #main function defining basic parameterised TP-profile from 
    #PARAMS and INDEX. INDEX = [Chi, T,P]
    
        INDEX = self.TPindex
        COUNT = self.TPcount
        
        X_params = PARAMS[:INDEX[0]]
        T_params = PARAMS[INDEX[0]:INDEX[1]]
        P_params = PARAMS[INDEX[1]:]

        #convert X into 2d arrays (previously done in fitting module but seems more appropriate here)
#         X   = np.zeros((self.ngas,self.nlayers))
        self.X[:] = 0.0
        for i in xrange(self.ngas):
            self.X[i,:] += X_params[i]
        
        
#         print INDEX, COUNT
#         print X_params, T_params, P_params
    
        if P is None:
            P = self.P
            
        if T is not None: 
            return T, P, self.X

        if COUNT[1] > 1:
            P_params =  [self.params.tp_max_pres] + list(P_params) + [np.min(P)]
            T_params = list(T_params) + [T_params[-1]]
#             print P_params
#             print T_params
            #creating linear T-P profile
            T = np.interp(np.log(P[::-1]), np.log(P_params[::-1]), T_params[::-1])
            return T[::-1], P, self.X
        
        if COUNT[1] == 1:
#             T = np.zeros_like(P)
            T = T_params
#             T += T_params
            return T, P, self.X
        
        
        
       # @profile #line-by-line profiling decorator
    def get_rho(self,T=None,P=None):
        #calculate atmospheric densities for given temperature and pressure
        
        if P is None:
            P = self.P
        if T is None:
            T = self.params.planet_temp
            
        return  (P)/(self.boltzmann*T)   
        
    # def cast_FIT_array(self,FIT,):     
        
        
        
        
        
            
    
    
    
    
    
        
