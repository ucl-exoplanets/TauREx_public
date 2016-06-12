#! /usr/bin/python 


###########################################################
# taurex_mlehess - running L-BFGS downhill with TauRex on
# previously best-fit parameters (from multinest) to calculate
# the hessian of the maximum likelihood 
# 
#
# Requirements: -python libraries: pylab, numpy, ConfigParser 
#             [these are the minimum requirements]
#
#
# Inputs: minimal or normal parameter file
#
# Outputs: spectrum.dat
#
# To Run: 
# ./taurex_mlehess.py [-p exonest.par] 
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Feb 2016      
#       
###########################################################


#loading libraries     
import sys, os, optparse, time, scipy
import numpy as np #nummerical array library 
import pylab as pl#science and plotting library for python
from ConfigParser import SafeConfigParser
import cPickle as pkl
import numdifftools as nd
from mpltools import special
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from algopy.utils import piv2det

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,fitting,atmosphere,data,preselector
from parameters import *
from emission import *
from transmission import *
from fitting import *
from atmosphere import *
from data import *


#loading libraries
import library_emission
import library_general
from library_emission import *
from library_general import *


class statistics(object):
    '''
    Class that will in the future calculate statistical quantities such as the 
    Hessian and Fisher information matrix of the Maximum Likelihood function 
    It may also contain some trace analysis tools (e.g. eigenvectors etc)
    '''

    def __init__(self,options=None, params=None):
        
        self.options = options
        if options is None and params is None:
            self.params = None
        else:
            if params is None:
                #initialising parameters object
                self.params = parameters(options.param_filename)
            else:
                self.params = params
        
        
    def init_stats(self):
        '''
        manually initialises the TauREx output analysis
        when not loaded distance functinos are still available
        ''' 
        
        self.dir = self.options.dir
        self.params.gen_manual_waverange = False
        self.params.nest_run = False
        self.params.mcmc_run = False
        self.params.downhill_run = True

        #setting up output storage 
        self.stats = {}
        
        #loading data from TauREx NEST output 
        self.load_traces_likelihood('nest_out.db')
        
        #loading parameter list 
        self.parameters = np.loadtxt(os.path.join(self.dir,'parameters.txt'),dtype='str')
        self.stats['parameters'] = self.parameters
        
        # initialising data object
        self.dataob = data(self.params)

        # initialising atmosphere object
        self.atmosphereob = atmosphere(self.dataob)

        # set forward model
        if self.params.gen_type == 'transmission':
            self.fmob = transmission(self.atmosphereob)
        elif self.params.gen_type == 'emission':
            self.fmob = emission(self.atmosphereob)
        
        #initialising fitting object 
        self.fitting = fitting(self.fmob)
        
    def load_traces_likelihood(self,nest_db_fname,solution_idx=0):
        '''
        loading following files from TauREx output: 
          NEST_tracedata.txt
          NEST_likelihood.txt
          NEST_out.db
        '''
        
#         self.NEST_trace = np.loadtxt(os.path.join(self.dir,trace_fname))
#         self.NEST_like  = np.loadtxt(os.path.join(self.dir,like_fname))
        with open(os.path.join(self.dir,nest_db_fname),'r') as file:
            self.NEST_db = pkl.load(file)
            
        self.NEST_trace = self.NEST_db['solutions'][solution_idx]['tracedata']
        self.NEST_like = self.NEST_db['solutions'][solution_idx]['weights']
            
        #getting maximum likelihood values 
        self.mle_values = self.NEST_trace[-1,:]
        
        print np.shape(self.NEST_trace)
        
        #calculating mean and error from traces 
        trace_std = []; trace_mean = []
        for i in range(len(self.NEST_trace[0,:])):
            tmp  = self.weighted_avg_and_std(self.NEST_trace[:,i], self.NEST_like[:,0])
            trace_mean.append(tmp[0])
            trace_std.append(tmp[1])

        
        #storing mle values 
        self.stats['mle_values'] = self.mle_values
        self.stats['trace_mean'] = trace_mean
        self.stats['mle_std']    = trace_std
        
    
    
    def compute_hessian_mle(self):
        '''
        wrapper that takes the MLE from the multinest trace 
        and starts the taurex native L-BFGS-B minimisation to 
        compute the Hessian matrix around MLE. 
        Returns: 
            Jacobian,
            Hessian, inverse Hessian, 
            Fisher information matrix (negative Hessian),
            Asymptotic Standard Error on parameters (CRLB)
        '''
        
        #setting starting parameters to MLE 
        self.fitting.fit_params = self.mle_values
        
        data = self.fitting.data.obs_spectrum[:,1] # observed data
        datastd = self.fitting.data.obs_spectrum[:,2] # data error
        
        def loglike(fit_params):
            # log-likelihood function called by multinest
            chi_t = self.fitting.chisq_trans(fit_params, data, datastd)
            loglike = np.sum(np.log(datastd*np.sqrt(2*np.pi))) - 0.5 * chi_t
            return loglike
        
        
        print 'Fine-tuning MLE'
        bounds =[]
        for bound in self.mle_values:
            tmp = bound*0.2
            bounds.append((bound-tmp,bound+tmp))
            
        MLE_fine = scipy.optimize.minimize(self.fitting.chisq_trans, 
                                           x0=self.mle_values,  
                                           method='TNC',
                                           args=(data,datastd),
                                           bounds=bounds)

        self.mle_values = MLE_fine['x']
        self.stats['mle_values'] = self.mle_values
        
        
        print 'Calculating Jacobian'
        fd = nd.Jacobian(loglike,step=1e-4)
        jacob = fd(self.mle_values)
        
        print 'Calculating Hessian'
        fdd = nd.Hessian(loglike,step=1e-4)
        hess = fdd(self.mle_values)
        
        hess_inv = np.linalg.inv(hess) #inverse hessian
        fisher = -1.0 * hess #the observed fisher information matrix is the negative hessian at MLE     
        
        #calculating assymptotic standard error of data 
        stderr_asymp = np.sqrt(np.diag(hess_inv))
          
        #store results
        self.stats['jacobian']     = jacob
        self.stats['hessian']      = hess
        self.stats['hessian_inv']  = hess_inv
        self.stats['cov_asymp']    = -1.0*hess_inv
        self.stats['stderr_asymp'] = stderr_asymp 
        self.stats['fisher']       = fisher
        
#         return MLE_fit['hess:'], MLE_fit['hess_inv']
        
        
    def compute_trace_eigenvectors(self):
        '''
        - computes the covariance from NEST_tracedata for parameters
          - computes covariance of total sample space
          - computes covariance of highest 10% likelihood space (i.e. the best fits)
        - does eigenvalue decomposition to compute eigenvector/value sets for covariances
        Return: 
            Covariance matrix, Covariance matrix of top 10% likelihood results 
            Eigenvalue/Eigenvector sets
            Rotational angles of eigenvectors @todo define basis 
        '''
        #subtracting MLE from traces
        trace_norm = self.NEST_trace - self.mle_values
#         trace_norm = trace_norm[:,:3]
        
        #getting index for top 10% 
        top_idx = int(len(trace_norm[:,0])*0.02)
        
        #getting covariance of centred traces
        postcov = np.cov(np.transpose(trace_norm))

#         print 'tracenorm ', np.shape(trace_norm)
        print 'Calculating covariance matrix'
        postcov_small = np.cov(np.transpose(trace_norm[-top_idx:,:])) 

        #LAPAC decomposition 
        print 'Calculating eigenvalue/eigenvector pairs'
        [eigval_small, eigvec_small] = np.linalg.eig(postcov_small)
        [eigval,eigvec] = np.linalg.eig(postcov)
   
        #need to figure out if angles should be w.r.t. each other 
        #or w.r.t. a basis vector
#         postangle = np.arctan2(eigvec[1,:],eigvec[0,:])
        
        #storing results
        self.stats['trace_norm']     = trace_norm
        self.stats['cov']            = postcov
        self.stats['cov_small']      = postcov_small
        self.stats['eigenval']       = eigval
        self.stats['eigenvec']       = eigvec
        self.stats['eigenval_small'] = eigval_small
        self.stats['eigenvec_small'] = eigvec_small
            

        
    #setting up storage object 
#     def setup_storage(self):
#         '''
#         small wrapper setting up storage, now just returns 
#         dictionary but may be more advanced later
#         '''
#         self.stats = {}  
    
#     def store_stats(self,type,value):
#         #small wrapper storing data in object but may be more sophisticated later         
#         
#         self.stats[type] = value
    def weighted_avg_and_std(self,values, weights):
        '''
        Function calculating weighted averages and standard-deviation
        Same function in taurex/library_general
        '''
        average = np.average(values, weights=weights)
        variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
        return (average, np.sqrt(variance))

        
    def get_stats(self,type):
        #wrapper to storage
        return self.stats[type]
    
    def dump(self,filename):
        '''
        wrapper dumping self.stats dict into pickle
        '''
        pkl.dump(self.stats,open(filename+'.db','wb'))
        
    #DISTANCE FUNCTIONS
    def distance_kullback(self,A, B):
        """Kullback leibler divergence between two covariance matrices A and B.
    
        :param A: First covariance matrix
        :param B: Second covariance matrix
        :returns: Kullback leibler divergence between A and B
    
        """
        dim = A.shape[0]
        logdet = np.log(np.linalg.det(B) / np.linalg.det(A))
        kl = np.trace(np.dot(np.linalg.inv(B), A)) - dim + logdet
        return 0.5 * kl
    
    
    def distance_kullback_right(self,A, B):
        """wrapper for right kullblack leibler div."""
        return self.distance_kullback(B, A)
    
    
    def distance_kullback_sym(self,A, B):
        """Symetrized kullback leibler divergence."""
        return self.distance_kullback(A, B) + self.distance_kullback_right(A, B)
    
    
    def distance_euclid(self,A, B):
        """Euclidean distance between two covariance matrices A and B.
    
        The Euclidean distance is defined by the Froebenius norm between the two
        matrices.
    
        .. math::
                d = \Vert \mathbf{A} - \mathbf{B} \Vert_F
    
        :param A: First covariance matrix
        :param B: Second covariance matrix
        :returns: Eclidean distance between A and B
    
        """
        return np.linalg.norm(A - B, ord='fro')
    
    
    def distance_logeuclid(self,A, B):
        """Log Euclidean distance between two covariance matrices A and B.
    
        .. math::
                d = \Vert \log(\mathbf{A}) - \log(\mathbf{B}) \Vert_F
    
        :param A: First covariance matrix
        :param B: Second covariance matrix
        :returns: Log-Eclidean distance between A and B
    
        """
        return self.distance_euclid(self.logm(A), self.logm(B))
    
    
    def distance_riemann(self,A, B):
        """Riemannian distance between two covariance matrices A and B.
    
        .. math::
                d = {\left( \sum_i \log(\lambda_i)^2 \\right)}^{-1/2}
    
        where :math:`\lambda_i` are the joint eigenvalues of A and B
    
        :param A: First covariance matrix
        :param B: Second covariance matrix
        :returns: Riemannian distance between A and B
    
        """
        return np.sqrt((np.log(scipy.linalg.eigvalsh(A, B))**2).sum())
        
    #GEODESIC FUNCTIONS
    def geodesic_riemann(self,A, B, alpha=0.5):
        """Return the matrix at the position alpha on the riemannian geodesic between A and B  :
    
        .. math::
                \mathbf{C} = \mathbf{A}^{1/2} \left( \mathbf{A}^{-1/2} \mathbf{B} \mathbf{A}^{-1/2} \\right)^\\alpha \mathbf{A}^{1/2}
    
        C is equal to A if alpha = 0 and B if alpha = 1
    
        :param A: the first coavriance matrix
        :param B: the second coavriance matrix
        :param alpha: the position on the geodesic
        :returns: the covariance matrix
    
        """
        sA = self.sqrtm(A)
        isA = self.invsqrtm(A)
        C = isA * B * isA
        D = self.powm(C, alpha)
        E = np.matrix(sA * D * sA)
        return E
        
    #TANGENT SPACE PROJECTIONS
    def tangent_space(self,covmats, Cref):
        """Project a set of covariance matrices in the tangent space according to the given reference point Cref
    
        :param covmats: Covariance matrices set, Ntrials X Nchannels X Nchannels
        :param Cref: The reference covariance matrix
        :returns: the Tangent space , a matrix of Ntrials X (Nchannels*(Nchannels+1)/2)
    
        """
        Nt, Ne, Ne = covmats.shape
        Cm12 = self.invsqrtm(Cref)
        idx = np.triu_indices_from(Cref)
        T = np.empty((Nt, Ne * (Ne + 1) / 2))
        coeffs = (
            np.sqrt(2) *
            np.triu(
                np.ones(
                    (Ne,
                     Ne)),
                1) +
            np.eye(Ne))[idx]
        for index in range(Nt):
            tmp = np.dot(np.dot(Cm12, covmats[index, :, :]), Cm12)
            tmp = self.logm(tmp)
            T[index, :] = np.multiply(coeffs, tmp[idx])
        return T
    
    
    def untangent_space(self,T, Cref):
        """Project a set of Tangent space vectors in the manifold according to the given reference point Cref
    
        :param T: the Tangent space , a matrix of Ntrials X (Nchannels*(Nchannels+1)/2)
        :param Cref: The reference covariance matrix
        :returns: A set of Covariance matrix, Ntrials X Nchannels X Nchannels
    
        """
        Nt, Nd = T.shape
        Ne = int((np.sqrt(1 + 8 * Nd) - 1) / 2)
        C12 = self.sqrtm(Cref)
    
        idx = np.triu_indices_from(Cref)
        covmats = np.empty((Nt, Ne, Ne))
        covmats[:, idx[0], idx[1]] = T
        for i in range(Nt):
            covmats[i] = np.diag(np.diag(covmats[i])) + np.triu(
                covmats[i], 1) / np.sqrt(2) + np.triu(covmats[i], 1).T / np.sqrt(2)
            covmats[i] = self.expm(covmats[i])
            covmats[i] = np.dot(np.dot(C12, covmats[i]), C12)
    
        return covmats
        
    #BASIC MATRIX MATH FUNCTIONS STOLEN FROM PyRIEMANN package 
    def sqrtm(self,Ci):
        """Return the matrix square root of a covariance matrix defined by :
    
        .. math::
                \mathbf{C} = \mathbf{V} \left( \mathbf{\Lambda} \\right)^{1/2} \mathbf{V}^T
    
        where :math:`\mathbf{\Lambda}` is the diagonal matrix of eigenvalues
        and :math:`\mathbf{V}` the eigenvectors of :math:`\mathbf{Ci}`
    
        :param Ci: the coavriance matrix
        :returns: the matrix square root
    
        """
        D, V = scipy.linalg.eigh(Ci)
        D = np.matrix(np.diag(np.sqrt(D)))
        V = np.matrix(V)
        Out = np.matrix(V * D * V.T)
        return Out
    
    
    def logm(self,Ci):
        """Return the matrix logarithm of a covariance matrix defined by :
    
        .. math::
                \mathbf{C} = \mathbf{V} \log{(\mathbf{\Lambda})} \mathbf{V}^T
    
        where :math:`\mathbf{\Lambda}` is the diagonal matrix of eigenvalues
        and :math:`\mathbf{V}` the eigenvectors of :math:`\mathbf{Ci}`
    
        :param Ci: the coavriance matrix
        :returns: the matrix logarithm
    
        """
        D, V = scipy.linalg.eigh(Ci)
        Out = np.dot(np.multiply(V, np.log(D)), V.T)
        return Out
    
    
    def expm(self,Ci):
        """Return the matrix exponential of a covariance matrix defined by :
    
        .. math::
                \mathbf{C} = \mathbf{V} \exp{(\mathbf{\Lambda})} \mathbf{V}^T
    
        where :math:`\mathbf{\Lambda}` is the diagonal matrix of eigenvalues
        and :math:`\mathbf{V}` the eigenvectors of :math:`\mathbf{Ci}`
    
        :param Ci: the coavriance matrix
        :returns: the matrix exponential
    
        """
        D, V = scipy.linalg.eigh(Ci)
        D = np.matrix(np.diag(np.exp(D)))
        V = np.matrix(V)
        Out = np.matrix(V * D * V.T)
        return Out
    
    
    def invsqrtm(self,Ci):
        """Return the inverse matrix square root of a covariance matrix defined by :
    
        .. math::
                \mathbf{C} = \mathbf{V} \left( \mathbf{\Lambda} \\right)^{-1/2} \mathbf{V}^T
    
        where :math:`\mathbf{\Lambda}` is the diagonal matrix of eigenvalues
        and :math:`\mathbf{V}` the eigenvectors of :math:`\mathbf{Ci}`
    
        :param Ci: the coavriance matrix
        :returns: the inverse matrix square root
    
        """
        D, V = scipy.linalg.eigh(Ci)
        D = np.matrix(np.diag(1.0 / np.sqrt(D)))
        V = np.matrix(V)
        Out = np.matrix(V * D * V.T)
        return Out
    
    
    def powm(self,Ci, alpha):
        """Return the matrix power :math:`\\alpha` of a covariance matrix defined by :
    
        .. math::
                \mathbf{C} = \mathbf{V} \left( \mathbf{\Lambda} \\right)^{\\alpha} \mathbf{V}^T
    
        where :math:`\mathbf{\Lambda}` is the diagonal matrix of eigenvalues
        and :math:`\mathbf{V}` the eigenvectors of :math:`\mathbf{Ci}`
    
        :param Ci: the coavriance matrix
        :param alpha: the power to apply
        :returns: the matrix power
    
        """
        D, V = scipy.linalg.eigh(Ci)
        D = np.matrix(np.diag(D**alpha))
        V = np.matrix(V)
        Out = np.matrix(V * D * V.T)
        return Out
    
    #OTHER VECTOR GEOMETRY FUNCTIONS
    def dotproduct(self,v1, v2):
      return sum((a*b) for a, b in zip(v1, v2))
    
    def length(self,v):
      return np.sqrt(self.dotproduct(v, v))
    
    def angle(self,v1, v2):
      return np.arccos(self.dotproduct(v1, v2) / (self.length(v1) * self.length(v2)))

    
    
        
        
        
#gets called when running from command line 
if __name__ == '__main__':
    
    parser = optparse.OptionParser()
    parser.add_option('-p', '--parfile',
                      dest="param_filename",
                      default="Parfiles/default.par",
    )
    parser.add_option('-d', '--dir',
                      dest="dir",
                      default="./",
    )
    parser.add_option('-o', '--outfname',
                      dest="outfname",
                      default="statistics",
    )
    parser.add_option('-c', '--cluster',
                  dest="cluster_dictionary",
                  default="None",
                  type = "string",
                  action="store"         
                  )
    parser.add_option('-i', '--cluster_index',
                  dest="cluster_procid",
                  default="None",
                  type = "string",
                  action="store"         
                  )

    options, remainder = parser.parse_args()
    #loading object
    if options.cluster_dictionary is not "None":
        params = parameters(options.param_filename)
        from cluster import cluster
        c = cluster()
        c.read_dict(dict_name=options.cluster_dictionary)
        params = c.modify_params(params,options.cluster_procid)
        statsob = statistics(options=options,params=params)
    else:
        statsob = statistics(options=options)
    statsob.init_stats()
    
    #compute hessian 
    statsob.compute_hessian_mle()
    
#     print statsob.stats['hessian']
#     print '-----------'
#     print statsob.stats['fisher']

    statsob.compute_trace_eigenvectors()
#     print statsob.stats.keys()
    
    statsob.dump(os.path.join(options.dir,options.outfname))


#Some plotting stuff 
############################
#     pl.figure(1)
# #         pl.imshow(postcov_small,interpolation='nearest')
# #         special.hinton(postcov_small,max_value=3.0)
#     special.hinton(statsob.stats['cov_small'])
#     pl.gca().invert_yaxis()
#           
#     pl.figure(2)
# #         special.hinton(postcov,max_value=3.0)
#     special.hinton(statsob.stats['cov'])
#     pl.gca().invert_yaxis()
#       
#     pl.figure(4)
#     special.hinton(statsob.stats['hessian'])
#     pl.gca().invert_yaxis()
#           
#     pl.figure(5)
#     special.hinton(statsob.stats['cov_asymp'])
#     pl.gca().invert_yaxis()    
 
#     eigenval = statsob.stats['eigenval']
#     eigenvec = statsob.stats['eigenvec']
#     
#     cov = statsob.stats['cov_small']
#     
# #     print cov
#     print cov[:3,:3]
#     
#     [eigval2,eigvec2] = np.linalg.eig(cov[:3,:3])
#     
# #     print np.shape(eigenval)
# #     print np.shape(eigenvec)
#     
# #     print eigenval
#     
#     
# #     print eigvec2
#     pivot = np.zeros((3,3))
# 
#     
#     
#     
# #     eigvec3 = eigval2*eigvec2
# #     print (eigenval[:,3]*eigenvec)[:3,:3]
#     
# #     
# #     soa =np.array( [ [0,0,3,2], [0,0,1,1],[0,0,9,9]]) 
# #     X,Y,U,V = zip(*soa)
# 
#     print eigval2
#     eigvec2[0,:] *= eigval2[0]
#     eigvec2[1,:] *= eigval2[1]
#     eigvec2[2,:] *= eigval2[2]
#     test = np.concatenate((pivot,eigvec2),axis=1)
#     X,Y,Z,U,V,W = zip(*test)
#     
#     print U 
#     print V 
#     print W
#     
# #     pl.figure()
# #     ax = plt.gca()
# #     ax.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=1)
# #     ax.set_xlim([-0.05,0.05])
# #     ax.set_xlim([-0.05,0.05])
# 
# 
# #     print X 
# #     print Y
# #     print Z
# #     print U
# #     print V 
# #     print W
# #     
# #     
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     ax.quiver(X,Y,Z,U,V,W, pivot='tail',color='red')
# #     ax.quiver(pivot,pivot,pivot,eigvec2[0,:],eigvec2[1,:],eigvec2[2,:], pivot='tail')
#     ax.set_xlim([-1,1])
#     ax.set_ylim([-1,1])
#     ax.set_zlim([-1,1])
# #  
# #      
#     pl.show()