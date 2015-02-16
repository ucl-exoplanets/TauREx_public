# TauREx (Tau Retrieval for Exoplanets)
## A next generation retrieval code for exoplanetary atmospheres

Spectroscopy of exoplanetary atmospheres has become a well established method for the characterisation of extrasolar planets. We here present a novel inverse retrieval code for exoplanetary atmospheres. TauRex (Tau Retrieval for Exoplanets) is a line-by-line radiative transfer fully Bayesian retrieval framework. 

TauRex includes the following features: 

1. The optimised use of molecular line-lists from the Exomol project
2. An unbiased atmospheric composition prior selection, through custom built pattern recognition software
3. The use of two independent algorithms to fully sample the Bayesian likelihood space: nested sampling as well as a more classical Markov Chain Monte Carlo approach
4. Iterative Bayesian parameter and model selection using the full Bayesian Evidence as well as the Savage-Dickey Ratio for nested models
5. The ability to fully map very large parameter spaces through optimal code parallelisation and scalability to cluster computing. In this publication we outline the TauRex framework and demonstrate, using a theoretical hot-Jupiter transmission spectrum, the parameter retrieval and model selection. 

## Installation:

run `./setup.py` before running `exonest.py` 
as this will *probably* take care of the c++ compilations required 

otherwise, read the code yourself I cannot be bothered writing a readme file right now. 

Good bye and good luck
