import numpy as np
from scipy.optimize import shgo, dual_annealing, differential_evolution, basinhopping, minimize

def Hamiltonian(state,E,xsize,ysize):
    Sum = 0.0;
    systemsize = xsize*ysize;
    for x in range(0,xsize):
        for y in range(0,ysize):
            site = x*ysize+y            
            xneighbour = (x+1)*ysize+y
            yneighbour = x*ysize+(y+1)
            if (x<xsize-1 and y<ysize-1):
                Sum += 0.5*(np.cos(state[site+systemsize]-state[xneighbour+systemsize]) \
                        *(1.0+np.cos(state[site]-state[xneighbour])) \
                        +np.cos(state[site+systemsize]+state[xneighbour+systemsize]) \
                        *(1.0-np.cos(state[site]-state[xneighbour])) \
                        +np.cos(state[site+systemsize]-state[yneighbour+systemsize]) \
                        *(1.0+np.cos(state[site]-state[yneighbour])) \
                        +np.cos(state[site+systemsize]+state[yneighbour+systemsize]) \
                        *(1.0-np.cos(state[site]-state[yneighbour]))) \
                        -E[site]*(np.sin(state[site+systemsize]+state[xneighbour+systemsize]) \
                        *np.sin(0.5*(state[site]-state[xneighbour])) \
                        *np.sin(0.5*(state[site]+state[xneighbour])) \
                        -np.sin(state[site+systemsize]-state[xneighbour+systemsize]) \
                        *np.cos(0.5*(state[site]-state[xneighbour])) \
                        *np.cos(0.5*(state[site]+state[xneighbour]))) \
                        +E[site+systemsize]*(np.sin(state[site+systemsize]+state[yneighbour+systemsize]) \
                        *np.sin(0.5*(state[site]-state[yneighbour])) \
                        *np.cos(0.5*(state[site]+state[yneighbour])) \
                        +np.sin(state[site+systemsize]-state[yneighbour+systemsize]) \
                        *np.cos(0.5*(state[site]-state[yneighbour])) \
                        *np.sin(0.5*(state[site]+state[yneighbour])))
            elif (x==xsize-1 and y!=ysize-1):
                Sum += 0.5*(np.cos(state[site+systemsize]-state[yneighbour+systemsize]) \
                *(1.0+np.cos(state[site]-state[yneighbour])) \
                +np.cos(state[site+systemsize]+state[yneighbour+systemsize]) \
                *(1.0-np.cos(state[site]-state[yneighbour]))) \
                +E[site+systemsize]*(np.sin(state[site+systemsize]+state[yneighbour+systemsize]) \
                *np.sin(0.5*(state[site]-state[yneighbour])) \
                *np.cos(0.5*(state[site]+state[yneighbour])) \
                +np.sin(state[site+systemsize]-state[yneighbour+systemsize]) \
                *np.cos(0.5*(state[site]-state[yneighbour])) \
                *np.sin(0.5*(state[site]+state[yneighbour])))

            elif (y==ysize-1 and x!=xsize-1):
                Sum += 0.5*(np.cos(state[site+systemsize]-state[xneighbour+systemsize]) \
                *(1.0+np.cos(state[site]-state[xneighbour])) \
                +np.cos(state[site+systemsize]+state[xneighbour+systemsize]) \
                *(1.0-np.cos(state[site]-state[xneighbour]))) \
                -E[site]*(np.sin(state[site+systemsize]+state[xneighbour+systemsize]) \
                *np.sin(0.5*(state[site]-state[xneighbour])) \
                *np.sin(0.5*(state[site]+state[xneighbour])) \
                -np.sin(state[site+systemsize]-state[xneighbour+systemsize]) \
                *np.cos(0.5*(state[site]-state[xneighbour])) \
                *np.cos(0.5*(state[site]+state[xneighbour])))
    return Sum

def dHamiltonian(state,E,xsize,ysize):
    systemsize = xsize*ysize
    result = np.zeros(2*systemsize)
    for x in range(0,xsize):
        for y in range(0,ysize):            
            site = x*ysize+y
            xposneighbour = (x+1)*ysize+y
            xnegneighbour = (x-1)*ysize+y
            yposneighbour = x*ysize+(y+1)
            ynegneighbour = x*ysize+(y-1)
            if ((x>0 and x<xsize-1) and (y>0 and y<ysize-1)):
                result[site] = 0.5*(np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[ynegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[ynegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[yposneighbour])) \
                            *np.sin(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[xnegneighbour]*(-np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])) \
                            -np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site]))) \
                            -0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *(-np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            +np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            *np.cos(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[xnegneighbour]*(np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(-np.sin(0.5*(state[xnegneighbour]+state[site]))) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) \
                            +0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *(-np.sin(0.5*(state[ynegneighbour]+state[site]))) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))

                result[site+systemsize] = 0.5*(-np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[xposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[xposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[yposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[yposneighbour])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[xnegneighbour]-state[site])+1.0) \
                            -np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[xnegneighbour]-state[site])) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[ynegneighbour]-state[site])+1.0) \
                            -np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[ynegneighbour]-state[site]))) \
                            -E[site]*(np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            -np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +E[site+systemsize]*(np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            -E[xnegneighbour]*(np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) \
                            +E[ynegneighbour+systemsize]*(np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))                            

            elif (x==0 and (y!=0 and y!=ysize-1)):
                result[site] = 0.5*(np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[ynegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[ynegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[yposneighbour])) \
                            *np.sin(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *(-np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            +np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            *np.cos(0.5*(state[site]-state[yposneighbour]))) \
                            +0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *(-np.sin(0.5*(state[ynegneighbour]+state[site]))) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))

                result[site+systemsize] = 0.5*(-np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[xposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[xposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[yposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[yposneighbour])) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize])  \
                            *(np.cos(state[ynegneighbour]-state[site])+1.0) \
                            -np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[ynegneighbour]-state[site]))) \
                            -E[site]*(np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            -np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +E[site+systemsize]*(np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            +E[ynegneighbour+systemsize]*(np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))

            elif (x==xsize-1 and (y!=0 and y!=ysize-1)):
                result[site] = 0.5*(+np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[ynegneighbour]-state[site])) \
                            +np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[ynegneighbour]-state[site]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[yposneighbour])) \
                            *np.sin(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[xnegneighbour]*(-np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])) \
                            -np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site]))) \
                            -0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])))   \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *(-np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            +np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            *np.cos(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[xnegneighbour]*(np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.cos(0.5*(state[xnegneighbour]+state[site]))) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            -np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) \
                            +0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *(np.sin(0.5*(state[ynegneighbour]+state[site]))) \
                            -np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))

                result[site+systemsize] = 0.5*(-np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[yposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[yposneighbour])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[xnegneighbour]-state[site])+1.0) \
                            -np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[xnegneighbour]-state[site])) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[ynegneighbour]-state[site])+1.0) \
                            -np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[ynegneighbour]-state[site]))) \
                            +E[site+systemsize]*(np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            -E[xnegneighbour]*(np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) \
                            +E[ynegneighbour+systemsize]*(np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))

            elif (y==0 and (x!=0 and x!=xsize-1)):
                result[site] = 0.5*(np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[yposneighbour]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[yposneighbour]))\
                            *np.sin(0.5*(state[site]-state[yposneighbour])))\
                            -0.5*E[xnegneighbour]*(-np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])) \
                            -np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *(-np.sin(0.5*(state[site]+state[yposneighbour])))
                            +np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            *np.cos(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[xnegneighbour]*(np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(-np.sin(0.5*(state[xnegneighbour]+state[site]))) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) 

                result[site+systemsize] = 0.5*(-np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[xposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[xposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[yposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[yposneighbour])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[xnegneighbour]-state[site])+1.0) \
                            -np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[xnegneighbour]-state[site]))) \
                            -E[site]*(np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            -np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +E[site+systemsize]*(np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            -E[xnegneighbour]*(np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])))

            elif (y==ysize-1 and (x!=0 and x!=xsize-1)):
                result[site] = 0.5*(np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[ynegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[ynegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour]))) \
                            -0.5*E[xnegneighbour]*(-np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])) \
                            -np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site]))) \
                            -0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            -0.5*E[xnegneighbour]*(np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(-np.sin(0.5*(state[xnegneighbour]+state[site]))) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) \
                            +0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *(-np.sin(0.5*(state[ynegneighbour]+state[site]))) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site]))) 
                            
                result[site+systemsize] = 0.5*(-np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[xposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[xposneighbour])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[xnegneighbour]-state[site])+1.0) \
                            -np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[xnegneighbour]-state[site])) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[ynegneighbour]-state[site])+1.0) \
                            -np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[ynegneighbour]-state[site]))) \
                            -E[site]*(np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            -np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            -E[xnegneighbour]*(np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) \
                            +E[ynegneighbour+systemsize]*(np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site]))) 

            elif (x==0 and y==0):
                result[site] = 0.5*(np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[yposneighbour]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[yposneighbour])) \
                            *np.sin(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *(-np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            +np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])))

                result[site+systemsize] = 0.5*(-np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[xposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[xposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) 
                            *(np.cos(state[site]-state[yposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[yposneighbour]))) \
                            -E[site]*(np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            -np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +E[site+systemsize]*(np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.sin(0.5*(state[site]+state[yposneighbour])))

            elif (x==0 and y==ysize-1):
                result[site] = 0.5*(np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[ynegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[xposneighbour])) \
                            +np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[ynegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour]))) \
                            -0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site]))) \
                            -0.5*E[site]*(np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            +np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *(-np.sin(0.5*(state[ynegneighbour]+state[site]))) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))

                result[site+systemsize] = 0.5*(-np.sin(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[xposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[xposneighbour])) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[ynegneighbour]-state[site])+1.0) \
                            -np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[ynegneighbour]-state[site]))) \
                            -E[site]*(np.cos(state[site+systemsize]+state[xposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[xposneighbour])) \
                            *np.sin(0.5*(state[site]-state[xposneighbour])) \
                            -np.cos(state[site+systemsize]-state[xposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[xposneighbour])) \
                            *np.cos(0.5*(state[site]-state[xposneighbour]))) \
                            +E[ynegneighbour+systemsize]*(np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))

            elif (x==xsize-1 and y==0):
                result[site] = 0.5*(+np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(-np.sin(state[site]-state[yposneighbour])) \
                            +np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(np.sin(state[site]-state[yposneighbour]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            -np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]+state[yposneighbour])) \
                            *np.sin(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[xnegneighbour]*(-np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])) \
                            -np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site]))) \
                            +0.5*E[site+systemsize]*(np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *(-np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            +np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            *np.cos(0.5*(state[site]-state[yposneighbour]))) \
                            -0.5*E[xnegneighbour]*(np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(-np.sin(0.5*(state[xnegneighbour]+state[site]))) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])))

                result[site+systemsize] = 0.5*(-np.sin(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *(np.cos(state[site]-state[yposneighbour])+1.0) \
                            -np.sin(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *(1.0-np.cos(state[site]-state[yposneighbour])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[xnegneighbour]-state[site])+1.0) \
                            -np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[xnegneighbour]-state[site]))) \
                            +E[site+systemsize]*(np.cos(state[site+systemsize]+state[yposneighbour+systemsize]) \
                            *np.sin(0.5*(state[site]-state[yposneighbour])) \
                            *np.cos(0.5*(state[site]+state[yposneighbour])) \
                            +np.cos(state[site+systemsize]-state[yposneighbour+systemsize]) \
                            *np.cos(0.5*(state[site]-state[yposneighbour])) \
                            *np.sin(0.5*(state[site]+state[yposneighbour]))) \
                            -E[xnegneighbour]*(np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])))

            elif (x==xsize-1 and y==ysize-1):
                result[site] = 0.5*(+np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.sin(state[ynegneighbour]-state[site])) \
                            +np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[xnegneighbour]-state[site])) \
                            +np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(-np.sin(state[ynegneighbour]-state[site]))) \
                            -0.5*E[xnegneighbour]*(-np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site])) \
                            -np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site]))) \
                            -0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site]))) \
                            -0.5*E[xnegneighbour]*(np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(-np.sin(0.5*(state[xnegneighbour]+state[site]))) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) \
                            +0.5*E[ynegneighbour+systemsize]*(np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *(-np.sin(0.5*(state[ynegneighbour]+state[site]))) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))

                result[site+systemsize] = 0.5*(+np.sin(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[xnegneighbour]-state[site])+1.0) \
                            -np.sin(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[xnegneighbour]-state[site])) \
                            +np.sin(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *(np.cos(state[ynegneighbour]-state[site])+1.0) \
                            -np.sin(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *(1.0-np.cos(state[ynegneighbour]-state[site]))) \
                            -E[xnegneighbour]*(np.cos(state[xnegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[xnegneighbour]+state[site])) \
                            *np.sin(0.5*(state[xnegneighbour]-state[site])) \
                            +np.cos(state[xnegneighbour+systemsize]-state[site+systemsize]) \
                            *np.cos(0.5*(state[xnegneighbour]+state[site])) \
                            *np.cos(0.5*(state[xnegneighbour]-state[site]))) \
                            +E[ynegneighbour+systemsize]*(np.cos(state[ynegneighbour+systemsize]+state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]-state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]+state[site])) \
                            -np.cos(state[ynegneighbour+systemsize]-state[site+systemsize]) \
                            *np.sin(0.5*(state[ynegneighbour]+state[site])) \
                            *np.cos(0.5*(state[ynegneighbour]-state[site])))
    
    return result

def minFunc(state):
    return Hamiltonian(state,E,xsize,ysize)

def minFunc_jac(state):
    return Hamiltonian(state,E,xsize,ysize), dHamiltonian(state,E,xsize,ysize)

def minFunc_jacOnly(state):
    return dHamiltonian(state,E,xsize,ysize)

xsize = 4
ysize = 4   

boundaries = []
neelstate = np.zeros(2*xsize*ysize)
for x in range(0,xsize):
    for y in range(0,ysize):
        site = x*ysize+y
        neelstate[xsize*ysize+site] = np.pi*(x+y)%2
        
for i in range(0,xsize*ysize):
    boundaries.append((0,2*np.pi))    
    
for i in range(0,xsize*ysize):
    boundaries.append((0,np.pi))

phimax = 2*np.pi*np.ones(xsize*ysize)
thetamax = np.pi*np.ones(xsize*ysize)
xmax = np.concatenate((phimax,thetamax))

class MyBounds(object):
    def __init__(self, xmax=xmax, xmin=np.zeros(2*xsize*ysize)):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin

minimizer_kwargs = {"method":"L-BFGS-B", "jac":True}
mybounds = MyBounds()

tempstate=np.array([-2.1299629 , -2.71203839,  2.45804659,  2.01630002, -1.50166928,
       -1.14373804,  1.09110667,  1.46388693, -1.06786288, -0.57864832,
        0.54537084,  1.03434621, -0.74743426, -0.29577916,  0.27642135,
        0.72471315, -0.31638314,  2.97025166, -0.18293004,  2.80441955,
        2.86371589, -0.09028476,  3.03434224, -0.29805912, -0.34338475,
        2.93894156, -0.21038126,  2.78928848,  2.68025153, -0.36351256,
        2.77354467, -0.46914034])

dincrement = 0.01
state = neelstate
for e in np.linspace(0,1,(1/dincrement)+1):
    E = e*np.ones(2*xsize*ysize)
    
    #result = basinhopping(minFunc,neelstate)
    #result = basinhopping(minFunc_jac, tempstate, minimizer_kwargs=minimizer_kwargs,
    #               niter=1000,accept_test=mybounds)
    #result = minimize(Hamiltonian,neelstate,args=(E,xsize,ysize), method='COBYLA', tol=None, callback=None, options={'rhobeg': 1.0, 'maxiter': 1000, 'catol': 0.0002})
    print(result)    
