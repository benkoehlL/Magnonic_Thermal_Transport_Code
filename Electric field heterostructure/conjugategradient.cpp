#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include "parameter.h"

using namespace std;

double hamiltonian(const gsl_vector *v,void *param){
    double *E = (double *) param;
    double* state = new double[2*systemsize];//phi in first half, theta in second
    for(int x=0;x<xsize;x++){
        for(int y=0;y<ysize;y++){
            state[x*ysize+y] = gsl_vector_get(v,x*ysize+y);
            state[x*ysize+y+systemsize] = gsl_vector_get(v,x*ysize+y+systemsize);
        }
    }
    double sum = 0.0;
    int site,xneighbour,yneighbour;
    for(int x=0;x<xsize;x++){
        for(int y=0;y<ysize;y++){
            site = x*ysize+y;
            xneighbour = (x+1)*ysize+y;
            yneighbour = x*ysize+(y+1);
            if(x<xsize-1 && y<ysize-1){
                sum += 0.5*(cos(state[site+systemsize]-state[xneighbour+systemsize])
                                *(1.0+cos(state[site]-state[xneighbour]))
                            +cos(state[site+systemsize]+state[xneighbour+systemsize])
                                *(1.0-cos(state[site]-state[xneighbour]))
                            +cos(state[site+systemsize]-state[yneighbour+systemsize])
                                *(1.0+cos(state[site]-state[yneighbour]))
                            +cos(state[site+systemsize]+state[yneighbour+systemsize])
                                *(1.0-cos(state[site]-state[yneighbour])))
                       -E[site]*(sin(state[site+systemsize]+state[xneighbour+systemsize])
                                *sin(0.5*(state[site]-state[xneighbour]))
                                *sin(0.5*(state[site]+state[xneighbour]))
                            -sin(state[site+systemsize]-state[xneighbour+systemsize])
                                *cos(0.5*(state[site]-state[xneighbour]))
                                *cos(0.5*(state[site]+state[xneighbour])))
                       +E[site+systemsize]*(sin(state[site+systemsize]+state[yneighbour+systemsize])
                                *sin(0.5*(state[site]-state[yneighbour]))
                                *cos(0.5*(state[site]+state[yneighbour]))
                            +sin(state[site+systemsize]-state[yneighbour+systemsize])
                                *cos(0.5*(state[site]-state[yneighbour]))
                                *sin(0.5*(state[site]+state[yneighbour])));
            }
            else if(x==xsize-1 && y!=ysize-1){
                sum += 0.5*(cos(state[site+systemsize]-state[yneighbour+systemsize])
                                *(1.0+cos(state[site]-state[yneighbour]))
                            +cos(state[site+systemsize]+state[yneighbour+systemsize])
                                *(1.0-cos(state[site]-state[yneighbour])))
                       +E[site+systemsize]*(sin(state[site+systemsize]+state[yneighbour+systemsize])
                                *sin(0.5*(state[site]-state[yneighbour]))
                                *cos(0.5*(state[site]+state[yneighbour]))
                            +sin(state[site+systemsize]-state[yneighbour+systemsize])
                                *cos(0.5*(state[site]-state[yneighbour]))
                                *sin(0.5*(state[site]+state[yneighbour])));
            }
            else if(y==ysize-1 && x!=xsize-1){
                sum += 0.5*(cos(state[site+systemsize]-state[xneighbour+systemsize])
                                *(1.0+cos(state[site]-state[xneighbour]))
                            +cos(state[site+systemsize]+state[xneighbour+systemsize])
                                *(1.0-cos(state[site]-state[xneighbour])))
                       -E[site]*(sin(state[site+systemsize]+state[xneighbour+systemsize])
                                *sin(0.5*(state[site]-state[xneighbour]))
                                *sin(0.5*(state[site]+state[xneighbour]))
                            -sin(state[site+systemsize]-state[xneighbour+systemsize])
                                *cos(0.5*(state[site]-state[xneighbour]))
                                *cos(0.5*(state[site]+state[xneighbour])))
                       // for periodic boundaries in y-direction:
                       /*+0.5*(cos(state[site+systemsize]-state[(yneighbour%(ysize))+systemsize])
                                *(1.0+cos(state[site]-state[yneighbour%ysize]))
                            +cos(state[site+systemsize]+state[(yneighbour%ysize)+systemsize])
                                *(1.0-cos(state[site]-state[(yneighbour%ysize)])))
                       +E[site+systemsize]*(sin(state[site+systemsize]+state[(yneighbour%ysize)+systemsize])
                                *sin(0.5*(state[site]-state[(yneighbour%ysize)]))
                                *cos(0.5*(state[site]+state[(yneighbour%ysize)]))
                            +sin(state[site+systemsize]-state[(yneighbour%ysize)+systemsize])
                                *cos(0.5*(state[site]-state[(yneighbour%ysize)]))
                                *sin(0.5*(state[site]+state[(yneighbour%ysize)])))*/;
            }
        }
    }
    delete [] state;
    return sum;
}

void dhamiltonian(const gsl_vector *v,void *param, gsl_vector *df){
    double *E = (double *) param;
    double* state = new double[2*systemsize];//phi in first half, theta in second half
    int site,xposneighbour,yposneighbour,xnegneighbour,ynegneighbour;
    for(int x=0;x<xsize;x++){
        for(int y=0;y<ysize;y++){
            state[x*ysize+y] = gsl_vector_get(v,x*ysize+y);
            state[x*ysize+y+systemsize] = gsl_vector_get(v,x*ysize+y+systemsize);
        }
    }
    for(int x=0;x<xsize;x++){
        for(int y=0;y<ysize;y++){
            site = x*ysize+y;
            xposneighbour = (x+1)*ysize+y;
            xnegneighbour = (x-1)*ysize+y;
            yposneighbour = x*ysize+(y+1);
            ynegneighbour = x*ysize+(y-1);
            // for periodic boundaries in y-direction:
            if(y==ysize-1){
                yposneighbour = x*ysize;
            }
            if(y==0){
                ynegneighbour = x*ysize+ysize-1;
            }
            if((x>0 && x<xsize-1) && (y>0 && y<ysize-1)){
                //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(-sin(state[site]-state[xposneighbour]))
                            +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[xnegneighbour]-state[site]))
                            +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(-sin(state[site]-state[yposneighbour]))
                            +cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[ynegneighbour]-state[site]))
                            +cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(sin(state[site]-state[xposneighbour]))
                            +cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[xnegneighbour]-state[site]))
                            +cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(sin(state[site]-state[yposneighbour]))
                            +cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[ynegneighbour]-state[site])))
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]-state[yposneighbour]))
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                 -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[yposneighbour]))
                                    *sin(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[xnegneighbour]*(-sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[xnegneighbour]+state[site]))
                                    *cos(0.5*(state[xnegneighbour]-state[site]))
                                 -sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site])))
                        -0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]-state[site]))
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                 -sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]+state[site]))
                                    *sin(0.5*(state[ynegneighbour]-state[site])))
                        //beta = phi_i+ phi_j derivatives:
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]-state[yposneighbour]))
                                    *(-sin(0.5*(state[site]+state[yposneighbour])))
                                 +sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                    *cos(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[xnegneighbour]*(sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site]))
                                 +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *(-sin(0.5*(state[xnegneighbour]+state[site])))
                                    *cos(0.5*(state[xnegneighbour]-state[site])))
                        +0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]-state[site]))
                                    *(-sin(0.5*(state[ynegneighbour]+state[site])))
                                 +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                    *cos(0.5*(state[ynegneighbour]-state[site])))
                        );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(-sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(cos(state[site]-state[xposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[xposneighbour]))
                             -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(cos(state[site]-state[yposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[yposneighbour]))
                             +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[xnegneighbour]-state[site])+1.0)
                             -sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[xnegneighbour]-state[site]))
                             +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[ynegneighbour]-state[site])+1.0)
                             -sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[ynegneighbour]-state[site])))
                        -E[site]*(cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *sin(0.5*(state[site]+state[xposneighbour]))
                                *sin(0.5*(state[site]-state[xposneighbour]))
                             -cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *cos(0.5*(state[site]+state[xposneighbour]))
                                *cos(0.5*(state[site]-state[xposneighbour])))
                        +E[site+systemsize]*(cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *sin(0.5*(state[site]-state[yposneighbour]))
                                *cos(0.5*(state[site]+state[yposneighbour]))
                             +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *cos(0.5*(state[site]-state[yposneighbour]))
                                *sin(0.5*(state[site]+state[yposneighbour])))
                        -E[xnegneighbour]*(cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[xnegneighbour]+state[site]))
                                *sin(0.5*(state[xnegneighbour]-state[site]))
                             +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *cos(0.5*(state[xnegneighbour]+state[site]))
                                *cos(0.5*(state[xnegneighbour]-state[site])))
                        +E[ynegneighbour+systemsize]*(cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]-state[site]))
                                *cos(0.5*(state[ynegneighbour]+state[site]))
                             -cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]+state[site]))
                                *cos(0.5*(state[ynegneighbour]-state[site])))
                );
            }
            else if(x==0 && (y!=0 && y!=ysize-1)){
                //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(-sin(state[site]-state[xposneighbour]))
                            +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(-sin(state[site]-state[yposneighbour]))
                            +cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[ynegneighbour]-state[site]))
                            +cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(sin(state[site]-state[xposneighbour]))
                            +cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(sin(state[site]-state[yposneighbour]))
                            +cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[ynegneighbour]-state[site])))
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]-state[yposneighbour]))
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                 -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[yposneighbour]))
                                    *sin(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]-state[site]))
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                 -sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]+state[site]))
                                    *sin(0.5*(state[ynegneighbour]-state[site])))
                        //beta = phi_i+ phi_j derivatives:
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]-state[yposneighbour]))
                                    *(-sin(0.5*(state[site]+state[yposneighbour])))
                                 +sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                    *cos(0.5*(state[site]-state[yposneighbour])))
                        +0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]-state[site]))
                                    *(-sin(0.5*(state[ynegneighbour]+state[site])))
                                 +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                    *cos(0.5*(state[ynegneighbour]-state[site])))
                        );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(-sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(cos(state[site]-state[xposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[xposneighbour]))
                             -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(cos(state[site]-state[yposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[yposneighbour]))
                             +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[ynegneighbour]-state[site])+1.0)
                             -sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[ynegneighbour]-state[site])))
                        -E[site]*(cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *sin(0.5*(state[site]+state[xposneighbour]))
                                *sin(0.5*(state[site]-state[xposneighbour]))
                             -cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *cos(0.5*(state[site]+state[xposneighbour]))
                                *cos(0.5*(state[site]-state[xposneighbour])))
                        +E[site+systemsize]*(cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *sin(0.5*(state[site]-state[yposneighbour]))
                                *cos(0.5*(state[site]+state[yposneighbour]))
                             +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *cos(0.5*(state[site]-state[yposneighbour]))
                                *sin(0.5*(state[site]+state[yposneighbour])))
                        +E[ynegneighbour+systemsize]*(cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]-state[site]))
                                *cos(0.5*(state[ynegneighbour]+state[site]))
                             -cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]+state[site]))
                                *cos(0.5*(state[ynegneighbour]-state[site])))
                );
            }
            else if(x==xsize-1 && (y!=0 && y!=ysize-1)){
                //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(+cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[xnegneighbour]-state[site]))
                            +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(-sin(state[site]-state[yposneighbour]))
                            +cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[ynegneighbour]-state[site]))
                            +cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[xnegneighbour]-state[site]))
                            +cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(sin(state[site]-state[yposneighbour]))
                            +cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[ynegneighbour]-state[site])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]-state[yposneighbour]))
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                 -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[yposneighbour]))
                                    *sin(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[xnegneighbour]*(-sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[xnegneighbour]+state[site]))
                                    *cos(0.5*(state[xnegneighbour]-state[site]))
                                 -sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site])))
                        -0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]-state[site]))
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                 -sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]+state[site]))
                                    *sin(0.5*(state[ynegneighbour]-state[site])))
                        //beta = phi_i+ phi_j derivatives:
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]-state[yposneighbour]))
                                    *(-sin(0.5*(state[site]+state[yposneighbour])))
                                 +sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                    *cos(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[xnegneighbour]*(sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *(-cos(0.5*(state[xnegneighbour]+state[site])))
                                    *sin(0.5*(state[xnegneighbour]-state[site]))
                                 -sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *sin(0.5*(state[xnegneighbour]+state[site]))
                                    *cos(0.5*(state[xnegneighbour]-state[site])))
                        +0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]-state[site]))
                                    *(sin(0.5*(state[ynegneighbour]+state[site])))
                                 -sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                    *cos(0.5*(state[ynegneighbour]-state[site])))
                        );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(-sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(cos(state[site]-state[yposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[yposneighbour]))
                             +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[xnegneighbour]-state[site])+1.0)
                             -sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[xnegneighbour]-state[site]))
                             +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[ynegneighbour]-state[site])+1.0)
                             -sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[ynegneighbour]-state[site])))
                        +E[site+systemsize]*(cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *sin(0.5*(state[site]-state[yposneighbour]))
                                *cos(0.5*(state[site]+state[yposneighbour]))
                             +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *cos(0.5*(state[site]-state[yposneighbour]))
                                *sin(0.5*(state[site]+state[yposneighbour])))
                        -E[xnegneighbour]*(cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[xnegneighbour]+state[site]))
                                *sin(0.5*(state[xnegneighbour]-state[site]))
                             +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *cos(0.5*(state[xnegneighbour]+state[site]))
                                *cos(0.5*(state[xnegneighbour]-state[site])))
                        +E[ynegneighbour+systemsize]*(cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]-state[site]))
                                *cos(0.5*(state[ynegneighbour]+state[site]))
                             -cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]+state[site]))
                                *cos(0.5*(state[ynegneighbour]-state[site])))
                );
            }
            else if(y==0 && (x!=0 && x!=xsize-1)){
                //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(-sin(state[site]-state[xposneighbour]))
                            +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[xnegneighbour]-state[site]))
                            +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(-sin(state[site]-state[yposneighbour]))
                            +cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(sin(state[site]-state[xposneighbour]))
                            +cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[xnegneighbour]-state[site]))
                            +cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(sin(state[site]-state[yposneighbour])))
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]-state[yposneighbour]))
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                 -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[yposneighbour]))
                                    *sin(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[xnegneighbour]*(-sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[xnegneighbour]+state[site]))
                                    *cos(0.5*(state[xnegneighbour]-state[site]))
                                 -sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site])))
                        //beta = phi_i+ phi_j derivatives:
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]-state[yposneighbour]))
                                    *(-sin(0.5*(state[site]+state[yposneighbour])))
                                 +sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                    *cos(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[xnegneighbour]*(sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site]))
                                 +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *(-sin(0.5*(state[xnegneighbour]+state[site])))
                                    *cos(0.5*(state[xnegneighbour]-state[site])))
                        );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(-sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(cos(state[site]-state[xposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[xposneighbour]))
                             -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(cos(state[site]-state[yposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[yposneighbour]))
                             +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[xnegneighbour]-state[site])+1.0)
                             -sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[xnegneighbour]-state[site])))
                        -E[site]*(cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *sin(0.5*(state[site]+state[xposneighbour]))
                                *sin(0.5*(state[site]-state[xposneighbour]))
                             -cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *cos(0.5*(state[site]+state[xposneighbour]))
                                *cos(0.5*(state[site]-state[xposneighbour])))
                        +E[site+systemsize]*(cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *sin(0.5*(state[site]-state[yposneighbour]))
                                *cos(0.5*(state[site]+state[yposneighbour]))
                             +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *cos(0.5*(state[site]-state[yposneighbour]))
                                *sin(0.5*(state[site]+state[yposneighbour])))
                        -E[xnegneighbour]*(cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[xnegneighbour]+state[site]))
                                *sin(0.5*(state[xnegneighbour]-state[site]))
                             +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *cos(0.5*(state[xnegneighbour]+state[site]))
                                *cos(0.5*(state[xnegneighbour]-state[site])))
                );
            }
            else if(y==ysize-1 && (x!=0 && x!=xsize-1)){
                //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(-sin(state[site]-state[xposneighbour]))
                            +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[xnegneighbour]-state[site]))
                            +cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[ynegneighbour]-state[site]))
                            +cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(sin(state[site]-state[xposneighbour]))
                            +cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[xnegneighbour]-state[site]))
                            +cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[ynegneighbour]-state[site])))
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour])))
                        -0.5*E[xnegneighbour]*(-sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[xnegneighbour]+state[site]))
                                    *cos(0.5*(state[xnegneighbour]-state[site]))
                                 -sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site])))
                        -0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]-state[site]))
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                 -sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]+state[site]))
                                    *sin(0.5*(state[ynegneighbour]-state[site])))
                        //beta = phi_i+ phi_j derivatives:
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour])))
                        -0.5*E[xnegneighbour]*(sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site]))
                                 +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *(-sin(0.5*(state[xnegneighbour]+state[site])))
                                    *cos(0.5*(state[xnegneighbour]-state[site])))
                        +0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]-state[site]))
                                    *(-sin(0.5*(state[ynegneighbour]+state[site])))
                                 +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                    *cos(0.5*(state[ynegneighbour]-state[site])))
                        );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(-sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(cos(state[site]-state[xposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[xposneighbour]))
                             +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[xnegneighbour]-state[site])+1.0)
                             -sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[xnegneighbour]-state[site]))
                             +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[ynegneighbour]-state[site])+1.0)
                             -sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[ynegneighbour]-state[site])))
                        -E[site]*(cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *sin(0.5*(state[site]+state[xposneighbour]))
                                *sin(0.5*(state[site]-state[xposneighbour]))
                             -cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *cos(0.5*(state[site]+state[xposneighbour]))
                                *cos(0.5*(state[site]-state[xposneighbour])))
                        -E[xnegneighbour]*(cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[xnegneighbour]+state[site]))
                                *sin(0.5*(state[xnegneighbour]-state[site]))
                             +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *cos(0.5*(state[xnegneighbour]+state[site]))
                                *cos(0.5*(state[xnegneighbour]-state[site])))
                        +E[ynegneighbour+systemsize]*(cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]-state[site]))
                                *cos(0.5*(state[ynegneighbour]+state[site]))
                             -cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]+state[site]))
                                *cos(0.5*(state[ynegneighbour]-state[site])))
                );
            }
            else if(x==0 && y==0){
                //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(-sin(state[site]-state[xposneighbour]))
                            +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(-sin(state[site]-state[yposneighbour]))
                            +cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(sin(state[site]-state[xposneighbour]))
                            +cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(sin(state[site]-state[yposneighbour])))
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]-state[yposneighbour]))
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                 -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[yposneighbour]))
                                    *sin(0.5*(state[site]-state[yposneighbour])))
                        //beta = phi_i+ phi_j derivatives:
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]-state[yposneighbour]))
                                    *(-sin(0.5*(state[site]+state[yposneighbour])))
                                 +sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                    *cos(0.5*(state[site]-state[yposneighbour])))
                );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(-sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(cos(state[site]-state[xposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[xposneighbour]))
                             -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(cos(state[site]-state[yposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[yposneighbour])))
                        -E[site]*(cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *sin(0.5*(state[site]+state[xposneighbour]))
                                *sin(0.5*(state[site]-state[xposneighbour]))
                             -cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *cos(0.5*(state[site]+state[xposneighbour]))
                                *cos(0.5*(state[site]-state[xposneighbour])))
                        +E[site+systemsize]*(cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *sin(0.5*(state[site]-state[yposneighbour]))
                                *cos(0.5*(state[site]+state[yposneighbour]))
                             +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *cos(0.5*(state[site]-state[yposneighbour]))
                                *sin(0.5*(state[site]+state[yposneighbour])))
                );
            }
            else if(x==0 && y==ysize-1){
                //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(-sin(state[site]-state[xposneighbour]))
                            +cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[ynegneighbour]-state[site]))
                            +cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(sin(state[site]-state[xposneighbour]))
                            +cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[ynegneighbour]-state[site])))
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour])))
                        -0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]-state[site]))
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                 -sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]+state[site]))
                                    *sin(0.5*(state[ynegneighbour]-state[site])))
                        //beta = phi_i+ phi_j derivatives:
                        -0.5*E[site]*(sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[xposneighbour]))
                                    *sin(0.5*(state[site]-state[xposneighbour]))
                                 +sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[xposneighbour]))
                                    *cos(0.5*(state[site]-state[xposneighbour])))
                        +0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]-state[site]))
                                    *(-sin(0.5*(state[ynegneighbour]+state[site])))
                                 +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                    *cos(0.5*(state[ynegneighbour]-state[site])))
                        );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(-sin(state[site+systemsize]-state[xposneighbour+systemsize])
                                *(cos(state[site]-state[xposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[xposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[xposneighbour]))
                             +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[ynegneighbour]-state[site])+1.0)
                             -sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[ynegneighbour]-state[site])))
                        -E[site]*(cos(state[site+systemsize]+state[xposneighbour+systemsize])
                                *sin(0.5*(state[site]+state[xposneighbour]))
                                *sin(0.5*(state[site]-state[xposneighbour]))
                             -cos(state[site+systemsize]-state[xposneighbour+systemsize])
                                *cos(0.5*(state[site]+state[xposneighbour]))
                                *cos(0.5*(state[site]-state[xposneighbour])))
                        +E[ynegneighbour+systemsize]*(cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]-state[site]))
                                *cos(0.5*(state[ynegneighbour]+state[site]))
                             -cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]+state[site]))
                                *cos(0.5*(state[ynegneighbour]-state[site])))
                );
            }
            else if(x==xsize-1 && y==0){
            //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(+cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[xnegneighbour]-state[site]))
                            +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(-sin(state[site]-state[yposneighbour]))
                            +cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[xnegneighbour]-state[site]))
                            +cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(sin(state[site]-state[yposneighbour])))
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]-state[yposneighbour]))
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                 -sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]+state[yposneighbour]))
                                    *sin(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[xnegneighbour]*(-sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[xnegneighbour]+state[site]))
                                    *cos(0.5*(state[xnegneighbour]-state[site]))
                                 -sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site])))
                        //beta = phi_i+ phi_j derivatives:
                        +0.5*E[site+systemsize]*(sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                    *sin(0.5*(state[site]-state[yposneighbour]))
                                    *(-sin(0.5*(state[site]+state[yposneighbour])))
                                 +sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                    *cos(0.5*(state[site]+state[yposneighbour]))
                                    *cos(0.5*(state[site]-state[yposneighbour])))
                        -0.5*E[xnegneighbour]*(sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site]))
                                 +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *(-sin(0.5*(state[xnegneighbour]+state[site])))
                                    *cos(0.5*(state[xnegneighbour]-state[site])))
                );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(-sin(state[site+systemsize]-state[yposneighbour+systemsize])
                                *(cos(state[site]-state[yposneighbour])+1.0)
                             -sin(state[site+systemsize]+state[yposneighbour+systemsize])
                                *(1.0-cos(state[site]-state[yposneighbour]))
                             +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[xnegneighbour]-state[site])+1.0)
                             -sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[xnegneighbour]-state[site])))
                        +E[site+systemsize]*(cos(state[site+systemsize]+state[yposneighbour+systemsize])
                                *sin(0.5*(state[site]-state[yposneighbour]))
                                *cos(0.5*(state[site]+state[yposneighbour]))
                             +cos(state[site+systemsize]-state[yposneighbour+systemsize])
                                *cos(0.5*(state[site]-state[yposneighbour]))
                                *sin(0.5*(state[site]+state[yposneighbour])))
                        -E[xnegneighbour]*(cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[xnegneighbour]+state[site]))
                                *sin(0.5*(state[xnegneighbour]-state[site]))
                             +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *cos(0.5*(state[xnegneighbour]+state[site]))
                                *cos(0.5*(state[xnegneighbour]-state[site])))
                );
            }
            else if(x==xsize-1 && y==ysize-1){
                //phi gradient
                gsl_vector_set(df, site, 
                        //alpha=phi_i-phi_j derivatives:
                        0.5*(+cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[xnegneighbour]-state[site]))
                            +cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(sin(state[ynegneighbour]-state[site]))
                            +cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[xnegneighbour]-state[site]))
                            +cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(-sin(state[ynegneighbour]-state[site])))
                        -0.5*E[xnegneighbour]*(-sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[xnegneighbour]+state[site]))
                                    *cos(0.5*(state[xnegneighbour]-state[site]))
                                 -sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site])))
                        -0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]-state[site]))
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                 -sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]+state[site]))
                                    *sin(0.5*(state[ynegneighbour]-state[site])))
                        //beta = phi_i+ phi_j derivatives:
                        -0.5*E[xnegneighbour]*(sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                    *cos(0.5*(state[xnegneighbour]+state[site]))
                                    *sin(0.5*(state[xnegneighbour]-state[site]))
                                 +sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                    *(-sin(0.5*(state[xnegneighbour]+state[site])))
                                    *cos(0.5*(state[xnegneighbour]-state[site])))
                        +0.5*E[ynegneighbour+systemsize]*(sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                    *sin(0.5*(state[ynegneighbour]-state[site]))
                                    *(-sin(0.5*(state[ynegneighbour]+state[site])))
                                 +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                    *cos(0.5*(state[ynegneighbour]+state[site]))
                                    *cos(0.5*(state[ynegneighbour]-state[site])))
                );
                //theta gradient
                gsl_vector_set(df, site+systemsize,
                        0.5*(+sin(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[xnegneighbour]-state[site])+1.0)
                             -sin(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[xnegneighbour]-state[site]))
                             +sin(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *(cos(state[ynegneighbour]-state[site])+1.0)
                             -sin(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *(1.0-cos(state[ynegneighbour]-state[site])))
                        -E[xnegneighbour]*(cos(state[xnegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[xnegneighbour]+state[site]))
                                *sin(0.5*(state[xnegneighbour]-state[site]))
                             +cos(state[xnegneighbour+systemsize]-state[site+systemsize])
                                *cos(0.5*(state[xnegneighbour]+state[site]))
                                *cos(0.5*(state[xnegneighbour]-state[site])))
                        +E[ynegneighbour+systemsize]*(cos(state[ynegneighbour+systemsize]+state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]-state[site]))
                                *cos(0.5*(state[ynegneighbour]+state[site]))
                             -cos(state[ynegneighbour+systemsize]-state[site+systemsize])
                                *sin(0.5*(state[ynegneighbour]+state[site]))
                                *cos(0.5*(state[ynegneighbour]-state[site])))
                );
            }
        }
    }
    delete [] state;
}

void pathdirection(const gsl_vector *v, void *param, double *f, gsl_vector *df){
    *f = hamiltonian(v, param);
    dhamiltonian(v,param, df);
}

int main(){
    double E;
    double help;
    size_t iter;
    int status;
    double size;
    
    double* parameter = new double[2*systemsize];
    for(int ix=0;ix<xsize;ix++){
        for(int y=0;y<ysize;y++){
                    parameter[ix*ysize+y] = 0.0;
                    parameter[ix*ysize+y+systemsize] = 0.0;
        }
    }
    
    gsl_vector *state;//phi in first half, theta in second half

    gsl_vector *steps;    
    steps = gsl_vector_alloc(2*systemsize);
    gsl_vector_set_all(steps,1e5);
    
    state = gsl_vector_alloc (2*systemsize);
    for(int ix=0;ix<xsize;ix++){
        for(int y=0;y<ysize;y++){
            gsl_vector_set (state,ix*ysize+y,0.0*M_PI);
            gsl_vector_set (state,ix*ysize+y+systemsize,M_PI*((y+ix)%2));//
            //initialise Neel state
        }
    }   
    ostringstream fm;
    fm << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << "/textureD" << 0 << "cJxsize_simplex.dat";
    ofstream writezero(fm.str().c_str(),ios_base::binary);
    for(int ix=0; ix<xsize;ix++){
        for(int y=0;y<ysize;y++){
            help = gsl_vector_get(state,ix*ysize+y);
            writezero.write((char*)& help, sizeof(double));
        }
    }
    for(int ix=0; ix<xsize;ix++){
        for(int y=0;y<ysize;y++){
            help = gsl_vector_get(state,ix*ysize+y+systemsize);
            writezero.write((char*)& help, sizeof(double));
        }
    }
    writezero.close();
    
    gsl_multimin_function_fdf heterostruct;
    //gsl_multimin_function simplex_min;
    for(int count=1;count<=(int)(Dmax/Dincrement)+1;count++){
        E = count*Dincrement;
        
        const gsl_multimin_fdfminimizer_type *T;
        gsl_multimin_fdfminimizer *s;
        //const gsl_multimin_fminimizer_type *T  = 
        //    gsl_multimin_fminimizer_nmsimplex2;
        //gsl_multimin_fminimizer *s = NULL;
        iter = 0;
        for(int ix=0;ix<xsize;ix++){
            for(int y=0;y<ysize;y++){
                if(ix>=freefieldsize && ix<(freefieldsize+fieldsize)){
                    parameter[ix*ysize+y] = E;
                    parameter[ix*ysize+y+systemsize] = E;
                }
                else{
                    parameter[ix*ysize+y+systemsize] = 0.0;
                }
            }
        }
        //initialise x-spiral state in field region and Neel state in the rest
        for(int ix=0;ix<xsize;ix++){
            for(int y=0;y<ysize;y++){
                if(ix>=freefieldsize && ix<(freefieldsize+fieldsize)){
                    gsl_vector_set (state,ix*ysize+y,0.0*M_PI);
                    gsl_vector_set (state,ix*ysize+y+systemsize,M_PI*((y+ix)%2)-(ix-freefieldsize)*atan(E));
                }
                else{
                    gsl_vector_set (state,ix*ysize+y,0.0*M_PI);
                    if(ix<freefieldsize){
                        gsl_vector_set (state,ix*ysize+y+systemsize,M_PI*((y+ix)%2));
                    }
                    else{
                        gsl_vector_set (state,ix*ysize+y+systemsize,
                                        M_PI*((y+ix)%2)-(fieldsize)*atan(E));
                    }
                }
            }
        }
        //initialise y-spiral state in field region and Neel state in the rest
        /*for(int ix=0;ix<xsize;ix++){
            for(int y=0;y<ysize;y++){
                if(ix>=freefieldsize && ix<(freefieldsize+fieldsize)){
                    gsl_vector_set (state,ix*ysize+y,-0.5*M_PI);
                    gsl_vector_set (state,ix*ysize+y+systemsize,M_PI*((y+ix)%2)+y*atan(E));
                }
                else{
                    gsl_vector_set (state,ix*ysize+y,0.5*M_PI);
                    gsl_vector_set (state,ix*ysize+y+systemsize,M_PI*((y+ix)%2));
                }
            }
        }*/
        
        heterostruct.n = 2*systemsize;
        heterostruct.f = hamiltonian;
        heterostruct.df = dhamiltonian;
        heterostruct.fdf = pathdirection;
        heterostruct.params = parameter;
        
        /*simplex_min.n = 2*systemsize;
        simplex_min.f = hamiltonian;
        simplex_min.params = parameter;*/
    
        //T = gsl_multimin_fdfminimizer_vector_bfgs2;
        //T = gsl_multimin_fminimizer_nmsimplex2
        T = gsl_multimin_fdfminimizer_conjugate_fr;        
        s = gsl_multimin_fdfminimizer_alloc (T,2*xsize*ysize);
        //s = gsl_multimin_fminimizer_alloc(T,2*xsize*ysize);
    
        gsl_multimin_fdfminimizer_set (s, &heterostruct, state, 1e-12, 1e-12);
    
        //gsl_multimin_fminimizer_set(s,&simplex_min,state,steps);
        
        do{
            iter++;
            status = gsl_multimin_fdfminimizer_iterate (s);
            //status = gsl_multimin_fminimizer_iterate(s);

            if(status)
                break;
            
            status = gsl_multimin_test_gradient (s->gradient, 1e-8);
            //size = gsl_multimin_fminimizer_size(s);
            //status = gsl_multimin_test_size(size,1e2);
        }           
        while (status == GSL_CONTINUE && iter <1000);
        
        //raise interaction in x- and y-direction one after an other-- start
        /*
        iter = 0;
        
        for(int ix=0;ix<xsize;ix++){
            for(int y=0;y<ysize;y++){
                if(ix>=freefieldsize && ix<(freefieldsize+fieldsize)){
                    parameter[ix*ysize+y] = E;
                }
                else{
                    parameter[ix*ysize+y] = 0.0;
                }
            }
        }
        heterostruct.n = 2*systemsize;
        heterostruct.f = hamiltonian;
        heterostruct.df = dhamiltonian;
        heterostruct.fdf = pathdirection;
        heterostruct.params = parameter;
        
        //T = gsl_multimin_fdfminimizer_vector_bfgs2;
        //T = gsl_multimin_fdfminimizer_conjugate_fr;
        //s = gsl_multimin_fdfminimizer_alloc (T,2*xsize*ysize);
    
        gsl_multimin_fdfminimizer_set (s, &heterostruct, state, 0.000000001, 1e-6);
    
        do{
            iter++;
            status = gsl_multimin_fdfminimizer_iterate (s);

            if(status)
                break;
            
            status = gsl_multimin_test_gradient (s->gradient, 1e-5);
        }           
        while (status == GSL_CONTINUE && iter <1000);
        //raise interaction in x- and y-direction one after an other-- end
              
        for(int ix=0; ix<xsize;ix++){
            for(int y=0;y<ysize;y++){
                gsl_vector_set(state,ix*ysize+y,gsl_vector_get(s->x,ix*ysize+y));
                gsl_vector_set(state,ix*ysize+y+systemsize,gsl_vector_get(s->x,ix*ysize+y+systemsize));
            }
        }
        */
        
        if(count % D_OUTPUT_COUNTER ==0){
            ostringstream fn;
            fn << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << "/textureD" << abs(E)*100 << "cJxsize_simplex.dat";
            ofstream write(fn.str().c_str(),ios_base::binary);
            
            for(int ix=0; ix<xsize;ix++){
                for(int y=0;y<ysize;y++){
                    help = gsl_vector_get(s->x,ix*ysize+y);
                    write.write((char*)& help, sizeof(double));
                }
            }
            for(int ix=0; ix<xsize;ix++){
                for(int y=0;y<ysize;y++){
                    help = gsl_vector_get(s->x,ix*ysize+y+systemsize);
                    write.write((char*)& help, sizeof(double));
                }
            }
            
            write.close();
        }
        gsl_multimin_fdfminimizer_free (s);
        //gsl_multimin_fminimizer_free(s);
    }
    gsl_vector_free (state);
    gsl_vector_free(steps);
    return 0;
}
