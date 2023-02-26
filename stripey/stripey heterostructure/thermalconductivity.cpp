#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameter.h"

using namespace std;

double bosefunction(int,double,double);
void intrastripe(double[xsize],gsl_matrix_complex*,double[xsize]);
void interstripe(double[xsize],gsl_matrix_complex*,double,double[xsize]);

typedef complex<double> dcomp;

int main(){
    dcomp I;
    I = -1;
    I = sqrt(I);
    int helpint;
    double T,interaction,fourierexponent;
    gsl_complex help;
    double* texture = new double[xsize];
    double* parameter = new double[xsize];
    double* energy = new double[2*xsize];
    double* zerodrude = new double[4];
    gsl_matrix_complex* helpmatrixX = gsl_matrix_complex_alloc(2*xsize,2*xsize);
    gsl_matrix_complex* helpmatrixY = gsl_matrix_complex_alloc(2*xsize,2*xsize);
    gsl_matrix_complex* helpmatrixH = gsl_matrix_complex_alloc(2*xsize,2*xsize);
    for(int t=1;t<4;t++){
        T = t*0.1;
        ostringstream of;
        of << "Drudeweightoffield/Tdepscat/T" << T*100 << "cJ/drudeweightxsize" 
            << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        ofstream writefield(of.str().c_str());//add when using binary format: ,ios_base::binary
        
        for(int count=0;count<=(int)(Dmax/Dincrement)+1;count++){
            int site,xnegneighbour,xposneighbour;
            double xprefactora,xprefactorb,yprefactora,yprefactorb;
            complex<double> drudesumxx,drudesumxy,drudesumyx,drudesumyy;
            complex<double> cutoffepsilon;
            cutoffepsilon = {0.0,cutoffdelta};
            double zalpha,zbeta;
            complex <double> prod1xx,prod1xy,prod1yx,prod1yy,prod2xx,prod2xy,prod2yx,prod2yy;
            double energydiff;
            double omega;
            double bosealpha;
            double bosebeta;
            drudesumxx = {0.0,0.0};
            drudesumxy = {0.0,0.0};
            drudesumyx = {0.0,0.0};
            drudesumyy = {0.0,0.0};
            interaction = count*Dincrement;
            for(int x=0;x<xsize;x++){
                if(x>=freefieldsize && x<(freefieldsize+fieldsize)){
                    parameter[x] = interaction;
                }
                else{
                    parameter[x] = 0.0;
                }
            }
            ostringstream fin1;
            fin1 << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << "/textureD" << count << "cJ.dat";
            ifstream read1(fin1.str().c_str(),ios_base::binary);
            
            //read texture
            for(int x=0;x<xsize;x++){
                read1.seekg((x+xsize)*sizeof(double));
                read1.read((char*)& texture[x], sizeof(double));                   
            }
            read1.close();
            for(int k=0;k<ysize;k++){
                //set up Hamiltonian from interstripe and intrastripe interactions
            gsl_matrix_complex* hamiltonian = gsl_matrix_complex_alloc(2*xsize,2*xsize);
            gsl_matrix_complex_set_zero(hamiltonian);
            
            intrastripe(texture,helpmatrixH,parameter);
            gsl_matrix_complex_add(hamiltonian,helpmatrixH);
            interstripe(texture,helpmatrixH,fourierexponent,parameter);
            gsl_matrix_complex_add(hamiltonian,helpmatrixH);
            interstripe(texture,helpmatrixH,-fourierexponent,parameter);
            gsl_matrix_complex_add(hamiltonian,helpmatrixH);
            
            //set up current operators from analytic results
            gsl_matrix_complex* jx = gsl_matrix_complex_alloc(2*xsize,2*xsize);
            gsl_matrix_complex* jy = gsl_matrix_complex_alloc(2*xsize,2*xsize);
            gsl_matrix_complex_set_zero(jx);
            gsl_matrix_complex_set_zero(jy);
            fourierexponent = 2*M_PI*(k*1.0/ysize);
            for(int x=0;x<xsize;x++){
                site = x;
                xnegneighbour = x-1;
                xposneighbour = x+1;
                if(x<xsize-1){
                    xprefactora = cos(texture[site]-texture[xposneighbour]+M_PI) 
                            + parameter[site]*sin(texture[site]-texture[xposneighbour]+M_PI);
                }
                
                if(x>0){
                    xprefactorb = cos(texture[xnegneighbour]-texture[site]+M_PI)
                            + parameter[xnegneighbour]*sin(texture[xnegneighbour]-texture[site]+M_PI);
                }
                yprefactora = cos(0.0);
                yprefactorb = cos(0.0);
            
                //\delta=-\Delta = \pm e_x
                if(x>0 && x<xsize-1){
                    gsl_matrix_complex_set(jx,site,xposneighbour,
                            gsl_complex_add(gsl_complex_rect(-xprefactorb*(xprefactora-1.0),0.0),
                            gsl_matrix_complex_get(jx,site,xposneighbour)));//b_{l}^\dagger b_{l+e_x}
                    
                    gsl_matrix_complex_set(jx,site,xnegneighbour,
                            gsl_complex_add(gsl_complex_rect(xprefactora*(xprefactorb-1.0),0.0),
                            gsl_matrix_complex_get(jx,site,xnegneighbour)));//b_{l}^\dagger b_{l-e_x}
                    
                    gsl_matrix_complex_set(jx,xposneighbour,site,
                            gsl_complex_add(gsl_complex_rect(xprefactorb*(xprefactora-1.0),0.0),
                            gsl_matrix_complex_get(jx,xposneighbour,site)));//b_{l+e_x}^\dagger b_{l}
                    
                    gsl_matrix_complex_set(jx,xnegneighbour,site,
                            gsl_complex_add(gsl_complex_rect(-xprefactora*(xprefactorb-1.0),0.0),
                            gsl_matrix_complex_get(jx,xnegneighbour,site)));//b_{l-e_x}^\dagger b_{l}
                    
                    gsl_matrix_complex_set(jx,site,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(-1.0*xprefactorb*(xprefactora+1.0),0.0),
                            gsl_matrix_complex_get(jx,site,xposneighbour+xsize)));//b_{l}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_complex_set(jx,site,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(1.0*xprefactora*(xprefactorb+1.0),0.0)
                            ,gsl_matrix_complex_get(jx,site,xnegneighbour+xsize)));//b_{l}^\dagger b_{l-e_x}^\dagger
                    
                    gsl_matrix_complex_set(jx,xposneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(-1.0*xprefactorb*(xprefactora+1.0),0.0)
                            ,gsl_matrix_complex_get(jx,xposneighbour,site+xsize)));//b_{l+e_x}^\dagger b_{l}\dagger
                    
                    gsl_matrix_complex_set(jx,xnegneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(1.0*xprefactora*(xprefactorb+1.0),0.0)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,site+xsize)));//b_{l-e_x}^\dagger b_{l}\dagger
                    
                    gsl_matrix_complex_set(jx,xnegneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_rect(0.5*((xprefactorb-1.)*(xprefactora-1.0)-(xprefactorb+1.0)*(xprefactora+1.0)),0.0)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xposneighbour)));//b_{l-e_x}^\dagger b_{l+e_x}
                        
                    gsl_matrix_complex_set(jx,xposneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_rect(-0.5*((xprefactora-1.)*(xprefactorb-1.)-(xprefactora+1.)*(xprefactorb+1.)),0.0)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xnegneighbour)));//b_{l+e_x}^\dagger b_{l-e_x}
                    
                    gsl_matrix_complex_set(jx,xposneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*(xprefactora-1.)*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xnegneighbour+xsize)));//b_{l+e_x}^\dagger b_{l-e_x}^\dagger
                    
                    gsl_matrix_complex_set(jx,xnegneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*(xprefactorb-1.)*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xposneighbour+xsize)));//b_{l-e_x}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_complex_set(jx,xnegneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*(xprefactora-1.)*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xposneighbour+xsize)));//b_{l-e_x}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_complex_set(jx,xposneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*(xprefactorb-1.)*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xnegneighbour+xsize)));//b_{l+e_x}^\dagger b_{l-e_x}^\dagger
                    }
                    
                    //\delta=-\Delta = \pm e_y
                    gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(-yprefactorb*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l}^\dagger b_{l+e_y}
                
                    gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(yprefactora*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l}^\dagger b_{l-e_y}
                    
                    gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(yprefactorb*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l+e_y}^\dagger b_{l}
                
                    gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(-yprefactora*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l-e_y}^\dagger b_{l}
                    
                    gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-1.0*yprefactorb*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l}^\dagger b_{l+e_y}^\dagger
                    
                    gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(1.0*yprefactora*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l}^\dagger b_{l-e_y}^\dagger
                    
                    gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-1.0*yprefactorb*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l+e_y}^\dagger b_{l}\dagger
                    
                    gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(1.0*yprefactora*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l-e_y}^\dagger b_{l}\dagger
                    
                    gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*((yprefactorb-1.)*(yprefactora-1.)-(yprefactorb+1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l-e_y}^\dagger b_{l+e_y}
                    
                    gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*((yprefactora-1.)*(yprefactorb-1.)-(yprefactora+1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l+e_y}^\dagger b_{l-e_y}
                        
                    gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*(yprefactora-1.)*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l+e_y}^\dagger b_{l-e_y}^\dagger
                    
                    gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*(yprefactorb-1.)*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l-e_y}^\dagger b_{l+e_y}^\dagger
                    
                    gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*(yprefactora-1.)*(yprefactorb+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l-e_y}^\dagger b_{l+e_y}^\dagger
                    
                    gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*(yprefactorb-1.)*(yprefactora+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l+e_y}^\dagger b_{l-e_y}^\dagger
                    
                    //\delta not (-\Delta)
                    //\Delta=\pm e_x
                    if(x<xsize-1){
                        gsl_matrix_complex_set(jx,site,xposneighbour,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactora*(xprefactora-1.),0.0)
                            ,gsl_matrix_complex_get(jx,site,xposneighbour)));//b_{l}^\dagger b_{l+e_x}  \delta = +e_y
                        
                        gsl_matrix_complex_set(jy,site,xposneighbour,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactora-1.),0.0)
                            ,gsl_matrix_complex_get(jy,site,xposneighbour)));
                    
                        gsl_matrix_complex_set(jx,site,xposneighbour,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactora-1.),0.0)
                            ,gsl_matrix_complex_get(jx,site,xposneighbour)));//b_{l}^\dagger b_{l+e_x} \delta = -ey
                        
                        gsl_matrix_complex_set(jy,site,xposneighbour,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactora-1.),0.0)
                            ,gsl_matrix_complex_get(jy,site,xposneighbour)));    
                        
                        gsl_matrix_complex_set(jx,xposneighbour,site,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactora-1.),0.0)
                            ,gsl_matrix_complex_get(jx,xposneighbour,site)));//b_{l+e_x}^\dagger b_{l}  \delta = +e_y 
                        
                        gsl_matrix_complex_set(jy,xposneighbour,site,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactora*(xprefactora-1.),0.0)
                            ,gsl_matrix_complex_get(jy,xposneighbour,site)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,site,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactorb*(xprefactora-1.),0.0)
                            ,gsl_matrix_complex_get(jx,xposneighbour,site)));//b_{l+e_x}^\dagger b_{l} \delta = -ey
                        
                        gsl_matrix_complex_set(jy,xposneighbour,site,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactorb*(xprefactora-1.),0.0)
                            ,gsl_matrix_complex_get(jy,xposneighbour,site)));
                        
                        gsl_matrix_complex_set(jx,site,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactora*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jx,site,xposneighbour+xsize)));//b_{l}^\dagger b_{l+e_x}^\dagger  \delta = +e_y
                    
                        gsl_matrix_complex_set(jy,site,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jy,site,xposneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,site,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jx,site,xposneighbour+xsize)));//b_{l}^\dagger b_{l+e_x}^\dagger \delta = -ey
                        
                        gsl_matrix_complex_set(jy,site,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jy,site,xposneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactora*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jx,xposneighbour,site+xsize)));//b_{l+e_x}^\dagger b_{l}^\dagger  \delta = +e_y
                        
                        gsl_matrix_complex_set(jy,xposneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jy,xposneighbour,site+xsize)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jx,xposneighbour,site+xsize)));//b_{l+e_x}^\dagger b_{l}^\dagger \delta = -ey
                        
                        gsl_matrix_complex_set(jy,xposneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactora+1.),0.0)
                            ,gsl_matrix_complex_get(jy,xposneighbour,site+xsize)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactora-1.)*(xprefactora-1.)-(yprefactora+1.)*(xprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour)));//b_{l+e_y}^\dagger b_{l+e_x}
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactora-1.)*(xprefactora-1.)-(yprefactora+1.)*(xprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactorb-1.)*(xprefactora-1.)-(yprefactorb+1.)*(xprefactora+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour)));//b_{l-e_y}^\dagger b_{l+e_x}
                    
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactorb-1.)*(xprefactora-1.)-(yprefactorb+1.)*(xprefactora+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactora-1.)*(xprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour+xsize)));//b_{l+e_y}^\dagger b_{l+e_x}^\dagger
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactora-1.)*(xprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactorb-1.)*(xprefactora+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour+xsize)));//b_{l-e_y}^\dagger b_{l+e_x}^\dagger 
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactorb-1.)*(xprefactora+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactora-1.)*(xprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour+xsize)));//b_{l+e_x}^\dagger b_{l+e_y}^\dagger
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactora-1.)*(xprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactorb-1.)*(xprefactora+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour+xsize)));//b_{l+e_x}^\dagger b_{l-e_y}^\dagger
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactorb-1.)*(xprefactora+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour+xsize)));
                        
                        //\Delta=\pm e_y
                        gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactora*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l}^\dagger b_{l+e_y}  \delta = +e_x
                        
                        gsl_matrix_complex_set(jx,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site)));
                        
                        gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l}^\dagger b_{l-e_y} \delta = +e_x
                        
                        gsl_matrix_complex_set(jx,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site)));
                        
                        gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l+e_y}^\dagger b_{l}  \delta = +e_x
                        
                        gsl_matrix_complex_set(jx,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactora*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site)));
                        
                        gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactora*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l-e_y}^\dagger b_{l} \delta = +e_x
                        
                        gsl_matrix_complex_set(jx,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactora*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site))); 
                        
                        gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactora*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l}^\dagger b_{l+e_y}^\dagger  \delta = +e_x
                        
                        gsl_matrix_complex_set(jx,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site+xsize)));
                        
                        gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l}^\dagger b_{l-e_y}^\dagger \delta = +e_x
                        
                        gsl_matrix_complex_set(jx,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site+xsize))); 
                        
                        gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactora*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l+e_y}^\dagger b_{l}^\dagger  \delta = +e_x
                        
                        gsl_matrix_complex_set(jx,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site+xsize)));
                        
                        gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l-e_y}^\dagger b_{l}^\dagger \delta = +e_x
                        
                        gsl_matrix_complex_set(jx,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactora*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site+xsize)));
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactora-1.)*(yprefactora-1.)-(xprefactora+1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour)));//b_{l+e_x}^\dagger b_{l+e_y}
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactora-1.)*(yprefactora-1.)-(xprefactora+1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour)));
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactora-1.)*(yprefactorb-1.)-(xprefactora+1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour)));//b_{l+e_x}^\dagger b_{l-e_y}
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactora-1.)*(yprefactorb-1.)-(xprefactora+1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour)));
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactora-1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour+xsize)));//b_{l+e_x}^\dagger b_{l+e_y}^\dagger
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactora-1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactora-1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour+xsize)));//b_{l+e_x}^\dagger b_{l-e_y}^\dagger
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactora-1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactora-1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour+xsize)));//b_{l+e_y}^\dagger b_{l+e_x}^\dagger
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactora-1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jy,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactora-1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xposneighbour,xposneighbour+xsize)));//b_{l-e_y}^\dagger b_{l+e_x}^\dagger
                        
                        gsl_matrix_complex_set(jx,xposneighbour,xposneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactora-1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xposneighbour,xposneighbour+xsize)));
                    }
                    
                    if(x>0){
                        gsl_matrix_complex_set(jx,site,xnegneighbour,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactorb-1.),0.0)
                            ,gsl_matrix_complex_get(jx,site,xnegneighbour)));//b_{l}^\dagger b_{l-e_x} \delta = +e_y
                        
                        gsl_matrix_complex_set(jy,site,xnegneighbour,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactorb-1.),0.0)
                            ,gsl_matrix_complex_get(jy,site,xnegneighbour)));
                        
                        gsl_matrix_complex_set(jx,site,xnegneighbour,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactorb*(xprefactorb-1.),0.0)
                            ,gsl_matrix_complex_get(jx,site,xnegneighbour)));//b_{l}^\dagger b_{l-e_x} \delta = -e_y
                        
                        gsl_matrix_complex_set(jy,site,xnegneighbour,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactorb-1.),0.0)
                            ,gsl_matrix_complex_get(jy,site,xnegneighbour)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,site,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactora*(xprefactorb-1.),0.0)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,site)));//b_{l-e_x}^\dagger b_{l} \delta = +e_y
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,site,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactora*(xprefactorb-1.),0.0)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,site)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,site,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactorb-1.),0.0)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,site)));//b_{l-e_x}^\dagger b_{l} \delta = -e_y
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,site,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactorb*(xprefactorb-1.),0.0)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,site)));
                        
                        gsl_matrix_complex_set(jx,site,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jx,site,xnegneighbour+xsize)));//b_{l}^\dagger b_{l-e_x}^\dagger \delta = +e_y
                        
                        gsl_matrix_complex_set(jy,site,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jy,site,xnegneighbour+xsize)));      
                        
                        gsl_matrix_complex_set(jx,site,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactorb*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jx,site,xnegneighbour+xsize)));//b_{l}^\dagger b_{l-e_x}^\dagger \delta = -e_y
                        
                        gsl_matrix_complex_set(jy,site,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jy,site,xnegneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,site+xsize)));//b_{l-e_x}^\dagger b_{l}^\dagger \delta = +e_y
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactora*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,site+xsize)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(0.5*yprefactorb*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,site+xsize)));//b_{l-e_x}^\dagger b_{l}^\dagger \delta = -e_y
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,site+xsize,
                            gsl_complex_add(gsl_complex_rect(-0.5*yprefactorb*(xprefactorb+1.),0.0)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,site+xsize)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactora-1.)*(xprefactorb-1.)-(yprefactora+1.)*(xprefactorb+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour)));//b_{l+e_y}^\dagger b_{l-e_x}
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactora-1.)*(xprefactorb-1.)-(yprefactora+1.)*(xprefactorb+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactorb-1.)*(xprefactorb-1.)-(yprefactorb+1.)*(xprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour)));//b_{l-e_y}^\dagger b_{l-e_x}
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactorb-1.)*(xprefactorb-1.)-(yprefactorb+1.)*(xprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactora-1.)*(xprefactorb+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour+xsize)));//b_{l+e_y}^\dagger b_{l-e_x}^\dagger
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactora-1.)*(xprefactorb+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactorb-1.)*(xprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour+xsize)));//b_{l-e_y}^\dagger b_{l-e_x}^\dagger
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactorb-1.)*(xprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactora-1.)*(xprefactorb+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour+xsize)));//b_{l-e_x}^\dagger b_{l+e_y}^\dagger
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactora-1.)*(xprefactorb+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((yprefactorb-1.)*(xprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour+xsize)));//b_{l-e_x}^\dagger b_{l-e_y}^\dagger
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((yprefactorb-1.)*(xprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour+xsize)));
                        
                        //\Delta=\pm e_y
                        gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l}^\dagger b_{l+e_y} \delta = -e_x
                        
                        gsl_matrix_complex_set(jx,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site)));
                        
                        gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactorb*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l}^\dagger b_{l-e_y} \delta = -e_x
                        
                        gsl_matrix_complex_set(jx,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site)));
                        
                        gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactorb*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l+e_y}^\dagger b_{l} \delta = -e_x
                        
                        gsl_matrix_complex_set(jx,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactorb*(yprefactora-1.),fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site)));
                        
                        gsl_matrix_complex_set(jy,site,site,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site)));//b_{l-e_y}^\dagger b_{l} \delta = -e_x
                        
                        gsl_matrix_complex_set(jx,site,site,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactorb*(yprefactorb-1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site)));
                        
                        gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l}^\dagger b_{l+e_y}^\dagger \delta = -e_x
                    
                        gsl_matrix_complex_set(jx,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site+xsize)));
                        
                        gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactorb*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l}^\dagger b_{l-e_y}^\dagger \delta = -e_x
                        
                        gsl_matrix_complex_set(jx,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site+xsize)));
                        
                        gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l+e_y}^\dagger b_{l}^\dagger \delta = -e_x
                    
                        gsl_matrix_complex_set(jx,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactora+1.),fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site+xsize)));
                        
                        gsl_matrix_complex_set(jy,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(0.5*xprefactorb*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,site,site+xsize)));//b_{l-e_y}^\dagger b_{l}^\dagger \delta = -e_x
                        
                        gsl_matrix_complex_set(jx,site,site+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.5*xprefactorb*(yprefactorb+1.),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,site,site+xsize)));
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactorb-1.)*(yprefactora-1.)-(xprefactorb+1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour)));//b_{l-e_x}^\dagger b_{l+e_y}
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactorb-1.)*(yprefactora-1.)-(xprefactorb+1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour)));
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactorb-1.)*(yprefactorb-1.)-(xprefactorb+1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour)));//b_{l-e_x}^\dagger b_{l-e_y}
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactorb-1.)*(yprefactorb-1.)-(xprefactorb+1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour)));
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactorb-1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour+xsize)));//b_{l-e_x}^\dagger b_{l+e_y}^\dagger
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactorb-1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactorb-1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour+xsize)));//b_{l-e_x}^\dagger b_{l-e_y}^\dagger
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactorb-1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactorb-1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour+xsize)));//b_{l+e_y}^\dagger b_{l-e_x}^\dagger
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactorb-1.)*(yprefactora+1.)),fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour+xsize)));
                        
                        gsl_matrix_complex_set(jy,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(0.25*((xprefactorb-1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jy,xnegneighbour,xnegneighbour+xsize)));//b_{l-e_y}^\dagger b_{l-e_x}^\dagger
                        
                        gsl_matrix_complex_set(jx,xnegneighbour,xnegneighbour+xsize,
                            gsl_complex_add(gsl_complex_polar(-0.25*((xprefactorb-1.)*(yprefactorb+1.)),-fourierexponent)
                            ,gsl_matrix_complex_get(jx,xnegneighbour,xnegneighbour+xsize)));
                    }
                }

                //add complex conjugated elements
                for(int x=0;x<xsize;x++){
                    for(int y=0;y<xsize;y++){
                        gsl_matrix_complex_set(jx,y+xsize,x+xsize,gsl_complex_conjugate(gsl_matrix_complex_get(jx,x,y)));
                        gsl_matrix_complex_set(jy,y+xsize,x+xsize,gsl_complex_conjugate(gsl_matrix_complex_get(jy,x,y)));
                    }
                }

                for(int x=0;x<xsize;x++){
                    for(int y=0;y<xsize;y++){
                        gsl_matrix_complex_set(jx,y+xsize,x,gsl_complex_negative(gsl_complex_conjugate(gsl_matrix_complex_get(jx,x,y+xsize))));
                        gsl_matrix_complex_set(jy,y+xsize,x,gsl_complex_negative(gsl_complex_conjugate(gsl_matrix_complex_get(jy,x,y+xsize))));
                    }
                }   
                    
                //multiply by prefactor S^2J^2/8
                gsl_matrix_complex_scale(jx,gsl_complex_rect(0.125*S*S*J*J,0.0));
                gsl_matrix_complex_scale(jy,gsl_complex_rect(0.125*S*S*J*J,0.0));

                //read Bogoliubov matrix
                ostringstream fin2;
                fin2  << "Bogoliubov/Eigenvectors/k" << k << "D" << count << "cJxsize" << xsize << "ysize" << ysize 
                        << "fieldsize" << fieldsize << ".dat";
                ifstream read2(fin2.str().c_str(),ios_base::binary);
                gsl_matrix_complex* bogo = gsl_matrix_complex_alloc(2*xsize,2*xsize);
                for(int x=0;x<2*xsize;x++){
                    for(int y=0;y<2*xsize;y++){
                        read2.seekg((2*x*xsize+y)*sizeof(gsl_complex));
                        read2.read((char*)& help, sizeof(gsl_complex));                     
                        gsl_matrix_complex_set(bogo,x,y,help);
                    }
                }
                read2.close();
                    
                //Check diagonalsiation property
                gsl_matrix_complex* bogohamilton = gsl_matrix_complex_alloc(2*xsize,2*xsize);
                gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),bogo,hamiltonian,gsl_complex_rect(0.0,0.0),helpmatrixH);
                gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),helpmatrixH,bogo,gsl_complex_rect(0.0,0.0),bogohamilton);
                
                gsl_matrix_complex_free (hamiltonian);
                
                //transform current operators
                gsl_matrix_complex* bogojx = gsl_matrix_complex_alloc(2*xsize,2*xsize);
                gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),bogo,jx,gsl_complex_rect(0.0,0.0),helpmatrixX);
                gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),helpmatrixX,bogo,gsl_complex_rect(0.0,0.0),bogojx);
                
                gsl_matrix_complex* bogojy = gsl_matrix_complex_alloc(2*xsize,2*xsize);
                gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),bogo,jy,gsl_complex_rect(0.0,0.0),helpmatrixY);
                gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),helpmatrixY,bogo,gsl_complex_rect(0.0,0.0),bogojy);
                
                gsl_matrix_complex_free (jx);
                gsl_matrix_complex_free (jy);
                gsl_matrix_complex_free (bogo);
                        
                //calculate Drude weight
                //read energy-eigenvalues
                for(int i=0;i<2*xsize;i++){
                    energy[i] = 2*GSL_REAL(gsl_matrix_complex_get(bogohamilton,i,i));
                }
                
                gsl_matrix_complex_free (bogohamilton);
                
                //ostringstream of1;
                //of1 << "Drudeweight/drudeweightnumD" << D*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
                //ofstream writedrude(of1.str().c_str());//add when using binary format: ,ios_base::binary
                for(int o=0;o<0.2/omegastep;o++){
                    omega = o*omegastep;
                    complex<double>* omegasum = new complex<double>[4];
                    for(int i=0; i<4;i++){
                        omegasum[i] = {0.0,0.0};
                    }
                    for(int alpha=0;alpha<2*xsize;alpha++){
                        for(int beta=0;beta<2*xsize;beta++){
                            if(alpha < xsize){
                                zalpha = 1.0;
                            }
                            else{
                                zalpha = -1.0;
                            }
                            if(beta < xsize){
                                zbeta = 1.0;
                            }
                            else{
                                zbeta = -1.0;
                            }
                            energydiff = zalpha*energy[alpha]-zbeta*energy[beta];
                            
                            help = gsl_complex_mul(gsl_matrix_complex_get(bogojx,alpha,beta),gsl_matrix_complex_get(bogojx,
                                        (alpha+xsize)%(2*xsize),(beta+xsize)%(2*xsize)));
                            prod1xx = gsl_complex_abs(help)*exp(I*gsl_complex_arg(help));
                            help = gsl_complex_mul(gsl_matrix_complex_get(bogojx,alpha,beta),gsl_matrix_complex_get(bogojx,beta,alpha));
                            prod2xx = gsl_complex_abs(help)*exp(I*gsl_complex_arg(help));
                            help = gsl_complex_mul(gsl_matrix_complex_get(bogojy,alpha,beta),gsl_matrix_complex_get(bogojy,
                                        (alpha+xsize)%(2*xsize),(beta+xsize)%(2*xsize)));
                            prod1yy = gsl_complex_abs(help)*exp(I*gsl_complex_arg(help));
                            help = gsl_complex_mul(gsl_matrix_complex_get(bogojy,alpha,beta),gsl_matrix_complex_get(bogojy,beta,alpha));
                            prod2yy = gsl_complex_abs(help)*exp(I*gsl_complex_arg(help));
                            help =  gsl_complex_mul(gsl_matrix_complex_get(bogojx,alpha,beta),gsl_matrix_complex_get(bogojy,
                                        (alpha+xsize)%(2*xsize),(beta+xsize)%(2*xsize)));
                            prod1xy = gsl_complex_abs(help)*exp(I*gsl_complex_arg(help));
                            help = gsl_complex_mul(gsl_matrix_complex_get(bogojx,alpha,beta),gsl_matrix_complex_get(bogojy,beta,alpha));
                            prod2xy = gsl_complex_abs(help)*exp(I*gsl_complex_arg(help));
                            help = gsl_complex_mul(gsl_matrix_complex_get(bogojy,alpha,beta),gsl_matrix_complex_get(bogojx,
                                        (alpha+xsize)%(2*xsize),(beta+xsize)%(2*xsize)));
                            prod1yx = gsl_complex_abs(help)*exp(I*gsl_complex_arg(help));
                            help = gsl_complex_mul(gsl_matrix_complex_get(bogojy,alpha,beta),gsl_matrix_complex_get(bogojx,beta,alpha));
                            prod2yx = gsl_complex_abs(help)*exp(I*gsl_complex_arg(help));
                            bosealpha = bosefunction(alpha,energy[alpha],T);
                            bosebeta = bosefunction((beta+xsize)%(2*xsize),energy[beta],T);
                            omegasum[0] += bosealpha*bosebeta*(prod1xx+prod2xx)/(omega - energydiff + T*T*cutoffepsilon);
                            omegasum[1] += bosealpha*bosebeta*(prod1xy+prod2xy)/(omega - energydiff + T*T*cutoffepsilon);
                            omegasum[2] += bosealpha*bosebeta*(prod1yx+prod2yx)/(omega - energydiff + T*T*cutoffepsilon);
                            omegasum[3] += bosealpha*bosebeta*(prod1yy+prod2yy)/(omega - energydiff + T*T*cutoffepsilon);
                            //drudesumxx += bosealpha*bosebeta*(prod1xx+prod2xx)/(omega - energydiff + T*T*cutoffepsilon);
                            //drudesumxy += bosealpha*bosebeta*(prod1xy+prod2xy)/(omega - energydiff + T*T*cutoffepsilon);
                            //drudesumyx += bosealpha*bosebeta*(prod1yx+prod2yx)/(omega - energydiff + T*T*cutoffepsilon);
                            //drudesumyy += bosealpha*bosebeta*(prod1yy+prod2yy)/(omega - energydiff + T*T*cutoffepsilon);
                        }
                    }
                    if(o!=0){
                        drudesumxx += omegasum[0]*(1.0-exp(-omega/T))/omega;
                        drudesumxy += omegasum[1]*(1.0-exp(-omega/T))/omega;
                        drudesumyx += omegasum[2]*(1.0-exp(-omega/T))/omega;
                        drudesumyy += omegasum[3]*(1.0-exp(-omega/T))/omega;
                    }
                    else{
                        drudesumxx += omegasum[0]*T;
                        drudesumxy += omegasum[1]*T;
                        drudesumyx += omegasum[2]*T;
                        drudesumyy += omegasum[3]*T;
                    }
                }
                if(count==0){
                    zerodrude[0] = imag(drudesumxx/((complex<double>)(xsize*ysize)));
                    zerodrude[1] = imag(drudesumxy/((complex<double>)(xsize*ysize)));
                    zerodrude[2] = imag(drudesumyx/((complex<double>)(xsize*ysize)));
                    zerodrude[3] = imag(drudesumyy/((complex<double>)(xsize*ysize)));
                }
                gsl_matrix_complex_free(bogojx);
                gsl_matrix_complex_free(bogojy);
            }
            cout << interaction   << '\t' << imag(drudesumxx/((complex<double>)(xsize*ysize))) 
                 << '\t' << imag(drudesumxy/((complex<double>)(xsize*ysize))) 
                 << '\t' << imag(drudesumyx/((complex<double>)(xsize*ysize))) 
                 << '\t' << imag(drudesumyy/((complex<double>)(xsize*ysize))) 
                 << '\n';
               
            writefield  << interaction << '\t' << imag(drudesumxx/((complex<double>)(xsize*ysize)))/zerodrude[0] 
                 << '\t' << imag(drudesumxy/((complex<double>)(xsize*ysize)))/zerodrude[1] 
                 << '\t' << imag(drudesumyx/((complex<double>)(xsize*ysize)))/zerodrude[2] 
                 << '\t' << imag(drudesumyy/((complex<double>)(xsize*ysize)))/zerodrude[3] 
                 << '\n';
        }
        writefield.close();
    }
    gsl_matrix_complex_free (helpmatrixH);
    gsl_matrix_complex_free (helpmatrixX);
    gsl_matrix_complex_free (helpmatrixY);
    delete [] texture;
    delete [] parameter;
    delete [] energy;
    delete [] zerodrude;    
    return 0;
}

double bosefunction(int index, double energy, double T){
    if(index<xsize){
        return(1.0/(exp(energy/T)-1));
    }
    else if(index>=xsize){
        return(1.0/(exp(energy/T)-1)+1);
    }
}

void intrastripe(double *texture, gsl_matrix_complex *edmatrix, double *D){
   double xprefactora,xprefactorb;
   int site,xnegneighbour,xposneighbour;
   gsl_matrix_complex_set_zero(edmatrix);
   for(int x=0;x<xsize;x++){
      site = x;
      xposneighbour = (x+1);
      xnegneighbour = (x-1);
           if(x==0){
               xprefactora = cos(texture[site]-texture[xposneighbour]+M_PI) 
                            +D[site]*sin(texture[site]-texture[xposneighbour]+M_PI);
               
               gsl_matrix_complex_set(edmatrix,x,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                                  gsl_complex_rect(-0.5*(xprefactora)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x+1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x+1),
                                  gsl_complex_rect(-0.5*(xprefactora)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x,x+1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+1),
                                  gsl_complex_rect(0.25*(xprefactora-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x),
                                  gsl_complex_rect(0.25*(xprefactora-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x+xsize),
                                  gsl_complex_rect(0.25*(xprefactora+1.),0.0)));
               gsl_matrix_complex_set(edmatrix,x,x+1+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+1+xsize),
                                  gsl_complex_rect(0.25*(xprefactora+1.),0.0)));
           }
           else if(x==xsize-1){
               xprefactorb = cos(texture[xnegneighbour]-texture[site]+M_PI) 
                            +D[xnegneighbour]*sin(texture[xnegneighbour]-texture[site]+M_PI);
                
               gsl_matrix_complex_set(edmatrix,x,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                                  gsl_complex_rect(-0.5*(xprefactorb)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x-1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x-1),
                                  gsl_complex_rect(-0.5*(xprefactorb)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x,x-1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x-1),
                                  gsl_complex_rect(0.25*(xprefactorb-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x),
                                  gsl_complex_rect(0.25*(xprefactorb-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x+xsize),
                                  gsl_complex_rect(0.25*(xprefactorb+1.),0.0)));
               gsl_matrix_complex_set(edmatrix,x,x-1+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x-1+xsize),
                                  gsl_complex_rect(0.25*(xprefactorb+1.),0.0)));
           }
           else{
               xprefactora = cos(texture[site]-texture[xposneighbour]+M_PI) 
                            +D[site]*sin(texture[site]-texture[xposneighbour]+M_PI);
               
               gsl_matrix_complex_set(edmatrix,x,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                                  gsl_complex_rect(-0.5*(xprefactora)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x+1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x+1),
                                  gsl_complex_rect(-0.5*(xprefactora)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x,x+1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+1),
                                  gsl_complex_rect(0.25*(xprefactora-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x),
                                  gsl_complex_rect(0.25*(xprefactora-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x+xsize),
                                  gsl_complex_rect(0.25*(xprefactora+1.),0.0)));
               gsl_matrix_complex_set(edmatrix,x,x+1+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+1+xsize),
                                  gsl_complex_rect(0.25*(xprefactora+1.),0.0)));
               
               xprefactorb = cos(texture[xnegneighbour]-texture[site]+M_PI) 
                            +D[xnegneighbour]*sin(texture[xnegneighbour]-texture[site]+M_PI);
                
               gsl_matrix_complex_set(edmatrix,x,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                                  gsl_complex_rect(-0.5*(xprefactorb)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x-1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x-1),
                                  gsl_complex_rect(-0.5*(xprefactorb)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x,x-1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x-1),
                                  gsl_complex_rect(0.25*(xprefactorb-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x),
                                  gsl_complex_rect(0.25*(xprefactorb-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x+xsize),
                                  gsl_complex_rect(0.25*(xprefactorb+1.),0.0)));
               gsl_matrix_complex_set(edmatrix,x,x-1+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x-1+xsize),
                                  gsl_complex_rect(0.25*(xprefactorb+1.),0.0)));
           }
   }
   //add complex conjugated elements
    for(int i=0;i<xsize;i++){
        for(int j=0;j<xsize;j++){
            gsl_matrix_complex_set(edmatrix,i+xsize,j+xsize, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j)));
        }
    }
    for(int i=0;i<xsize;i++){
        for(int j=0;j<xsize;j++){
            gsl_matrix_complex_set(edmatrix,j+xsize,i, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j+xsize)));
        }
    }
    gsl_matrix_complex_scale(edmatrix,gsl_complex_rect(-0.5*J*S,0.0));
}

void interstripe(double *texture, gsl_matrix_complex *edmatrix, double fourierexponent, double *D){
   double yprefactora,yprefactorb;
   gsl_matrix_complex_set_zero(edmatrix);
   for(int x=0;x<xsize;x++){
      yprefactora = cos(0.0);
      yprefactorb = cos(0.0);
      gsl_matrix_complex_set(edmatrix,x,x,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                    gsl_complex_rect(-0.5*(yprefactora)-epsilon,0.0)));
      gsl_matrix_complex_set(edmatrix,x,x,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                    gsl_complex_rect(-0.5*(yprefactorb)-epsilon,0.0)));
      gsl_matrix_complex_set(edmatrix,x,x,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                    gsl_complex_polar(0.25*(yprefactora-1.),fourierexponent)));
      gsl_matrix_complex_set(edmatrix,x,x,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                    gsl_complex_polar(0.25*(yprefactorb-1.),fourierexponent)));
      gsl_matrix_complex_set(edmatrix,x,x+xsize,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+xsize),
                    gsl_complex_polar(0.25*(yprefactora+1.),fourierexponent)));
      gsl_matrix_complex_set(edmatrix,x,x+xsize,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+xsize),
                    gsl_complex_polar(0.25*(yprefactorb+1.),fourierexponent)));
    
    }
    //add complex conjugated elements
    for(int i=0;i<xsize;i++){
        for(int j=0;j<xsize;j++){
            gsl_matrix_complex_set(edmatrix,i+xsize,j+xsize, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j)));
        }
    }
    for(int i=0;i<xsize;i++){
        for(int j=0;j<xsize;j++){
            gsl_matrix_complex_set(edmatrix,j+xsize,i, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j+xsize)));
        }
    }
    gsl_matrix_complex_scale(edmatrix,gsl_complex_rect(-0.5*J*S,0.0));
}