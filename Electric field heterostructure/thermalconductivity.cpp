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
double bosefunction1(int,double,double);
void commutator(gsl_matrix*,gsl_matrix*,gsl_matrix*);
void energydensity(double[systemsize],gsl_matrix*,int,int,double[systemsize]);

int main(){
    int helpint;
    double D,help, helpx, helpy;
    double* texture = new double[systemsize];
    double* parameter = new double[systemsize];
    double* energy = new double[2*systemsize];
    
    gsl_matrix* helpmatrixX = gsl_matrix_alloc(2*systemsize,2*systemsize);
    gsl_matrix* helpmatrixY = gsl_matrix_alloc(2*systemsize,2*systemsize);
    gsl_matrix* helpmatrixH = gsl_matrix_alloc(2*systemsize,2*systemsize);
    //gsl_matrix* paraunity = gsl_matrix_alloc(2*systemsize,2*systemsize);
    //gsl_matrix_set_identity(paraunity);
    /*for(int x=0;x<systemsize;x++){
        gsl_matrix_set(paraunity,x+systemsize,x+systemsize,-1.0);
    }*/
    double T;
    for(int t=1;t<4;t++){
    T = t*0.1;
    ostringstream of;
    of << "Drudeweightoffield/Tdepscat/T" << T*100 << "cJ/drudeweightxsize" 
        << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
    ofstream writefield(of.str().c_str());//add when using binary format: ,ios_base::binary
    
    for(int count=0;count<=(int)(Dmax/Dincrement)+1;count++){
        D = count*Dincrement;
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                if(x>=freefieldsize && x<(freefieldsize+fieldsize)){
                    parameter[x*ysize+y] = D;
                }
                else{
                    parameter[x*ysize+y] = 0.0;
                }
            }
        }
        ostringstream fin1;
        fin1 << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" 
            << fieldsize << "/textureD" << D*100 << "cJxsize_simplex.dat";
        ifstream read1(fin1.str().c_str(),ios_base::binary);
        
        //read texture
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                read1.seekg((x*ysize+y+systemsize)*sizeof(double));
                read1.read((char*)& texture[x*ysize+y], sizeof(double));                   
            }
        }
        read1.close();
        
        //set up hamiltonian from energydensity
        gsl_matrix* energydensityX = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix* energydensityY = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix* hamiltonian = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix_set_zero(energydensityX);
        gsl_matrix_set_zero(energydensityY);
        gsl_matrix_set_zero(hamiltonian);
        
        for(int i=0;i<xsize;i++){
            for(int j=0;j<ysize;j++){
                energydensity(texture,helpmatrixH,i,j,parameter);
                
                gsl_matrix_add(hamiltonian,helpmatrixH);
                gsl_matrix_memcpy(helpmatrixX,helpmatrixH);
                gsl_matrix_scale(helpmatrixX,i);
                gsl_matrix_memcpy(helpmatrixY,helpmatrixH);
                gsl_matrix_scale(helpmatrixY,j);
                gsl_matrix_add(energydensityX,helpmatrixX);
                gsl_matrix_add(energydensityY,helpmatrixY);
            }
        }
        /*if(count==Dmax*100){
        //if(count==1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(hamiltonian,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //set up current operators numerically
        gsl_matrix* njx = gsl_matrix_alloc(2*systemsize,2*systemsize);
        commutator(hamiltonian,energydensityX,njx);
        
        gsl_matrix* njy = gsl_matrix_alloc(2*systemsize,2*systemsize);
        commutator(hamiltonian,energydensityY,njy);
        
        gsl_matrix_free (energydensityX);
        gsl_matrix_free (energydensityY);
        
        //set up current operators from analytic results
        gsl_matrix* jx = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix* jy = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix_set_zero(jx);
        gsl_matrix_set_zero(jy);
        
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                int site,xnegneighbour,xposneighbour,ynegneighbour,yposneighbour;
                double xprefactora,xprefactorb,yprefactora,yprefactorb;
                
                site = x*ysize+y;
                xnegneighbour = (x-1)*ysize+y;
                xposneighbour = (x+1)*ysize+y;
                ynegneighbour = x*ysize+y-1;
                yposneighbour = x*ysize+y+1;
                
                if(y==0){
                    ynegneighbour = x*ysize+ysize-1;
                }
                
                if(y==ysize-1){
                    yposneighbour = x*ysize;
                }
                
                if(x<xsize-1){
                    xprefactora = cos(texture[site]-texture[xposneighbour]+M_PI) 
                            + parameter[site]*sin(texture[site]-texture[xposneighbour]+M_PI);
                }
                
                if(x>0){
                    xprefactorb = cos(texture[xnegneighbour]-texture[site]+M_PI)
                            + parameter[xnegneighbour]*sin(texture[xnegneighbour]-texture[site]+M_PI);
                }
                yprefactora = cos(texture[site]-texture[yposneighbour]+M_PI);
                yprefactorb = cos(texture[ynegneighbour]-texture[site]+M_PI);
            
                //\delta=-\Delta = \pm e_x
                if(x>0 && x<xsize-1){
                    gsl_matrix_set(jx,site,xposneighbour,
                               -xprefactorb*(xprefactora-1)
                               +gsl_matrix_get(jx,site,xposneighbour));//b_{l}^\dagger b_{l+e_x}
                    
                    gsl_matrix_set(jx,site,xnegneighbour,
                               xprefactora*(xprefactorb-1)
                               +gsl_matrix_get(jx,site,xnegneighbour));//b_{l}^\dagger b_{l-e_x}
                    
                    gsl_matrix_set(jx,xposneighbour,site,
                               xprefactorb*(xprefactora-1)
                               +gsl_matrix_get(jx,xposneighbour,site));//b_{l+e_x}^\dagger b_{l}
                    
                    gsl_matrix_set(jx,xnegneighbour,site,
                               -xprefactora*(xprefactorb-1)
                               +gsl_matrix_get(jx,xnegneighbour,site));//b_{l-e_x}^\dagger b_{l}
                    
                    gsl_matrix_set(jx,site,xposneighbour+systemsize,
                               -2*xprefactorb*(xprefactora+1)
                               +gsl_matrix_get(jx,site,xposneighbour+systemsize));//b_{l}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_set(jx,site,xnegneighbour+systemsize,
                               2*xprefactora*(xprefactorb+1)
                               +gsl_matrix_get(jx,site,xnegneighbour+systemsize));//b_{l}^\dagger b_{l-e_x}^\dagger
                    
                    gsl_matrix_set(jx,xposneighbour,site+systemsize,
                               -2*xprefactorb*(xprefactora+1)
                               +gsl_matrix_get(jx,xposneighbour,site+systemsize));//b_{l+e_x}^\dagger b_{l}\dagger
                    
                    gsl_matrix_set(jx,xnegneighbour,site+systemsize,
                               2*xprefactora*(xprefactorb+1)
                               +gsl_matrix_get(jx,xnegneighbour,site+systemsize));//b_{l-e_x}^\dagger b_{l}\dagger
                    
                    gsl_matrix_set(jx,xnegneighbour,xposneighbour,
                               0.5*((xprefactorb-1)*(xprefactora-1)-(xprefactorb+1)*(xprefactora+1))
                               +gsl_matrix_get(jx,xnegneighbour,xposneighbour));//b_{l-e_x}^\dagger b_{l+e_x}
                    
                    gsl_matrix_set(jx,xposneighbour,xnegneighbour,
                               -0.5*((xprefactora-1)*(xprefactorb-1)-(xprefactora+1)*(xprefactorb+1))
                               +gsl_matrix_get(jx,xposneighbour,xnegneighbour));//b_{l+e_x}^\dagger b_{l-e_x}
                    
                    gsl_matrix_set(jx,xposneighbour,xnegneighbour+systemsize,
                               (xprefactora-1)*(xprefactorb+1)
                               +gsl_matrix_get(jx,xposneighbour,xnegneighbour+systemsize));//b_{l+e_x}^\dagger b_{l-e_x}^\dagger
                    
                    gsl_matrix_set(jx,xnegneighbour,xposneighbour+systemsize,
                               -(xprefactorb-1)*(xprefactora+1)
                               +gsl_matrix_get(jx,xnegneighbour,xposneighbour+systemsize));//b_{l-e_x}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_set(jx,xnegneighbour,xposneighbour+systemsize,
                               (xprefactora-1)*(xprefactorb+1)
                               +gsl_matrix_get(jx,xnegneighbour,xposneighbour+systemsize));//b_{l-e_x}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_set(jx,xposneighbour,xnegneighbour+systemsize,
                               -(xprefactorb-1)*(xprefactora+1)
                               +gsl_matrix_get(jx,xposneighbour,xnegneighbour+systemsize));//b_{l+e_x}^\dagger b_{l-e_x}^\dagger
                }
                
                //\delta=-\Delta = \pm e_y
                gsl_matrix_set(jy,site,yposneighbour,
                               -yprefactorb*(yprefactora-1)
                               +gsl_matrix_get(jy,site,yposneighbour));//b_{l}^\dagger b_{l+e_y}
                
                gsl_matrix_set(jy,site,ynegneighbour,
                               yprefactora*(yprefactorb-1)
                               +gsl_matrix_get(jy,site,ynegneighbour));//b_{l}^\dagger b_{l-e_y}
                
                gsl_matrix_set(jy,yposneighbour,site,
                               yprefactorb*(yprefactora-1)
                               +gsl_matrix_get(jy,yposneighbour,site));//b_{l+e_y}^\dagger b_{l}
                
                gsl_matrix_set(jy,ynegneighbour,site,
                               -yprefactora*(yprefactorb-1)
                               +gsl_matrix_get(jy,ynegneighbour,site));//b_{l-e_y}^\dagger b_{l}
                
                gsl_matrix_set(jy,site,yposneighbour+systemsize,
                               -2*yprefactorb*(yprefactora+1)
                               +gsl_matrix_get(jy,site,yposneighbour+systemsize));//b_{l}^\dagger b_{l+e_y}^\dagger
                
                gsl_matrix_set(jy,site,ynegneighbour+systemsize,
                               2*yprefactora*(yprefactorb+1)
                               +gsl_matrix_get(jy,site,ynegneighbour+systemsize));//b_{l}^\dagger b_{l-e_y}^\dagger
                
                gsl_matrix_set(jy,yposneighbour,site+systemsize,
                               -2*yprefactorb*(yprefactora+1)
                               +gsl_matrix_get(jy,yposneighbour,site+systemsize));//b_{l+e_y}^\dagger b_{l}\dagger
                
                gsl_matrix_set(jy,ynegneighbour,site+systemsize,
                               2*yprefactora*(yprefactorb+1)
                               +gsl_matrix_get(jy,ynegneighbour,site+systemsize));//b_{l-e_y}^\dagger b_{l}\dagger
                
                gsl_matrix_set(jy,ynegneighbour,yposneighbour,
                               0.5*((yprefactorb-1)*(yprefactora-1)-(yprefactorb+1)*(yprefactora+1))
                               +gsl_matrix_get(jy,ynegneighbour,yposneighbour));//b_{l-e_y}^\dagger b_{l+e_y}
                
                gsl_matrix_set(jy,yposneighbour,ynegneighbour,
                               -0.5*((yprefactora-1)*(yprefactorb-1)-(yprefactora+1)*(yprefactorb+1))
                               +gsl_matrix_get(jy,yposneighbour,ynegneighbour));//b_{l+e_y}^\dagger b_{l-e_y}
                    
                gsl_matrix_set(jy,yposneighbour,ynegneighbour+systemsize,
                               (yprefactora-1)*(yprefactorb+1)
                               +gsl_matrix_get(jy,yposneighbour,ynegneighbour+systemsize));//b_{l+e_y}^\dagger b_{l-e_y}^\dagger
                
                gsl_matrix_set(jy,ynegneighbour,yposneighbour+systemsize,
                               -(yprefactorb-1)*(yprefactora+1)
                               +gsl_matrix_get(jy,ynegneighbour,yposneighbour+systemsize));//b_{l-e_y}^\dagger b_{l+e_y}^\dagger
                
                gsl_matrix_set(jy,ynegneighbour,yposneighbour+systemsize,
                               (yprefactora-1)*(yprefactorb+1)
                               +gsl_matrix_get(jy,ynegneighbour,yposneighbour+systemsize));//b_{l-e_y}^\dagger b_{l+e_y}^\dagger
                
                gsl_matrix_set(jy,yposneighbour,ynegneighbour+systemsize,
                               -(yprefactorb-1)*(yprefactora+1)
                               +gsl_matrix_get(jy,yposneighbour,ynegneighbour+systemsize));//b_{l+e_y}^\dagger b_{l-e_y}^\dagger
                
                //\delta not (-\Delta)
                //\Delta=\pm e_x
                if(x<xsize-1){
                    gsl_matrix_set(jx,site,xposneighbour,
                               -0.5*yprefactora*(xprefactora-1)
                               +gsl_matrix_get(jx,site,xposneighbour));//b_{l}^\dagger b_{l+e_x}  \delta = +e_y
                    
                    gsl_matrix_set(jy,site,xposneighbour,
                               0.5*yprefactora*(xprefactora-1)
                               +gsl_matrix_get(jy,site,xposneighbour));
                    
                    gsl_matrix_set(jx,site,xposneighbour,
                               -0.5*yprefactorb*(xprefactora-1)
                               +gsl_matrix_get(jx,site,xposneighbour));//b_{l}^\dagger b_{l+e_x} \delta = -ey
                    
                    gsl_matrix_set(jy,site,xposneighbour,
                               -0.5*yprefactorb*(xprefactora-1)
                               +gsl_matrix_get(jy,site,xposneighbour));    
                    
                    gsl_matrix_set(jx,xposneighbour,site,
                               0.5*yprefactora*(xprefactora-1)
                               +gsl_matrix_get(jx,xposneighbour,site));//b_{l+e_x}^\dagger b_{l}  \delta = +e_y 
                    
                    gsl_matrix_set(jy,xposneighbour,site,
                               -0.5*yprefactora*(xprefactora-1)
                               +gsl_matrix_get(jy,xposneighbour,site));
                    
                    gsl_matrix_set(jx,xposneighbour,site,
                               0.5*yprefactorb*(xprefactora-1)
                               +gsl_matrix_get(jx,xposneighbour,site));//b_{l+e_x}^\dagger b_{l} \delta = -ey
                    
                    gsl_matrix_set(jy,xposneighbour,site,
                               0.5*yprefactorb*(xprefactora-1)
                               +gsl_matrix_get(jy,xposneighbour,site));
                    
                    gsl_matrix_set(jx,site,xposneighbour+systemsize,
                               -yprefactora*(xprefactora+1)
                               +gsl_matrix_get(jx,site,xposneighbour+systemsize));//b_{l}^\dagger b_{l+e_x}^\dagger  \delta = +e_y
                    
                    gsl_matrix_set(jy,site,xposneighbour+systemsize,
                               yprefactora*(xprefactora+1)
                               +gsl_matrix_get(jy,site,xposneighbour+systemsize));
                    
                    gsl_matrix_set(jx,site,xposneighbour+systemsize,
                               -yprefactorb*(xprefactora+1)
                               +gsl_matrix_get(jx,site,xposneighbour+systemsize));//b_{l}^\dagger b_{l+e_x}^\dagger \delta = -ey
                    
                    gsl_matrix_set(jy,site,xposneighbour+systemsize,
                               -yprefactorb*(xprefactora+1)
                               +gsl_matrix_get(jy,site,xposneighbour+systemsize));
                    
                    gsl_matrix_set(jx,xposneighbour,site+systemsize,
                               -yprefactora*(xprefactora+1)
                               +gsl_matrix_get(jx,xposneighbour,site+systemsize));//b_{l+e_x}^\dagger b_{l}^\dagger  \delta = +e_y
                    
                    gsl_matrix_set(jy,xposneighbour,site+systemsize,
                               yprefactora*(xprefactora+1)
                               +gsl_matrix_get(jy,xposneighbour,site+systemsize));
                    
                    gsl_matrix_set(jx,xposneighbour,site+systemsize,
                               -yprefactorb*(xprefactora+1)
                               +gsl_matrix_get(jx,xposneighbour,site+systemsize));//b_{l+e_x}^\dagger b_{l}^\dagger \delta = -ey
                    
                    gsl_matrix_set(jy,xposneighbour,site+systemsize,
                               -yprefactorb*(xprefactora+1)
                               +gsl_matrix_get(jy,xposneighbour,site+systemsize));
                    
                    gsl_matrix_set(jx,yposneighbour,xposneighbour,
                               0.25*((yprefactora-1)*(xprefactora-1)-(yprefactora+1)*(xprefactora+1))
                               +gsl_matrix_get(jx,yposneighbour,xposneighbour));//b_{l+e_y}^\dagger b_{l+e_x}
                    
                    gsl_matrix_set(jy,yposneighbour,xposneighbour,
                               -0.25*((yprefactora-1)*(xprefactora-1)-(yprefactora+1)*(xprefactora+1))
                               +gsl_matrix_get(jy,yposneighbour,xposneighbour));
                    
                    gsl_matrix_set(jx,ynegneighbour,xposneighbour,
                               0.25*((yprefactorb-1)*(xprefactora-1)-(yprefactorb+1)*(xprefactora+1))
                               +gsl_matrix_get(jx,ynegneighbour,xposneighbour));//b_{l-e_y}^\dagger b_{l+e_x}
                    
                    gsl_matrix_set(jy,ynegneighbour,xposneighbour,
                               0.25*((yprefactorb-1)*(xprefactora-1)-(yprefactorb+1)*(xprefactora+1))
                               +gsl_matrix_get(jy,ynegneighbour,xposneighbour));
                    
                    gsl_matrix_set(jx,yposneighbour,xposneighbour+systemsize,
                               -0.5*((yprefactora-1)*(xprefactora+1))
                               +gsl_matrix_get(jx,yposneighbour,xposneighbour+systemsize));//b_{l+e_y}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_set(jy,yposneighbour,xposneighbour+systemsize,
                               0.5*((yprefactora-1)*(xprefactora+1))
                               +gsl_matrix_get(jy,yposneighbour,xposneighbour+systemsize));
                    
                    gsl_matrix_set(jx,ynegneighbour,xposneighbour+systemsize,
                               -0.5*((yprefactorb-1)*(xprefactora+1))
                               +gsl_matrix_get(jx,ynegneighbour,xposneighbour+systemsize));//b_{l-e_y}^\dagger b_{l+e_x}^\dagger 
                    
                    gsl_matrix_set(jy,ynegneighbour,xposneighbour+systemsize,
                               -0.5*((yprefactorb-1)*(xprefactora+1))
                               +gsl_matrix_get(jy,ynegneighbour,xposneighbour+systemsize));
                    
                    gsl_matrix_set(jx,xposneighbour,yposneighbour+systemsize,
                               -0.5*((yprefactora-1)*(xprefactora+1))
                               +gsl_matrix_get(jx,xposneighbour,yposneighbour+systemsize));//b_{l+e_x}^\dagger b_{l+e_y}^\dagger
                    
                    gsl_matrix_set(jy,xposneighbour,yposneighbour+systemsize,
                               0.5*((yprefactora-1)*(xprefactora+1))
                               +gsl_matrix_get(jy,xposneighbour,yposneighbour+systemsize));
                    
                    gsl_matrix_set(jx,xposneighbour,ynegneighbour+systemsize,
                               -0.5*((yprefactorb-1)*(xprefactora+1))
                               +gsl_matrix_get(jx,xposneighbour,ynegneighbour+systemsize));//b_{l+e_x}^\dagger b_{l-e_y}^\dagger
                    
                    gsl_matrix_set(jy,xposneighbour,ynegneighbour+systemsize,
                               -0.5*((yprefactorb-1)*(xprefactora+1))
                               +gsl_matrix_get(jy,xposneighbour,ynegneighbour+systemsize));
                    
                    //\Delta=\pm e_y
                    gsl_matrix_set(jy,site,yposneighbour,
                               -0.5*xprefactora*(yprefactora-1)
                               +gsl_matrix_get(jy,site,yposneighbour));//b_{l}^\dagger b_{l+e_y}  \delta = +e_x
                    gsl_matrix_set(jx,site,yposneighbour,
                               0.5*xprefactora*(yprefactora-1)
                               +gsl_matrix_get(jx,site,yposneighbour));
                    
                    gsl_matrix_set(jy,site,ynegneighbour,
                               0.5*xprefactora*(yprefactorb-1)
                               +gsl_matrix_get(jy,site,ynegneighbour));//b_{l}^\dagger b_{l-e_y} \delta = +e_x
                    
                    gsl_matrix_set(jx,site,ynegneighbour,
                               0.5*xprefactora*(yprefactorb-1)
                               +gsl_matrix_get(jx,site,ynegneighbour));
                    
                    gsl_matrix_set(jy,yposneighbour,site,
                               0.5*xprefactora*(yprefactora-1)
                               +gsl_matrix_get(jy,yposneighbour,site));//b_{l+e_y}^\dagger b_{l}  \delta = +e_x
                    
                    gsl_matrix_set(jx,yposneighbour,site,
                               -0.5*xprefactora*(yprefactora-1)
                               +gsl_matrix_get(jx,yposneighbour,site));
                    
                    gsl_matrix_set(jy,ynegneighbour,site,
                               -0.5*xprefactora*(yprefactorb-1)
                               +gsl_matrix_get(jy,ynegneighbour,site));//b_{l-e_y}^\dagger b_{l} \delta = +e_x
                    
                    gsl_matrix_set(jx,ynegneighbour,site,
                               -0.5*xprefactora*(yprefactorb-1)
                               +gsl_matrix_get(jx,ynegneighbour,site)); 
                    
                    gsl_matrix_set(jy,site,yposneighbour+systemsize,
                               -xprefactora*(yprefactora+1)
                               +gsl_matrix_get(jy,site,yposneighbour+systemsize));//b_{l}^\dagger b_{l+e_y}^\dagger  \delta = +e_x
                    
                    gsl_matrix_set(jx,site,yposneighbour+systemsize,
                               xprefactora*(yprefactora+1)
                               +gsl_matrix_get(jx,site,yposneighbour+systemsize));
                    
                    gsl_matrix_set(jy,site,ynegneighbour+systemsize,
                               xprefactora*(yprefactorb+1)
                               +gsl_matrix_get(jy,site,ynegneighbour+systemsize));//b_{l}^\dagger b_{l-e_y}^\dagger \delta = +e_x
                    
                    gsl_matrix_set(jx,site,ynegneighbour+systemsize,
                               xprefactora*(yprefactorb+1)
                               +gsl_matrix_get(jx,site,ynegneighbour+systemsize)); 
                    
                    gsl_matrix_set(jy,yposneighbour,site+systemsize,
                               -xprefactora*(yprefactora+1)
                               +gsl_matrix_get(jy,yposneighbour,site+systemsize));//b_{l+e_y}^\dagger b_{l}^\dagger  \delta = +e_x
                    
                    gsl_matrix_set(jx,yposneighbour,site+systemsize,
                               xprefactora*(yprefactora+1)
                               +gsl_matrix_get(jx,yposneighbour,site+systemsize));
                    
                    gsl_matrix_set(jy,ynegneighbour,site+systemsize,
                               xprefactora*(yprefactorb+1)
                               +gsl_matrix_get(jy,ynegneighbour,site+systemsize));//b_{l-e_y}^\dagger b_{l}^\dagger \delta = +e_x
                    
                    gsl_matrix_set(jx,ynegneighbour,site+systemsize,
                               xprefactora*(yprefactorb+1)
                               +gsl_matrix_get(jx,ynegneighbour,site+systemsize));
                    
                    gsl_matrix_set(jy,xposneighbour,yposneighbour,
                               0.25*((xprefactora-1)*(yprefactora-1)-(xprefactora+1)*(yprefactora+1))
                               +gsl_matrix_get(jy,xposneighbour,yposneighbour));//b_{l+e_x}^\dagger b_{l+e_y}
                    
                    gsl_matrix_set(jx,xposneighbour,yposneighbour,
                               -0.25*((xprefactora-1)*(yprefactora-1)-(xprefactora+1)*(yprefactora+1))
                               +gsl_matrix_get(jx,xposneighbour,yposneighbour));
                    
                    gsl_matrix_set(jy,xposneighbour,ynegneighbour,
                               -0.25*((xprefactora-1)*(yprefactorb-1)-(xprefactora+1)*(yprefactorb+1))
                               +gsl_matrix_get(jy,xposneighbour,ynegneighbour));//b_{l+e_x}^\dagger b_{l-e_y}
                    
                    gsl_matrix_set(jx,xposneighbour,ynegneighbour,
                               -0.25*((xprefactora-1)*(yprefactorb-1)-(xprefactora+1)*(yprefactorb+1))
                               +gsl_matrix_get(jx,xposneighbour,ynegneighbour));
                    
                    gsl_matrix_set(jy,xposneighbour,yposneighbour+systemsize,
                               -0.5*((xprefactora-1)*(yprefactora+1))
                               +gsl_matrix_get(jy,xposneighbour,yposneighbour+systemsize));//b_{l+e_x}^\dagger b_{l+e_y}^\dagger
                    
                    gsl_matrix_set(jx,xposneighbour,yposneighbour+systemsize,
                               0.5*((xprefactora-1)*(yprefactora+1))
                               +gsl_matrix_get(jx,xposneighbour,yposneighbour+systemsize));
                    
                    gsl_matrix_set(jy,xposneighbour,ynegneighbour+systemsize,
                               0.5*((xprefactora-1)*(yprefactorb+1))
                               +gsl_matrix_get(jy,xposneighbour,ynegneighbour+systemsize));//b_{l+e_x}^\dagger b_{l-e_y}^\dagger
                    
                    gsl_matrix_set(jx,xposneighbour,ynegneighbour+systemsize,
                               0.5*((xprefactora-1)*(yprefactorb+1))
                               +gsl_matrix_get(jx,xposneighbour,ynegneighbour+systemsize));
                    
                    gsl_matrix_set(jy,yposneighbour,xposneighbour+systemsize,
                               -0.5*((xprefactora-1)*(yprefactora+1))
                               +gsl_matrix_get(jy,yposneighbour,xposneighbour+systemsize));//b_{l+e_y}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_set(jx,yposneighbour,xposneighbour+systemsize,
                               0.5*((xprefactora-1)*(yprefactora+1))
                               +gsl_matrix_get(jx,yposneighbour,xposneighbour+systemsize));
                    
                    gsl_matrix_set(jy,ynegneighbour,xposneighbour+systemsize,
                               0.5*((xprefactora-1)*(yprefactorb+1))
                               +gsl_matrix_get(jy,ynegneighbour,xposneighbour+systemsize));//b_{l-e_y}^\dagger b_{l+e_x}^\dagger
                    
                    gsl_matrix_set(jx,ynegneighbour,xposneighbour+systemsize,
                               0.5*((xprefactora-1)*(yprefactorb+1))
                               +gsl_matrix_get(jx,ynegneighbour,xposneighbour+systemsize));
                }
                
                if(x>0){
                    gsl_matrix_set(jx,site,xnegneighbour,
                               0.5*yprefactora*(xprefactorb-1)
                               +gsl_matrix_get(jx,site,xnegneighbour));//b_{l}^\dagger b_{l-e_x} \delta = +e_y
                    
                    gsl_matrix_set(jy,site,xnegneighbour,
                               0.5*yprefactora*(xprefactorb-1)
                               +gsl_matrix_get(jy,site,xnegneighbour));
                    
                    gsl_matrix_set(jx,site,xnegneighbour,
                               0.5*yprefactorb*(xprefactorb-1)
                               +gsl_matrix_get(jx,site,xnegneighbour));//b_{l}^\dagger b_{l-e_x} \delta = -e_y
                    
                    gsl_matrix_set(jy,site,xnegneighbour,
                               -0.5*yprefactorb*(xprefactorb-1)
                               +gsl_matrix_get(jy,site,xnegneighbour));
                    
                    gsl_matrix_set(jx,xnegneighbour,site,
                               -0.5*yprefactora*(xprefactorb-1)
                               +gsl_matrix_get(jx,xnegneighbour,site));//b_{l-e_x}^\dagger b_{l} \delta = +e_y
                    
                    gsl_matrix_set(jy,xnegneighbour,site,
                               -0.5*yprefactora*(xprefactorb-1)
                               +gsl_matrix_get(jy,xnegneighbour,site));
                    
                    gsl_matrix_set(jx,xnegneighbour,site,
                               -0.5*yprefactorb*(xprefactorb-1)
                               +gsl_matrix_get(jx,xnegneighbour,site));//b_{l-e_x}^\dagger b_{l} \delta = -e_y
                    
                    gsl_matrix_set(jy,xnegneighbour,site,
                               0.5*yprefactorb*(xprefactorb-1)
                               +gsl_matrix_get(jy,xnegneighbour,site));
                    
                    gsl_matrix_set(jx,site,xnegneighbour+systemsize,
                               yprefactora*(xprefactorb+1)
                               +gsl_matrix_get(jx,site,xnegneighbour+systemsize));//b_{l}^\dagger b_{l-e_x}^\dagger \delta = +e_y
                    
                    gsl_matrix_set(jy,site,xnegneighbour+systemsize,
                               yprefactora*(xprefactorb+1)
                               +gsl_matrix_get(jy,site,xnegneighbour+systemsize));      
                    
                    gsl_matrix_set(jx,site,xnegneighbour+systemsize,
                               yprefactorb*(xprefactorb+1)
                               +gsl_matrix_get(jx,site,xnegneighbour+systemsize));//b_{l}^\dagger b_{l-e_x}^\dagger \delta = -e_y
                    
                    gsl_matrix_set(jy,site,xnegneighbour+systemsize,
                               -yprefactorb*(xprefactorb+1)
                               +gsl_matrix_get(jy,site,xnegneighbour+systemsize));
                    
                    gsl_matrix_set(jx,xnegneighbour,site+systemsize,
                               yprefactora*(xprefactorb+1)
                               +gsl_matrix_get(jx,xnegneighbour,site+systemsize));//b_{l-e_x}^\dagger b_{l}^\dagger \delta = +e_y
                    
                    gsl_matrix_set(jy,xnegneighbour,site+systemsize,
                               yprefactora*(xprefactorb+1)
                               +gsl_matrix_get(jy,xnegneighbour,site+systemsize));
                    
                    gsl_matrix_set(jx,xnegneighbour,site+systemsize,
                               yprefactorb*(xprefactorb+1)
                               +gsl_matrix_get(jx,xnegneighbour,site+systemsize));//b_{l-e_x}^\dagger b_{l}^\dagger \delta = -e_y
                    
                    gsl_matrix_set(jy,xnegneighbour,site+systemsize,
                               -yprefactorb*(xprefactorb+1)
                               +gsl_matrix_get(jy,xnegneighbour,site+systemsize));
                    
                    gsl_matrix_set(jx,yposneighbour,xnegneighbour,
                               -0.25*((yprefactora-1)*(xprefactorb-1)-(yprefactora+1)*(xprefactorb+1))
                               +gsl_matrix_get(jx,yposneighbour,xnegneighbour));//b_{l+e_y}^\dagger b_{l-e_x}
                    
                    gsl_matrix_set(jy,yposneighbour,xnegneighbour,
                               -0.25*((yprefactora-1)*(xprefactorb-1)-(yprefactora+1)*(xprefactorb+1))
                               +gsl_matrix_get(jy,yposneighbour,xnegneighbour));
                    
                    gsl_matrix_set(jx,ynegneighbour,xnegneighbour,
                               -0.25*((yprefactorb-1)*(xprefactorb-1)-(yprefactorb+1)*(xprefactorb+1))
                               +gsl_matrix_get(jx,ynegneighbour,xnegneighbour));//b_{l-e_y}^\dagger b_{l-e_x}
                    
                    gsl_matrix_set(jy,ynegneighbour,xnegneighbour,
                               0.25*((yprefactorb-1)*(xprefactorb-1)-(yprefactorb+1)*(xprefactorb+1))
                               +gsl_matrix_get(jy,ynegneighbour,xnegneighbour));
                    
                    gsl_matrix_set(jx,yposneighbour,xnegneighbour+systemsize,
                               0.5*((yprefactora-1)*(xprefactorb+1))
                               +gsl_matrix_get(jx,yposneighbour,xnegneighbour+systemsize));//b_{l+e_y}^\dagger b_{l-e_x}^\dagger
                    
                    gsl_matrix_set(jy,yposneighbour,xnegneighbour+systemsize,
                               0.5*((yprefactora-1)*(xprefactorb+1))
                               +gsl_matrix_get(jy,yposneighbour,xnegneighbour+systemsize));
                    
                    gsl_matrix_set(jx,ynegneighbour,xnegneighbour+systemsize,
                               0.5*((yprefactorb-1)*(xprefactorb+1))
                               +gsl_matrix_get(jx,ynegneighbour,xnegneighbour+systemsize));//b_{l-e_y}^\dagger b_{l-e_x}^\dagger
                    
                    gsl_matrix_set(jy,ynegneighbour,xnegneighbour+systemsize,
                               -0.5*((yprefactorb-1)*(xprefactorb+1))
                               +gsl_matrix_get(jy,ynegneighbour,xnegneighbour+systemsize));
                    
                    gsl_matrix_set(jx,xnegneighbour,yposneighbour+systemsize,
                               0.5*((yprefactora-1)*(xprefactorb+1))
                               +gsl_matrix_get(jx,xnegneighbour,yposneighbour+systemsize));//b_{l-e_x}^\dagger b_{l+e_y}^\dagger
                    
                    gsl_matrix_set(jy,xnegneighbour,yposneighbour+systemsize,
                               0.5*((yprefactora-1)*(xprefactorb+1))
                               +gsl_matrix_get(jy,xnegneighbour,yposneighbour+systemsize));
                    
                    gsl_matrix_set(jx,xnegneighbour,ynegneighbour+systemsize,
                               0.5*((yprefactorb-1)*(xprefactorb+1))
                               +gsl_matrix_get(jx,xnegneighbour,ynegneighbour+systemsize));//b_{l-e_x}^\dagger b_{l-e_y}^\dagger
                    
                    gsl_matrix_set(jy,xnegneighbour,ynegneighbour+systemsize,
                               -0.5*((yprefactorb-1)*(xprefactorb+1))
                               +gsl_matrix_get(jy,xnegneighbour,ynegneighbour+systemsize));
                    
                    //\Delta=\pm e_y
                    gsl_matrix_set(jy,site,yposneighbour,
                               -0.5*xprefactorb*(yprefactora-1)
                               +gsl_matrix_get(jy,site,yposneighbour));//b_{l}^\dagger b_{l+e_y} \delta = -e_x
                    
                    gsl_matrix_set(jx,site,yposneighbour,
                               -0.5*xprefactorb*(yprefactora-1)
                               +gsl_matrix_get(jx,site,yposneighbour));
                    
                    gsl_matrix_set(jy,site,ynegneighbour,
                               0.5*xprefactorb*(yprefactorb-1)
                               +gsl_matrix_get(jy,site,ynegneighbour));//b_{l}^\dagger b_{l-e_y} \delta = -e_x
                    
                    gsl_matrix_set(jx,site,ynegneighbour,
                               -0.5*xprefactorb*(yprefactorb-1)
                               +gsl_matrix_get(jx,site,ynegneighbour));
                    
                    gsl_matrix_set(jy,yposneighbour,site,
                               0.5*xprefactorb*(yprefactora-1)
                               +gsl_matrix_get(jy,yposneighbour,site));//b_{l+e_y}^\dagger b_{l} \delta = -e_x
                    
                    gsl_matrix_set(jx,yposneighbour,site,
                               0.5*xprefactorb*(yprefactora-1)
                               +gsl_matrix_get(jx,yposneighbour,site));
                    
                    gsl_matrix_set(jy,ynegneighbour,site,
                               -0.5*xprefactorb*(yprefactorb-1)
                               +gsl_matrix_get(jy,ynegneighbour,site));//b_{l-e_y}^\dagger b_{l} \delta = -e_x
                    
                    gsl_matrix_set(jx,ynegneighbour,site,
                               0.5*xprefactorb*(yprefactorb-1)
                               +gsl_matrix_get(jx,ynegneighbour,site));
                    
                    gsl_matrix_set(jy,site,yposneighbour+systemsize,
                               -xprefactorb*(yprefactora+1)
                               +gsl_matrix_get(jy,site,yposneighbour+systemsize));//b_{l}^\dagger b_{l+e_y}^\dagger \delta = -e_x
                
                    gsl_matrix_set(jx,site,yposneighbour+systemsize,
                               -xprefactorb*(yprefactora+1)
                               +gsl_matrix_get(jx,site,yposneighbour+systemsize));
                    
                    gsl_matrix_set(jy,site,ynegneighbour+systemsize,
                               xprefactorb*(yprefactorb+1)
                               +gsl_matrix_get(jy,site,ynegneighbour+systemsize));//b_{l}^\dagger b_{l-e_y}^\dagger \delta = -e_x
                    
                    gsl_matrix_set(jx,site,ynegneighbour+systemsize,
                               -xprefactorb*(yprefactorb+1)
                               +gsl_matrix_get(jx,site,ynegneighbour+systemsize));
                    
                    gsl_matrix_set(jy,yposneighbour,site+systemsize,
                               -xprefactorb*(yprefactora+1)
                               +gsl_matrix_get(jy,yposneighbour,site+systemsize));//b_{l+e_y}^\dagger b_{l}^\dagger \delta = -e_x
                
                    gsl_matrix_set(jx,yposneighbour,site+systemsize,
                               -xprefactorb*(yprefactora+1)
                               +gsl_matrix_get(jx,yposneighbour,site+systemsize));
                    
                    gsl_matrix_set(jy,ynegneighbour,site+systemsize,
                               xprefactorb*(yprefactorb+1)
                               +gsl_matrix_get(jy,ynegneighbour,site+systemsize));//b_{l-e_y}^\dagger b_{l}^\dagger \delta = -e_x
                    
                    gsl_matrix_set(jx,ynegneighbour,site+systemsize,
                               -xprefactorb*(yprefactorb+1)
                               +gsl_matrix_get(jx,ynegneighbour,site+systemsize));
                    
                    gsl_matrix_set(jy,xnegneighbour,yposneighbour,
                               0.25*((xprefactorb-1)*(yprefactora-1)-(xprefactorb+1)*(yprefactora+1))
                               +gsl_matrix_get(jy,xnegneighbour,yposneighbour));//b_{l-e_x}^\dagger b_{l+e_y}
                    
                    gsl_matrix_set(jx,xnegneighbour,yposneighbour,
                               0.25*((xprefactorb-1)*(yprefactora-1)-(xprefactorb+1)*(yprefactora+1))
                               +gsl_matrix_get(jx,xnegneighbour,yposneighbour));
                    
                    gsl_matrix_set(jy,xnegneighbour,ynegneighbour,
                               -0.25*((xprefactorb-1)*(yprefactorb-1)-(xprefactorb+1)*(yprefactorb+1))
                               +gsl_matrix_get(jy,xnegneighbour,ynegneighbour));//b_{l-e_x}^\dagger b_{l-e_y}
                    
                    gsl_matrix_set(jx,xnegneighbour,ynegneighbour,
                               0.25*((xprefactorb-1)*(yprefactorb-1)-(xprefactorb+1)*(yprefactorb+1))
                               +gsl_matrix_get(jx,xnegneighbour,ynegneighbour));
                    
                    gsl_matrix_set(jy,xnegneighbour,yposneighbour+systemsize,
                               -0.5*((xprefactorb-1)*(yprefactora+1))
                               +gsl_matrix_get(jy,xnegneighbour,yposneighbour+systemsize));//b_{l-e_x}^\dagger b_{l+e_y}^\dagger
                    
                    gsl_matrix_set(jx,xnegneighbour,yposneighbour+systemsize,
                               -0.5*((xprefactorb-1)*(yprefactora+1))
                               +gsl_matrix_get(jx,xnegneighbour,yposneighbour+systemsize));
                    
                    gsl_matrix_set(jy,xnegneighbour,ynegneighbour+systemsize,
                               0.5*((xprefactorb-1)*(yprefactorb+1))
                               +gsl_matrix_get(jy,xnegneighbour,ynegneighbour+systemsize));//b_{l-e_x}^\dagger b_{l-e_y}^\dagger
                    
                    gsl_matrix_set(jx,xnegneighbour,ynegneighbour+systemsize,
                               -0.5*((xprefactorb-1)*(yprefactorb+1))
                               +gsl_matrix_get(jx,xnegneighbour,ynegneighbour+systemsize));
                    
                    gsl_matrix_set(jy,yposneighbour,xnegneighbour+systemsize,
                               -0.5*((xprefactorb-1)*(yprefactora+1))
                               +gsl_matrix_get(jy,yposneighbour,xnegneighbour+systemsize));//b_{l+e_y}^\dagger b_{l-e_x}^\dagger
                    
                    gsl_matrix_set(jx,yposneighbour,xnegneighbour+systemsize,
                               -0.5*((xprefactorb-1)*(yprefactora+1))
                               +gsl_matrix_get(jx,yposneighbour,xnegneighbour+systemsize));
                    
                    gsl_matrix_set(jy,ynegneighbour,xnegneighbour+systemsize,
                               0.5*((xprefactorb-1)*(yprefactorb+1))
                               +gsl_matrix_get(jy,ynegneighbour,xnegneighbour+systemsize));//b_{l-e_y}^\dagger b_{l-e_x}^\dagger
                    
                    gsl_matrix_set(jx,ynegneighbour,xnegneighbour+systemsize,
                               -0.5*((xprefactorb-1)*(yprefactorb+1))
                               +gsl_matrix_get(jx,ynegneighbour,xnegneighbour+systemsize));
                }
            }
        }

        //add complex conjugated elements
        for(int x=0;x<systemsize;x++){
            for(int y=0;y<systemsize;y++){
                gsl_matrix_set(jx,y+systemsize,x+systemsize,gsl_matrix_get(jx,x,y));
                gsl_matrix_set(jy,y+systemsize,x+systemsize,gsl_matrix_get(jy,x,y));
            }
        }

        for(int x=0;x<systemsize;x++){
            for(int y=0;y<systemsize;y++){
                gsl_matrix_set(jx,y+systemsize,x,-gsl_matrix_get(jx,x,y+systemsize));
                gsl_matrix_set(jy,y+systemsize,x,-gsl_matrix_get(jy,x,y+systemsize));
            }
        }   
        
        //multiply by prefactor S^2J^2/8
        gsl_matrix_scale(jx,0.125*S*S*J*J);
        gsl_matrix_scale(jy,0.125*S*S*J*J);

        //if(count==1){
        if(count==Dmax*100){
            ostringstream fx;
            fx << "jxa.dat";
            ofstream writex(fx.str().c_str());
            ostringstream fy;
            fy << "jya.dat";
            ofstream writey(fy.str().c_str());
            ostringstream fnx;
            fnx << "jxn.dat";
            ofstream writenx(fnx.str().c_str());
            ostringstream fny;
            fny << "jyn.dat";
            ofstream writeny(fny.str().c_str());
            for(int i=0;i<2*systemsize;i++){
                for(int j=0;j<2*systemsize;j++){
                    writex << gsl_matrix_get(jx,i,j) << '\t';
                    writey << gsl_matrix_get(jy,i,j) << '\t';
                    writenx << gsl_matrix_get(njx,i,j) << '\t';
                    writeny << gsl_matrix_get(njy,i,j) << '\t';  
                }
                writex << '\n';
                writey << '\n';
                writenx << '\n';
                writeny << '\n';
            }
            writex.close();
            writey.close();
            writenx.close();
            writeny.close();
        }
        
        gsl_matrix_free(njx);
        gsl_matrix_free(njy);
        
        //read Bogoliubov matrix
        ostringstream fin2;
        fin2 << "Bogoliubov/eigenvectorD" << D*100 << "cJxsize" << xsize 
            << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";;
        ifstream read2(fin2.str().c_str(),ios_base::binary);
        gsl_matrix* bogo = gsl_matrix_alloc(2*systemsize,2*systemsize);
        for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                read2.seekg((2*x*systemsize+y)*sizeof(double));
                read2.read((char*)& help, sizeof(double));                     
                gsl_matrix_set(bogo,x,y,help);
            }
        }
        read2.close();
        
        //Check diagonalsiation property
        gsl_matrix* bogohamilton = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,bogo,hamiltonian,0.0,helpmatrixH);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,helpmatrixH,bogo,0.0,bogohamilton);
        
        gsl_matrix_free (hamiltonian);
        /*if(count == 1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(bogohamilton,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //check paraunitary property
        /*
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,bogo,paraunity,0.0,helpmatrixH);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,helpmatrixH,bogo,0.0,helpmatrixX);
        
        
        if(count == 1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(helpmatrixX,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //transform current operators
        gsl_matrix* bogojx = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,bogo,jx,0.0,helpmatrixX);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,helpmatrixX,bogo,0.0,bogojx);
        
        gsl_matrix* bogojy = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,bogo,jy,0.0,helpmatrixY);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,helpmatrixY,bogo,0.0,bogojy);
        
        gsl_matrix_free (jx);
        gsl_matrix_free (jy);
        gsl_matrix_free (bogo);
        
        /*if(count == 1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(bogojx,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //calculate Drude weight
        //read energy-eigenvalues
        for(int i=0;i<2*systemsize;i++){
            energy[i] = 2*gsl_matrix_get(bogohamilton,i,i);
        }
        
        gsl_matrix_free (bogohamilton);
        
        //ostringstream of1;
        //of1 << "Drudeweight/drudeweightnumD" << D*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        //ofstream writedrude(of1.str().c_str());//add when using binary format: ,ios_base::binary
        complex<double> drudesumxx,drudesumxy,drudesumyx,drudesumyy;
        complex<double> cutoffepsilon;
        cutoffepsilon = {0.0,cutoffdelta};
        double zalpha,zbeta;
        double prod1xx,prod1xy,prod1yx,prod1yy,prod2xx,prod2xy,prod2yx,prod2yy;
        double energydiff;
        double omega = 0.0;
        double bosealpha;
        double bosebeta;
        
        //for(double T=0.001;T<0.5;T+=0.01){
        
        drudesumxx = {0.0,0.0};
        drudesumxy = {0.0,0.0};
        drudesumyx = {0.0,0.0};
        drudesumyy = {0.0,0.0};
            
        for(int alpha=0;alpha<2*systemsize;alpha++){
            for(int beta=0;beta<2*systemsize;beta++){
                if(alpha < systemsize){
                    zalpha = 1.0;
                }
                else{
                    zalpha = -1.0;
                }
                if(beta < systemsize){
                    zbeta = 1.0;
                }
                else{
                    zbeta = -1.0;
                }
                energydiff = zalpha*energy[alpha]-zbeta*energy[beta];
                prod1xx = gsl_matrix_get(bogojx,alpha,beta)*gsl_matrix_get(bogojx,(alpha+systemsize)%(2*systemsize),(beta+systemsize)%(2*systemsize));
                prod2xx = gsl_matrix_get(bogojx,alpha,beta)*gsl_matrix_get(bogojx,beta,alpha);
                prod1yy = gsl_matrix_get(bogojy,alpha,beta)*gsl_matrix_get(bogojy,(alpha+systemsize)%(2*systemsize),(beta+systemsize)%(2*systemsize));
                prod2yy = gsl_matrix_get(bogojy,alpha,beta)*gsl_matrix_get(bogojy,beta,alpha);
                prod1xy = gsl_matrix_get(bogojx,alpha,beta)*gsl_matrix_get(bogojy,(alpha+systemsize)%(2*systemsize),(beta+systemsize)%(2*systemsize));
                prod2xy = gsl_matrix_get(bogojx,alpha,beta)*gsl_matrix_get(bogojy,beta,alpha);
                prod1yx = gsl_matrix_get(bogojy,alpha,beta)*gsl_matrix_get(bogojx,(alpha+systemsize)%(2*systemsize),(beta+systemsize)%(2*systemsize));
                prod2yx = gsl_matrix_get(bogojy,alpha,beta)*gsl_matrix_get(bogojx,beta,alpha);
                bosealpha = bosefunction(alpha,energy[alpha],T);
                bosebeta = bosefunction((beta+systemsize)%(2*systemsize),energy[beta],T);                        
                drudesumxx += bosealpha*bosebeta*(prod1xx+prod2xx)/(omega - energydiff + T*T*cutoffepsilon);
                drudesumxy += bosealpha*bosebeta*(prod1xy+prod2xy)/(omega - energydiff + T*T*cutoffepsilon);
                drudesumyx += bosealpha*bosebeta*(prod1yx+prod2yx)/(omega - energydiff + T*T*cutoffepsilon);
                drudesumyy += bosealpha*bosebeta*(prod1yy+prod2yy)/(omega - energydiff + T*T*cutoffepsilon);
            }
        }
    
    writefield << D << '\t' << imag(drudesumxx/((complex<double>)(systemsize))) << '\t' << imag(drudesumxy/((complex<double>)(systemsize)))<< '\t' << imag(drudesumyx/((complex<double>)(2*systemsize)))<< '\t' << imag(drudesumyy/((complex<double>)(systemsize)))<< '\n';
    
        //writedrude << T << '\t' << imag(drudesumxx/((complex<double>)(systemsize))) << '\t' << imag(drudesumxy/((complex<double>)(systemsize)))<< '\t' << imag(drudesumyx/((complex<double>)(systemsize)))<< '\t' << imag(drudesumyy/((complex<double>)(systemsize)))<< '\n';
        //}
        //writedrude.close();
        
        
        //current-current-correlationsfunction of omega for T==0.9
        //if(D>0.8999 && D<0.9001){
        /*complex<double> omegasumxx,omegasumxy,omegasumyx,omegasumyy;
        //double T = 1.0;
        ostringstream of2;
        of2 << "Cofomega/cofomegaD" << D*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        ofstream writeC(of2.str().c_str());//add when using binary format: ,ios_base::binary
        for(omega=-0.01;omega<2.501;omega+=0.01){
            omegasumxx = {0.0,0.0};
            omegasumxy = {0.0,0.0};
            omegasumyx = {0.0,0.0};
            omegasumyy = {0.0,0.0};
        
            for(int alpha=0;alpha<2*systemsize;alpha++){
                for(int beta=0;beta<2*systemsize;beta++){
                    if(alpha < systemsize){
                        zalpha = 1.0;
                    }
                    else{
                        zalpha = -1.0;
                    }
                    if(beta < systemsize){
                        zbeta = 1.0;
                    }
                    else{
                        zbeta = -1.0;
                    }
                    energydiff = zalpha*energy[alpha]-zbeta*energy[beta];
                    prod1xx = gsl_matrix_get(bogojx,alpha,beta)*gsl_matrix_get(bogojx,(alpha+systemsize)%(2*systemsize),(beta+systemsize)%(2*systemsize));
                    prod2xx = gsl_matrix_get(bogojx,alpha,beta)*gsl_matrix_get(bogojx,beta,alpha);
                    prod1yy = gsl_matrix_get(bogojy,alpha,beta)*gsl_matrix_get(bogojy,(alpha+systemsize)%(2*systemsize),(beta+systemsize)%(2*systemsize));
                    prod2yy = gsl_matrix_get(bogojy,alpha,beta)*gsl_matrix_get(bogojy,beta,alpha);
                    prod1xy = gsl_matrix_get(bogojx,alpha,beta)*gsl_matrix_get(bogojy,(alpha+systemsize)%(2*systemsize),(beta+systemsize)%(2*systemsize));
                    prod2xy = gsl_matrix_get(bogojx,alpha,beta)*gsl_matrix_get(bogojy,beta,alpha);
                    prod1yx = gsl_matrix_get(bogojy,alpha,beta)*gsl_matrix_get(bogojx,(alpha+systemsize)%(2*systemsize),(beta+systemsize)%(2*systemsize));
                    prod2yx = gsl_matrix_get(bogojy,alpha,beta)*gsl_matrix_get(bogojx,beta,alpha);
                    bosealpha = bosefunction(alpha,energy[alpha],T);
                    bosebeta = bosefunction((beta+systemsize)%(2*systemsize),energy[beta],T);                        
                    omegasumxx += bosealpha*bosebeta*(prod1xx+prod2xx)/(omega - energydiff + cutoffepsilon);
                    omegasumxy += bosealpha*bosebeta*(prod1xy+prod2xy)/(omega - energydiff + cutoffepsilon);
                    omegasumyx += bosealpha*bosebeta*(prod1yx+prod2yx)/(omega - energydiff + cutoffepsilon);
                    omegasumyy += bosealpha*bosebeta*(prod1yy+prod2yy)/(omega - energydiff + cutoffepsilon);
                }
            }
            writeC << omega << '\t' << imag(omegasumxx/((complex<double>)(systemsize))) << '\t' << imag(omegasumxy/((complex<double>)(systemsize)))<< '\t' << imag(omegasumyx/((complex<double>)(systemsize)))<< '\t' << imag(omegasumyy/((complex<double>)(systemsize)))<< '\n';
        }
        writeC.close();
        */
        //}
        gsl_matrix_free(bogojx);
        gsl_matrix_free(bogojy);
    
        cerr << 0.01*count << "\r";
    }
    writefield.close();
    }
    gsl_matrix_free (helpmatrixH);
    gsl_matrix_free (helpmatrixX);
    gsl_matrix_free (helpmatrixY);
    delete [] texture;
    delete [] parameter;
    delete [] energy;
    
    return 0;
}

double bosefunction1(int index, double energy, double T){
    return(1.0/(exp(energy/T)-1));
}

double bosefunction(int index, double energy, double T){
    if(index<systemsize){
        return(1.0/(exp(energy/T)-1));
    }
    else if(index>=systemsize){
        return(1.0/(exp(energy/T)-1)+1);
    }
}

void commutator(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C){
    
    gsl_matrix* AB = gsl_matrix_alloc(2*systemsize,2*systemsize);
    gsl_matrix* BA = gsl_matrix_alloc(2*systemsize,2*systemsize);
    
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,B,0,AB);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B,A,0,BA);
    gsl_matrix_sub(AB,BA);
    for(int i=0;i<2*systemsize;i++){
            for(int j=0;j<2*systemsize;j++){
                gsl_matrix_set(C,i,j,gsl_matrix_get(AB,i,j));
            }
    }
    gsl_matrix_free (AB);
    gsl_matrix_free (BA);
}

void energydensity(double *texture, gsl_matrix *edmatrix, int x, int y, double *D){
    double xprefactora,xprefactorb,yprefactora,yprefactorb;
    int site,xnegneighbour,xposneighbour,ynegneighbour,yposneighbour;
    
    gsl_matrix_set_zero(edmatrix);
    site = x*ysize + y;
    xnegneighbour = (x-1)*ysize+y;
    xposneighbour = (x+1)*ysize+y;
    ynegneighbour = x*ysize+y-1;
    yposneighbour = x*ysize+y+1;
    if(y==0){
        ynegneighbour = x*ysize+ysize-1;
    }
    if(y==ysize-1){
        yposneighbour = x*ysize;
    }
    if(x<xsize-1){
        xprefactora = cos(texture[site]-texture[xposneighbour]+M_PI) 
                      +D[site] * sin(texture[site]-texture[xposneighbour]+M_PI);
    }
    if(x>0){
        xprefactorb = cos(texture[xnegneighbour]-texture[site]+M_PI)
                      +D[xnegneighbour] * sin(texture[xnegneighbour]-texture[site]+M_PI);
    }
    yprefactora = cos(texture[site]-texture[yposneighbour]+M_PI);
    yprefactorb = cos(texture[ynegneighbour]-texture[site]+M_PI);
    
    //ATTENTION diagonal terms have epsilon addend!!!
    if(x<xsize-1){
        gsl_matrix_set(edmatrix,site,xposneighbour,
                       -0.25*(xprefactora-1)
                       +gsl_matrix_get(edmatrix,site,xposneighbour));//b_l^\dagger b_{l+1}
        
        gsl_matrix_set(edmatrix,xposneighbour,site,
                       -0.25*(xprefactora-1)
                       +gsl_matrix_get(edmatrix,xposneighbour,site));//b_{l+1}^\dagger b_{l}
        
        gsl_matrix_set(edmatrix,xposneighbour,xposneighbour,
                       0.5*xprefactora+epsilon
                       +gsl_matrix_get(edmatrix,xposneighbour,xposneighbour));//b_{l+1}^\dagger b_{l+1}
        
        gsl_matrix_set(edmatrix,site,site,
                       0.5*xprefactora+epsilon
                       +gsl_matrix_get(edmatrix,site,site));//b_l^\dagger b_{l}
        
        gsl_matrix_set(edmatrix,xposneighbour,site+systemsize,
                       -0.25*(xprefactora+1)
                       +gsl_matrix_get(edmatrix,xposneighbour,site+systemsize));//b_{l+1}^\dagger b_l^\dagger
        
        gsl_matrix_set(edmatrix,site,xposneighbour+systemsize,
                       -0.25*(xprefactora+1)
                       +gsl_matrix_get(edmatrix,site,xposneighbour+systemsize));//b_l^\dagger b_{l+1}^\dagger
    }
    
    if(x>0){
        gsl_matrix_set(edmatrix,site,xnegneighbour,
                       -0.25*(xprefactorb-1)
                       +gsl_matrix_get(edmatrix,site,xnegneighbour));//b_l^\dagger b_{l-1}
        
        gsl_matrix_set(edmatrix,xnegneighbour,site,
                       -0.25*(xprefactorb-1)
                       +gsl_matrix_get(edmatrix,xnegneighbour,site));//b_{l-1}^\dagger b_{l}
        
        gsl_matrix_set(edmatrix,xnegneighbour,xnegneighbour,
                       0.5*xprefactorb+epsilon
                       +gsl_matrix_get(edmatrix,xnegneighbour,xnegneighbour));//b_{l-1}^\dagger b_{l-1}
        
        gsl_matrix_set(edmatrix,site,site,
                       0.5*xprefactorb+epsilon
                       +gsl_matrix_get(edmatrix,site,site));//b_l^\dagger b_{l}
        
        gsl_matrix_set(edmatrix,site,xnegneighbour+systemsize,
                       -0.25*(xprefactorb+1)
                       +gsl_matrix_get(edmatrix,site,xnegneighbour+systemsize));//b_l^\dagger b_{l-1}^\dagger
        
        gsl_matrix_set(edmatrix,xnegneighbour,site+systemsize,
                       -0.25*(xprefactorb+1)
                       +gsl_matrix_get(edmatrix,xnegneighbour,site+systemsize));//b_{l-1}^\dagger b_l^\dagger
    }
    
    gsl_matrix_set(edmatrix,site,yposneighbour,
                       -0.25*(yprefactora-1)
                       +gsl_matrix_get(edmatrix,site,yposneighbour));//b_l^\dagger b_{l+1}
        
    gsl_matrix_set(edmatrix,site,ynegneighbour,
                       -0.25*(yprefactorb-1)
                       +gsl_matrix_get(edmatrix,site,ynegneighbour));//b_l^\dagger b_{l-1}
    
    gsl_matrix_set(edmatrix,yposneighbour,site,
                       -0.25*(yprefactora-1)
                       +gsl_matrix_get(edmatrix,yposneighbour,site));//b_{l+1}^\dagger b_{l}
        
    gsl_matrix_set(edmatrix,ynegneighbour,site,
                       -0.25*(yprefactorb-1)
                       +gsl_matrix_get(edmatrix,ynegneighbour,site));//b_{l-1}^\dagger b_{l}

    gsl_matrix_set(edmatrix,yposneighbour,yposneighbour,
                       0.5*yprefactora+epsilon
                       +gsl_matrix_get(edmatrix,yposneighbour,yposneighbour));//b_{l+1}^\dagger b_{l+1}
        
    gsl_matrix_set(edmatrix,ynegneighbour,ynegneighbour,
                       0.5*yprefactorb+epsilon
                       +gsl_matrix_get(edmatrix,ynegneighbour,ynegneighbour));//b_{l-1}^\dagger b_{l-1}
    
    gsl_matrix_set(edmatrix,site,site,
                       0.5*yprefactora+epsilon
                       +gsl_matrix_get(edmatrix,site,site));//b_l^\dagger b_{l}
        
    gsl_matrix_set(edmatrix,site,site,
                       0.5*yprefactorb+epsilon
                       +gsl_matrix_get(edmatrix,site,site));//b_l^\dagger b_{l}
    
    gsl_matrix_set(edmatrix,site,yposneighbour+systemsize,
                       -0.25*(yprefactora+1)
                       +gsl_matrix_get(edmatrix,site,yposneighbour+systemsize));//b_l^\dagger b_{l+1}^\dagger
        
    gsl_matrix_set(edmatrix,site,ynegneighbour+systemsize,
                       -0.25*(yprefactorb+1)
                       +gsl_matrix_get(edmatrix,site,ynegneighbour+systemsize));//b_l^\dagger b_{l-1}^\dagger
    
    gsl_matrix_set(edmatrix,yposneighbour,site+systemsize,
                       -0.25*(yprefactora+1)
                       +gsl_matrix_get(edmatrix,yposneighbour,site+systemsize));//b_{l+1}^\dagger b_l^\dagger
        
    gsl_matrix_set(edmatrix,ynegneighbour,site+systemsize,
                       -0.25*(yprefactorb+1)
                       +gsl_matrix_get(edmatrix,ynegneighbour,site+systemsize));//b_{l-1}^\dagger b_l^\dagger
    
    //add complex conjugated elements
    for(int i=0;i<systemsize;i++){
        for(int j=0;j<systemsize;j++){
            gsl_matrix_set(edmatrix,i+systemsize,j+systemsize, gsl_matrix_get(edmatrix,i,j));
        }
    }
    for(int i=0;i<systemsize;i++){
        for(int j=0;j<systemsize;j++){
            gsl_matrix_set(edmatrix,j+systemsize,i, gsl_matrix_get(edmatrix,i,j+systemsize));
        }
    }
    gsl_matrix_scale(edmatrix,0.5*S*J);
}
