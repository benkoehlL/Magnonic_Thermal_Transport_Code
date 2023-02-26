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
    of << "Drudeweightoffield/Tdepscat/T" << T*100 << "cJ/drudeweightxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
    ofstream writefield(of.str().c_str());//add when using binary format: ,ios_base::binary
    
    for(int count=1;count<=Dmax*100/*+1*/;count++){
        D = count*0.01;
        for(int i=0;i<xsize;i++){
            for(int j=0;j<ysize;j++){
                if(i>=freefieldsize && i<(freefieldsize+fieldsize)){
                    if(i==freefieldsize && j==0){
                        parameter[i*ysize+j] = -D;
                    }
                    else if(i!=freefieldsize && j==0){
                        parameter[i*ysize+j] = (-1)*parameter[(i-1)*ysize+j];
                    }
                    else{
                        parameter[i*ysize+j] = (-1)*parameter[i*ysize+j-1];
                    }
                }
                else{
                    parameter[i*ysize+j] = 0.0;
                }
            }
        }
        ostringstream fin1;
        fin1 << "Texture/textureD" << D*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        ifstream read1(fin1.str().c_str(),ios_base::binary);
        
        //read texture
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                read1.seekg((x*ysize+y)*sizeof(double));
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
        //if(count==(int)(Dmax/Dincrement)/1000){
        /*if(count==1){
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
                int xbefore,xafter,ybefore,yafter;
                double xprefactora,xprefactorb,yprefactora,yprefactorb;

                if(x==0){
                    xbefore = xsize-1;
                    xafter = x+1;
                }
                else{
                    xbefore = x-1;
                    xafter = (x+1)%xsize;
                }
                    
                if(y==0){
                    ybefore = ysize-1;
                    yafter = y+1;
                }
                else{
                    ybefore = y-1;
                    yafter = (y+1)%ysize;
                }
                
                xprefactora = cos(texture[x*ysize+y]-texture[xafter*ysize+y])- parameter[x*ysize+y]*sin(texture[x*ysize+y]-texture[xafter*ysize+y]);
                xprefactorb = cos(texture[x*ysize+y]-texture[xbefore*ysize+y])-parameter[x*ysize+y]*sin(texture[x*ysize+y]-texture[xbefore*ysize+y]);
                yprefactora = cos(texture[x*ysize+y]-texture[x*ysize+yafter])- parameter[x*ysize+y]*sin(texture[x*ysize+y]-texture[x*ysize+yafter]);
                yprefactorb = cos(texture[x*ysize+y]-texture[x*ysize+ybefore])-parameter[x*ysize+y]*sin(texture[x*ysize+y]-texture[x*ysize+ybefore]);
            
                //\delta=-\Delta = \pm e_x
                gsl_matrix_set(jx,x*ysize+y,xafter*ysize+y , xprefactorb*(xprefactora-1)+gsl_matrix_get(jx,x*ysize+y,xafter*ysize+y));//b_{l}^\dagger b_{l+e_x}
                gsl_matrix_set(jx,x*ysize+y,xbefore*ysize+y,-xprefactora*(xprefactorb-1)+gsl_matrix_get(jx,x*ysize+y,xbefore*ysize+y));//b_{l}^\dagger b_{l-e_x}
                
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+y,-xprefactorb*(xprefactora-1)+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+y));//b_{l+e_x}^\dagger b_{l}
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+y,xprefactora*(xprefactorb-1)+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+y));//b_{l-e_x}^\dagger b_{l}
                
                gsl_matrix_set(jx,x*ysize+y,xafter*ysize+y+systemsize,2*xprefactorb*(xprefactora+1)+gsl_matrix_get(jx,x*ysize+y,xafter*ysize+y+systemsize));//b_{l}^\dagger b_{l+e_x}^\dagger
                gsl_matrix_set(jx,x*ysize+y,xbefore*ysize+y+systemsize,-2*xprefactora*(xprefactorb+1)+gsl_matrix_get(jx,x*ysize+y,xbefore*ysize+y+systemsize));//b_{l}^\dagger b_{l-e_x}^\dagger
                
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+y+systemsize,2*xprefactorb*(xprefactora+1)+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+y+systemsize));//b_{l+e_x}^\dagger b_{l}\dagger
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+y+systemsize,-2*xprefactora*(xprefactorb+1)+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+y+systemsize));//b_{l-e_x}^\dagger b_{l}\dagger
                
                gsl_matrix_set(jx,xbefore*ysize+y,xafter*ysize+y,0.5*((xprefactorb-1)*(xprefactora-1)-(xprefactorb+1)*(xprefactora+1))+gsl_matrix_get(jx,xbefore*ysize+y,xafter*ysize+y));//b_{l-e_x}^\dagger b_{l+e_x}
                gsl_matrix_set(jx,xafter*ysize+y,xbefore*ysize+y,-0.5*((xprefactora-1)*(xprefactorb-1)-(xprefactora+1)*(xprefactorb+1))+gsl_matrix_get(jx,xafter*ysize+y,xbefore*ysize+y));//b_{l+e_x}^\dagger b_{l-e_x}
                    
                gsl_matrix_set(jx,xafter*ysize+y,xbefore*ysize+y+systemsize,(xprefactora-1)*(xprefactorb+1)+gsl_matrix_get(jx,xafter*ysize+y,xbefore*ysize+y+systemsize));//b_{l+e_x}^\dagger b_{l-e_x}^\dagger
                gsl_matrix_set(jx,xbefore*ysize+y,xafter*ysize+y+systemsize,-(xprefactorb-1)*(xprefactora+1)+gsl_matrix_get(jx,xbefore*ysize+y,xafter*ysize+y+systemsize));//b_{l-e_x}^\dagger b_{l+e_x}^\dagger
                
                gsl_matrix_set(jx,xbefore*ysize+y,xafter*ysize+y+systemsize,(xprefactora-1)*(xprefactorb+1)+gsl_matrix_get(jx,xbefore*ysize+y,xafter*ysize+y+systemsize));//b_{l-e_x}^\dagger b_{l+e_x}^\dagger
                gsl_matrix_set(jx,xafter*ysize+y,xbefore*ysize+y+systemsize,-(xprefactorb-1)*(xprefactora+1)+gsl_matrix_get(jx,xafter*ysize+y,xbefore*ysize+y+systemsize));//b_{l+e_x}^\dagger b_{l-e_x}^\dagger
                
                //\delta=-\Delta = \pm e_y
                gsl_matrix_set(jy,x*ysize+y,x*ysize+yafter,yprefactorb*(yprefactora-1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+yafter));//b_{l}^\dagger b_{l+e_y}
                gsl_matrix_set(jy,x*ysize+y,x*ysize+ybefore,-yprefactora*(yprefactorb-1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+ybefore));//b_{l}^\dagger b_{l-e_y}
                
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+y,-yprefactorb*(yprefactora-1)+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+y));//b_{l+e_y}^\dagger b_{l}
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+y,yprefactora*(yprefactorb-1)+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+y));//b_{l-e_y}^\dagger b_{l}
                
                gsl_matrix_set(jy,x*ysize+y,x*ysize+yafter+systemsize,2*yprefactorb*(yprefactora+1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+yafter+systemsize));//b_{l}^\dagger b_{l+e_y}^\dagger
                gsl_matrix_set(jy,x*ysize+y,x*ysize+ybefore+systemsize,-2*yprefactora*(yprefactorb+1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+ybefore+systemsize));//b_{l}^\dagger b_{l-e_y}^\dagger
                
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+y+systemsize,2*yprefactorb*(yprefactora+1)+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+y+systemsize));//b_{l+e_y}^\dagger b_{l}\dagger
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+y+systemsize,-2*yprefactora*(yprefactorb+1)+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+y+systemsize));//b_{l-e_y}^\dagger b_{l}\dagger
                
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+yafter,0.5*((yprefactorb-1)*(yprefactora-1)-(yprefactorb+1)*(yprefactora+1))+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+yafter));//b_{l-e_y}^\dagger b_{l+e_y}
                
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+ybefore,-0.5*((yprefactora-1)*(yprefactorb-1)-(yprefactora+1)*(yprefactorb+1))+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+ybefore));//b_{l+e_y}^\dagger b_{l-e_y}
                    
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+ybefore+systemsize,(yprefactora-1)*(yprefactorb+1)+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+ybefore+systemsize));//b_{l+e_y}^\dagger b_{l-e_y}^\dagger
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+yafter+systemsize,-(yprefactorb-1)*(yprefactora+1)+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+yafter+systemsize));//b_{l-e_y}^\dagger b_{l+e_y}^\dagger
                
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+yafter+systemsize,(yprefactora-1)*(yprefactorb+1)+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+yafter+systemsize));//b_{l-e_y}^\dagger b_{l+e_y}^\dagger
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+ybefore+systemsize,-(yprefactorb-1)*(yprefactora+1)+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+ybefore+systemsize));//b_{l+e_y}^\dagger b_{l-e_y}^\dagger
                
                
                //\delta not (-\Delta)
                //\Delta=\pm e_x
                
                gsl_matrix_set(jx,x*ysize+y,xafter*ysize+y,0.5*yprefactora*(xprefactora-1)+gsl_matrix_get(jx,x*ysize+y,xafter*ysize+y));//b_{l}^\dagger b_{l+e_x}  \delta = +e_y 
                gsl_matrix_set(jy,x*ysize+y,xafter*ysize+y,-0.5*yprefactora*(xprefactora-1)+gsl_matrix_get(jy,x*ysize+y,xafter*ysize+y));
                gsl_matrix_set(jx,x*ysize+y,xafter*ysize+y,0.5*yprefactorb*(xprefactora-1)+gsl_matrix_get(jx,x*ysize+y,xafter*ysize+y));//b_{l}^\dagger b_{l+e_x} \delta = -ey
                gsl_matrix_set(jy,x*ysize+y,xafter*ysize+y,0.5*yprefactorb*(xprefactora-1)+gsl_matrix_get(jy,x*ysize+y,xafter*ysize+y));
                gsl_matrix_set(jx,x*ysize+y,xbefore*ysize+y,-0.5*yprefactora*(xprefactorb-1)+gsl_matrix_get(jx,x*ysize+y,xbefore*ysize+y));//b_{l}^\dagger b_{l-e_x} \delta = +e_y
                gsl_matrix_set(jy,x*ysize+y,xbefore*ysize+y,-0.5*yprefactora*(xprefactorb-1)+gsl_matrix_get(jy,x*ysize+y,xbefore*ysize+y));                        
                gsl_matrix_set(jx,x*ysize+y,xbefore*ysize+y,-0.5*yprefactorb*(xprefactorb-1)+gsl_matrix_get(jx,x*ysize+y,xbefore*ysize+y));//b_{l}^\dagger b_{l-e_x} \delta = -e_y
                gsl_matrix_set(jy,x*ysize+y,xbefore*ysize+y,0.5*yprefactorb*(xprefactorb-1)+gsl_matrix_get(jy,x*ysize+y,xbefore*ysize+y));
                        
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+y,-0.5*yprefactora*(xprefactora-1)+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+y));//b_{l+e_x}^\dagger b_{l}  \delta = +e_y 
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+y,0.5*yprefactora*(xprefactora-1)+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+y));
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+y,-0.5*yprefactorb*(xprefactora-1)+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+y));//b_{l+e_x}^\dagger b_{l} \delta = -ey
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+y,-0.5*yprefactorb*(xprefactora-1)+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+y));
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+y,0.5*yprefactora*(xprefactorb-1)+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+y));//b_{l-e_x}^\dagger b_{l} \delta = +e_y
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+y,0.5*yprefactora*(xprefactorb-1)+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+y));                        
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+y,0.5*yprefactorb*(xprefactorb-1)+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+y));//b_{l-e_x}^\dagger b_{l} \delta = -e_y
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+y,-0.5*yprefactorb*(xprefactorb-1)+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+y));
                        
                gsl_matrix_set(jx,x*ysize+y,xafter*ysize+y+systemsize,yprefactora*(xprefactora+1)+gsl_matrix_get(jx,x*ysize+y,xafter*ysize+y+systemsize));//b_{l}^\dagger b_{l+e_x}^\dagger  \delta = +e_y 
                gsl_matrix_set(jy,x*ysize+y,xafter*ysize+y+systemsize,-yprefactora*(xprefactora+1)+gsl_matrix_get(jy,x*ysize+y,xafter*ysize+y+systemsize));
                gsl_matrix_set(jx,x*ysize+y,xafter*ysize+y+systemsize,yprefactorb*(xprefactora+1)+gsl_matrix_get(jx,x*ysize+y,xafter*ysize+y+systemsize));//b_{l}^\dagger b_{l+e_x}^\dagger \delta = -ey
                gsl_matrix_set(jy,x*ysize+y,xafter*ysize+y+systemsize,yprefactorb*(xprefactora+1)+gsl_matrix_get(jy,x*ysize+y,xafter*ysize+y+systemsize));
                gsl_matrix_set(jx,x*ysize+y,xbefore*ysize+y+systemsize,-yprefactora*(xprefactorb+1)+gsl_matrix_get(jx,x*ysize+y,xbefore*ysize+y+systemsize));//b_{l}^\dagger b_{l-e_x}^\dagger \delta = +e_y
                gsl_matrix_set(jy,x*ysize+y,xbefore*ysize+y+systemsize,-yprefactora*(xprefactorb+1)+gsl_matrix_get(jy,x*ysize+y,xbefore*ysize+y+systemsize));                        
                gsl_matrix_set(jx,x*ysize+y,xbefore*ysize+y+systemsize,-yprefactorb*(xprefactorb+1)+gsl_matrix_get(jx,x*ysize+y,xbefore*ysize+y+systemsize));//b_{l}^\dagger b_{l-e_x}^\dagger \delta = -e_y
                gsl_matrix_set(jy,x*ysize+y,xbefore*ysize+y+systemsize,yprefactorb*(xprefactorb+1)+gsl_matrix_get(jy,x*ysize+y,xbefore*ysize+y+systemsize));
                        
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+y+systemsize,yprefactora*(xprefactora+1)+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+y+systemsize));//b_{l+e_x}^\dagger b_{l}^\dagger  \delta = +e_y 
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+y+systemsize,-yprefactora*(xprefactora+1)+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+y+systemsize));
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+y+systemsize,yprefactorb*(xprefactora+1)+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+y+systemsize));//b_{l+e_x}^\dagger b_{l}^\dagger \delta = -ey
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+y+systemsize,yprefactorb*(xprefactora+1)+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+y+systemsize));
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+y+systemsize,-yprefactora*(xprefactorb+1)+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+y+systemsize));//b_{l-e_x}^\dagger b_{l}^\dagger \delta = +e_y
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+y+systemsize,-yprefactora*(xprefactorb+1)+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+y+systemsize));                        
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+y+systemsize,-yprefactorb*(xprefactorb+1)+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+y+systemsize));//b_{l-e_x}^\dagger b_{l}^\dagger \delta = -e_y
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+y+systemsize,yprefactorb*(xprefactorb+1)+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+y+systemsize));
                
                gsl_matrix_set(jx,x*ysize+yafter,xafter*ysize+y,0.25*((yprefactora-1)*(xprefactora-1)-(yprefactora+1)*(xprefactora+1))+gsl_matrix_get(jx,x*ysize+yafter,xafter*ysize+y));//b_{l+e_y}^\dagger b_{l+e_x}
                gsl_matrix_set(jy,x*ysize+yafter,xafter*ysize+y,-0.25*((yprefactora-1)*(xprefactora-1)-(yprefactora+1)*(xprefactora+1))+gsl_matrix_get(jy,x*ysize+yafter,xafter*ysize+y));
                gsl_matrix_set(jx,x*ysize+ybefore,xafter*ysize+y,0.25*((yprefactorb-1)*(xprefactora-1)-(yprefactorb+1)*(xprefactora+1))+gsl_matrix_get(jx,x*ysize+ybefore,xafter*ysize+y));//b_{l-e_y}^\dagger b_{l+e_x} 
                gsl_matrix_set(jy,x*ysize+ybefore,xafter*ysize+y,0.25*((yprefactorb-1)*(xprefactora-1)-(yprefactorb+1)*(xprefactora+1))+gsl_matrix_get(jy,x*ysize+ybefore,xafter*ysize+y));
                gsl_matrix_set(jx,x*ysize+yafter,xbefore*ysize+y,-0.25*((yprefactora-1)*(xprefactorb-1)-(yprefactora+1)*(xprefactorb+1))+gsl_matrix_get(jx,x*ysize+yafter,xbefore*ysize+y));//b_{l+e_y}^\dagger b_{l-e_x} 
                gsl_matrix_set(jy,x*ysize+yafter,xbefore*ysize+y,-0.25*((yprefactora-1)*(xprefactorb-1)-(yprefactora+1)*(xprefactorb+1))+gsl_matrix_get(jy,x*ysize+yafter,xbefore*ysize+y));
                gsl_matrix_set(jx,x*ysize+ybefore,xbefore*ysize+y,-0.25*((yprefactorb-1)*(xprefactorb-1)-(yprefactorb+1)*(xprefactorb+1))+gsl_matrix_get(jx,x*ysize+ybefore,xbefore*ysize+y));//b_{l-e_y}^\dagger b_{l-e_x}
                gsl_matrix_set(jy,x*ysize+ybefore,xbefore*ysize+y,0.25*((yprefactorb-1)*(xprefactorb-1)-(yprefactorb+1)*(xprefactorb+1))+gsl_matrix_get(jy,x*ysize+ybefore,xbefore*ysize+y));
                
                gsl_matrix_set(jx,x*ysize+yafter,xafter*ysize+y+systemsize,-0.5*((yprefactora-1)*(xprefactora+1))+gsl_matrix_get(jx,x*ysize+yafter,xafter*ysize+y+systemsize));//b_{l+e_y}^\dagger b_{l+e_x}^\dagger
                gsl_matrix_set(jy,x*ysize+yafter,xafter*ysize+y+systemsize,0.5*((yprefactora-1)*(xprefactora+1))+gsl_matrix_get(jy,x*ysize+yafter,xafter*ysize+y+systemsize));
                gsl_matrix_set(jx,x*ysize+ybefore,xafter*ysize+y+systemsize,-0.5*((yprefactorb-1)*(xprefactora+1))+gsl_matrix_get(jx,x*ysize+ybefore,xafter*ysize+y+systemsize));//b_{l-e_y}^\dagger b_{l+e_x}^\dagger 
                gsl_matrix_set(jy,x*ysize+ybefore,xafter*ysize+y+systemsize,-0.5*((yprefactorb-1)*(xprefactora+1))+gsl_matrix_get(jy,x*ysize+ybefore,xafter*ysize+y+systemsize));
                gsl_matrix_set(jx,x*ysize+yafter,xbefore*ysize+y+systemsize,0.5*((yprefactora-1)*(xprefactorb+1))+gsl_matrix_get(jx,x*ysize+yafter,xbefore*ysize+y+systemsize));//b_{l+e_y}^\dagger b_{l-e_x}^\dagger
                gsl_matrix_set(jy,x*ysize+yafter,xbefore*ysize+y+systemsize,0.5*((yprefactora-1)*(xprefactorb+1))+gsl_matrix_get(jy,x*ysize+yafter,xbefore*ysize+y+systemsize));
                gsl_matrix_set(jx,x*ysize+ybefore,xbefore*ysize+y+systemsize,0.5*((yprefactorb-1)*(xprefactorb+1))+gsl_matrix_get(jx,x*ysize+ybefore,xbefore*ysize+y+systemsize));//b_{l-e_y}^\dagger b_{l-e_x}^\dagger
                gsl_matrix_set(jy,x*ysize+ybefore,xbefore*ysize+y+systemsize,-0.5*((yprefactorb-1)*(xprefactorb+1))+gsl_matrix_get(jy,x*ysize+ybefore,xbefore*ysize+y+systemsize));
                        
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+yafter+systemsize,-0.5*((yprefactora-1)*(xprefactora+1))+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+yafter+systemsize));//b_{l+e_x}^\dagger b_{l+e_y}^\dagger
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+yafter+systemsize,0.5*((yprefactora-1)*(xprefactora+1))+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+yafter+systemsize));
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+ybefore+systemsize,-0.5*((yprefactorb-1)*(xprefactora+1))+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+ybefore+systemsize));//b_{l+e_x}^\dagger b_{l-e_y}^\dagger 
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+ybefore+systemsize,-0.5*((yprefactorb-1)*(xprefactora+1))+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+ybefore+systemsize));
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+yafter+systemsize,0.5*((yprefactora-1)*(xprefactorb+1))+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+yafter+systemsize));//b_{l-e_x}^\dagger b_{l+e_y}^\dagger
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+yafter+systemsize,0.5*((yprefactora-1)*(xprefactorb+1))+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+yafter+systemsize));
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+ybefore+systemsize,0.5*((yprefactorb-1)*(xprefactorb+1))+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+ybefore+systemsize));//b_{l-e_x}^\dagger b_{l-e_y}^\dagger
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+ybefore+systemsize,-0.5*((yprefactorb-1)*(xprefactorb+1))+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+ybefore+systemsize));
                
                //\Delta=\pm e_y
                
                gsl_matrix_set(jy,x*ysize+y,x*ysize+yafter,0.5*xprefactora*(yprefactora-1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+yafter));//b_{l}^\dagger b_{l+e_y}  \delta = +e_x
                gsl_matrix_set(jx,x*ysize+y,x*ysize+yafter,-0.5*xprefactora*(yprefactora-1)+gsl_matrix_get(jx,x*ysize+y,x*ysize+yafter));
                gsl_matrix_set(jy,x*ysize+y,x*ysize+yafter,0.5*xprefactorb*(yprefactora-1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+yafter));//b_{l}^\dagger b_{l+e_y} \delta = -e_x
                gsl_matrix_set(jx,x*ysize+y,x*ysize+yafter,0.5*xprefactorb*(yprefactora-1)+gsl_matrix_get(jx,x*ysize+y,x*ysize+yafter));
                gsl_matrix_set(jy,x*ysize+y,x*ysize+ybefore,-0.5*xprefactora*(yprefactorb-1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+ybefore));//b_{l}^\dagger b_{l-e_y} \delta = +e_x
                gsl_matrix_set(jx,x*ysize+y,x*ysize+ybefore,-0.5*xprefactora*(yprefactorb-1)+gsl_matrix_get(jx,x*ysize+y,x*ysize+ybefore));                        
                gsl_matrix_set(jy,x*ysize+y,x*ysize+ybefore,-0.5*xprefactorb*(yprefactorb-1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+ybefore));//b_{l}^\dagger b_{l-e_y} \delta = -e_x
                gsl_matrix_set(jx,x*ysize+y,x*ysize+ybefore,0.5*xprefactorb*(yprefactorb-1)+gsl_matrix_get(jx,x*ysize+y,x*ysize+ybefore));
                        
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+y,-0.5*xprefactora*(yprefactora-1)+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+y));//b_{l+e_y}^\dagger b_{l}  \delta = +e_x
                gsl_matrix_set(jx,x*ysize+yafter,x*ysize+y,0.5*xprefactora*(yprefactora-1)+gsl_matrix_get(jx,x*ysize+yafter,x*ysize+y));
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+y,-0.5*xprefactorb*(yprefactora-1)+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+y));//b_{l+e_y}^\dagger b_{l} \delta = -e_x
                gsl_matrix_set(jx,x*ysize+yafter,x*ysize+y,-0.5*xprefactorb*(yprefactora-1)+gsl_matrix_get(jx,x*ysize+yafter,x*ysize+y));
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+y,0.5*xprefactora*(yprefactorb-1)+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+y));//b_{l-e_y}^\dagger b_{l} \delta = +e_x
                gsl_matrix_set(jx,x*ysize+ybefore,x*ysize+y,0.5*xprefactora*(yprefactorb-1)+gsl_matrix_get(jx,x*ysize+ybefore,x*ysize+y));                        
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+y,0.5*xprefactorb*(yprefactorb-1)+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+y));//b_{l-e_y}^\dagger b_{l} \delta = -e_x
                gsl_matrix_set(jx,x*ysize+ybefore,x*ysize+y,-0.5*xprefactorb*(yprefactorb-1)+gsl_matrix_get(jx,x*ysize+ybefore,x*ysize+y));
                        
                gsl_matrix_set(jy,x*ysize+y,x*ysize+yafter+systemsize,xprefactora*(yprefactora+1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+yafter+systemsize));//b_{l}^\dagger b_{l+e_y}^\dagger  \delta = +e_x 
                gsl_matrix_set(jx,x*ysize+y,x*ysize+yafter+systemsize,-xprefactora*(yprefactora+1)+gsl_matrix_get(jx,x*ysize+y,x*ysize+yafter+systemsize));
                gsl_matrix_set(jy,x*ysize+y,x*ysize+yafter+systemsize,xprefactorb*(yprefactora+1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+yafter+systemsize));//b_{l}^\dagger b_{l+e_y}^\dagger \delta = -e_x
                gsl_matrix_set(jx,x*ysize+y,x*ysize+yafter+systemsize,xprefactorb*(yprefactora+1)+gsl_matrix_get(jx,x*ysize+y,x*ysize+yafter+systemsize));
                gsl_matrix_set(jy,x*ysize+y,x*ysize+ybefore+systemsize,-xprefactora*(yprefactorb+1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+ybefore+systemsize));//b_{l}^\dagger b_{l-e_y}^\dagger \delta = +e_x
                gsl_matrix_set(jx,x*ysize+y,x*ysize+ybefore+systemsize,-xprefactora*(yprefactorb+1)+gsl_matrix_get(jx,x*ysize+y,x*ysize+ybefore+systemsize));                        
                gsl_matrix_set(jy,x*ysize+y,x*ysize+ybefore+systemsize,-xprefactorb*(yprefactorb+1)+gsl_matrix_get(jy,x*ysize+y,x*ysize+ybefore+systemsize));//b_{l}^\dagger b_{l-e_y}^\dagger \delta = -e_x
                gsl_matrix_set(jx,x*ysize+y,x*ysize+ybefore+systemsize,xprefactorb*(yprefactorb+1)+gsl_matrix_get(jx,x*ysize+y,x*ysize+ybefore+systemsize));
                        
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+y+systemsize,xprefactora*(yprefactora+1)+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+y+systemsize));//b_{l+e_y}^\dagger b_{l}^\dagger  \delta = +e_x
                gsl_matrix_set(jx,x*ysize+yafter,x*ysize+y+systemsize,-xprefactora*(yprefactora+1)+gsl_matrix_get(jx,x*ysize+yafter,x*ysize+y+systemsize));
                gsl_matrix_set(jy,x*ysize+yafter,x*ysize+y+systemsize,xprefactorb*(yprefactora+1)+gsl_matrix_get(jy,x*ysize+yafter,x*ysize+y+systemsize));//b_{l+e_y}^\dagger b_{l}^\dagger \delta = -e_x
                gsl_matrix_set(jx,x*ysize+yafter,x*ysize+y+systemsize,xprefactorb*(yprefactora+1)+gsl_matrix_get(jx,x*ysize+yafter,x*ysize+y+systemsize));
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+y+systemsize,-xprefactora*(yprefactorb+1)+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+y+systemsize));//b_{l-e_y}^\dagger b_{l}^\dagger \delta = +e_x
                gsl_matrix_set(jx,x*ysize+ybefore,x*ysize+y+systemsize,-xprefactora*(yprefactorb+1)+gsl_matrix_get(jx,x*ysize+ybefore,x*ysize+y+systemsize));                        
                gsl_matrix_set(jy,x*ysize+ybefore,x*ysize+y+systemsize,-xprefactorb*(yprefactorb+1)+gsl_matrix_get(jy,x*ysize+ybefore,x*ysize+y+systemsize));//b_{l-e_y}^\dagger b_{l}^\dagger \delta = -e_x
                gsl_matrix_set(jx,x*ysize+ybefore,x*ysize+y+systemsize,xprefactorb*(yprefactorb+1)+gsl_matrix_get(jx,x*ysize+ybefore,x*ysize+y+systemsize));
                
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+yafter,0.25*((xprefactora-1)*(yprefactora-1)-(xprefactora+1)*(yprefactora+1))+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+yafter));//b_{l+e_x}^\dagger b_{l+e_y}
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+yafter,-0.25*((xprefactora-1)*(yprefactora-1)-(xprefactora+1)*(yprefactora+1))+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+yafter));
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+yafter,0.25*((xprefactorb-1)*(yprefactora-1)-(xprefactorb+1)*(yprefactora+1))+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+yafter));//b_{l-e_x}^\dagger b_{l+e_y} 
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+yafter,0.25*((xprefactorb-1)*(yprefactora-1)-(xprefactorb+1)*(yprefactora+1))+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+yafter));
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+ybefore,-0.25*((xprefactora-1)*(yprefactorb-1)-(xprefactora+1)*(yprefactorb+1))+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+ybefore));//b_{l+e_x}^\dagger b_{l-e_y} 
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+ybefore,-0.25*((xprefactora-1)*(yprefactorb-1)-(xprefactora+1)*(yprefactorb+1))+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+ybefore));
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+ybefore,-0.25*((xprefactorb-1)*(yprefactorb-1)-(xprefactorb+1)*(yprefactorb+1))+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+ybefore));//b_{l-e_x}^\dagger b_{l-e_y}
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+ybefore,0.25*((xprefactorb-1)*(yprefactorb-1)-(xprefactorb+1)*(yprefactorb+1))+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+ybefore));
                        
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+yafter+systemsize,-0.5*((xprefactora-1)*(yprefactora+1))+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+yafter+systemsize));//b_{l+e_x}^\dagger b_{l+e_y}^\dagger
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+yafter+systemsize,0.5*((xprefactora-1)*(yprefactora+1))+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+yafter+systemsize));
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+yafter+systemsize,-0.5*((xprefactorb-1)*(yprefactora+1))+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+yafter+systemsize));//b_{l-e_x}^\dagger b_{l+e_y}^\dagger 
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+yafter+systemsize,-0.5*((xprefactorb-1)*(yprefactora+1))+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+yafter+systemsize));
                gsl_matrix_set(jy,xafter*ysize+y,x*ysize+ybefore+systemsize,0.5*((xprefactora-1)*(yprefactorb+1))+gsl_matrix_get(jy,xafter*ysize+y,x*ysize+ybefore+systemsize));//b_{l+e_x}^\dagger b_{l-e_y}^\dagger
                gsl_matrix_set(jx,xafter*ysize+y,x*ysize+ybefore+systemsize,0.5*((xprefactora-1)*(yprefactorb+1))+gsl_matrix_get(jx,xafter*ysize+y,x*ysize+ybefore+systemsize));
                gsl_matrix_set(jy,xbefore*ysize+y,x*ysize+ybefore+systemsize,0.5*((xprefactorb-1)*(yprefactorb+1))+gsl_matrix_get(jy,xbefore*ysize+y,x*ysize+ybefore+systemsize));//b_{l-e_x}^\dagger b_{l-e_y}^\dagger
                gsl_matrix_set(jx,xbefore*ysize+y,x*ysize+ybefore+systemsize,-0.5*((xprefactorb-1)*(yprefactorb+1))+gsl_matrix_get(jx,xbefore*ysize+y,x*ysize+ybefore+systemsize));
                        
                gsl_matrix_set(jy,x*ysize+yafter,xafter*ysize+y+systemsize,-0.5*((xprefactora-1)*(yprefactora+1))+gsl_matrix_get(jy,x*ysize+yafter,xafter*ysize+y+systemsize));//b_{l+e_y}^\dagger b_{l+e_x}^\dagger
                gsl_matrix_set(jx,x*ysize+yafter,xafter*ysize+y+systemsize,0.5*((xprefactora-1)*(yprefactora+1))+gsl_matrix_get(jx,x*ysize+yafter,xafter*ysize+y+systemsize));
                gsl_matrix_set(jy,x*ysize+yafter,xbefore*ysize+y+systemsize,-0.5*((xprefactorb-1)*(yprefactora+1))+gsl_matrix_get(jy,x*ysize+yafter,xbefore*ysize+y+systemsize));//b_{l+e_y}^\dagger b_{l-e_x}^\dagger 
                gsl_matrix_set(jx,x*ysize+yafter,xbefore*ysize+y+systemsize,-0.5*((xprefactorb-1)*(yprefactora+1))+gsl_matrix_get(jx,x*ysize+yafter,xbefore*ysize+y+systemsize));
                gsl_matrix_set(jy,x*ysize+ybefore,xafter*ysize+y+systemsize,0.5*((xprefactora-1)*(yprefactorb+1))+gsl_matrix_get(jy,x*ysize+ybefore,xafter*ysize+y+systemsize));//b_{l-e_y}^\dagger b_{l+e_x}^\dagger
                gsl_matrix_set(jx,x*ysize+ybefore,xafter*ysize+y+systemsize,0.5*((xprefactora-1)*(yprefactorb+1))+gsl_matrix_get(jx,x*ysize+ybefore,xafter*ysize+y+systemsize));
                gsl_matrix_set(jy,x*ysize+ybefore,xbefore*ysize+y+systemsize,0.5*((xprefactorb-1)*(yprefactorb+1))+gsl_matrix_get(jy,x*ysize+ybefore,xbefore*ysize+y+systemsize));//b_{l-e_y}^\dagger b_{l-e_x}^\dagger
                gsl_matrix_set(jx,x*ysize+ybefore,xbefore*ysize+y+systemsize,-0.5*((xprefactorb-1)*(yprefactorb+1))+gsl_matrix_get(jx,x*ysize+ybefore,xbefore*ysize+y+systemsize));

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
        if(count==(int)(Dmax/Dincrement)/1000){
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
        
        
        //Fourier-transform current operators
        /*
        gsl_matrix_complex* zjx = gsl_matrix_complex_alloc(2*systemsize,2*systemsize);
        gsl_matrix_complex* zjy = gsl_matrix_complex_alloc(2*systemsize,2*systemsize);
        for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                gsl_matrix_complex_set(zjx,x,y,gsl_complex_rect(gsl_matrix_get(jx,x,y),0.0));
                gsl_matrix_complex_set(zjy,x,y,gsl_complex_rect(gsl_matrix_get(jy,x,y),0.0));
            }
        }
        gsl_matrix_complex* fourier = gsl_matrix_complex_alloc(2*systemsize,2*systemsize);
        double fourierexponent;
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                for(int kx=0;kx<xsize;kx++){
                    for(int ky=0;ky<ysize;ky++){
                        fourierexponent = 2*M_PI*(kx*x*1.0/xsize+ky*y*1.0/ysize);
                        gsl_matrix_complex_set(fourier,x*ysize+y,kx*ysize+ky,gsl_complex_polar(sqrt(1.0/(xsize*ysize)),fourierexponent));
                        gsl_matrix_complex_set(fourier,x*ysize+y+systemsize,kx*ysize+ky,gsl_complex_rect(0.0,0.0));
                        gsl_matrix_complex_set(fourier,x*ysize+y,kx*ysize+ky+systemsize,gsl_complex_rect(0.0,0.0));
                        gsl_matrix_complex_set(fourier,x*ysize+y+systemsize,kx*ysize+ky+systemsize,gsl_complex_polar(sqrt(1.0/(xsize*ysize)),fourierexponent));
                    }
                }
            }
        }
        
        gsl_matrix_complex* helpzmatrix = gsl_matrix_complex_alloc(2*systemsize,2*systemsize);
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1.0,0.0), fourier, zjx, gsl_complex_rect(0.0,0.0), helpzmatrix);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.0,0.0), helpzmatrix, fourier, gsl_complex_rect(0.0,0.0), zjx);
        for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                gsl_matrix_set(jx,x,y,GSL_IMAG(gsl_matrix_complex_get(zjx,x,y)));
            }
        }

        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1.0,0.0), fourier, zjy, gsl_complex_rect(0.0,0.0), helpzmatrix);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.0,0.0), helpzmatrix, fourier, gsl_complex_rect(0.0,0.0), zjy);
        for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                gsl_matrix_set(jy,x,y,GSL_IMAG(gsl_matrix_complex_get(zjy,x,y)));
            }
        }
    
        gsl_matrix_complex_free(zjx);
        gsl_matrix_complex_free(zjy);
        gsl_matrix_complex_free(helpzmatrix);
        gsl_matrix_complex_free(fourier);
        
        /*if(count==(int)(Dmax/Dincrement)/1000){
        //if(count==1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(jx,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //read Bogoliubov matrix
        ostringstream fin2;
        fin2 << "Bogoliubov/eigenvectorD" << D*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
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
    
        cerr << count/((Dmax/Dincrement)/1000) << "\r";
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
    int xbefore,xafter,ybefore,yafter;
    gsl_matrix_set_zero(edmatrix);
    
    if(x==0){
        xbefore = xsize-1;
        xafter = x+1;
    }
    else{
        xbefore = x-1;
        xafter = (x+1)%xsize;
    }
                    
    if(y==0){
        ybefore = ysize-1;
        yafter = y+1;
    }
    else{
        ybefore = y-1;
        yafter = (y+1)%ysize;
    }
    
    xprefactora = cos(texture[x*ysize+y]-texture[xafter*ysize+y])- D[x*ysize+y]*sin(texture[x*ysize+y]-texture[xafter*ysize+y]);
    xprefactorb = cos(texture[x*ysize+y]-texture[xbefore*ysize+y])-D[x*ysize+y]*sin(texture[x*ysize+y]-texture[xbefore*ysize+y]);
    yprefactora = cos(texture[x*ysize+y]-texture[x*ysize+yafter])- D[x*ysize+y]*sin(texture[x*ysize+y]-texture[x*ysize+yafter]);
    yprefactorb = cos(texture[x*ysize+y]-texture[x*ysize+ybefore])-D[x*ysize+y]*sin(texture[x*ysize+y]-texture[x*ysize+ybefore]);
    
    //ATTENTION diagonal terms have epsilon addent!!!
    gsl_matrix_set(edmatrix,x*ysize+y,xafter*ysize+y,0.25*(xprefactora-1)+gsl_matrix_get(edmatrix,x*ysize+y,xafter*ysize+y));//b_l^\dagger b_{l+1}
    gsl_matrix_set(edmatrix,x*ysize+y,xbefore*ysize+y,0.25*(xprefactorb-1)+gsl_matrix_get(edmatrix,x*ysize+y,xbefore*ysize+y));//b_l^\dagger b_{l-1}
    
    gsl_matrix_set(edmatrix,xafter*ysize+y,x*ysize+y,0.25*(xprefactora-1)+gsl_matrix_get(edmatrix,xafter*ysize+y,x*ysize+y));//b_{l+1}^\dagger b_{l}
    gsl_matrix_set(edmatrix,xbefore*ysize+y,x*ysize+y,0.25*(xprefactorb-1)+gsl_matrix_get(edmatrix,xbefore*ysize+y,x*ysize+y));//b_{l-1}^\dagger b_{l}
        
    gsl_matrix_set(edmatrix,x*ysize+y,x*ysize+yafter,0.25*(yprefactora-1)+gsl_matrix_get(edmatrix,x*ysize+y,x*ysize+yafter));//b_l^\dagger b_{l+1}
    gsl_matrix_set(edmatrix,x*ysize+y,x*ysize+ybefore,0.25*(yprefactorb-1)+gsl_matrix_get(edmatrix,x*ysize+y,x*ysize+ybefore));//b_l^\dagger b_{l-1}
    
    gsl_matrix_set(edmatrix,x*ysize+yafter,x*ysize+y,0.25*(yprefactora-1)+gsl_matrix_get(edmatrix,x*ysize+yafter,x*ysize+y));//b_{l+1}^\dagger b_{l}
    gsl_matrix_set(edmatrix,x*ysize+ybefore,x*ysize+y,0.25*(yprefactorb-1)+gsl_matrix_get(edmatrix,x*ysize+ybefore,x*ysize+y));//b_{l-1}^\dagger b_{l}

    gsl_matrix_set(edmatrix,x*ysize+yafter,x*ysize+yafter,0.5*yprefactora+epsilon+gsl_matrix_get(edmatrix,x*ysize+yafter,x*ysize+yafter));//b_{l+1}^\dagger b_{l+1}
    gsl_matrix_set(edmatrix,x*ysize+ybefore,x*ysize+ybefore,0.5*yprefactorb+epsilon+gsl_matrix_get(edmatrix,x*ysize+ybefore,x*ysize+ybefore));//b_{l-1}^\dagger b_{l-1}
    
    gsl_matrix_set(edmatrix,xafter*ysize+y,xafter*ysize+y,0.5*xprefactora+epsilon+gsl_matrix_get(edmatrix,xafter*ysize+y,xafter*ysize+y));//b_{l+1}^\dagger b_{l+1}
    gsl_matrix_set(edmatrix,xbefore*ysize+y,xbefore*ysize+y,0.5*xprefactorb+epsilon+gsl_matrix_get(edmatrix,xbefore*ysize+y,xbefore*ysize+y));//b_{l-1}^\dagger b_{l-1}
    
    gsl_matrix_set(edmatrix,x*ysize+y,x*ysize+y,0.5*xprefactora+epsilon+gsl_matrix_get(edmatrix,x*ysize+y,x*ysize+y));//b_l^\dagger b_{l}
    gsl_matrix_set(edmatrix,x*ysize+y,x*ysize+y,0.5*xprefactorb+epsilon+gsl_matrix_get(edmatrix,x*ysize+y,x*ysize+y));//b_l^\dagger b_{l}
    
    gsl_matrix_set(edmatrix,x*ysize+y,x*ysize+y,0.5*yprefactora+epsilon+gsl_matrix_get(edmatrix,x*ysize+y,x*ysize+y));//b_l^\dagger b_{l}
    gsl_matrix_set(edmatrix,x*ysize+y,x*ysize+y,0.5*yprefactorb+epsilon+gsl_matrix_get(edmatrix,x*ysize+y,x*ysize+y));//b_l^\dagger b_{l}
    
    gsl_matrix_set(edmatrix,x*ysize+y,xafter*ysize+y+systemsize,0.25*(xprefactora+1)+gsl_matrix_get(edmatrix,x*ysize+y,xafter*ysize+y+systemsize));//b_l^\dagger b_{l+1}^\dagger
    gsl_matrix_set(edmatrix,x*ysize+y,xbefore*ysize+y+systemsize,0.25*(xprefactorb+1)+gsl_matrix_get(edmatrix,x*ysize+y,xbefore*ysize+y+systemsize));//b_l^\dagger b_{l-1}^\dagger
    
    gsl_matrix_set(edmatrix,xafter*ysize+y,x*ysize+y+systemsize,0.25*(xprefactora+1)+gsl_matrix_get(edmatrix,xafter*ysize+y,x*ysize+y+systemsize));//b_{l+1}^\dagger b_l^\dagger
    gsl_matrix_set(edmatrix,xbefore*ysize+y,x*ysize+y+systemsize,0.25*(xprefactorb+1)+gsl_matrix_get(edmatrix,xbefore*ysize+y,x*ysize+y+systemsize));//b_{l-1}^\dagger b_l^\dagger
    
    gsl_matrix_set(edmatrix,x*ysize+y,x*ysize+yafter+systemsize,0.25*(yprefactora+1)+gsl_matrix_get(edmatrix,x*ysize+y,x*ysize+yafter+systemsize));//b_l^\dagger b_{l+1}^\dagger
    gsl_matrix_set(edmatrix,x*ysize+y,x*ysize+ybefore+systemsize,0.25*(yprefactorb+1)+gsl_matrix_get(edmatrix,x*ysize+y,x*ysize+ybefore+systemsize));//b_l^\dagger b_{l-1}^\dagger
    
    gsl_matrix_set(edmatrix,x*ysize+yafter,x*ysize+y+systemsize,0.25*(yprefactora+1)+gsl_matrix_get(edmatrix,x*ysize+yafter,x*ysize+y+systemsize));//b_{l+1}^\dagger b_l^\dagger
    gsl_matrix_set(edmatrix,x*ysize+ybefore,x*ysize+y+systemsize,0.25*(yprefactorb+1)+gsl_matrix_get(edmatrix,x*ysize+ybefore,x*ysize+y+systemsize));//b_{l-1}^\dagger b_l^\dagger
    
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
