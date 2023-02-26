#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include "parameter.h"

using namespace std;

int main(){
    double xprefactor,yprefactor,help;
    double interaction;
    double* state = new double[systemsize];
    double* hamiltonian = new double[4*systemsize*systemsize];
    double* D = new double[systemsize];
    for(int count=1;count<=(int)(Dmax/Dincrement)/1000+1;count++){
        interaction = count * Dincrement *1000;
        ostringstream fin;
        fin << "Texture/textureD" << interaction*100 << "cJ.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                read.seekg((x*ysize+y)*sizeof(double));
                read.read((char*)& state[x*ysize+y], sizeof(double));
                
                for(int xhelp=0;xhelp<xsize;xhelp++){
                    for(int yhelp=0;yhelp<ysize;yhelp++){
                        hamiltonian[2*(x*ysize+y)*systemsize+(xhelp*ysize+yhelp)] = 0.0;
                        hamiltonian[2*(x*ysize+y)*systemsize+(xhelp*ysize+yhelp)+systemsize] = 0.0;
                        hamiltonian[2*((x*ysize+y)+systemsize)*systemsize+(xhelp*ysize+yhelp)] = 0.0;
                        hamiltonian[2*((x*ysize+y)+systemsize)*systemsize+(xhelp*ysize+yhelp)+systemsize] = 0.0;
                    }
                }
            }
        }
       
        read.close();
        for(int i=0;i<xsize;i++){
            for(int j=0;j<ysize;j++){
                if(i>=freefieldsize && i<(freefieldsize+fieldsize)){
                    if(i==freefieldsize && j==0){
                        D[i*ysize+j] = interaction;
                    }
                    else if(i!=freefieldsize && j==0){
                        D[i*ysize+j] = (-1)*D[(i-1)*ysize+j];
                    }
                    else{
                        D[i*ysize+j] = (-1)*D[i*ysize+j-1];
                    }
                }
                else{
                    D[i*ysize+j] = 0.0;
                }
            }
        }
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                //y-neighbour
                yprefactor = S*(J*cos(state[x*ysize+y]-state[x*ysize+(y+1)%ysize])-D[x*ysize+y]*sin(state[x*ysize+y]-state[x*ysize+(y+1)%ysize]));
                
                hamiltonian[(x*ysize+y)*2*systemsize+((x*ysize+(y+1)%ysize))] += 0.5*(yprefactor-J*S);//a_i^\dagger*a_{i+1}
                hamiltonian[(x*ysize+(y+1)%ysize)*2*systemsize+x*ysize+y] += 0.5*(yprefactor-J*S);//a_{i-1}^\dagger*a_{i}
                
                hamiltonian[((x*ysize+y))*2*systemsize+(x*ysize+y)] += yprefactor;//a_i^\dagger*a_{i}
                
                hamiltonian[((x*ysize+(y+1)%ysize))*2*systemsize+(x*ysize+(y+1)%ysize)] += yprefactor;//a_{i+1)^\dagger*a_{i+1}
                
                hamiltonian[(x*ysize+y)*2*systemsize+((x*ysize+(y+1)%ysize)+systemsize)] += 0.5*(yprefactor+J*S);//a_i^\dagger*a_{i+1}^\dagger
                hamiltonian[(x*ysize+(y+1)%ysize)*2*systemsize+((x*ysize+y)+systemsize)] += 0.5*(yprefactor+J*S);//a_i^\dagger*a_{i-1}^\dagger
                
                //x-neighbour (not periodic)
                /*if(x<(xsize-1)){
                    xprefactor = S*(J*cos(state[x*ysize+y]-state[(x+1)*ysize+y])-D[x*ysize+y]*sin(state[x*ysize+y]-state[(x+1)*ysize+y]));
                
                    hamiltonian[(x*ysize+y)*2*systemsize+((x+1)*ysize+y)] += 0.5*(xprefactor-J*S);//a_i^\dagger*a_{i+1}
                    hamiltonian[((x+1)*ysize+y)*2*systemsize+x*ysize+y] += 0.5*(xprefactor-J*S);//a_{i-1}^\dagger*a_{i}
                    
                    hamiltonian[(x*ysize+y)*2*systemsize+((x)*ysize+y)] += xprefactor;//a_i^\dagger*a_{i}
                    
                    hamiltonian[((x+1)*ysize+y)*2*systemsize+((x+1)*ysize+y)] += xprefactor;//a_{i+1}^\dagger*a_{i+1}

                    hamiltonian[(x*ysize+y)*2*systemsize+((x+1)*ysize+y+systemsize)] += 0.5*(xprefactor+J*S);//a_i^\dagger*a_{i+1}^\dagger
                    hamiltonian[((x+1)*ysize+y)*2*systemsize+(x*ysize+y+systemsize)] += 0.5*(xprefactor+J*S);//a_i^\dagger*a_{i-1}^\dagger
                }*/
                //x-neighbour (periodic x)
                xprefactor = S*(J*cos(state[x*ysize+y]-state[((x+1)%xsize)*ysize+y])-D[x*ysize+y]*sin(state[x*ysize+y]-state[((x+1)%xsize)*ysize+y]));
                
                hamiltonian[(x*ysize+y)*2*systemsize+((((x+1)%xsize)*ysize+y))] += 0.5*(xprefactor-J*S);//a_i^\dagger*a_{i+1}
                hamiltonian[(((x+1)%xsize)*ysize+y)*2*systemsize+x*ysize+y] += 0.5*(xprefactor-J*S);//a_{i-1}^\dagger*a_{i}
                
                hamiltonian[((x*ysize+y))*2*systemsize+(x*ysize+y)] += xprefactor;//a_i^\dagger*a_{i}
                
                hamiltonian[((((x+1)%xsize)*ysize+y))*2*systemsize+(((x+1)%xsize)*ysize+y)] += xprefactor;//a_{i+1)^\dagger*a_{i+1}
                
                hamiltonian[(x*ysize+y)*2*systemsize+((x+1)%xsize)*ysize+y+systemsize] += 0.5*(xprefactor+J*S);//a_i^\dagger*a_{i+1}^\dagger
                hamiltonian[(((x+1)%xsize)*ysize+y)*2*systemsize+(x*ysize+y)+systemsize] += 0.5*(xprefactor+J*S);//a_i^\dagger*a_{i-1}^\dagger
            }
            
        }
        
        //add complex conjugated elements
        for(int x=0;x<systemsize;x++){
            for(int y=0;y<systemsize;y++){
                hamiltonian[(x+systemsize)*2*systemsize+y+systemsize] = hamiltonian[x*2*systemsize+y];
            }
        }
        for(int x=0;x<systemsize;x++){
            for(int y=0;y<systemsize;y++){
                hamiltonian[(y+systemsize)*2*systemsize+x] = hamiltonian[x*2*systemsize+y+systemsize];
            }
        }
    
        if(count==(int)(Dmax/Dincrement)/1000){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << hamiltonian[2*x*systemsize+y] << '\t';
                }
                cout << '\n';
            }
        }
        ostringstream fn;
        fn << "Hamiltonian/hamiltonianD" << interaction*100 << "cJ.dat";
        ofstream write(fn.str().c_str(),ios_base::binary);
        for(int i=0;i<4*systemsize*systemsize;i++){
            help = hamiltonian[i];
            write.write((char*)& help, sizeof(double));
        }
    }
    return 0;
}
