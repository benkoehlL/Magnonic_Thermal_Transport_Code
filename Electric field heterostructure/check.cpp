#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "parameter.h"

using namespace std;

int main(){
    double interaction;
    double* state = new double[systemsize];
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
                //state[x*ysize+y] = 0.0;
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
        double sum = 0.0;
        for(int i=0;i<xsize-1;i++){
            for(int j=0;j<ysize;j++){
                sum += -J*(sin(state[i*(ysize)+j]-state[(i+1)*ysize+j])+sin(state[i*(ysize)+j]-state[i*ysize+(j+1)%(ysize)]));
                if(i!=freefieldsize+fieldsize-1){
                    sum += D[i*ysize+j]*(cos(state[i*(ysize)+(j)]-state[(i+1)*ysize+j])+cos(state[i*(ysize)+j]-state[i*ysize+(j+1)%(ysize)]));
                }
                else{
                    sum += D[i*ysize+j]*(cos(state[i*(ysize)+j]-state[i*ysize+(j+1)%(ysize)]));
                }
            }
        }
        cout << "D = " << interaction << '\t' << "H1 = " << sum << '\n';
    }
    for(int i=0;i<xsize;i++){
        for(int j=0;j<ysize;j++){
            if(state[i*ysize+j]!=0){
            //cout << state[i*ysize+j] << '\t';
            }
        }
        //cout << '\n';
    }
    return 0;
}
