#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include "parameter.h"

using namespace std;


int main(){
    double E;
    double help;
    
    double* parameter = new double[xsize];
    for(int ix=0;ix<xsize;ix++){
        parameter[ix] = 0.0;
    }
    gsl_vector *state = gsl_vector_alloc(2*xsize);//phi in first half, theta in second half
    for(int ix=0;ix<xsize;ix++){
        gsl_vector_set (state,ix,0.0*M_PI);
        gsl_vector_set (state,ix+xsize,M_PI*((ix)%2));//
        //initialise Neel state
    }
    ostringstream fm;
    fm << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << "/textureD" << 0 << "cJ.dat";
    ofstream writezero(fm.str().c_str(),ios_base::binary);
    for(int ix=0; ix<xsize;ix++){
        help = gsl_vector_get(state,ix);
        writezero.write((char*)& help, sizeof(double));
    }
    for(int ix=0; ix<xsize;ix++){
        help = gsl_vector_get(state,ix+xsize);
        writezero.write((char*)& help, sizeof(double));
    }
    writezero.close();
    for(int count=0;count<=(int)(Dmax/Dincrement)+1;count++){
        E = count*Dincrement;
        for(int ix=0;ix<xsize;ix++){
            if(ix>=freefieldsize && ix<(freefieldsize+fieldsize)){
                parameter[ix] = E;
            }
            else{
                parameter[ix] = 0.0;
            }
        }
        //initialise x-spiral state in field region and Neel state in the rest
        for(int ix=0;ix<xsize;ix++){
            if(ix>=freefieldsize && ix<(freefieldsize+fieldsize)){
                gsl_vector_set (state,ix,0.0*M_PI);
                gsl_vector_set (state,ix+xsize,M_PI*((ix)%2)-(ix-freefieldsize)*atan(E));
                
            }
            else{
                gsl_vector_set (state,ix,0.0*M_PI);
                if(ix<freefieldsize){
                    gsl_vector_set (state,ix+xsize,M_PI*((ix)%2));
                }
                else{
                    gsl_vector_set (state,ix+xsize,
                                    M_PI*((ix)%2)-(fieldsize)*atan(E));
                }
            }
        }
        if(count % D_OUTPUT_COUNTER ==0){
            ostringstream fn;
            fn << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << "/textureD" << count << "cJ.dat";
            ofstream write(fn.str().c_str(),ios_base::binary);
            
            for(int ix=0; ix<xsize;ix++){
                help = gsl_vector_get(state,ix);
                write.write((char*)& help, sizeof(double));
            }
            for(int ix=0; ix<xsize;ix++){
                help = gsl_vector_get(state,ix+xsize);
                write.write((char*)& help, sizeof(double));
            }
            write.close();
        }
    }
    gsl_vector_free (state);
    return 0;
}
