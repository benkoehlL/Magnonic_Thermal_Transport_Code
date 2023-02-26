#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include "parameter.h"

using namespace std;

double cubicinterpolate(double[4],double);
double bicubicinterpolate(double[4][4],double,double);

int main(){
    int dklimit = dksize+1;
    double helppoint[4][4];
    double* helpsort = new double[ksize*(dklimit)*(dklimit)];
    double help;
    double T,b;
    double blist[6] = {0.7,0.75,0.78,0.8,0.85,0.9};
    complex<double>* point = new complex<double>[ksize*ksize];
    //while(b<1.0){
    for(int i=0;i<6;i++){
        b = 100*blist[i];
        T = Tstart;
        while(T<Tend){
            ostringstream fin;
            //fin << "Tempselfenergy/tempselfenergyB" << b*100 << "cBs.dat";
            //fin << "SelfenergyT10cJvariousBfields/tempselfenergyB" << b*100 << "cBs.dat";
            fin << "SelfenergyBconst/" << b << "BsvariousTemperatures/tempselfenergyT" << T*100 << "cJ.dat";
            //fin << "Tempselfenergy/tempselfenergyT" << T*100 << "cJ.dat";
            ifstream read(fin.str().c_str(),ios_base::binary);
            ostringstream fout;
            //fout << "Tempinterpolateddecayrate/interpolateddecayrateT" << T*100 << "cJ.dat";
            //fout << "SelfenergyT10cJvariousBfields/interpolateddecayrateB" << b*100 << "cBs.dat";
            fout << "SelfenergyBconst/" << b << "BsvariousTemperatures/interpolateddecayrateT" << T*100 << "cJ.dat";
            //fout << "Tempselfenergy/interpolateddecayrateT" << T*100 << "cJ.dat";
            ofstream write(fout.str().c_str(),ios_base::binary);
            for(int kx=0; kx<ksize; kx++){
                for(int ky=0; ky<ksize; ky++){
                    read.seekg((kx*ksize+ky)*sizeof(complex<double>));
                    read.read((char*)& point[kx*ksize+ky], sizeof(complex<double>));
                }
            }
            for(int kx=0; kx<ksize; kx++){
                for(int ky=0; ky<ksize;ky++){
                    for(int x=0;x<4;x++){
                        for(int y=0;y<4;y++){
                            if((kx==0 && x==0)){
                                if((ky==0 && y==0)){
                                    helppoint[x][y]= -imag(point[(kx+x+1)*ksize+(ky+y+1)]);
                                }
                                else if(ky==ksize-2 && y==3){
                                    helppoint[x][y]= -imag(point[(kx+x+1)*ksize+(ky+y-3)]);
                                }
                                else if(ky==ksize-1 && y==2){
                                    helppoint[x][y] = -imag(point[(kx+x+1)*ksize+(ky+y-2)]);
                                }
                                else{
                                    helppoint[x][y] = -imag(point[(kx+x+1)*ksize+(ky+y-1)]);
                                }
                            }
                            else if((kx==ksize-2 && x==3)){
                                if((ky==0 && y==0)){
                                    helppoint[x][y]= -imag(point[(kx+x-3)*ksize+(ky+y+1)]);
                                }
                                else if(ky==ksize-2 && y==3){
                                    helppoint[x][y]= -imag(point[(kx+x-3)*ksize+(ky+y-3)]);
                                }
                                else if(ky==ksize-1 && y==2){
                                    helppoint[x][y] = -imag(point[(kx+x-3)*ksize+(ky+y-2)]);
                                }
                                else{
                                    helppoint[x][y] = -imag(point[(kx+x-3)*ksize+(ky+y-1)]);
                                }
                            }
                            else if((kx==ksize-1) && x==2){
                                if(ky==0 && y==0){
                                    helppoint[x][y] = -imag(point[(kx+x-2)*ksize+(ky+y+1)]);
                                }
                                else if(ky==ksize-2 && y==3){
                                    helppoint[x][y] = -imag(point[(kx+x-2)*ksize+(ky+y-3)]);
                                }
                                else if(ky==ksize-1 && y==2){
                                    helppoint[x][y] = -imag(point[(kx+x-2)*ksize+(ky+y-2)]);
                                }
                                else{
                                    helppoint[x][y] = -imag(point[(kx+x-2)*ksize+(ky+y-1)]);
                                }
                            }                            
                            else{
                                if((ky==0 && y==0)){
                                    helppoint[x][y]= -imag(point[(kx+x-1)*ksize+(ky+y+1)]);
                                }
                                else if(ky==ksize-2 && y==3){
                                    helppoint[x][y]= -imag(point[(kx+x-1)*ksize+(ky+y-3)]);
                                }
                                else if(ky==ksize-1 && y==2){
                                    helppoint[x][y] = -imag(point[(kx+x-1)*ksize+(ky+y-2)]);
                                }
                                else{
                                    helppoint[x][y] = -imag(point[(kx+x-1)*ksize+(ky+y-1)]);
                                }
                            }
                        }
                    }
                    if((kx!=ksize-1) && (ky!=ksize-1)){
                        for(int dkx=0; dkx<dklimit; dkx++){
                            for(int dky=0; dky<dklimit;dky++){
                                helpsort[ky*(dklimit*dklimit)+(dkx*dklimit+dky)] = bicubicinterpolate(helppoint,dkx*1.0/(dklimit),dky*1.0/(dklimit));
                            }
                        }
                    }
                    else if((kx==ksize-1) && (ky!=ksize-1)){
                        for(int dky=0; dky<dklimit; dky++){
                            int dkx = 0;
                            helpsort[ky*(dklimit*dklimit)+(dkx*dklimit+dky)] = bicubicinterpolate(helppoint,dkx*1.0/(dklimit),dky*1.0/(dklimit));
                        }
                    }
                    else if((ky==ksize-1) && (kx!=ksize-1)){
                        for(int dkx=0; dkx<dklimit; dkx++){
                            int dky = 0;
                            helpsort[ky*(dklimit*dklimit)+(dkx*dklimit+dky)] = bicubicinterpolate(helppoint,dkx*1.0/(dklimit),dky*1.0/(dklimit));
                        }
                    }
                    else{
                        int dkx = 0;
                        int dky = 0;
                        helpsort[ky*(dklimit*dklimit)+(dkx*dklimit+dky)] = bicubicinterpolate(helppoint,dkx*1.0/(dklimit),dky*1.0/(dklimit));
                    }
                }
                for(int dkx=0;dkx<dklimit;dkx++){
                    for(int ky=0; ky<ksize;ky++){
                        if((kx!=ksize-1) && (ky!=ksize-1)){
                            for(int dky=0; dky<dklimit;dky++){
                                help = helpsort[ky*(dklimit*dklimit)+(dkx*dklimit+dky)];
                                write.write((char*)& help, sizeof(double));
                            }
                        }
                        else if((kx==ksize-1) && (ky!=ksize-1)){
                            if(dkx==0){
                                for(int dky=0; dky<dklimit;dky++){
                                    help = helpsort[ky*(dklimit*dklimit)+(dkx*dklimit+dky)];
                                    write.write((char*)& help, sizeof(double));
                                }
                            }
                        }
                        else if((ky==ksize-1) && (kx!=ksize-1)){
                            int dky = 0;
                            help = helpsort[ky*(dklimit*dklimit)+(dkx*dklimit+dky)];
                            write.write((char*)& help, sizeof(double));
                        }
                        else{
                            if(dkx==0){
                                int dky = 0;
                                help = helpsort[ky*(dklimit*dklimit)+(dkx*dklimit+dky)];
                                write.write((char*)& help, sizeof(double));
                            }
                        }
                    }
                }
            }
            write.close();
            read.close();
            cerr << (T-Tstart)/Tend << "\r";        
            T = T + Tstep;
            //b = b + 0.05;
        }
    }
    return 0;
}

double cubicinterpolate(double x[4], double dx){
    double a = 0.5*(x[3]-3*x[2]+3*x[1]-x[0]);
    double b = 0.5*(-x[3]+4*x[2]-5*x[1]+2*x[0]);
    double c = 0.5*(x[2]-x[0]);
    double d = x[1];
    return(a*dx*dx*dx+b*dx*dx+c*dx+d);
}

double bicubicinterpolate(double points[4][4],double x,double y){
    double help[4];
    help[0] = cubicinterpolate(points[0],y);
    help[1] = cubicinterpolate(points[1],y);
    help[2] = cubicinterpolate(points[2],y);
    help[3] = cubicinterpolate(points[3],y);
    return cubicinterpolate(help,x);
}
