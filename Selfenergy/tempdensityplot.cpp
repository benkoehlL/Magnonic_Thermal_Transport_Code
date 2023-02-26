#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "parameter.h"

using namespace std;

int main(){
    double T=0.2;
    int limit = (ksize)*(dksize+1);
    int kspacesize = limit-dksize;
    double* help = new double[kspacesize*kspacesize];
    double tlist[3] = {0.1,0.2,0.3};
    for(int b=10;b<1000;b+=10){
    //for(int b=70;b<95;b=b+5){
        //T = Tstart;
	for(int i=0;i<3;i++){
        	T=tlist[i];
    	//while(T<Tend){
        	ostringstream fin;
		fin << "Compactinterpolateddecayrate/interpolateddecayrateB" << b << "mBsT" << T*100 << "cJ.dat";
		//fin << "SelfenergyBconst/" << b << "BsvariousTemperatures/interpolateddecayrateT" << T*100 << "cJ.dat";
        	ifstream read(fin.str().c_str(),ios_base::binary);
        	for(int kx=0;kx<kspacesize;kx++){
            		for(int ky=0;ky<kspacesize;ky++){
                		read.seekg((kx*(kspacesize)+ky)*sizeof(double));
                		read.read((char*)& help[kx*kspacesize+ky], sizeof(double));
                    }
            }
	        read.close();
        	ostringstream fout;
		fout << "test/tempdensityplotB" << b << "mBsT" << T*100 << "cJ.dat";        	
		//fout << "SelfenergyBconst/" << b << "BsvariousTemperatures/tempdensityplotT" << T*100 << "cJ.dat";
        	ofstream write(fout.str().c_str());
        	for(int kx=0;kx<kspacesize;kx++){
        	    for(int ky=0;ky<kspacesize;ky++){
                	write << kx*M_PI/(kspacesize-1) << '\t' << ky*M_PI/(kspacesize-1) << '\t' << help[kx*kspacesize+ky] << '\n';
                	//write << kx << '\t' << ky << '\t' << help[kx*kspacesize+ky] << '\n';
            	}
            }
        	write.close();
        	/*if(b<0.7){
            		b = b + 0.2;
        	}
        	else{
            		b = b + 0.025;
        	}*/
		//cerr << (T-Tstart)/Tend << "\r";
        	cerr << double(i/3) << '\t' << b*0.0001 << "\r";
		//T = T + Tstep;
    	}
    }
    return 0;
}
