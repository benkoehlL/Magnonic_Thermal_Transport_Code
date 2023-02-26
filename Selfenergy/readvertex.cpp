#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#define size 100

using namespace std;


int main()
{
    double* help = new double[size*size*size*size];
    double b=0.1;
    while(b<1.0){
        ostringstream fin;
        fin << "vertexB" << b*100 << "cBs.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        for(int kx=0; kx<size; kx++){
            for(int ky=0; ky<size; ky++){
                for(int qx=0; qx<size; qx++){
                    for(int qy=0; qy<size; qy++){
                        read.seekg(((size*size)*(kx*size+ky)+qx*size+qy)*sizeof(double));
                        read.read((char*)& help[(size*size)*(kx*size+ky)+size*qx+qy], sizeof(double));
                    }
                }
            }
        }
        read.close();
        cout << '\n';
        b = b + 0.1;
    }
    return 0;
}

 
