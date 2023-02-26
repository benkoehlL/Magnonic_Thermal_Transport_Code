#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

#define size 200

using namespace std;

int main(){
    int count, kx, ky;;
    double b = 0.750;
    double* help= new double[4*(size+1)];
    double klength;
    while(b<1.0){
        ostringstream fin;
        fin << "Dispersion/dispB" << b*1000 << "mBs.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        count = 0;
        kx = 0;
        ky = 0;
        while(count<4*(size+1)){
            if(count<(size+1)){
                read.seekg(((kx*(size+1)+ky))*sizeof(double));
                read.read((char*)& help[count], sizeof(double));
                //cout << help[count] << '\t' << kx << '\t' << ky << '\n';
                if(kx<size){
                kx++;
                ky++;}
                count++;
            }
            else if(count>=(size+1) && count<2*(size+1)){
                if(kx>0){
                    kx--;
                }
                read.seekg(((kx*(size+1)+ky))*sizeof(double));
                read.read((char*)& help[count], sizeof(double));
                //cout << help[count] << '\t' << kx << '\t' << ky << '\n';
                count++;
            }
            else if(count>=2*(size+1) && count <= 3*(size+1)){
                if(kx<size){
                    kx++;
                    ky--;
                }
                read.seekg(((kx*(size+1)+ky))*sizeof(double));
                read.read((char*)& help[count], sizeof(double));
                //cout << help[count] << '\t' << kx << '\t' << ky << '\n';
                count++;
            }
            else if(count>=3*(size+1)){
                if(kx>0){
                    kx--;
                }
                read.seekg(((kx*(size+1)+ky))*sizeof(double));
                read.read((char*)& help[count], sizeof(double));
               // cout << help[count] << '\t' << kx << '\t' << ky << '\n';
                count++;
            }
        }
        read.close();
        ostringstream fout;
        fout << "Dispersion/dispplotB" << b*1000 << "mBs.dat";
        ofstream write(fout.str().c_str());
        klength = 0.0;
        for(count=0; count<4*(size+1);count++){
            if(count<(size+1)){
                write << klength << '\t' << help[count] << '\n';
                klength = klength + sqrt(2.0)*M_PI/(size);
            }
            else if(count>=(size+1) && count<2*(size+1)){
                write << klength << '\t' << help[count] << '\n';
                klength = klength + M_PI/size;
            }
            else if(count>=2*(size+1) && count <= 3*(size+1)){
                write << klength << '\t' << help[count] << '\n';
                klength = klength + sqrt(2.0)*M_PI/size;
            }
            else if(count>=3*(size+1)){
                write << klength << '\t' << help[count] << '\n';
                klength = klength + M_PI/size;
            }
        }
        write.close();
        b = b + 0.025;
    }
    return 0;
}
