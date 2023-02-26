/*for details see Colpa '78: whole procedure is described there*/

#include <iostream>
#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <fstream>
#include <sstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "parameter.h"

using namespace std;

void energydensity(double[systemsize],gsl_matrix*,int,int,double[systemsize]);

int main(){
    double interaction;
    double help;
    double* texture = new double[systemsize];
    double* parameter = new double[systemsize];
    int helporder;
    gsl_vector *helpvector = gsl_vector_alloc(2*systemsize);
    gsl_matrix* paraunity = gsl_matrix_alloc(2*systemsize,2*systemsize);
    gsl_matrix* helpmatrixX = gsl_matrix_alloc(2*systemsize,2*systemsize);
    gsl_matrix* helpmatrixY = gsl_matrix_alloc(2*systemsize,2*systemsize);
    gsl_matrix* helpmatrixH = gsl_matrix_alloc(2*systemsize,2*systemsize);
    gsl_matrix_set_identity(paraunity);
    for(int x=0;x<systemsize;x++){
        gsl_matrix_set(paraunity,x+systemsize,x+systemsize,-1.0);
    }
    
    for(int count=1;count<=(int)(Dmax/Dincrement)/1000+1;count++){
        //gsl_matrix_complex* zhamiltonian = gsl_matrix_complex_alloc(2*systemsize,2*systemsize);
        gsl_matrix* diagonal = gsl_matrix_alloc(2*systemsize,2*systemsize);
        interaction = count * Dincrement *1000;
        
        //read hamiltonian from file
        /*ostringstream fin;
        fin << "Hamiltonian/hamiltonianD" << interaction*100 << "cJ.dat";
        ifstream read(fin.str().c_str(),ios_base::binary);
        
        //for k-space hamiltonian
        /*for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                read.seekg((2*x*systemsize+y)*sizeof(double));
                read.read((char*)& help, sizeof(double));                     
                gsl_matrix_complex_set(zhamiltonian,x,y,gsl_complex_rect(help,0.0));
            }
        }
        
        // for real-space Hamiltonian
        gsl_matrix* hamiltonian = gsl_matrix_alloc(2*systemsize,2*systemsize);
        for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                read.seekg((2*x*systemsize+y)*sizeof(double));
                read.read((char*)& help, sizeof(double));                     
                gsl_matrix_set(hamiltonian,x,y,help);
            }
        }
        read.close();
        */

        //Calculate Hamiltonian directly in programm
        for(int i=0;i<xsize;i++){
            for(int j=0;j<ysize;j++){
                if(i>=freefieldsize && i<(freefieldsize+fieldsize)){
                    if(i==freefieldsize && j==0){
                        parameter[i*ysize+j] = -interaction;
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
        fin1 << "Texture/textureD" << interaction*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        ifstream read1(fin1.str().c_str(),ios_base::binary);
        
        //read texture
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                read1.seekg((x*ysize+y)*sizeof(double));
                read1.read((char*)& texture[x*ysize+y], sizeof(double));                   
            }
        }
        read1.close();
        
        gsl_matrix* hamiltonian = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix_set_zero(hamiltonian);
        
        for(int i=0;i<xsize;i++){
            for(int j=0;j<ysize;j++){
                energydensity(texture,helpmatrixH,i,j,parameter);
                gsl_matrix_add(hamiltonian,helpmatrixH);
            }
        }
        
        
        //Fourier matrix
        /*
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
        
        //Fourier transformation of hamiltonian
        
        gsl_matrix_complex* helpzmatrix = gsl_matrix_complex_alloc(2*systemsize,2*systemsize);
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1.0,0.0), fourier, zhamiltonian, gsl_complex_rect(0.0,0.0), helpzmatrix);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.0,0.0), helpzmatrix, fourier, gsl_complex_rect(0.0,0.0), zhamiltonian);
        gsl_matrix* hamiltonian = gsl_matrix_alloc(2*systemsize,2*systemsize);
        for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                gsl_matrix_set(hamiltonian,x,y,GSL_REAL(gsl_matrix_complex_get(zhamiltonian,x,y)));
                /*if(GSL_IMAG(gsl_matrix_complex_get(zhamiltonian,x,y))>0.00001){
                    cout << "Something is fishy ;) " << '\n';
                }
            }
        }
        
        /*if(count == 1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(hamiltonian,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //bring hamiltonian in blockdiagonal form
        /*
        gsl_matrix* helpmatrix = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix_set_zero(helpmatrix);
        for(int x=0;x<systemsize;x++){
                gsl_matrix_set(helpmatrix,2*x,2*x,gsl_matrix_get(hamiltonian,x,x));
                gsl_matrix_set(helpmatrix,2*x+1,2*x+1,gsl_matrix_get(hamiltonian,x+systemsize,x+systemsize));
                gsl_matrix_set(helpmatrix,2*x+1,2*x,gsl_matrix_get(hamiltonian,x+systemsize,x));
                gsl_matrix_set(helpmatrix,2*x,2*x+1,gsl_matrix_get(hamiltonian,x,x+systemsize));
        }
        /*for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                gsl_matrix_set(hamiltonian,x,y,gsl_matrix_get(helpmatrix,x,y));
            }
        }
        
        ostringstream bdH;
        bdH << "Blockdiagonal Hamiltonian/blockhamiltonianD" << interaction*100 << "cJ.dat";
        ofstream writebdH(bdH.str().c_str());//add when using binary format: ,ios_base::binary
        for(int i=0;i<2*systemsize;i++){
            for(int j=0;j<2*systemsize;j++){
                //help = gsl_matrix_get(helpmatrix,i,j);
                //writebdH.write((char*)& help, sizeof(double));
                writebdH << gsl_matrix_get(helpmatrix,i,j) << '\t'; //without binary
            }
            writebdH << '\n';
        }
        writebdH.close();
        
        //gsl_matrix_free(helpmatrix);
        gsl_matrix_complex_free(zhamiltonian);
        gsl_matrix_complex_free(helpzmatrix);
        gsl_matrix_complex_free(fourier);
        */
        
//build matrix to be diagonalised:
        //1st: with cholesky decomposition:        
        gsl_linalg_cholesky_decomp1(hamiltonian);
        gsl_matrix* choleskymatrix = gsl_matrix_alloc(2*systemsize,2*systemsize);
        for(int i=0;i<2*systemsize;i++){
            for(int j=0;j<2*systemsize;j++){
                if(j<=i){
                    gsl_matrix_set(choleskymatrix,i,j,gsl_matrix_get(hamiltonian,i,j));
                }
                else{
                    gsl_matrix_set(choleskymatrix,i,j,0);
                }
            }
        }//this yields lower (!) triangular cholesky matrix
        
        /*if(count == 1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(choleskymatrix,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        gsl_matrix* paraunitcholsesky = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,paraunity,choleskymatrix,0.0, paraunitcholsesky);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans,1.0,choleskymatrix,paraunitcholsesky,0.0, diagonal);
        gsl_matrix_free (paraunitcholsesky);
        gsl_matrix_free (hamiltonian);
        
        /*if(count == 1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(diagonal,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //2nd: without cholesky decomposition
        //gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0,paraunity, hamiltonian,0.0, diagonal);
        
        //diagonalise
        gsl_vector *eval = gsl_vector_alloc (2*systemsize);
        gsl_matrix *evec = gsl_matrix_alloc (2*systemsize, 2*systemsize);
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (2*systemsize);
        gsl_eigen_symmv (diagonal, eval, evec, w);
        gsl_eigen_symmv_free (w);
        
        //Print k-dependent eigenvalues before ordering:
        /*if(count == 1){
            for(int x=0;x<2*systemsize;x++){
                cout << 2*gsl_vector_get(eval,x) << '\n';
            }
        }*/
        
        //sort eigensystem in descending order
        //gsl_eigen_symmv_sort (eval, evec,GSL_EIGEN_SORT_VAL_DESC);
        
        int index = 0;
        int maximumindex;
        int order[2*systemsize];
        for(int i=0;i<2*systemsize;i++){
            order[i] = i;
        }
        while(index<(2*systemsize-1)){
            maximumindex = index;
            for(int i=index+1;i<2*systemsize;i++){
                if(gsl_vector_get(eval,maximumindex)<gsl_vector_get(eval,i)){
                    maximumindex = i;
                }
            }
            help = gsl_vector_get(eval,index);
            helporder = order[index];
            gsl_vector_set(eval,index,gsl_vector_get(eval,maximumindex));
            order[index] = order[maximumindex];
            gsl_vector_set(eval,maximumindex,help);
            order[maximumindex] = helporder;
            index++;
        }

        //check eigenvalue symmetry (E(n)=-E(2N-n))
        /*for(int i=0;i<systemsize;i++){
            if(abs(gsl_vector_get(eval,i))-abs(gsl_vector_get(eval,2*systemsize-1-i))>=epsilon){
                cout << "Energies are not equal!!!" << '\t' << gsl_vector_get(eval,i) << '\t' << -gsl_vector_get(eval,2*systemsize-1-i) << '\n';
            }
        }*/

        //sort negative energies so that E(n)=-E(n+N)
        /*for(int i=0;i<systemsize;i++){
            gsl_vector_set(eval,i,-gsl_vector_get(eval,i+systemsize));
            gsl_vector_set(helpvector,i,order[i]);
        }
        
        for(int i=0;i<systemsize;i++){
            order[systemsize-1-i] = (int) gsl_vector_get(helpvector,i);
        }*/
        
        
        /*ostringstream of;
        of << "Order/orderD" << interaction*100 << "cJ.dat";
        ofstream writeorder(of.str().c_str(),ios_base::binary);//add when using binary format: ,ios_base::binary
        for(int i=0;i<2*systemsize;i++){
            writeorder.write((char*)& order[i], sizeof(int));
            //writeorder << order[i] << '\n'; //without binary
        }
        writeorder.close();
        */
        
        gsl_matrix *orderedevec = gsl_matrix_alloc(2*systemsize,2*systemsize);
        for(int i=0;i<2*systemsize;i++){
            gsl_matrix_get_col(helpvector,evec,order[i]);
            gsl_matrix_set_col(orderedevec,i,helpvector);
        }
        gsl_matrix_free (evec);
        
        /*if(count==1){
            for(int x=0;x<2*systemsize;x++){
                for(int y=0;y<2*systemsize;y++){
                    cout << gsl_matrix_get(orderedevec,x,y) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //calculate Bogoliubov-transformation
        gsl_matrix *Lmatrix = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix_set_identity(Lmatrix);
        for(int x=0;x<2*systemsize;x++){                        
            gsl_matrix_set(Lmatrix,x,x,sqrt(abs(1.0/gsl_vector_get(eval,x))));
        }
            
        gsl_matrix *invbogo = gsl_matrix_alloc(2*systemsize,2*systemsize);
        
        gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,Lmatrix,orderedevec,0.0, diagonal);
        gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,diagonal,choleskymatrix,0.0,invbogo);
        
        gsl_matrix_free (Lmatrix);
        
        //H=B^dagger.D.B -> the real transformation matrix is the inverse of B
        gsl_matrix* bogo = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,paraunity,invbogo,0.0,diagonal);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,diagonal,paraunity,0.0,bogo);
        
        gsl_matrix_free (diagonal);
        gsl_matrix_free (choleskymatrix);
        gsl_matrix_free (orderedevec);
        gsl_matrix_free (invbogo);

        //construct Bogoliubov transformation which ensures that if d_i^\dagger=v_{ij}^A c_j^\dagger+v_{ij}^B c_j then d_i=(v_{ij}^A)^\ast c_j+(v_{ij}^B)^\ast c_j^\dagger
        for(int x=0;x<systemsize;x++){
            for(int y=0;y<systemsize;y++){
                gsl_matrix_set(bogo,x,y+systemsize,gsl_matrix_get(bogo,x+systemsize,y));
                gsl_matrix_set(bogo,x+systemsize,y+systemsize,gsl_matrix_get(bogo,x,y));
            }
        }
        
        //write eigenvalues (diagonalelements of Lmatrix) and bogo matrix into a file
        
        ostringstream fn,fm;
        fn << "Bogoliubov/eigenvalueD" << interaction*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        ofstream writeeval(fn.str().c_str(),ios_base::binary);//add when using binary format: ,ios_base::binary
        fm << "Bogoliubov/eigenvectorD" << interaction*100 << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        ofstream writeevec(fm.str().c_str(),ios_base::binary);//add when using binary format: ,ios_base::binary
        for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                if(x==y){
                    if(x<systemsize){
                        help = 2*gsl_vector_get(eval,x);
                        writeeval.write((char*)& help, sizeof(double));
                        //writeeval << 2*gsl_vector_get(eval,x) << '\n'; //without binary
                    }
                    else{
                        help = 2*gsl_vector_get(eval,x);
                        writeeval.write((char*)& help, sizeof(double));
                        //writeeval << 2*gsl_vector_get(eval,x) << '\n'; //without binary
                    }
                        
                }
                    help = gsl_matrix_get(bogo,x,y);
                    writeevec.write((char*)& help, sizeof(double));
                    //writeevec << gsl_matrix_get(bogo,x,y) << '\t'; //without binary
            }
            //writeevec << '\n'; //without binary
        }
        
        /*if(count ==1){
            for(int i=0;i<2*systemsize;i++){
                for(int j=0;j<2*systemsize;j++){
                    cout << gsl_matrix_get(bogo,i,j) << '\t';
                }
                cout << '\n';
            }
        }*/
        
        //3rd: procedure for positive semidefinite hamiltonian
        /*
        gsl_vector *eval = gsl_vector_alloc (2*systemsize);
        gsl_matrix *evec = gsl_matrix_alloc (2*systemsize, 2*systemsize);
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (2*systemsize);
        gsl_eigen_symmv (hamiltonian, eval, evec, w);
        gsl_eigen_symmv_free (w);
        gsl_matrix_free (hamiltonian);
        
        gsl_matrix *Fmatrix = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix *Amatrix = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix *helpmatrix = gsl_matrix_alloc(2*systemsize,2*systemsize);
        
        gsl_matrix_set_identity(Fmatrix);
        for(int x=0;x<2*systemsize;x++){                        
            gsl_matrix_set(Fmatrix,x,x,sqrt(abs(gsl_vector_get(eval,x))));
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,evec, Fmatrix,0.0, helpmatrix);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0,helpmatrix, evec,0.0, Amatrix);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,Amatrix, paraunity,0.0, helpmatrix);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,helpmatrix, Amatrix,0.0, diagonal);
         
        
        gsl_eigen_symmv_workspace * v = gsl_eigen_symmv_alloc (2*systemsize);
        gsl_eigen_symmv (diagonal, eval, evec, v);
        gsl_eigen_symmv_free (v);
        gsl_matrix_free(diagonal);
        gsl_eigen_symmv_sort (eval, evec,GSL_EIGEN_SORT_VAL_DESC);
        
        
        //calculate Bogoliubov-transformation
        gsl_matrix *Lmatrix = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_matrix_set_identity(Lmatrix);
        for(int x=0;x<2*systemsize;x++){                        
            gsl_matrix_set(Lmatrix,x,x,sqrt(abs(1.0/gsl_vector_get(eval,x))));
        }
        gsl_matrix *bogo = gsl_matrix_alloc(2*systemsize,2*systemsize);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,evec,Amatrix,0.0, helpmatrix);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Lmatrix,helpmatrix,0.0,bogo);
               
        
        gsl_matrix_free (Amatrix);
        gsl_matrix_free (Lmatrix);
        gsl_matrix_free (Fmatrix);
        gsl_matrix_free (helpmatrix);
         //write eigenvalues (diagonalelements of Lmatrix) and bogo matrix into a file
        ostringstream fn,fm;
        fn << "Bogoliubov/eigenvalueD" << interaction*100 << "cJ.dat";
        ofstream writeeval(fn.str().c_str());
        fm << "Bogoliubov/eigenvectorD" << interaction*100 << "cJ.dat";
        ofstream writeevec(fm.str().c_str());
        for(int x=0;x<2*systemsize;x++){
            for(int y=0;y<2*systemsize;y++){
                if(x==y){
                    if(x<systemsize){
                        writeeval << abs(gsl_vector_get(eval,x)) << '\n';
                    }
                    else{
                        writeeval << abs(gsl_vector_get(eval,x)) << '\n';
                    }
                        
                }
                    writeevec << gsl_matrix_get(bogo,x,y) << '\t';
            }
            writeevec << '\n';
        }
        */
        writeevec.close();
        writeeval.close();
        
        gsl_vector_free (eval);
        gsl_matrix_free (bogo);
    }
    delete [] texture;
    delete [] parameter;
    gsl_vector_free (helpvector);
    gsl_matrix_free (paraunity); 
    
    return 0;
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
