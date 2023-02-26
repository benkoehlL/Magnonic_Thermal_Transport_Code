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
    
    for(int count=0;count<=(int)(Dmax/Dincrement)+1;count++){
        //gsl_matrix_complex* zhamiltonian = gsl_matrix_complex_alloc(2*systemsize,2*systemsize);
        gsl_matrix* diagonal = gsl_matrix_alloc(2*systemsize,2*systemsize);
        interaction = count * Dincrement;

        //Calculate Hamiltonian directly in programm
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                if(x>=freefieldsize && x<(freefieldsize+fieldsize)){
                    parameter[x*ysize+y] = interaction;
                }
                else{
                    parameter[x*ysize+y] = 0.0;
                }
            }
        }
        ostringstream fin1;
        fin1 << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" 
            << fieldsize << "/textureD" << interaction*100 << "cJxsize_simplex.dat";
        ifstream read1(fin1.str().c_str(),ios_base::binary);
        
        //read texture
        for(int x=0;x<xsize;x++){
            for(int y=0;y<ysize;y++){
                read1.seekg((x*ysize+y+systemsize)*sizeof(double));
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
        fn << "Bogoliubov/eigenvalueD" << interaction*100 << "cJxsize" << xsize 
            << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        ofstream writeeval(fn.str().c_str(),ios_base::binary);//add when using binary format: ,ios_base::binary
        fm << "Bogoliubov/eigenvectorD" << interaction*100 << "cJxsize" << xsize 
            << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
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
                       -0.25*(yprefactorb-1)+gsl_matrix_get(edmatrix,ynegneighbour,site));//b_{l-1}^\dagger b_{l}

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
