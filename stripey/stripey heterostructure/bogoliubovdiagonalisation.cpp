/*for details see Colpa '78: whole procedure is described there*/

#include <iostream>
#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <fstream>
#include <sstream>
#include <complex>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "parameter.h"

using namespace std;

void intrastripe(double[xsize],gsl_matrix_complex*,double[xsize]);
void interstripe(double[xsize],gsl_matrix_complex*,double,double[xsize]);

typedef complex<double> dcomp;

int main(){
    dcomp I;
    I = -1;
    I = sqrt(I);
    double fourierexponent,interaction,help;
    gsl_complex helpcomplex;
    int helporder;
    double* texture = new double[xsize];
    double* parameter = new double[xsize];
    
    gsl_vector_complex *helpvector = gsl_vector_complex_alloc(2*xsize);
    gsl_matrix_complex* paraunity = gsl_matrix_complex_alloc(2*xsize,2*xsize);
    gsl_matrix_complex* helpmatrixX = gsl_matrix_complex_alloc(2*xsize,2*xsize);
    gsl_matrix* helpmatrixY = gsl_matrix_alloc(2*xsize,2*xsize);
    gsl_matrix_complex* helpmatrixH = gsl_matrix_complex_alloc(2*xsize,2*xsize);
    gsl_matrix_complex_set_identity(paraunity);
    
    for(int x=0;x<xsize;x++){
        gsl_matrix_complex_set(paraunity,x+xsize,x+xsize,gsl_complex_rect(-1.0,0.0));
    }
    for(int count=0;count<=(int)(Dmax/Dincrement)+1;count++){
        gsl_vector *eval = gsl_vector_alloc (2*xsize);
        interaction = count * Dincrement;
        ostringstream fn;
        fn << "Bogoliubov/Eigenvalues/D" << count << "cJxsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << ".dat";
        ofstream writeeval(fn.str().c_str());//add when using binary format: ,ios_base::binary
        ostringstream fin1;
        fin1 << "Texture/xsize" << xsize << "ysize" << ysize << "fieldsize" << fieldsize << "/textureD" << count << "cJ.dat";
        ifstream read1(fin1.str().c_str(),ios_base::binary);
        
        //read texture
        for(int x=0;x<xsize;x++){
            read1.seekg((x+xsize)*sizeof(double));
            read1.read((char*)& texture[x], sizeof(double));                   
        }
        
        for(int x=0;x<xsize;x++){
            if(x>=freefieldsize && x<(freefieldsize+fieldsize)){
                parameter[x] = interaction;
            }
            else{
                parameter[x] = 0.0;
            }
        }

        read1.close();
        for(int k=0;k<ysize;k++){
            fourierexponent = 2*M_PI*(k*1.0/ysize);
            gsl_matrix* diagonal = gsl_matrix_alloc(2*xsize,2*xsize);
    
            //Calculate Hamiltonian directly in programm
            gsl_matrix_complex* hamiltonian = gsl_matrix_complex_alloc(2*xsize,2*xsize);
            gsl_matrix_complex_set_zero(hamiltonian);
   
            intrastripe(texture,helpmatrixH,parameter);
            gsl_matrix_complex_add(hamiltonian,helpmatrixH);
            interstripe(texture,helpmatrixH,fourierexponent,parameter);
            gsl_matrix_complex_add(hamiltonian,helpmatrixH);
            interstripe(texture,helpmatrixH,-fourierexponent,parameter);
            gsl_matrix_complex_add(hamiltonian,helpmatrixH);
        
            
            if(k==0 && count==100){
            for(int i=0;i<2*xsize;i++){
                for(int j=0;j<2*xsize;j++){
                    cout <<
                    gsl_complex_abs(gsl_matrix_complex_get(hamiltonian,i,j))
			*exp(I*gsl_complex_arg(gsl_matrix_complex_get(hamiltonian,i,j))) << '\t';
            	}
            	cout << '\n';
            }
            cout << '\n' << k << '\n';
            }
            //build matrix to be diagonalised:
            gsl_linalg_complex_cholesky_decomp(hamiltonian);
            gsl_matrix_complex* choleskymatrix = gsl_matrix_complex_alloc(2*xsize,2*xsize);
            for(int i=0;i<2*xsize;i++){
                for(int j=0;j<2*xsize;j++){
                    if(j<=i){
                        gsl_matrix_complex_set(choleskymatrix,i,j,gsl_matrix_complex_get(hamiltonian,i,j));
                    }
                    else{
                        gsl_matrix_complex_set(choleskymatrix,i,j,gsl_complex_rect(0.0,0.0));
                    }
                }
        }//this yields lower (!) triangular cholesky matrix
        
        gsl_matrix_complex* paraunitcholsesky = gsl_matrix_complex_alloc(2*xsize,2*xsize);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,gsl_complex_rect(1.0,0.0),paraunity,choleskymatrix,gsl_complex_rect(0.0,0.0), paraunitcholsesky);
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans,gsl_complex_rect(1.0,0.0),choleskymatrix,paraunitcholsesky,gsl_complex_rect(0.0,0.0), hamiltonian);
        gsl_matrix_complex_free (paraunitcholsesky);
        
        //diagonalise
        gsl_matrix_complex *evec = gsl_matrix_complex_alloc (2*xsize, 2*xsize);
        gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc (2*xsize);
        gsl_eigen_hermv(hamiltonian, eval, evec, w);
        gsl_eigen_hermv_free(w);
        int index = 0;
        int maximumindex;
        int order[2*xsize];
        for(int i=0;i<2*xsize;i++){
            order[i] = i;
        }
        while(index<(2*xsize-1)){
            maximumindex = index;
            for(int i=index+1;i<2*xsize;i++){
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

        gsl_matrix_complex *orderedevec = gsl_matrix_complex_alloc(2*xsize,2*xsize);
        for(int i=0;i<2*xsize;i++){
            gsl_matrix_complex_get_col(helpvector,evec,order[i]);
            gsl_matrix_complex_set_col(orderedevec,i,helpvector);
        }
        gsl_matrix_complex_free (evec);
                
        //calculate Bogoliubov-transformation
        gsl_matrix_complex *Lmatrix = gsl_matrix_complex_alloc(2*xsize,2*xsize);
        gsl_matrix_complex_set_identity(Lmatrix);
        for(int x=0;x<2*xsize;x++){                        
            gsl_matrix_complex_set(Lmatrix,x,x,gsl_complex_rect(sqrt(abs(1.0/gsl_vector_get(eval,x))),0.0));
        }
        
        gsl_matrix_complex *invbogo = gsl_matrix_complex_alloc(2*xsize,2*xsize);
        
        gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),Lmatrix,orderedevec,gsl_complex_rect(0.0,0.0), hamiltonian);
        gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),hamiltonian,choleskymatrix,gsl_complex_rect(0.0,0.0),invbogo);
        
        gsl_matrix_complex_free (Lmatrix);
        //H=B^dagger.D.B -> the real transformation matrix is the inverse of B
        gsl_matrix_complex* bogo = gsl_matrix_complex_alloc(2*xsize,2*xsize);
        gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),paraunity,invbogo,gsl_complex_rect(0.0,0.0),hamiltonian);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),hamiltonian,paraunity,gsl_complex_rect(0.0,0.0),bogo);
        
        gsl_matrix_complex_free (hamiltonian);
        gsl_matrix_complex_free (choleskymatrix);
        gsl_matrix_complex_free (orderedevec);
        gsl_matrix_complex_free (invbogo);

        //construct Bogoliubov transformation which ensures that if d_i^\dagger=v_{ij}^A c_j^\dagger+v_{ij}^B c_j then d_i=(v_{ij}^A)^\ast c_j+(v_{ij}^B)^\ast c_j^\dagger
        for(int x=0;x<xsize;x++){
            for(int y=0;y<xsize;y++){
                gsl_matrix_complex_set(bogo,x,y+xsize,gsl_complex_conjugate(gsl_matrix_complex_get(bogo,x+xsize,y)));
                gsl_matrix_complex_set(bogo,x+xsize,y+xsize,gsl_complex_conjugate(gsl_matrix_complex_get(bogo,x,y)));
            }
        }
        
        //check paraunity property
        /*gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),bogo,paraunity,gsl_complex_rect(0.0,0.0),helpmatrixX);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),helpmatrixX,bogo,gsl_complex_rect(0.0,0.0),helpmatrixH);
        
        for(int x=0;x<2*xsize;x++){
            for(int y=0;y<2*xsize;y++){
                if(gsl_complex_abs(gsl_matrix_complex_get(helpmatrixH,x,y)) >= epsilon && x!=y){
                    cout << "k= " << k << ": " << '\t' << x << '\t' << y << '\t' 
                    << gsl_complex_abs(gsl_matrix_complex_get(helpmatrixH,x,y))
                        *exp(I*gsl_complex_arg(gsl_matrix_complex_get(helpmatrixH,x,y))) << '\n';
                }
            }
        }*/
        /*double sum = 0.0;
        gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),bogo,paraunity,gsl_complex_rect(0.0,0.0),helpmatrixX);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),helpmatrixX,bogo,gsl_complex_rect(0.0,0.0),helpmatrixH);
        for(int i=0;i<xsize;i++){
            sum += gsl_complex_abs(gsl_matrix_complex_get(helpmatrixH,i,i));
        }
        cout <<  k << '\t' << sum/xsize << '\n';
        */
        //write eigenvalues (diagonalelements of Lmatrix) and bogo matrix into a file
        ostringstream fm;
        fm << "Bogoliubov/Eigenvectors/k" << k << "D" << count << "cJxsize" << xsize << "ysize" << ysize 
	   << "fieldsize" << fieldsize << ".dat";
        ofstream writeevec(fm.str().c_str(),ios_base::binary);//add when using binary format: ,ios_base::binary
        for(int x=0;x<2*xsize;x++){
            for(int y=0;y<2*xsize;y++){
                if(x==y){
                    if(x<xsize){
//                        help = 2*gsl_vector_get(eval,x);
//                        writeeval.write((char*)& help, sizeof(double));
                        writeeval << 2*gsl_vector_get(eval,x) << '\n'; //without binary
                    }
                    else{
//                        help = 2*gsl_vector_get(eval,x);
//                        writeeval.write((char*)& help, sizeof(double));
                        writeeval << 2*gsl_vector_get(eval,x) << '\n'; //without binary
                    }
                        
                }
                    helpcomplex = gsl_matrix_complex_get(bogo,x,y);
                    writeevec.write((char*)& helpcomplex, sizeof(gsl_complex));
                    //writeevec << gsl_matrix_get(bogo,x,y) << '\t'; //without binary
            }
            //writeevec << '\n'; //without binary
        }
        writeevec.close();
        gsl_matrix_complex_free (bogo);
    }
    gsl_vector_free (eval);
    writeeval.close();
    cerr << 0.01*count << '\r';
    }
    gsl_vector_complex_free (helpvector);
    gsl_matrix_complex_free (paraunity); 
    return 0;
}

void intrastripe(double *texture, gsl_matrix_complex *edmatrix, double *D){
   double xprefactora,xprefactorb;
   int site,xnegneighbour,xposneighbour;
   gsl_matrix_complex_set_zero(edmatrix);
   for(int x=0;x<xsize;x++){
      site = x;
      xposneighbour = (x+1);
      xnegneighbour = (x-1);
           if(x==0){
               xprefactora = cos(texture[site]-texture[xposneighbour]+M_PI) 
                            +D[site]*sin(texture[site]-texture[xposneighbour]+M_PI);
               
               gsl_matrix_complex_set(edmatrix,x,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                                  gsl_complex_rect(-0.5*(xprefactora)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x+1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x+1),
                                  gsl_complex_rect(-0.5*(xprefactora)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x,x+1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+1),
                                  gsl_complex_rect(0.25*(xprefactora-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x),
                                  gsl_complex_rect(0.25*(xprefactora-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x+xsize),
                                  gsl_complex_rect(0.25*(xprefactora+1.),0.0)));
               gsl_matrix_complex_set(edmatrix,x,x+1+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+1+xsize),
                                  gsl_complex_rect(0.25*(xprefactora+1.),0.0)));
           }
           else if(x==xsize-1){
               xprefactorb = cos(texture[xnegneighbour]-texture[site]+M_PI) 
                            +D[xnegneighbour]*sin(texture[xnegneighbour]-texture[site]+M_PI);
                
               gsl_matrix_complex_set(edmatrix,x,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                                  gsl_complex_rect(-0.5*(xprefactorb)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x-1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x-1),
                                  gsl_complex_rect(-0.5*(xprefactorb)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x,x-1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x-1),
                                  gsl_complex_rect(0.25*(xprefactorb-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x),
                                  gsl_complex_rect(0.25*(xprefactorb-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x+xsize),
                                  gsl_complex_rect(0.25*(xprefactorb+1.),0.0)));
               gsl_matrix_complex_set(edmatrix,x,x-1+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x-1+xsize),
                                  gsl_complex_rect(0.25*(xprefactorb+1.),0.0)));
           }
           else{
               xprefactora = cos(texture[site]-texture[xposneighbour]+M_PI) 
                            +D[site]*sin(texture[site]-texture[xposneighbour]+M_PI);
               
               gsl_matrix_complex_set(edmatrix,x,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                                  gsl_complex_rect(-0.5*(xprefactora)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x+1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x+1),
                                  gsl_complex_rect(-0.5*(xprefactora)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x,x+1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+1),
                                  gsl_complex_rect(0.25*(xprefactora-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x),
                                  gsl_complex_rect(0.25*(xprefactora-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x+1,x+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x+1,x+xsize),
                                  gsl_complex_rect(0.25*(xprefactora+1.),0.0)));
               gsl_matrix_complex_set(edmatrix,x,x+1+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+1+xsize),
                                  gsl_complex_rect(0.25*(xprefactora+1.),0.0)));
               
               xprefactorb = cos(texture[xnegneighbour]-texture[site]+M_PI) 
                            +D[xnegneighbour]*sin(texture[xnegneighbour]-texture[site]+M_PI);
                
               gsl_matrix_complex_set(edmatrix,x,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                                  gsl_complex_rect(-0.5*(xprefactorb)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x-1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x-1),
                                  gsl_complex_rect(-0.5*(xprefactorb)-epsilon,0.0)));
               gsl_matrix_complex_set(edmatrix,x,x-1,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x-1),
                                  gsl_complex_rect(0.25*(xprefactorb-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x),
                                  gsl_complex_rect(0.25*(xprefactorb-1.0),0.0)));
               gsl_matrix_complex_set(edmatrix,x-1,x+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x-1,x+xsize),
                                  gsl_complex_rect(0.25*(xprefactorb+1.),0.0)));
               gsl_matrix_complex_set(edmatrix,x,x-1+xsize,
                                  gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x-1+xsize),
                                  gsl_complex_rect(0.25*(xprefactorb+1.),0.0)));
           }
   }
   //add complex conjugated elements
    for(int i=0;i<xsize;i++){
        for(int j=0;j<xsize;j++){
            gsl_matrix_complex_set(edmatrix,i+xsize,j+xsize, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j)));
        }
    }
    for(int i=0;i<xsize;i++){
        for(int j=0;j<xsize;j++){
            gsl_matrix_complex_set(edmatrix,j+xsize,i, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j+xsize)));
        }
    }
    gsl_matrix_complex_scale(edmatrix,gsl_complex_rect(-0.5*J*S,0.0));
}

void interstripe(double *texture, gsl_matrix_complex *edmatrix, double fourierexponent, double *D){
   double yprefactora,yprefactorb;
   gsl_matrix_complex_set_zero(edmatrix);
   for(int x=0;x<xsize;x++){
      yprefactora = cos(0.0);
      yprefactorb = cos(0.0);
      gsl_matrix_complex_set(edmatrix,x,x,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                    gsl_complex_rect(-0.5*(yprefactora)-epsilon,0.0)));
      gsl_matrix_complex_set(edmatrix,x,x,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                    gsl_complex_rect(-0.5*(yprefactorb)-epsilon,0.0)));
      gsl_matrix_complex_set(edmatrix,x,x,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                    gsl_complex_polar(0.25*(yprefactora-1.),fourierexponent)));
      gsl_matrix_complex_set(edmatrix,x,x,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x),
                    gsl_complex_polar(0.25*(yprefactorb-1.),fourierexponent)));
      gsl_matrix_complex_set(edmatrix,x,x+xsize,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+xsize),
                    gsl_complex_polar(0.25*(yprefactora+1.),fourierexponent)));
      gsl_matrix_complex_set(edmatrix,x,x+xsize,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,x,x+xsize),
                    gsl_complex_polar(0.25*(yprefactorb+1.),fourierexponent)));
    
    }
    //add complex conjugated elements
    for(int i=0;i<xsize;i++){
        for(int j=0;j<xsize;j++){
            gsl_matrix_complex_set(edmatrix,i+xsize,j+xsize, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j)));
        }
    }
    for(int i=0;i<xsize;i++){
        for(int j=0;j<xsize;j++){
            gsl_matrix_complex_set(edmatrix,j+xsize,i, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j+xsize)));
        }
    }
    gsl_matrix_complex_scale(edmatrix,gsl_complex_rect(-0.5*J*S,0.0));
}
