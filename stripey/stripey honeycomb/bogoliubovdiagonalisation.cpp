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

void intrastripe(gsl_matrix_complex*);
void interstripe(gsl_matrix_complex*,double);
void DMI(gsl_matrix_complex*,double);

typedef complex<double> dcomp;

int main(){
    double help, fourierexponent;
    dcomp I;
    I = -1;
    I = sqrt(I);
    int helporder;
    gsl_vector_complex *helpvector = gsl_vector_complex_alloc(2*N);
    gsl_matrix_complex* paraunity = gsl_matrix_complex_alloc(2*N,2*N);
    gsl_matrix_complex* helpmatrixX = gsl_matrix_complex_alloc(2*N,2*N);
    gsl_matrix* helpmatrixY = gsl_matrix_alloc(2*N,2*N);
    gsl_matrix_complex* helpmatrixH = gsl_matrix_complex_alloc(2*N,2*N);
    gsl_matrix_complex_set_identity(paraunity);
    
    ostringstream fn;
    fn << "eigenvalues.dat";
    ofstream writeeval(fn.str().c_str());//add when using binary format: ,ios_base::binary
    for(int x=0;x<N;x++){
        gsl_matrix_complex_set(paraunity,x+N,x+N,gsl_complex_rect(-1.0,0.0));
    }
    for(int k=0;k<n;k++){
      fourierexponent = 2*M_PI*(k*1.0/(n));
      gsl_matrix* diagonal = gsl_matrix_alloc(2*N,2*N);
    
      //Calculate Hamiltonian directly in programm
      gsl_matrix_complex* hamiltonian = gsl_matrix_complex_alloc(2*N,2*N);
      gsl_matrix_complex_set_zero(hamiltonian);
   
      intrastripe(helpmatrixH);
      gsl_matrix_complex_add(hamiltonian,helpmatrixH);
      interstripe(helpmatrixH,fourierexponent);
      gsl_matrix_complex_add(hamiltonian,helpmatrixH);
      interstripe(helpmatrixH,-fourierexponent);
      gsl_matrix_complex_add(hamiltonian,helpmatrixH);
      DMI(helpmatrixH,fourierexponent);
      gsl_matrix_complex_add(hamiltonian,helpmatrixH);
      if(k==0){
  	for(int i=0;i<2*N;i++){
            for(int j=0;j<2*N;j++){
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
      gsl_matrix_complex* choleskymatrix = gsl_matrix_complex_alloc(2*N,2*N);
      for(int i=0;i<2*N;i++){
          for(int j=0;j<2*N;j++){
              if(j<=i){
                  gsl_matrix_complex_set(choleskymatrix,i,j,gsl_matrix_complex_get(hamiltonian,i,j));
              }
              else{
                gsl_matrix_complex_set(choleskymatrix,i,j,gsl_complex_rect(0.0,0.0));
              }
          }
        }//this yields lower (!) triangular cholesky matrix
        
        gsl_matrix_complex* paraunitcholsesky = gsl_matrix_complex_alloc(2*N,2*N);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,gsl_complex_rect(1.0,0.0),paraunity,choleskymatrix,gsl_complex_rect(0.0,0.0), paraunitcholsesky);
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans,gsl_complex_rect(1.0,0.0),choleskymatrix,paraunitcholsesky,gsl_complex_rect(0.0,0.0), hamiltonian);
        gsl_matrix_complex_free (paraunitcholsesky);
        
        //diagonalise
        gsl_vector *eval = gsl_vector_alloc (2*N);
        gsl_matrix_complex *evec = gsl_matrix_complex_alloc (2*N, 2*N);
        gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc (2*N);
        gsl_eigen_hermv(hamiltonian, eval, evec, w);
        gsl_eigen_hermv_free(w);
        int index = 0;
        int maximumindex;
        int order[2*N];
        for(int i=0;i<2*N;i++){
            order[i] = i;
        }
        while(index<(2*N-1)){
            maximumindex = index;
            for(int i=index+1;i<2*N;i++){
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

        gsl_matrix_complex *orderedevec = gsl_matrix_complex_alloc(2*N,2*N);
        for(int i=0;i<2*N;i++){
            gsl_matrix_complex_get_col(helpvector,evec,order[i]);
            gsl_matrix_complex_set_col(orderedevec,i,helpvector);
        }
        gsl_matrix_complex_free (evec);
                
        //calculate Bogoliubov-transformation
        gsl_matrix_complex *Lmatrix = gsl_matrix_complex_alloc(2*N,2*N);
        gsl_matrix_complex_set_identity(Lmatrix);
        for(int x=0;x<2*N;x++){                        
            gsl_matrix_complex_set(Lmatrix,x,x,gsl_complex_rect(sqrt(abs(1.0/gsl_vector_get(eval,x))),0.0));
        }
        
        gsl_matrix_complex *invbogo = gsl_matrix_complex_alloc(2*N,2*N);
        
        gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),Lmatrix,orderedevec,gsl_complex_rect(0.0,0.0), hamiltonian);
        gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),hamiltonian,choleskymatrix,gsl_complex_rect(0.0,0.0),invbogo);
        
        gsl_matrix_complex_free (Lmatrix);
        //H=B^dagger.D.B -> the real transformation matrix is the inverse of B
        gsl_matrix_complex* bogo = gsl_matrix_complex_alloc(2*N,2*N);
        gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),paraunity,invbogo,gsl_complex_rect(0.0,0.0),hamiltonian);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),hamiltonian,paraunity,gsl_complex_rect(0.0,0.0),bogo);
        
        gsl_matrix_complex_free (hamiltonian);
        gsl_matrix_complex_free (choleskymatrix);
        gsl_matrix_complex_free (orderedevec);
        gsl_matrix_complex_free (invbogo);

        //construct Bogoliubov transformation which ensures that if d_i^\dagger=v_{ij}^A c_j^\dagger+v_{ij}^B c_j then d_i=(v_{ij}^A)^\ast c_j+(v_{ij}^B)^\ast c_j^\dagger
        for(int x=0;x<N;x++){
            for(int y=0;y<N;y++){
                gsl_matrix_complex_set(bogo,x,y+N,gsl_complex_conjugate(gsl_matrix_complex_get(bogo,x+N,y)));
                gsl_matrix_complex_set(bogo,x+N,y+N,gsl_complex_conjugate(gsl_matrix_complex_get(bogo,x,y)));
            }
        }
        
        //check paraunity property
        /*gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),bogo,paraunity,gsl_complex_rect(0.0,0.0),helpmatrixX);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),helpmatrixX,bogo,gsl_complex_rect(0.0,0.0),helpmatrixH);
        
        for(int x=0;x<2*N;x++){
            for(int y=0;y<2*N;y++){
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
        for(int i=0;i<N;i++){
            sum += gsl_complex_abs(gsl_matrix_complex_get(helpmatrixH,i,i));
        }
        cout <<  k << '\t' << sum/N << '\n';
        */
        //write eigenvalues (diagonalelements of Lmatrix) and bogo matrix into a file
        ostringstream fm;
        fm << "eigenvector/k=" << k << ".dat";
        ofstream writeevec(fm.str().c_str());//add when using binary format: ,ios_base::binary
        for(int x=0;x<2*N;x++){
            for(int y=0;y<2*N;y++){
                if(x==y){
                    if(x<N){
                        writeeval << fourierexponent << '\t' << gsl_vector_get(eval,x) << '\n'; //without binary
                    }
                    else{
                        writeeval << fourierexponent << '\t' << gsl_vector_get(eval,x) << '\n'; //without binary
                    }
                }
                    writeevec << gsl_complex_abs(gsl_matrix_complex_get(bogo,x,y))*exp(I*gsl_complex_arg(gsl_matrix_complex_get(bogo,x,y))) << '\t'; //without binary
            }
            writeevec << '\n'; //without binary
        }
        writeevec.close();
        gsl_vector_free (eval);
        gsl_matrix_complex_free (bogo);
    }
    writeeval.close();
    
    
    gsl_vector_complex_free (helpvector);
    gsl_matrix_complex_free (paraunity); 
    
    return 0;
}

void intrastripe(gsl_matrix_complex *edmatrix){
    gsl_matrix_complex_set_zero(edmatrix);
    double xprefactora,xprefactorb;
    for(int i=0;i<N;i++){
        //ATTENTION diagonal terms have epsilon addent!!!
        if(i==0){
            xprefactora = sin(asin(B))*sin(asin(B));
            gsl_matrix_complex_set(edmatrix,i,i,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i),
                    gsl_complex_rect(0.5*t+epsilon,0.0)));
            gsl_matrix_complex_set(edmatrix,i+1,i+1,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i+1),
                    gsl_complex_rect(0.5*t+epsilon,0.0)));
            gsl_matrix_complex_set(edmatrix,i,i+1,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i+1),
                    gsl_complex_rect(-0.5*t*xprefactora,0.0)));
            gsl_matrix_complex_set(edmatrix,i+1,i,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i),
                    gsl_complex_rect(-0.5*t*xprefactora,0.0)));
            gsl_matrix_complex_set(edmatrix,i+1,i+N,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i+N),
                    gsl_complex_rect(0.5*t*(1.-xprefactora),0.0)));
            gsl_matrix_complex_set(edmatrix,i,i+1+N,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i+1+N),
                    gsl_complex_rect(0.5*t*(1.-xprefactora),0.0)));
        }
        
        else if(i==N-1){}
        
        else{
            xprefactora = sin(asin(B))*sin(asin(B));
            gsl_matrix_complex_set(edmatrix,i,i,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i),
                    gsl_complex_rect(0.5*t+epsilon,0.0)));
            gsl_matrix_complex_set(edmatrix,i+1,i+1,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i+1),
                    gsl_complex_rect(0.5*t+epsilon,0.0)));
            gsl_matrix_complex_set(edmatrix,i,i+1,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i+1),
                    gsl_complex_rect(-0.5*t*xprefactora,0.0)));
            gsl_matrix_complex_set(edmatrix,i+1,i,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i),
                    gsl_complex_rect(-0.5*t*xprefactora,0.0)));
            gsl_matrix_complex_set(edmatrix,i+1,i+N,
                    gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i+N),
                    gsl_complex_rect(0.5*t*(1.-xprefactora),0.0)));
            gsl_matrix_complex_set(edmatrix,i,i+1+N,
                gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i+1+N),
                gsl_complex_rect(0.5*t*(1.-xprefactora),0.0)));
        }
    }
    //add complex conjugated elements
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            gsl_matrix_complex_set(edmatrix,i+N,j+N,gsl_matrix_complex_get(edmatrix,i,j));
        }
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            gsl_matrix_complex_set(edmatrix,j+N,i, gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j+N)));
        }
    }
    gsl_matrix_complex_scale(edmatrix,gsl_complex_rect(S,0.0));
}

void interstripe(gsl_matrix_complex *edmatrix, double fourierexponent){
    gsl_matrix_complex_set_zero(edmatrix);
    double yprefactora;
    yprefactora = sin(asin(B))*sin(asin(B));
    for(int i=0;i<N;i++){
        //ATTENTION diagonal terms have epsilon addent!!!
        gsl_matrix_complex_set(edmatrix,i,i,
                        gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i),
                        gsl_complex_rect(0.25*tprime+epsilon,0.0)));
    }
    for(int i=0;i<N;i+=2){
        gsl_matrix_complex_set(edmatrix,i,i+1,
                        gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i+1),
                        gsl_complex_polar(-0.25*tprime*(yprefactora),sqrt(3.0)/(2.0)*fourierexponent)));
        gsl_matrix_complex_set(edmatrix,i+1,i,
                        gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i),
                        gsl_complex_polar(-0.25*tprime*(yprefactora),-sqrt(3.0)/(2.0)*fourierexponent)));
        gsl_matrix_complex_set(edmatrix,i,i+1+N,
                        gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i+1+N),
                        gsl_complex_polar(0.25*tprime*(1.-yprefactora),sqrt(3.0)/(2.0)*fourierexponent)));
        gsl_matrix_complex_set(edmatrix,i+1,i+N,
                        gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i+N),
                        gsl_complex_polar(0.25*tprime*(1.-yprefactora),-sqrt(3.0)/(2.0)*fourierexponent)));
    }
    //add complex conjugated elements
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            gsl_matrix_complex_set(edmatrix,i+N,j+N,
                        gsl_matrix_complex_get(edmatrix,i,j));
        }
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            gsl_matrix_complex_set(edmatrix,j+N,i,
                        gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j+N)));
        }
    }
    gsl_matrix_complex_scale(edmatrix,gsl_complex_rect(S,0.0));
}

void DMI(gsl_matrix_complex *edmatrix, double fourierexponent){
    gsl_matrix_complex_set_zero(edmatrix);
    double xprefactora;
    xprefactora = sin(asin(B))*sin(asin(B));
    for(int i=0;i<N;i++){
        gsl_matrix_complex_set(edmatrix,i,i,
            gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i),
            gsl_complex_polar(-0.5*D*sqrt(xprefactora),0.5*fourierexponent+0.5*M_PI)));
        gsl_matrix_complex_set(edmatrix,i+1,i+1,
            gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i+1),
            gsl_complex_polar(0.5*D*sqrt(xprefactora),-0.5*fourierexponent+0.5*M_PI)));
        gsl_matrix_complex_set(edmatrix,i,i,
            gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i),
            gsl_complex_polar(0.5*D*sqrt(xprefactora),-0.5*fourierexponent+0.5*M_PI)));
        gsl_matrix_complex_set(edmatrix,i+1,i+1,
            gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+1,i+1),
            gsl_complex_polar(-0.5*D*sqrt(xprefactora),0.5*fourierexponent+0.5*M_PI)));
        if(i<N-2){
            gsl_matrix_complex_set(edmatrix,i,i+2,
                gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i+2),
                gsl_complex_rect(0.0,-D*(0.5-(i%2))*sqrt(xprefactora))));
            gsl_matrix_complex_set(edmatrix,i+2,i,
                gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+2,i),
                gsl_complex_rect(0.0,D*(0.5-(i%2))*sqrt(xprefactora))));
            gsl_matrix_complex_set(edmatrix,i,i+2,
                gsl_complex_add(gsl_matrix_complex_get(edmatrix,i,i+2),
                gsl_complex_polar(-D*(0.5-(i%2))*sqrt(xprefactora),sqrt(3.0)/(2.0)*fourierexponent+0.5*M_PI)));
            gsl_matrix_complex_set(edmatrix,i+2,i,
                gsl_complex_add(gsl_matrix_complex_get(edmatrix,i+2,i),
                gsl_complex_polar(D*(0.5-(i%2))*sqrt(xprefactora),sqrt(3.0)/(2.0)*fourierexponent+0.5*M_PI)));
        }
    }
    //add complex conjugated elements
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            gsl_matrix_complex_set(edmatrix,i+N,j+N,
                gsl_matrix_complex_get(edmatrix,i,j));
        }
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            gsl_matrix_complex_set(edmatrix,j+N,i,
                gsl_complex_conjugate(gsl_matrix_complex_get(edmatrix,i,j+N)));
        }
    }
    gsl_matrix_complex_scale(edmatrix,gsl_complex_rect(S,0.0));
}
