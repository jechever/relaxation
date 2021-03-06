#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.14159265359
#define SIZE 20  //size of nxn matrix for point charge
//Shortcomming: SIZE N defines an NxN matrix, thus restricting us to work with square matrices for conveninence purposes 

//Declarations for different PDE solvers. Numbers 1-4 correspond to geometries
void relax_point(double M[SIZE][SIZE][SIZE],double q[SIZE][SIZE][SIZE], int size, double Mnew, double omega, float h);
void relax_dipole(double M[SIZE][SIZE][SIZE],double q[SIZE][SIZE][SIZE], int size, double Mnew, double omega, float h);
//Electric field function
void Efield(double V[SIZE][SIZE][SIZE], int size);

int main() {
    double   M[SIZE][SIZE][SIZE], q[SIZE][SIZE][SIZE],Mnew;
    float h=0.1;
    int ii,jj,kk, n = 8000; //n number of iterations, ii and jj are indices
    double omega = 2/(1+pi/n);
    
    //Initialize grid 
    for(ii=0;ii<SIZE;ii++){
        for(jj=0;jj<SIZE;jj++){
            for(kk=0;kk<SIZE;kk++){
                q[ii][jj][kk]=0.0;
                M[ii][jj][kk]=0.0;
            }
            
        }
    }
    
    //Main loop to test MONOPOLE AND DIPOLE. Uncomment appropriately
    for(ii=0;ii<n;ii++){
      relax_point(M, q, SIZE, Mnew, omega, h);
      //relax_dipole(M, q, SIZE, Mnew, omega, h);
    }

    for(ii=0; ii<SIZE; ii++){
        for(jj=0; jj<SIZE; jj++){
            printf("%lf ", M[ii][jj][SIZE/2]);
        }
        printf("\n");
    }
    printf("\n");

    Efield(M, SIZE);
}

void relax_point(double M[SIZE][SIZE][SIZE],double q[SIZE][SIZE][SIZE], int size, double Mnew, double omega,float h){
    int i, j, k;
    q[size/2][size/2][size/2]=1.0; 
    Mnew = 0;   
    for(i=1; i<size-1; i++){
        for(j=1; j<size-1; j++){
            for(k=1; k<size-1; k++){
	      //if(M[i][j][k] != M[size/2][size/2][size/2]){
	      Mnew = (1/6.0)*(M[i+1][j][k]+M[i-1][j][k]+M[i][j+1][k]+M[i][j-1][k]+M[i][j][k+1]+M[i][j][k-1]+h*h*q[i][j][k]);
	      M[i][j][k]=M[i][j][k]+omega*(Mnew-M[i][j][k]);
		//  }
            }
        }
    }
    
}

void relax_dipole(double M[SIZE][SIZE][SIZE],double q[SIZE][SIZE][SIZE], int size, double Mnew, double omega,float h){
    int i, j, k;
    q[size/4][size/2][size/2]=1.0; 
    q[3*size/4][size/2][size/2]=-1.0; 
    Mnew = 0;   
    for(i=1; i<size-1; i++){
        for(j=1; j<size-1; j++){
            for(k=1; k<size-1; k++){
	      //if(M[i][j][k] != M[size/2][size/2][size/2]){
	      Mnew = (1/6.0)*(M[i+1][j][k]+M[i-1][j][k]+M[i][j+1][k]+M[i][j-1][k]+M[i][j][k+1]+M[i][j][k-1]+h*h*q[i][j][k]);
	      M[i][j][k]=M[i][j][k]+omega*(Mnew-M[i][j][k]);
		//  }
            }
        }
    }
    
}

//ELECTRIC FIELD  
void Efield(double V[SIZE][SIZE][SIZE], int size){
    int i, j; 
    double Ex[size][size], Ey[size][size];
    printf("\n"); //several spaces to easily distinguish between data sets
    printf("\n");
    printf("\n");
    printf("E_x: \n");
    for (i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            Ex[i][j] = -(V[i+1][j][size/2] - V[i-1][j][size/2])/2;
            //Ey[i][j] = -(V[i][j+1] - V[i][j-1])/2;
            printf("%lf ",Ex[i][j]);
        }
        printf("\n");
    }
    printf("\n"); //several spaces to easily distinguish between data sets
    printf("\n");
    printf("\n");
    printf("E_y: \n");
    for (i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            //Ex[i][j] = -(V[i+1][j] - V[i-1][j])/2;
            Ey[i][j] = -(V[i][j+1][size/2] - V[i][j-1][size/2])/2;
            printf("%lf ",Ey[i][j]);
        }
        printf("\n");
    }
}
