#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define COLS 20 //number of columns. Use multiples of 3
#define ROWS 20 //number of rows. Use multiples of 3
#define pi 3.14159265359

//Function declarations to establish boundary conditions and geometry
void bound1(double V[ROWS][COLS],int rows, int cols);//Box
void bound2(double V[ROWS][COLS],int rows, int cols);//Parallel lines
void bound3(double V[ROWS][COLS],int rows, int cols);//Box + line
void bound4(double V[ROWS][COLS],int rows, int cols);//Lightning rod

//Declarations for different PDE solvers. Numbers 1-4 correspond to geometries 
void relax1_j(double V[ROWS][COLS],double VNEW[ROWS][COLS],int rows,int cols);
void relax1_g(double V[ROWS][COLS],int rows, int cols);
void relax1_s(double V[ROWS][COLS],int rows, int cols, int n);
void relax2_j(double V[ROWS][COLS],double VNEW[ROWS][COLS],int rows,int cols);
void relax2_g(double V[ROWS][COLS],int rows, int cols);
void relax2_s(double V[ROWS][COLS],int rows, int cols, int n);
void relax3_j(double V[ROWS][COLS],double VNEW[ROWS][COLS],int rows,int cols);
void relax3_g(double V[ROWS][COLS],int rows, int cols);
void relax3_s(double V[ROWS][COLS],int rows, int cols, int n);
void relax4_j(double V[ROWS][COLS],double VNEW[ROWS][COLS],int rows,int cols);
void relax4_g(double V[ROWS][COLS],int rows, int cols);
void relax4_s(double V[ROWS][COLS],int rows, int cols, int n);

//Declarations for function to compute electric field 
void Efield(double V[ROWS][COLS], int rows, int cols);

//Function for printing results 
void print_stuff(double V[ROWS][COLS], int rows, int cols);

int main() {
    double V[ROWS][COLS];
    double VNEW[ROWS][COLS];
    int ii,jj, n = 10; //n is number of iterations, ii and jj are indices

    //Grid initialization -set all points to zero
    for(ii=0;ii<ROWS;ii++){
        for(jj=0;jj<COLS;jj++){
            V[ii][jj]=0.0;
        }
    }
    //Grid initialization -set all points to zero
    for(ii=0;ii<ROWS;ii++){
        for(jj=0;jj<COLS;jj++){
            VNEW[ii][jj]=0.0;
        }
    }
    
    //Set-up GEOMETRY of problem. Uncomment to test different geometries 
    //bound1(V,ROWS,COLS); //Uncomment to test geometry 1 
    bound2(V,ROWS,COLS); //Uncomment to test geometry 2
    //bound3(V,ROWS,COLS); //Uncomment to test geometry 3
    //bound4(V,ROWS,COLS); //Uncomment to test geometry 4 

    //Initialize VNEW with BC for JACOBI. Uncomment for different geometries
    //bound1(VNEW,ROWS,COLS);
    bound2(VNEW,ROWS,COLS);
    //bound3(VNEW,ROWS,COLS);
    //bound4(VNEW,ROWS,COLS);


    //print_stuff(V,ROWS,COLS); //check geometries are set properly

    //Start timer to test computational time 
    clock_t start, end; 
    double cpu_time_used; 
    start = clock();

    //MAIN LOOP. Change upper bound to test N-iterations
    for(ii=0; ii<n; ii++){
      //relax1_g(V, ROWS, COLS); //Gauss Method
      relax2_g(V, ROWS, COLS); //Gauss Method
      //relax3_g(V, ROWS, COLS); //Gauss Method
      //relax4_g(V, ROWS, COLS); //Gauss Method
      //relax1_j(V, VNEW, ROWS, COLS); //Jacobi Method
      //relax2_j(V, VNEW, ROWS, COLS); //Jacobi Method
      //relax3_j(V, VNEW, ROWS, COLS); //Jacobi Method
      //relax4_j(V, VNEW, ROWS, COLS); //Jacobi Method
      //relax1_s(V, ROWS, COLS, n); //SOR method
      //relax2_s(V, ROWS, COLS, n); //SOR method
      //relax3_s(V, ROWS, COLS, n); //SOR method
      //relax4_s(V, ROWS, COLS, n); //SOR method
    }

    print_stuff(V,ROWS,COLS); //For printing results 
    //Efield(V, ROWS, COLS);

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("time elapsed is %f\n", cpu_time_used);
}

//Functions to set geometry
void bound1(double V[ROWS][COLS],int rows, int cols){
  int i, j;
    //Fig. 1: BC for box, just change the -100 above to +100
    for(j = cols/3; j <= 2*cols/3; j++){
      V[rows/3][j] = 100;
      V[2*rows/3][j] = 100;
    }
    for(i = rows/3; i <= 2*rows/3; i++){
      V[i][cols/3]=100.0;
      V[i][2*cols/3]=100.0;
  }
}
void bound2(double V[ROWS][COLS],int rows, int cols){
  int i;
 //Fig. 2: BC for parallel lines 
  for(i = rows/3; i <= 2*rows/3; i++){
    V[i][cols/3]=100.0;
    V[i][2*cols/3]=-100.0;
  }
}
void bound3(double V[ROWS][COLS],int rows, int cols){
  int i, j;
  //Fig. 3: BC line plus box.
  for(j = cols/2; j <= 5*cols/6; j++){
      V[rows/3][j] = 100;
      V[2*rows/3][j] = 100;
    }
  for(i = rows/3; i<=2*rows/3; i++){
    V[i][cols/2]=100.0;
    V[i][5*cols/6]=100.0;
  }
  for(i = rows/4; i<=3*rows/4; i++){
    V[i][cols/4] = -100;
  }
}
void bound4(double V[ROWS][COLS],int rows, int cols){
  int i, j;
  //Fig. 4: BC for lightning rod 
  for(i = rows/3; i <= 9*rows/10; i++){
    V[i][cols/2] = -100;
  }
  for(j = cols/4; j <= 3*rows/4; j++){
    V[rows/4][j] = 100;
  }
}

//Function for printing potential values 
void print_stuff(double V[ROWS][COLS], int rows, int cols){
  int i, j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            printf("%lf ",V[i][j]);
        }
        printf("\n");
    }
}

//GAUSS RELAXATION METHOD
void relax2_g(double V[ROWS][COLS],int rows, int cols){
    int i, j;   
    //Geometry 2: Parallel lines  
    for(i=1;i<=rows/3;i++){
        for(j=1;j<cols-1;j++){
            V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
        }
    }    
    for(i=2*rows/3;i<rows-1;i++){
        for(j=1;j<cols-1;j++){
            V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
        }
    }
    for(i=rows/3+1;i<2*rows/3;i++){
        for(j=1;j<cols-1;j++){
            if(j!=cols/3){
                if(j!=2*cols/3) {
                    V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
                }
            }
        }
    }
}
void relax1_g(double V[ROWS][COLS],int rows, int cols){
    int i, j;
    //Geometry 1: box
    for(i=1; i<rows/3; i++){
      for(j=1; j<cols-1; j++){
	V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
      }
    }
    for(i=2*rows/3+1; i<rows-1; i++){
        for(j=1;j<cols-1;j++){
            V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
        }
    }    
    for(i=rows/3; i<2*rows/3+1; i++){
        for(j=1;j < cols/3;j++){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
        }
	for(j=2*cols/3+1;j < cols-1;j++){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
        }
    }
}
void relax3_g(double V[ROWS][COLS],int rows, int cols){
    int i, j;
    //Geometry 3: box + line 
    for(i=1; i<rows/4; i++){
      for(j=1; j<=cols/2; j++){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=rows/4; i<= 3*rows/4; i++){
      for(j=1; j<cols/2; j++){
	if(j!=cols/4){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
	}
      }
    }
    for(i=3*rows/4+1; i< rows-1; i++){
      for(j=1; j<=cols/2; j++){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=1; i<rows/3; i++){
      for(j=cols/2; j<cols-1; j++){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=2*rows/3+1; i<rows-1;i++){
      for(j=cols/2; j<cols-1; j++){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=rows/3; i<2*rows/3+1; i++){
      for(j=5*cols/6+1; j<cols-1; j++){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
}
void relax4_g(double V[ROWS][COLS],int rows, int cols){
    int i, j;  
    //Geometry 4: lightning rod 
    for(i=1; i<rows/3; i++){
      for(j=cols/4; j<3*cols/4+1; j++){
	if(i != rows/4){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	}
      }
    }
    for(i=rows/3; i<9*rows/10+1; i++){
        for(j=1;j<cols-1;j++){
	  if(j != cols/2){
            V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  }
        }
    }    
    for(i=9*rows/10+1; i<rows-1; i++){
        for(j=1;j < cols-1;j++){
	  V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
        }
    }
    for(i=1; i<rows/3; i++){
      for(j=1; j < cols/4; j++){
	V[i][j] = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
      }
      for(j=3*cols/4+1; j<cols-1; j++){
	V[i][j] = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
      }
    }
}

//JACOBI RELAXATION METHOD
void relax2_j(double V[ROWS][COLS],double VNEW[ROWS][COLS], int rows, int cols){
  int i, j;
    for(i=1;i<rows/3;i++){
        for(j=1;j<cols-1;j++){
            VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
        }
    }    
    for(i=2*rows/3;i<rows-1;i++){
        for(j=1;j<cols-1;j++){
            VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
        }
    }  
    for(i=rows/3;i<2*rows/3;i++){
        for(j=1;j<cols-1;j++){
            if(j!=cols/3){
                if(j!=2*cols/3) {
                    VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
                }
            }
        }
    }   
    for(i=1;i<rows-1;i++){
        for(j=1;j<cols-1;j++){
            V[i][j]=VNEW[i][j];
        }
    }
}
void relax1_j(double V[ROWS][COLS],double VNEW[ROWS][COLS],int rows, int cols){
    int i, j;
    //Geometry 1: box
    for(i=1; i<rows/3; i++){
      for(j=1; j<cols-1; j++){
	VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
      }
    }
    for(i=2*rows/3+1; i<rows-1; i++){
        for(j=1;j<cols-1;j++){
            VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
        }
    }    
    for(i=rows/3; i<2*rows/3+1; i++){
        for(j=1;j < cols/3;j++){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
        }
	for(j=2*cols/3+1;j < cols-1;j++){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
        }
    }
    for(i=1;i<rows-1;i++){
        for(j=1;j<cols-1;j++){
            V[i][j]=VNEW[i][j];
        }
    }
}
void relax3_j(double V[ROWS][COLS],double VNEW[ROWS][COLS],int rows, int cols){
    int i, j;
    //Geometry 3: box + line 
    for(i=1; i<rows/4; i++){
      for(j=1; j<=cols/2; j++){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=rows/4; i<= 3*rows/4; i++){
      for(j=1; j<cols/2; j++){
	if(j!=cols/4){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
	}
      }
    }
    for(i=3*rows/4+1; i< rows-1; i++){
      for(j=1; j<=cols/2; j++){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=1; i<rows/3; i++){
      for(j=cols/2; j<cols-1; j++){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=2*rows/3+1; i<rows-1;i++){
      for(j=cols/2; j<cols-1; j++){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=rows/3; i<2*rows/3+1; i++){
      for(j=5*cols/6+1; j<cols-1; j++){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
      }
    }
    for(i=1;i<rows-1;i++){
        for(j=1;j<cols-1;j++){
            V[i][j]=VNEW[i][j];
        }
    }
}
void relax4_j(double V[ROWS][COLS],double VNEW[ROWS][COLS],int rows, int cols){
    int i, j;  
    //Geometry 4: lightning rod 
    for(i=1; i<rows/3; i++){
      for(j=cols/4; j<3*cols/4+1; j++){
	if(i != rows/4){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	}
      }
    }
    for(i=rows/3; i<9*rows/10+1; i++){
        for(j=1;j<cols-1;j++){
	  if(j != cols/2){
            VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  }
        }
    }    
    for(i=9*rows/10+1; i<rows-1; i++){
        for(j=1;j < cols-1;j++){
	  VNEW[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]); 
        }
    }
    for(i=1; i<rows/3; i++){
      for(j=1; j < cols/4; j++){
	VNEW[i][j] = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
      }
      for(j=3*cols/4+1; j<cols-1; j++){
	VNEW[i][j] = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
      }
    }
    for(i=1;i<rows-1;i++){
        for(j=1;j<cols-1;j++){
            V[i][j]=VNEW[i][j];
        }
    }
}



//GAUSS WITH SUCCESSIVE OVER RELAXATION 
void relax2_s(double V[ROWS][COLS],int rows, int cols, int n){
    int i, j;   
    double Vnew, omega = 2/(1+pi/n);
    //Geometry 2: Parallel lines  
    for(i=1;i<=rows/3;i++){
        for(j=1;j<cols-1;j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
        }
    }    
    for(i=2*rows/3;i<rows-1;i++){
        for(j=1;j<cols-1;j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
        }
    }
    for(i=rows/3+1;i<2*rows/3;i++){
        for(j=1;j<cols-1;j++){
            if(j!=cols/3){
                if(j!=2*cols/3) {
               	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
		  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
                }
            }
        }
    }
}
void relax1_s(double V[ROWS][COLS],int rows, int cols, int n){
    int i, j;
    double Vnew, omega = 2/(1+pi/n);
    //Geometry 1: box
    for(i=1; i<rows/3; i++){
      for(j=1; j<cols-1; j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
      }
    }
    for(i=2*rows/3+1; i<rows-1; i++){
        for(j=1;j<cols-1;j++){
 	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
        }
    }    
    for(i=rows/3; i<2*rows/3+1; i++){
        for(j=1;j < cols/3;j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]); 
        }
	for(j=2*cols/3+1;j < cols-1;j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
        }
    }
}
void relax3_s(double V[ROWS][COLS],int rows, int cols, int n){
    int i, j;
    double Vnew, omega = 2/(1+pi/n);
    //Geometry 3: box + line 
    for(i=1; i<rows/4; i++){
      for(j=1; j<=cols/2; j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
      }
    }
    for(i=rows/4; i<= 3*rows/4; i++){
      for(j=1; j<cols/2; j++){
	if(j!=cols/4){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
	}
      }
    }
    for(i=3*rows/4+1; i< rows-1; i++){
      for(j=1; j<=cols/2; j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]); 
      }
    }
    for(i=1; i<rows/3; i++){
      for(j=cols/2; j<cols-1; j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
      }
    }
    for(i=2*rows/3+1; i<rows-1;i++){
      for(j=cols/2; j<cols-1; j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
      }
    }
    for(i=rows/3; i<2*rows/3+1; i++){
      for(j=5*cols/6+1; j<cols-1; j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]); 
      }
    }
}
void relax4_s(double V[ROWS][COLS],int rows, int cols, int n){
    int i, j;  
    double Vnew, omega = 2/(1+pi/n);
    //Geometry 4: lightning rod 
    for(i=1; i<rows/3; i++){
      for(j=cols/4; j<3*cols/4+1; j++){
	if(i != rows/4){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
	}
      }
    }
    for(i=rows/3; i<9*rows/10+1; i++){
        for(j=1;j<cols-1;j++){
	  if(j != cols/2){
	    Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	    V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
	  }
        }
    }    
    for(i=9*rows/10+1; i<rows-1; i++){
        for(j=1;j < cols-1;j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
        }
    }
    for(i=1; i<rows/3; i++){
      for(j=1; j < cols/4; j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
      }
      for(j=3*cols/4+1; j<cols-1; j++){
	  Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
	  V[i][j]=V[i][j]+omega*(Vnew-V[i][j]);
      }
    }
}


//ELECTRIC FIELD  
void Efield(double V[ROWS][COLS], int rows, int cols){
  int i, j; 
  double Ex[rows][cols], Ey[rows][cols];
  printf("\n"); //several spaces to easily distinguish between data sets
  printf("\n");
  printf("\n");
  printf("E_x: \n");
  for (i = 1; i < rows-1; i++){
    for(j = 1; j < cols-1; j++){
      Ex[i][j] = -(V[i+1][j] - V[i-1][j])/2;
      //Ey[i][j] = -(V[i][j+1] - V[i][j-1])/2;
      printf("%lf ",Ex[i][j]);
    }
    printf("\n");
  }
  printf("\n"); //several spaces to easily distinguish between data sets
  printf("\n");
  printf("\n");
  printf("E_y: \n");
  for (i = 1; i < rows-1; i++){
    for(j = 1; j < cols-1; j++){
      //Ex[i][j] = -(V[i+1][j] - V[i-1][j])/2;
      Ey[i][j] = -(V[i][j+1] - V[i][j-1])/2;
      printf("%lf ",Ey[i][j]);
    }
    printf("\n");
  }
}
