/*
 * PDEsolve.c - solve the input 2D PDE
 *
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

void PDEsolve(double* x, double* y, double* zin, double* ztype,double* zrate, double* zout, double deltat, double* constants, mwSize length){/*x is vector of x positions, y is vector of y positions, z is a 3d matrix of heat values, tspan contains t0, deltat, and tf*/
    //if solving wave eqn, set rho==0 and use ztype as historical vector
    mwSize sideLength = (mwSize)sqrt((double)length);
    double north;
    double south;
    double east;
    double west;
    double k = constants[0];
    double c = constants[1];
    double rho = constants[2];
    double deltax = x[1]-x[0];//set the distance between points based on the distance between the first x points
    double deltay = y[sideLength]-y[0];
    int ns;
    int ew;
    //printf("1|%f|",deltay);
    //printf("2|%f\n",deltax);
    double maxzout = 0;
    if(rho==0){
        for(mwSize j=1;j<sideLength-1;j++){// loop through each y value, outermost border should be the boundary conditions, or atleast should not matter
        
        for(mwSize h=1;h<sideLength-1;h++){// loop through each x value
            if(ztype[h+sideLength*j]>0&&ztype[h+sideLength*j]<2){//only perform calculation if the POI is identified as an interior point
            //printf("1loop|%f\n",zin[h+(sideLength*j)]);
            //printf("3x-1xx+1|%f|%f|%f\n",x[h-1+(sideLength*j)],x[h+(sideLength*j)],x[h+1+(sideLength*(j+1))]);
            //printf("4y-1yy+1|%f|%f|%f\n",y[h+(sideLength*(j-1))],y[h+(sideLength*j)],y[h+(sideLength*(j+1))]);
            //printf("5z-1zz+1|%f|%f|%f\n",zin[h+(sideLength*(j-1))],zin[h+(sideLength*j)],zin[h+(sideLength*(j+1))]);
            //printf("6z-1zz+1|%f|%f|%f\n",zin[h-1+(sideLength*(j))],zin[h+(sideLength*j)],zin[h+1+(sideLength*(j))]);
                north = (zin[(h+1)+(sideLength*j)]-zin[h+(sideLength*j)])/deltay;
                south = (-zin[h+(sideLength*j)]+zin[(h-1)+(sideLength*j)])/deltay;
                west = (-zin[h+(sideLength*j)]+zin[h+(sideLength*(j-1))])/deltax;
                east = (zin[h+(sideLength*(j+1))]-zin[h+(sideLength*j)])/deltax;
                zout[h+(sideLength*j)] = 2*zin[h+(sideLength*j)] -zrate[h+(sideLength*j)] + ((deltat*deltat*k)/c)*( ( (north+south) / (deltay) ) + ( (east+west) / (deltax) ) );//combine values
            //printf("2loop|%f|%f|%f|%f|\n",zrate[h+(sideLength*j)],zin[h+1+(sideLength*j)],zin[h+(sideLength*j)],zout[h+(sideLength*j)]);
            //printf("4loop nsew |%f|%f|%f|%f|\n",north,south,east,west);
            //printf("3loop|%f|%f|\n",((deltat*deltat*k)/(c)),( ( (north+south) / (deltay) ) + ( (east+west) / (deltax) ) ));
            
            if(zout[h+sideLength*j]>maxzout){
            maxzout=zout[h+sideLength*j];
            }
          }
        }
    }
    }
    else{//set last constant to 0 to use wave equation
    for(mwSize j=1;j<sideLength-1;j++){// loop through each y value, outermost border should be the boundary conditions, or atleast should not matter
        
        for(mwSize h=1;h<sideLength-1;h++){// loop through each x value
            if(ztype[h+sideLength*j]>0&&ztype[h+sideLength*j]<2){//only perform calculation if the POI is identified as an interior point
            //printf("1loop|%f\n",zin[h+(sideLength*j)]);
            //printf("3x-1xx+1|%f|%f|%f\n",x[h-1+(sideLength*j)],x[h+(sideLength*j)],x[h+1+(sideLength*(j+1))]);
            //printf("4y-1yy+1|%f|%f|%f\n",y[h+(sideLength*(j-1))],y[h+(sideLength*j)],y[h+(sideLength*(j+1))]);
            //printf("5z-1zz+1|%f|%f|%f\n",zin[h+(sideLength*(j-1))],zin[h+(sideLength*j)],zin[h+(sideLength*(j+1))]);
            //printf("6z-1zz+1|%f|%f|%f\n",zin[h-1+(sideLength*(j))],zin[h+(sideLength*j)],zin[h+1+(sideLength*(j))]);
                ns = 1;//set markers for neumann conditions
                ew = 1;
                if(ztype[h+1+sideLength*j]>2){//if the boundary is neumann
                    north = zrate[h+1+sideLength*j];
                    ns = 0;//marker that means to ignore the south calculation as per neumann implementation
                }
                else{//if the boundary is dirichlet
                    north = (k/(c*rho))*(zin[(h+1)+(sideLength*j)]-zin[h+(sideLength*j)])/(deltax*deltax);
                }
                if(ztype[h-1+sideLength*j]>2){
                    south = zrate[h-1+sideLength*j];
                    north = 0;//enforces neumann condition
                }
                else if(ns==0){
                    south = 0;//enforces neumann condition
                }
                else{
                    south = (k/(c*rho))*(-zin[h+(sideLength*j)]+zin[(h-1)+(sideLength*j)])/(deltax*deltax);
                }
                if(ztype[h+sideLength*(j-1)]>2){
                    west = zrate[h+sideLength*(j-1)];
                    ew = 0;
                }
                else{
                    west = (k/(c*rho))*(-zin[h+(sideLength*j)]+zin[h+(sideLength*(j-1))])/(deltax*deltax);
                }
                if(ztype[h+sideLength*(j+1)]>2){
                    east = zrate[h+sideLength*(j+1)];
                    west = 0;
                }
                else if(ew==0){
                    east=0;
                }
                else{
                    east = (k/(c*rho))*(zin[h+(sideLength*(j+1))]-zin[h+(sideLength*j)])/(deltax*deltax);
                }
                zout[h+(sideLength*j)] = zin[h+(sideLength*j)] + (deltat*(north+south+east+west));//combine values
             //printf("2loop|%f\n",zout[h+(sideLength*j)]);
             //   printf("4loop nsew |%f|%f|%f|%f|\n",north,south,east,west);
            //printf("3loop|%f|%f|\n",((deltat*k)/(c*rho)),( ( (north+south) / (deltay*deltay) ) + ( (east+west) / (deltax*deltax) ) ));
            
            if(zout[h+sideLength*j]>maxzout){
            maxzout=zout[h+sideLength*j];
            }
          }
        }
    }

    
    }
    //printf("7maxzout|%f\n",maxzout);
    //printf("8|%f\n",deltat);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) /* cite this as a butchered MATLAB tutorial file */
{
    double deltat;              /* input scalar */
    double *zin;               /* Nx1 input matrix */
    double *xin;               /* Nx1 input matrix */
    double *yin;               /* Nx1 input matrix */
    double *ztype;
    double *zrate;
    double *constants;               /* 1x3 input matrix */
    size_t ncols;                   /* size of matrix */
    double *zout;              /* Nx1 output matrix */

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[2]);
    //separate and set the input matricies
    xin = mxGetPr(prhs[0]);
    yin = mxGetPr(prhs[1]);
    zin = mxGetPr(prhs[2]);
    ztype = mxGetPr(prhs[3]);
    zrate = mxGetPr(prhs[4]);//rate condition for BC
    deltat = mxGetScalar(prhs[5]);
    constants = mxGetPr(prhs[6]);
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);
    zout = mxGetPr(plhs[0]);//set pointer of output to zout so it can be passed to the function
    /* call the computational routine */
    PDEsolve(xin,yin,zin,ztype,zrate,zout,deltat,constants,(mwSize)ncols);

}