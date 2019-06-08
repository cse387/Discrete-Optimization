#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern int dim,a,b,Fev;//se mia synarthsh
struct Point{
  double *value;//[dim];
  double fvalue;
};
double ObjFunc(double x[]);
double rand_float(double a, double b);
int rand_int(int a, int b);
void randomPick1(struct Point *A,struct Point *B,int n,int N);
double intializationOfTemperature(double ratio,int m,double* x_start);
void intializationArray(int N,struct Point  *D);
void generateTrial(struct Point *G,int size,double *y);
void setMinMAX(struct Point *D,int size);
double std(double *M,int size,double mean);
void setRandDirect(double *P,int dim);
double power(int*x,int*y);
double quadInter(double len,double *x,double *P);
double localObjectFun(double m,double *x,double *P);
/* data */
/*template <int n>;
void size(double (&x)[n]);
*/
