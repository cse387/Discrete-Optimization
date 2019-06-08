#include <math.h>
#include <stdio.h>
#include <stdlib.h>
extern int dim;
double costFun(int *x,double dim[dim][dim]);
double rand_float(double a, double b);
void optchange2(int *x,int *y);
void optchange3(int *x,int *y);
int rand_int(int a, int b);
void set(int t,int *x,int *y,int stop);
int countDist(int t,int count,int*y,int stop);
int findNode(int t,int* x,int d,int stop);
void makePermutation(int d1,int d2,int d3,int *x,int*y);
double intializationOfTemperature(int *x,int *r,int L,double dist[dim][dim]);
int functorial(int x);
int *copy(int *r,int dim);
