#include "ContinuousOpt.h"
#define M_PI 3.14159265358979323846
int dim=6;
int a=-5;//[a,b] search interval
int b=5;
int Fev;
int main(){
  int N=600;
  int n=3;
  int stop=0;
  int iter=0;
  double rate=2.9;
  int L=10*dim;
  double *x=(double*)malloc(sizeof(double)*dim);//[dim];
  double len=1;
  for(int i=0; i<dim; i++){
    int random=rand_float(a,b);
    x[i]=random;
    len=len/(b-a);
  }
  //double *T=(double*)malloc(1);
  double *Fmean=(double*)malloc(sizeof(double)*1);

  double T=intializationOfTemperature(0.9,5,x);
  double Tprev=0;
  struct Point *A=(struct Point*)malloc((N)*sizeof(struct Point));
  struct Point *G=(struct Point*)malloc((n+1)*sizeof(struct Point));
  double *MC=(double*)malloc(sizeof(double)*L);
  double *y=(double*)malloc(sizeof(double)*dim);
  double *P=(double*)malloc(sizeof(double)*dim);
  intializationArray(N,A);
  while(stop!=1){
    double mean=0;
    for(int i=0; i<L; i++){
      int accept=0;
      double q=rand_float(0,1);///buug
      if(q>=0.99){
      randomPick1(A,G,n+1,N);
      generateTrial(G,n+1,y);
      }
      else{
      setRandDirect(P,dim);
      len=quadInter(len,x,P);
    //  printf(" ");
    //  printf("len==%f\n",len );
        for(int i=0; i<dim; i++){
        //  P[i]=rand_float(a,b);
          y[i]=x[i]+len*P[i];
        }
      }

      double fy=ObjFunc(y);//
      double prob=rand_float(0,1);
      double D=fy-A[N-1].fvalue;
      if(D<=0){
        accept=1;
      }
      else if(exp(-D/T)>prob){
        accept=1;
      }
      if(accept=1){
        A[N-1].fvalue=fy;//difff wayyy;;; pick randim in 1 N-1
        for(int i=0; i<dim; i++){
          A[N-1].value[i]=y[i];
          x[i]=y[i];
        }
        setMinMAX(A,N);//paralliaze it
      }
      double fx=ObjFunc(x);
      MC[i]=fx;
      mean+=fx;
    }
    Fmean=realloc(Fmean,sizeof(double)*(iter+1));
    Fmean[iter]=mean/L;
    if(iter>100){
    double stoping=(Fmean[iter]-Fmean[iter-100])*(T)/((T-Tprev)*Fmean[0]);
    if((fabs(stoping)<0.0004) || (T==0)){
      stop=1;
      }
    }
    //T=realloc(T,sizeof(double)*(iter+2));
    Tprev=T;
    T =T*pow(1+(T*log(1+iter))/(3*std(MC,L,mean)),-1);
    printf("fvalue==%f\n iter==%i\n T=%f\n ",A[0].fvalue,iter,T);
    iter++;
  }

  return 0;
}
