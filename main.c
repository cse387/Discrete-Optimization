#include "discreteOpt.h"
#include <omp.h>
int dim=8;
int main(){
  int *x=(int*)malloc(sizeof(int)*dim);
  int *r=(int*)malloc(sizeof(int)*dim);
  
  for(int i=0; i<dim; i++){
    x[i]=(i+1)%dim;
    r[x[i]]=i;
  }
  double dist[dim][dim];
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      if(i==j){
        dist[i][j]==0;
      }
      else{
        if(i>j){
          dist[i][j]=rand_float(50,1000);
          dist[j][i]=dist[i][j];
        }
      }
    }
  }
  
  int L=functorial(dim-1);
  int Lmin=L/(double)2;
  double F=costFun(x,dist);
  double min=F;
  
  double T=0;
  int iter=0;
  
  T=intializationOfTemperature(x,r,5,dist);
  printf("Temp=%f\n",T );
  while (1) {
    int accept=0;
    double Fnew=0;
    int *state;
  	int *revstate;
    for(int i=0; i<L; i++){
      state=copy(x,dim);
      revstate=copy(r,dim);
      optchange2(x,r);
      Fnew=costFun(x,dist);
      if(Fnew<F){
        F=Fnew; 
        accept++;
      }
      else if(exp(F-Fnew)/T>rand_float(0,1)){
        F=Fnew;
        accept++;
      }
      else{
        x=copy(state,dim);
        r=copy(revstate,dim);
      }
    }
    if(F<min){
        min=F;
      }
    L=Lmin-accept;
    T=0.99*T;
    iter++;
    printf("Temperature==%f\n",T );
    printf("Cost==%f\n",min );
    if(T==0.0 || iter==3500){
      break;
    }
  }
  printf("initial=%f\n",initial );
  printf("min=%f\n",min );

  return 0;
}
