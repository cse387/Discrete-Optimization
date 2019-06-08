#include "discreteOpt.h"
double costFun(int *x,double d[dim][dim]){
  double sum=0;
  for(int i=0; i<dim-1; i++){
    sum+=d[x[i]][x[i+1]];
  }
  return sum+d[x[dim-1]][x[0]];
}
double rand_float(double a, double b) {
	return ((double)rand() / RAND_MAX) * (b - a) + a;
}
int rand_int(int a, int b) {
	return rand() % (b - a + 1) + a;
}
int functorial(int x){
  if(x==0) return 1;
  else return functorial(x-1)*x;
}
int *copy(int *x,int dim){
  int *ret=(int*)malloc(sizeof(int)*dim);
  for(int i=0; i<dim; i++){
    ret[i]=x[i];
  }
  return ret;
}
double intializationOfTemperature(int *x,int *r,int L,double dist[dim][dim]){
  double ratio=0;
  double c=0.1;
  double F=costFun(x,dist);
  double Fnew=0;
  int *state=(int*)(malloc(sizeof(int)*dim));
  int *revstate=(int*)(malloc(sizeof(int)*dim));
  int count=0;
  while(ratio<0.8){
      printf("%s\n","STEEEEEEP" );
      printf("ratio==%f\n",ratio );
      printf("count==%i\n",count );
      printf("temp=%f\n",c );
      state=copy(x,dim);
      revstate=copy(r,dim);
      optchange2(x,r);
      Fnew=costFun(x,dist);
      if(Fnew<F){
        F=Fnew; count++;
      }
      else if((exp(F-Fnew)/c)>rand_float(0,1)){
        F=Fnew; count++;
      }
      else{
        x=copy(state,dim);
        r=copy(revstate,dim);
      }
      ratio=count/(double)L;
    c=c*1.1;
  }
  return c;
}
void optchange2(int *x,int *y){
  int d1=rand_int(0,dim-1);
  //printf("%i\n",d1 );
  int stop=rand_int(2,dim-2);
  int d2=findNode(d1,x,0,stop);
//  printf("%i\n",d2 );
  //printf("distance==%d\n",countDist(d2,0,y,d1) );
  int distd1d2=countDist(d1,0,y,d2);
  int distd2d1=countDist(d2,0,y,d1);
//  if(distd1d2<distd2d1){//peritoo
    x[x[d1]]=x[d2];
    set(d2,x,y,x[d1]);
    x[d1]=d2;
  /*}
 else{
    x[x[d2]]=x[d1];
    set(d1,x,y,x[d2]);
    x[d2]=d1;
  }*/
  //printf("distance==%d\n",countDist(d1,0,y,d2) );
  for(int i=0; i<dim; i++){
    y[x[i]]=i;
  }
}
void optchange3(int *x,int *y){
  int d1=rand_int(0,dim-1);
  int stop=rand_int(2,dim-2);
  //printf("stop==%i\n",stop );
  int d2=findNode(d1,x,0,stop);
  int d3=0;//rand_int(0,dim-1);
//  int distd1d2=dim-countDist(d1,0,y,d2);
//  int distd2d1=dim-countDist(d2,0,y,d1);
  while(d3==d1 || d3==y[d1] || d3==x[d1] || d3==d2 ||d3==x[d2]|| d3==y[d2]){
    d3=rand_int(0,dim-1);
  }
  //printf(" d1=%i\n d2=%i\n d3=%i\n ",d1,d2,d3 );

  int dist1=dim-countDist(d1,0,y,d2);
  int dist2=dim-countDist(d1,0,y,d3);
  //printf("%i\n %i\n",dist1,dist2 );

  if(dist1<=dist2){
    makePermutation(d1,d2,d3,x,y);
  }
  else{
    makePermutation(d1,d3,d2,x,y);
  }
    for(int i=0; i<dim; i++){
      y[x[i]]=i;
    }
}
void makePermutation(int d1,int d2,int d3,int *x,int*y){
  int c=rand_int(0,3);
  int t=0;
  //printf("case=%i\n",c );
  switch (c) {
    case  0:
            x[x[d1]]=x[d2];
            set(d2,x,y,x[d1]);
            x[x[d3]]=d2;
            set(d1,x,y,x[d3]);
            x[d3]=d1;
            break;
    case  1:
            x[x[d1]]=d3;
            x[x[d2]]=x[d3];
            set(d3,x,y,x[d2]);
            set(d2,x,y,x[d1]);
            x[d1]=d2;
            break;
    case  2:
            x[x[d1]]=d3;
            x[x[d2]]=d1;
            x[x[d3]]=d2;
            t=x[d3];
            set(d3,x,y,x[d2]);
            set(d2,x,y,x[d1]);
            set(d1,x,y,t);
            break;
    case  3:
            x[x[d1]]=x[d3];
            t=x[d2];
            set(d2,x,y,x[d1]);
            x[d1]=t;
            x[d3]=d2;
            break;
  }
}
int findNode(int t,int* x,int d,int stop){
  if(d==stop) return t;
  else{
    d++;
    return findNode(x[t],x,d,stop);
  }
}
void set(int t,int*x ,int*y,int stop){
  if(t==stop) return;
  //printf("x=%i y=%i t=%i ",x[t],y[t] ,t);
  x[t]=y[t];
  set(y[t],x,y,stop);
}
int countDist(int t,int count,int*y,int stop){
  if(t==stop) return count;
  else{
    count++;
    return countDist(y[t],count,y,stop);
  }
}
