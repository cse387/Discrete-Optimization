#include "ContiniousOpt.h"
/*double ObjFunc(double x[]){
 double p=1;
  double sum=0;
  for(int i=0; i<dim; i++){
    p=p*cos(x[i])/sqrt(i+1);
    sum=sum+ pow(x[i],2);
  }
	sum=sum/200;
	return 1+sum-p;
/*  double ret=4*pow(x[0],2)-2.1*pow(x[0],4)+pow(x[0],6)/3.0
          +x[0]*x[1]-4*pow(x[1],2)+4*pow(x[1],4);*/
    /*double ret=0;
    for(int i=0; i<dim; i++){
      ret+=pow(x[i],2)-cos(18*x[i]);
    }
    return ret;  */
/*    double p1=1+pow(x[0]+x[1]+1,2)*
    (19-14*x[0]+3*pow(x[0],2)-14*x[1]+6*x[0]*x[1]+3*pow(x[1],2));
    double p2=30+pow(2*x[0]-3*x[1],2)*
    (18-32*x[0]+12*pow(x[0],2)+48*x[1]-36*x[0]*x[1]+27*pow(x[1],2));
    return p1*p2;*/
  /*  double p1=pow(sin(3*M_PI*x[0]),2)/10;
    double p2=0;
    for(int i=1; i<dim-1; i++){
      p2+=pow(x[i]-1,2)*(1+pow(sin(3*M_PI*x[i+1]),2));
    }
    double p3=pow(x[dim-1]-1,2)*(1+pow(sin(2*M_PI*x[dim-1]),2));
    return p1*p2+p3;
    double sum=0;
    for(int i=0; i<dim; i++){
      sum+=pow(x[i],4)-16*pow(x[i],2)+5*x[i];
    }
    return sum/2;*/
//}
double ObjFunc(double x[]){
  int m=31;
  double sum=0;
	for(int i=0; i<m-1; i++){
		double a1=(double)i/29.0;//+1;///ATTTeenntio
		double sum1=0;
		for(int j=0; j<=dim-2; j++){
			sum1+=x[j]*(j-1)*pow(a1,j);
		}
		double sum2=0;
		for(int j=0; j<=dim-1; j++){
			sum2+=pow(a1,j)*x[j];
		}
		sum+=pow(sum1-pow(sum2,2)-1,2);
	}
  sum+=pow(x[0],2);
  Fev++;
	return sum;
}
double localObjectFun(double len,double *x,double *P){
  double xnew[dim];
  for(int i=0; i<dim; i++){
    xnew[i]=x[i]+len*P[i];
  }
  return ObjFunc(xnew);
}
double quadInter(double xm,double *x,double *P){
  double a=xm-1,b=xm+1;
  double h=(b-a)/2.0;
  double fa,fm,fb,x1,x2,f1,f2;
  fa=localObjectFun(a,x,P);
  fm=localObjectFun(xm,x,P);
  fb=localObjectFun(b,x,P);
  while(1){
      if(fa<fb){
        xm=a+h; fm=localObjectFun(xm,x,P);
        if(fm>fa){
          h=h/2.0; b=xm; fb=fm;
        }
        else{break;  }
      }
      else{
        xm=b-h; fm=localObjectFun(xm,x,P);
        if(fm>fb){
          h=h/2.0; a=xm; fa=fm;
        }
        else{break;  }
      }
      //printf("%s\n","looping" );
  }
  //double min=(a+xm)/2.0+((fa-fm)*(fm-fb)*(fb-fa))/(2*((xm-b)*fa+(b-a)*fm+(a-xm)*fb));
  double min=(a+b)/2.0-((xm-b)*(fb-fa)*(xm-a))/2.0*((fm-fa)*(b-a)-(fb-fa)*(xm-a));
  double fmi=localObjectFun(min,x,P);
  while(b-a>=0.0004){
    if(xm<min){
      x1=xm; x2=min; f1=fm; f2=fmi;
    }
    else{
      x1=min; x2=xm; f1=fmi; f2=fm;
    }
    if(f1<f2){
      b=x2; xm=x1; fb=f2; fm=f1;
    }
    else{
      a=x1; xm=x2; fa=f1; fm=f2;
    }
    //min=(a+xm)/2+((fa-fm)*(fm-fb)*(fb-fa))/(2*((xm-b)*fa+(b-a)*fm+(a-xm)*fb));
    min=(a+b)/2.0-((xm-b)*(fb-fa)*(xm-a))/2.0*((fm-fa)*(b-a)-(fb-fa)*(xm-a));
    if(min-a>2*(b-min)){
      min=(a+2*b)/3.0;
    }
    else if(b-min>2*(min-a)){
      min=(2*a+b)/3.0;
    }
    fm=localObjectFun(min,x,P);
  }
  double ret1=(b+a)/2.0+(b-a)/2.0;
  double ret2=(b+a)/2.0-(b-a)/2.0;
  if(localObjectFun(ret1,x,P)>localObjectFun(ret2,x,P)){
    return ret2;
  }
  return ret1;
}
double rand_float(double a, double b) {
	return ((double)rand() / RAND_MAX) * (b - a) + a;
}
int rand_int(int a, int b) {
	return rand() % (b - a + 1) + a;
}
void randomPick1(struct Point *A,struct Point *B,int n,int N){
  B[0].fvalue=A[0].fvalue;
  B[0].value=(double*)malloc((dim)*sizeof(double));
  for(int j=0; j<dim; j++){
    B[0].value[j]=A[0].value[j];
  }
  for(int i=1; i<n; i++){
    int random_index=rand_int(1,N-1);
    B[i].fvalue=A[random_index].fvalue;
    B[i].value=(double*)malloc((dim)*sizeof(double));
    for(int j=0; j<dim; j++){
      B[i].value[j]=A[random_index].value[j];
    }
  }
}
void intializationArray(int N,struct Point  *D){
  for(int i=0; i<N; i++){
    D[i].value=(double*)malloc((dim)*sizeof(double));
    for(int j=0; j<dim; j++){
      D[i].value[j]=rand_float(a,b);
    }
    D[i].fvalue=ObjFunc(D[i].value);
  }
  setMinMAX(D,N);
}
void generateTrial(struct Point *G,int size,double *y){
// double *trial=(double*)malloc(dim);
  double C[dim];//={ 0 } ;
  for(int i=0; i<size-1; i++){
    for(int j=0; j<dim; j++){
      C[j]+=G[i].value[j];
    }
  }
  for(int i=0; i<dim; i++){
    C[i]=C[i]/(size-1);
    y[i]=2*C[i]-G[size-1].value[i];
  }
//  return trial;
}
void setRandDirect(double *P,int dim){
  double dist=0;
  for(int i=0; i<dim; i++){
  double rand_val=rand_float(a,b);
  P[i]=rand_val;
  dist+=pow(P[i],2);
  }
  dist=sqrt(dist);
  for(int i=0; i<dim; i++){
    P[i]=P[i]/dist;
  }
}
void setMinMAX(struct Point *D,int size){
  struct Point *min,*max={NULL,NULL};//=NULL,NULL;
  min=D;
  max=D;
  for(int i=0; i<size; i++){
    if((*min).fvalue>D[i].fvalue){
      min=D+i;
    }
    if((*max).fvalue<D[i].fvalue){
      max=D+i;
    }
  }
  int max_index=(max-D);
  int min_index=(min-D);
  struct Point t=D[max_index];
  struct Point m=D[min_index];
  D[max_index]=D[size-1];
  D[min_index]=D[0];
  D[size-1]=t;
  D[0]=m;
}
double intializationOfTemperature(double ratio,int m,double *x_start){
    int m1,m2=0;
    double T=0;
    double fy=0;
    double fx=ObjFunc(x_start);
    double x[dim];//=(double*)malloc((dim)*sizeof(double));
    double avg=0;
    for(int i=0; i<m; i++){
      for(int j=0; j<dim; j++){
      int y=rand_float(a,b);
        x[j]=y;
      }
      fy=ObjFunc(x);
      if(fx-fy<=0){
        m1++;
      }
      else{
        m2++;
        avg+=fx-fy;
      }
    }
    avg=avg/m2;
    double par=m2/(m2*ratio+(1-ratio)*m1);
    T=avg/(log(par));
    return T;
}
double std(double *M,int size,double mean){
  double ret=0;
  for(int i=0; i<size; i++){
    ret+=pow(M[i]-mean,2);
  }
  ret=ret/(size-1);
  ret=sqrt(ret);
  return ret;
}
double power(int*x,int*y){
  return pow(*x,*y);
}
/*template <int n>
void size(double (&x)[n]){
  printf("%d\n",sizeof(double)*n) );
}*/
