#include <stdio.h>
#include <math.h>


void leadfrog_k(double t, double  *y1,  double  *y3,  double  delta_t);
void leadfrog_U(double t, double  *y1, double *y2, double  *y3, double *y4 , double  delta_t );

double f1(double t, double y1 ,double y3 );
double f2(double t, double y1, double y2 ,double y3 , double y4 );
double f3(double t, double y1 ,double y3 );
double f4(double t, double y1, double y2 ,double y3 , double y4 );

int main(){

  double T=2800;
  double delta_t=0.006;
  int n_step = (int)(T/delta_t);
  double q1;
  double q3;
  double p1;
  double p3;
   
    double c1;
    double c2;
    double t;
    
    int j;
    
    for (j=0; j<100; j++) {

  

    t=0.0;
    q1= 0.4325;//0422181052424519642417968844237018294;
    q3= 0.04 + 0.05*j;
    p1= 0.0;
    p3= 0.0;

  

  int i=0;
        n_step=300000;
  for(i=0;i<=n_step;i++){
    if(p1 < 0.0005 && p1> -0.0005){
    printf(" %.10f %.10f %.10f %.10f %.10f \n", t, q1, q3 , p1, p3);
    }

    
      leadfrog_k(t, &q1, &p1, delta_t);
      
      leadfrog_U(t, &q1, &q3, &p1, &p3, delta_t);
      
      t += delta_t;
      
     }
    }
  return 0;
}

/*
y1=q1
y2=q3


y3=p1
y4=p3

 */

/* p.1 */
double f1(double t, double y1,double y3){
  return y3;
}
/* p.3 */
double f2(double t, double y1, double y2 ,double y3 , double y4){
  return y4;
}
/* q.1 */
double f3(double t, double y1 ,double y3){

  return -2.0*y1/pow((4.0*y1*y1+1),1.5);
}
/* q.3 */
double f4(double t, double y1, double y2 ,double y3 , double y4){

  return ((y1-y2)/pow(((y1-y2)*(y1-y2)+0.25),1.5)) - ((y1+y2)/pow(((y1+y2)*(y1+y2)+0.25),1.5));
}


void leadfrog_k(double t, double  *y1,  double  *y3, double  delta_t ){


  double yn1, yn2, yn3, yn4;

  
  yn1=*y1;
  yn3=*y3;
  

    //c1=0.675604;
    yn3 += 0.5*f3(t,yn1,yn3)*delta_t;
    //c2=1.35121;
    yn1+=1.0*f1(t,yn1,yn3)*delta_t;
    //c1=-0.175604;
    yn3 += 0.5*f3(t,yn1,yn3)*delta_t;
//    //c2=-1.70241;
//    yn1 += -1.70241*f1(t,yn1,yn3)*delta_t;
//    //c1=-0.175604;
//    yn3 += -0.175604*f3(t,yn1,yn3)*delta_t;
//    //c2=1.35121;
//    yn1 += 1.35121*f1(t,yn1,yn3)*delta_t;
//    //c1=0.675604;
//    yn3 += 0.675604*f3(t,yn1,yn3)*delta_t;
//    //c2=0.0;
//    yn1 += 0.0*f1(t,yn1,yn3)*delta_t;

    *y1=yn1;
    *y3=yn3;

}


void leadfrog_U(double t, double  *y1, double *y2, double  *y3, double *y4 , double  delta_t){


  double yn1, yn2, yn3, yn4;

  
  yn1=*y1;
  yn2=*y2;
  yn3=*y3;
  yn4=*y4;

    //c1=0.675604;
    yn4 += 0.5*f4(t,yn1,yn2,yn3,yn4)*delta_t;
    //c2=1.35121;
    yn2+= 1.0*f2(t,yn1,yn2,yn3,yn4)*delta_t;
    //c1=-0.175604;
    yn4 += 0.5*f4(t,yn1,yn2,yn3,yn4)*delta_t;
//    //c2=-1.70241;
//    yn2 += -1.70241*f2(t,yn1,yn2,yn3,yn4)*delta_t;
//    //c1=-0.175604;
//    yn4 += -0.175604*f4(t,yn1,yn2,yn3,yn4)*delta_t;
//    //c2=1.35121;
//    yn2 += 1.35121*f2(t,yn1,yn2,yn3,yn4)*delta_t;
//    //c1=0.675604;
//    yn4 += 0.675604*f4(t,yn1,yn2,yn3,yn4)*delta_t;
//    //c2=0.0;
//    yn2 += 0.0*f2(t,yn1,yn2,yn3,yn4)*delta_t;


    
    
    *y1=yn1;
    *y2=yn2;
    *y3=yn3;
    *y4=yn4;
    

}


