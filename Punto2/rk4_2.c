#include <stdio.h>
#include <math.h>


/* inicializando funciones*/
void runge(double t, double  *y1, double  *y3, double  delta_t);
void runge2(double t, double  *y1, double *y2, double  *y3, double *y4 , double  delta_t);
double f1(double t, double y1,double y3  );
double f2(double t, double y1, double y2 ,double y3 , double y4 );
double f3(double t, double y1,double y3  );
double f4(double t, double y1, double y2 ,double y3 , double y4 );



/*main*/
int main(){
    
  double T=2800;
  double delta_t=0.006;
  int n_step = (int)(T/delta_t);
  double qi1;
  double qi2;
  double qi3;
  double qi4;
    int j;
  
    for(j=0; j<=100;j++)
    {
    qi1=0.353553390593273762200422181052424519642417968844237018294;
    qi2=-0.04+0.05*j;
    qi3=0.0;
    qi4=0.0;

  
  double t=0.0;
  int i=0;
   /* n_step=15;*/
        
  for(i=0;i<=n_step;i++){
      if(qi3 < 0.0005 && qi3> -0.0005){
          printf(" %.10f %.10f %.10f %.10f %.10f \n", t, qi1, qi2 , qi3, qi4);
      }
    
    runge(t, &qi1, &qi3, delta_t);
    runge2(t, &qi1, &qi2, &qi3, &qi4, delta_t);
    t += delta_t;
    

  
  }
        
        
    }
  return 0;
}

/* Definiendo las funciones asociadas al hamiltoniano*/

/*
y1=q1
y2=q3


y3=p1
y4=p3

 */

/* q.1 */
double f1(double t, double y1, double y3 ){
  return y3;
}
/* q.3 */
double f2(double t, double y1, double y2 ,double y3 , double y4){
  return y4;
}
/* p.1 */
double f3(double t, double y1, double y3 ){
  double e=1.0;
  return -2.0*y1/pow((4.0*pow(y1,2)+1),1.5);
}
/* p.3 */
double f4(double t, double y1, double y2 ,double y3 , double y4){
  double e=1.0;
   
  return ((y1-y2)/pow(((y1-y2)*(y1-y2)+0.25),1.5)) - ((y1+y2)/pow(((y1+y2)*(y1+y2)+0.25),1.5));
}


/* Runge Kutta 4th */


void runge(double t, double  *y1, double  *y3 , double  delta_t){
  
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double yn1, yn2, yn3, yn4;
  double c1, c2, c3, c4;

  yn1=*y1;
  yn3=*y3;

    k1 = f1(t,yn1, yn3);
    l1 = f3(t,yn1, yn3);

    k2 = f1(t + delta_t*0.5, yn1  + k1*delta_t*0.5, yn3 + l1*delta_t*0.5 );
    l2 = f3(t + delta_t*0.5, yn1  + k1*delta_t*0.5, yn3 + l1*delta_t*0.5  );

    k3 = f1(t + delta_t*0.5, yn1 + k2*delta_t*0.5,  yn3 + l2*delta_t*0.5);
    l3 = f3(t + delta_t*0.5, yn1 + k2*delta_t*0.5,  yn3 + l2*delta_t*0.5);

    k4 = f1(t + delta_t, yn1 + k3*delta_t, yn3 + l3*delta_t);
    l4 = f3(t + delta_t, yn1 + k3*delta_t, yn3 + l3*delta_t);

    c1 = (1.0/6.0)*(k1+2*k2+2*k3+k4)*delta_t;
    c3 = (1.0/6.0)*(l1+2*l2+2*l3+l4)*delta_t;

    yn1 += c1;
    yn3 += c3;


    /* almacena la infromacion nueva*/
  
  *y1=yn1;
  *y3=yn3;
}

void runge2(double t, double  *y1, double *y2, double  *y3, double *y4 , double  delta_t){
    
    double k1, k2, k3, k4;
    double g1, g2, g3, g4;
    double l1, l2, l3, l4;
    double m1, m2, m3, m4;
    double yn1, yn2, yn3, yn4;
    double c1, c2, c3, c4;
    
    yn1=*y1;
    yn2=*y2;
    yn3=*y3;
    yn4=*y4;
    
 //   k1 = f1(t,yn1, yn3);
    g1 = f2(t,yn1,yn2, yn3, yn4);
 //   l1 = f3(t,yn1, yn3);
    m1 = f4(t,yn1,yn2, yn3, yn4);
    
 //   k2 = f1(t + delta_t*0.5, yn1  + k1*delta_t*0.5, yn3 + l1*delta_t*0.5 );
    g2 = f2(t + delta_t*0.5, yn1  + k1*delta_t*0.5, yn2 + g1*delta_t*0.5, yn3 + l1*delta_t*0.5 ,yn4 + m1*delta_t*0.5);
 //   l2 = f3(t + delta_t*0.5, yn1  + k1*delta_t*0.5, yn3 + l1*delta_t*0.5  );
    m2 = f4(t + delta_t*0.5, yn1  + k1*delta_t*0.5, yn2 + g1*delta_t*0.5, yn3 + l1*delta_t*0.5 ,yn4 + m1*delta_t*0.5);
    
 //   k3 = f1(t + delta_t*0.5, yn1 + k2*delta_t*0.5,  yn3 + l2*delta_t*0.5);
    g3 = f2(t + delta_t*0.5, yn1 + k2*delta_t*0.5, yn2 + g2*delta_t*0.5, yn3 + l2*delta_t*0.5, yn4 + m2*delta_t*0.5);
 //   l3 = f3(t + delta_t*0.5, yn1 + k2*delta_t*0.5,  yn3 + l2*delta_t*0.5);
    m3 = f4(t + delta_t*0.5, yn1 + k2*delta_t*0.5, yn2 + g2*delta_t*0.5, yn3 + l2*delta_t*0.5, yn4 + m2*delta_t*0.5);
    
 //   k4 = f1(t + delta_t, yn1 + k3*delta_t, yn3 + l3*delta_t);
    g4 = f2(t + delta_t, yn1 + k3*delta_t, yn2 + g3*delta_t, yn3 + l3*delta_t, yn4 + m3*delta_t);
 //   l4 = f3(t + delta_t, yn1 + k3*delta_t, yn3 + l3*delta_t);
    m4 = f4(t + delta_t, yn1 + k3*delta_t, yn2 + g3*delta_t, yn3 + l3*delta_t, yn4 + m3*delta_t);
    
 //   c1 = (1.0/6.0)*(k1+2*k2+2*k3+k4)*delta_t;
    c2 = (1.0/6.0)*(g1+2*g2+2*g3+g4)*delta_t;
 //   c3 = (1.0/6.0)*(l1+2*l2+2*l3+l4)*delta_t;
    c4 = (1.0/6.0)*(m1+2*m2+2*m3+m4)*delta_t;
    
    /* printf(" %.10f %.10f \n", yn1, yn2);
     printf("%.10f %.10f %.10f %.10f \n",c1 ,k2, k3 ,k4);
     printf("%.10f %.10f %.10f %.10f \n", g1, g2, g3 ,g4);*/
    
    
    //yn1 += c1;
    yn2 += c2;
    //yn3 += c3;
    yn4 += c4;
    
    /* almacena la infromacion nueva*/
    
    *y1=yn1;
    *y2=yn2;
    *y3=yn3;
    *y4=yn4;
}


