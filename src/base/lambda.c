/*------------------------------------------------------------------------------
* lambda.c:integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/
#include "base.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000*100           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info) trace(2,"LD factorization error\n");
        //fprintf(stderr,"%s:LD factorization error\n",__FILE__);
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;
    
    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1;          /* current level */
    dist[k]=0.0;    /* f(z) without level k */
    zb[k]=zs[k];    /* zs: float ambiguities, zb: defined ambiguities, ref.[2] (17-1)*/
    z[k]=ROUND(zb[k]); /* candidate integer ambiguity */
    y=zb[k]-z[k];   /* (z_n-\bar z_n) ref.[2] (18) molecule in last item */
    step[k]=SGN(y);
    for (c=0;c<LOOPMAX;c++) {
    //while (1) {
        newdist=dist[k]+y*y/D[k]; /* f(z) from last item to k item ref.[2] (18) */
        if (newdist<maxdist) {
            /* generate a set of candidate ambiguities */
            if (k!=0) {
                dist[--k]=newdist;
                /* ref.[2] (17) */
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                /* next candidate ambiguity */
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
            }
            /* after generate a set of candidate ambiguities */
            else {
                /* set m sets of candidate ambituities to zn and f(z) to s 
                *set first solution of zn as z directly, 
                *and set second solution of zn as z with z[0] as the next integer ...
                *record max f(z)
                 */
                if (nn<m) {
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                /* after set m sets amb to zn */
                else {
                    /* if f(z)<s[imax], replace imax solution with current solution */
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        /* find the max f(z) */
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    /* set max f(z) as maxdist, continue to find z[0] with f(z)<s[imax] */
                    maxdist=s[imax];
                }
                /* set z[0] to the next candidate integer, 
                *z[0] changed, y changed, newdist(f(z)) changed */
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        /* no z[0] could make f(z)<s[imax] */
        else {
            if (k==n-1) break;
            /* change z[1] to make f(z)<s[imax], ... , z[n-1] to make f(z)<s[imax] */
            else {
                k++;
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s, small to large */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        //fprintf(stderr,"%s:search loop count overflow\n",__FILE__);
        trace(2,"search loop count overflow\n");
        return -1;
    }
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args  :int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return:status (0:ok,other:error)
* notes :matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s)
{
    int info;
    double *L,*D,*Z,*z,*E;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);
    
    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {
        
        /* lambda reduction */
        reduction(n,L,D,Z);
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */
        
        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {
            
            info=solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args  :int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
*          double *L     O  L of Qz (n x n)
*          double *D     O  D of Qz (n x 1)
* return:status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambdaReduction(int n, const double *Q, double *Z, double *L, 
    double *D)
{
    int i,j,info;

    if (n<=0) return -1;

    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        Z[i+j*n]=i==j?1.0:0.0;
        L[i+j*n]=0.0;
    }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) return info;

    /* lambda reduction */
    reduction(n,L,D,Z);

    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args  :int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1) before Z transformation
*          double *L     I  (n x n) L of Qz
*          double *D     I  (n x 1) D of Qz
*          double *Z     I  (n x n) Z transformation matrix
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return:status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambdaSearch(int n, int m, const double *a, const double *L, 
    const double *D, const double *Z, double *F, double *s)
{
    double *z,*E;
    int info;
    
    if (n<=0||m<=0) return -1;

    z=mat(n,1); E=mat(n,m);
    matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */

    /* mlambda search */
    if (!(info=search(n,m,L,D,z,E,s))) {
        info=solve("T",Z,E,n,m,F); /* F=Z'\E */
    }

    free(z); free(E);
    return info;
}
static double normalCDF(double x)
{
    return erfc(-x/sqrt(2.0))/2.0;
}
/* success rate --------------------------------------------------------------*/
static double successRate(const double *D, int n)
{
    double ps=1.0;
    int i;

    for (i=0;i<n;i++) ps*=(2*normalCDF(0.5/sqrt(D[i]))-1.0);
    return ps;
}
/* ar by lambda --------------------------------------------------------------*/
extern int lambdaAR(int n, int m, const double *a, const double *Qa, double *F,
    double *test)
{
    double *Z,*L,*D,s[2],rate;
    int stat=1;

    Z=mat(n,n);L=mat(n,n);D=mat(n,1);
    if (lambdaReduction(n,Qa,Z,L,D)) {
        stat=0;
        test[0]=test[1]=0.0;
    }
    rate=successRate(D,n);
    if (rate<test[1]) {
        test[0]=0.0;
        test[1]=rate;
    }
    else {
        if (lambdaSearch(n,m,a,L,D,Z,F,s)) {
            stat=0;
            test[0]=test[1]=0.0;
        }
        else {
            /* ratio */
            if (s[0]<=0.0) test[0]=0.0;
            else test[0]=s[1]/s[0];
            test[1]=rate;
        }
    }

    free(Z); free(L); free(D);
    return stat;
}
#if 0
/* ar by lambda --------------------------------------------------------------*/
extern int lambdaAR(int n, int m, const double *a, const double *Qa, double *F,
    double *test)
{
    double *Z,*L,*D,s[2];
    int stat=1;

    Z=mat(n,n); L=mat(n,n); D=mat(n,1);
    if (lambdaReduction(n,Qa,Z,L,D)||lambdaSearch(n,m,a,L,D,Z,F,s)) {
        stat=0;
        test[0]=test[1]=0.0;
    }
    else {
        /* ratio */
        if (s[0]<=0.0) test[0]=0.0;
        else test[0]=s[1]/s[0];
        /* success rate */
        test[1]=successRate(D,n);
    }

    free(Z); free(L); free(D);
    return stat;
}
#endif
void main_lam()
{
    double Q[9]={6.29,5.97,1.54,5.97,7.29,4.34,1.54,4.34,9.28};
    double Z[9],L[9],D[3];
    //int i;

    lambdaReduction(3,Q,Z,L,D);
    printf("Q=\n");
    matprint(Q,3,3,6,2);
    printf("Z=\n");
    matprint(Z,3,3,6,2);
    //printf("L=\n");
    //matprint(L,3,3,6,2);
    //printf("D=\n");
    //matprint(D,1,3,6,2);

    double zq[9],Qz[9];
    matmul("TN",3,3,3,1.0,Z,Q,0.0,zq);
    matmul("NN",3,3,3,1.0,zq,Z,0.0,Qz);
    printf("Qz=\n");
    matprint(Qz,3,3,6,2);

    double N[3]={1,1,1},a[3];
    matmul("TN",3,1,3,1.0,Z,N,0.0,a);
    printf("a=\n");
    matprint(a,1,3,6,2);

    //double QD[9]={0},ld[9],LDL[9];
    //for (i=0;i<3;i++) QD[i+i*3]=D[i];
    //printf("QD=\n");
    //matprint(QD,3,3,6,2);
    //matmul("TN",3,3,3,1.0,L,QD,0.0,ld);
    //matmul("NN",3,3,3,1.0,ld,L,0.0,LDL);
    //printf("LDL=\n");
    //matprint(LDL,3,3,6,2);

    //QD[0]=0.0;
    //for (i=0;i<3;i++) L[i]=L[i*3]=0.0;
    //printf("L=\n");
    //matprint(L,3,3,6,2);
    //printf("QD=\n");
    //matprint(QD,3,3,6,2);
    //matmul("TN",3,3,3,1.0,L,QD,0.0,ld);
    //matmul("NN",3,3,3,1.0,ld,L,0.0,LDL);
    //printf("LDL=\n");
    //matprint(LDL,3,3,6,2);

    //for (i=0;i<3;i++) Qz[i]=Qz[i*3]=0.0;
    Qz[0]=0.0;
    printf("Qz=\n");
    matprint(Qz,3,3,6,2);
    if (matinv(Z,3)) return;
    matmul("TN",3,3,3,1.0,Z,Qz,0.0,zq);
    matmul("NN",3,3,3,1.0,zq,Z,0.0,Qz);
    printf("Qz=\n");
    matprint(Qz,3,3,6,2);

    a[0]=0;
    matmul("TN",3,1,3,1.0,Z,a,0.0,N);
    printf("N=\n");
    matprint(N,1,3,6,2);


    system("pause");

}