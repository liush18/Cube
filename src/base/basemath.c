
#include "base.h"

static fatalfunc_t *fatalfunc=NULL;/* fatal callback function */


#ifdef MKL
#define LAPACK
#define dgemm_      dgemm
#define dgetrf_     dgetrf
#define dgetri_     dgetri
#define dgetrs_     dgetrs
#define mkl_free_buffers_    mkl_free_buffers
#endif
#ifdef LAPACK
extern void dgemm_(char *, char *, int *, int *, int *, double *, double *,
    int *, double *, int *, double *, double *, int *);
extern void dgetrf_(int *, int *, double *, int *, int *, int *);
extern void dgetri_(int *, double *, int *, int *, double *, int *, int *);
extern void dgetrs_(char *, int *, int *, double *, int *, int *, double *,
    int *, int *);
extern void mkl_free_buffers_();
#endif

/* fatal error ---------------------------------------------------------------*/
static void fatalerr(const char* format,...)
{
    char msg[1024];
    va_list ap;
    va_start(ap,format);vsprintf(msg,format,ap);va_end(ap);
    if (fatalfunc) fatalfunc(msg);
    else fprintf(stderr,"%s",msg);
    exit(-9);
}
/* add fatal callback function -------------------------------------------------
* add fatal callback function for mat(),zeros(),imat()
* args  :fatalfunc_t *func I  callback function
* return:none
* notes :if malloc() failed in return:none
*-----------------------------------------------------------------------------*/
extern void add_fatal(fatalfunc_t *func)
{
    fatalfunc=func;
}
/* new matrix ------------------------------------------------------------------
* allocate memory of matrix
* args  :int    n,m       I   number of rows and columns of matrix
* return:matrix pointer (if n<=0 or m<=0,return NULL)
*-----------------------------------------------------------------------------*/
extern double *mat(int n,int m)
{
    double* p;

    if (n<= 0||m<= 0) return NULL;
    if (!(p=(double*)malloc(sizeof(double)*n*m))) {
        fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* new integer matrix ----------------------------------------------------------
* allocate memory of integer matrix
* args  :int    n,m       I   number of rows and columns of matrix
* return:matrix pointer (if n<=0 or m<=0,return NULL)
*-----------------------------------------------------------------------------*/
extern int *imat(int n, int m)
{
    int *p;

    if (n<=0||m<=0) return NULL;
    if (!(p=(int*)malloc(sizeof(int)*n*m))) {
        fatalerr("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* zero matrix -----------------------------------------------------------------
* generate new zero matrix
* args  :int    n,m       I   number of rows and columns of matrix
* return:matrix pointer (if n<=0 or m<=0,return NULL)
*-----------------------------------------------------------------------------*/
extern double *zeros(int n, int m)
{
    double *p;

    if (n<=0||m<=0) return NULL;
    if (!(p=(double*)calloc(sizeof(double),n*m))) {
        fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* zero integer matrix ---------------------------------------------------------*/
extern int *izeros(int n,int m)
{
    int *p;

    if (n<= 0|| m<= 0) return NULL;
    if (!(p=(int*)calloc(sizeof(int),n*m))) {
        fatalerr("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
    }

    return p;
}
/* identity matrix -------------------------------------------------------------
* generate new identity matrix
* args  :int    n         I   number of rows and columns of matrix
* return:matrix pointer (if n<=0,return NULL)
*-----------------------------------------------------------------------------*/
extern double *eye(int n)
{
    double *p;
    int i;

    if((p=zeros(n,n))) for(i=0;i<n;i++) p[i+i*n]=1.0;
    return p;
}
/* inner product ---------------------------------------------------------------
* inner product of vectors
* args  :double *a,*b     I   vector a,b (n xp 1)
*          int    n         I   size of vector a,b
* return:a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a,const double *b,int n)
{
    double c=0.0;

    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args  :double *a        I   vector a (n xp 1)
*          int    n         I   size of vector a
* return :||a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double *a,int n)
{
    return sqrt(dot(a,a,n));
}
/* outer product of 3d vectors -------------------------------------------------
* outer product of 3d vectors
* args  :double *a,*b     I   vector a,b (3 xp 1)
*          double *c        O   outer product (a xp b) (3 xp 1)
* return:none
*-----------------------------------------------------------------------------*/
extern void cross3(const double *a,const double *b,double *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}
/* normalize 3d vector ---------------------------------------------------------
* normalize 3d vector
* args  :double *a        I   vector a (3 xp 1)
*          double *b        O   normlized vector (3 xp 1)||b ||=1
* return:status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int normv3(const double *a,double *b)
{
    double r;
    if ((r=norm(a,3))<=0.0) return 0;
    b[0]=a[0]/r;
    b[1]=a[1]/r;
    b[2]=a[2]/r;
    return 1;
}
/* copy matrix -----------------------------------------------------------------
* copy matrix
* args  :double *A        O   destination matrix A (n xp m)
*          double *B        I   source matrix B (n xp m)
*          int    n,m       I   number of rows and columns of matrix
* return:none
*-----------------------------------------------------------------------------*/
extern void matcpy(double *A,const double *B,int n,int m)
{
    memcpy(A,B,sizeof(double)*n*m);
}
extern void imatcpy(int *A,const int *B,int n,int m)
{
    memcpy(A,B,sizeof(int)*n*m);
}

/* matrix routines -----------------------------------------------------------*/
#ifdef LAPACK /* with LAPACK/BLAS or MKL */

/* multiply matrix (wrapper of blas dgemm) -------------------------------------
* multiply matrix by matrix (C=alpha*A*B+beta*C)
* args  :char   *tr       I  transpose flags ("N":normal,"T":transpose)
*          int    n,k,m     I  size of (transposed) matrix A,B
*          double alpha     I  alpha
*          double *A,*B     I  (transposed) matrix A (n x m), B (m x k)
*          double beta      I  beta
*          double *C        IO matrix C (n x k)
* return:none
*-----------------------------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
    const double *A, const double *B, double beta, double *C)
{
    int lda=tr[0]=='T'?m:n,ldb=tr[1]=='T'?k:m;

    dgemm_((char *)tr,(char *)tr+1,&n,&k,&m,&alpha,(double *)A,&lda,(double *)B,
        &ldb,&beta,C,&n);
    mkl_free_buffers_();
}
/* inverse of matrix -----------------------------------------------------------
* inverse of matrix (A=A^-1)
* args  :double *A        IO  matrix (n x n)
*          int    n         I   size of matrix A
* return:status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double *work;
    int info,lwork=n*16,*ipiv=imat(n,1);

    work=mat(lwork,1);
    dgetrf_(&n,&n,A,&n,ipiv,&info);
    if (!info) dgetri_(&n,A,&n,ipiv,work,&lwork,&info);
    mkl_free_buffers_();
    free(ipiv); free(work);
    return info;
}
/* solve linear equation -------------------------------------------------------
* solve linear equation (X=A\Y or X=A'\Y)
* args  :char   *tr       I   transpose flag ("N":normal,"T":transpose)
*          double *A        I   input matrix A (n x n)
*          double *Y        I   input matrix Y (n x m)
*          int    n,m       I   size of matrix A,Y
*          double *X        O   X=A\Y or X=A'\Y (n x m)
* return:status (0:ok,0>:error)
* notes :matirix stored by column-major order (fortran convention)
*          X can be same as Y
*-----------------------------------------------------------------------------*/
extern int solve(const char *tr, const double *A, const double *Y, int n,
    int m, double *X)
{
    double *B=mat(n,n);
    int info,*ipiv=imat(n,1);

    matcpy(B,A,n,n);
    matcpy(X,Y,n,m);
    dgetrf_(&n,&n,B,&n,ipiv,&info);
    if (!info) dgetrs_((char *)tr,&n,&m,B,&n,ipiv,X,&n,&info);
    mkl_free_buffers_();
    free(ipiv); free(B); 
    return info;
}

#else /* without LAPACK/BLAS or MKL */
/* multiply matrix -----------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
    const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

    for (i=0;i<n;i++) for (j=0;j<k;j++) {
        d=0.0;
        switch (f) {
        case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m];break;
        case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k];break;
        case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m];break;
        case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k];break;
        }
        if (beta==0.0) C[i+j*n]=alpha*d;else C[i+j*n]=alpha*d+beta*C[i+j*n];
    }
}
/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double *A, int n, int *indx, double *d)
{
    double big,s,tmp,*vv=mat(n,1);
    int i,imax=0,j,k;

    *d=1.0;
    for (i=0;i<n;i++) {
        big=0.0;for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
        if (big>0.0) vv[i]=1.0/big;else { free(vv);return -1;}
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            s=A[i+j*n];for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n];A[i+j*n]=s;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            s=A[i+j*n];for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n];A[i+j*n]=s;
            if ((tmp=vv[i]*fabs(s))>=big) { big=tmp;imax=i;}
        }
        if (j!=imax) {
            for (k=0;k<n;k++) {
                tmp=A[imax+k*n];A[imax+k*n]=A[j+k*n];A[j+k*n]=tmp;
            }
            *d=-(*d);vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (A[j+j*n]==0.0) { free(vv);return -1;}
        if (j!=n-1) {
            tmp=1.0/A[j+j*n];for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
        }
    }
    free(vv);
    return 0;
}
/* LU back-substitution ------------------------------------------------------*/
static void lubksb(const double *A, int n, const int *indx, double *b)
{
    double s;
    int i,ii=-1,ip,j;

    for (i=0;i<n;i++) {
        ip=indx[i];s=b[ip];b[ip]=b[i];
        if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j];else if (s) ii=i;
        b[i]=s;
    }
    for (i=n-1;i>=0;i--) {
        s=b[i];for(j=i+1;j<n;j++) s-=A[i+j*n]*b[j];b[i]=s/A[i+i*n];
    }
}
/* inverse of matrix ---------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double d,*B;
    int i,j,*indx;

    indx=imat(n,1);B=mat(n,n);matcpy(B,A,n,n);
    if (ludcmp(B,n,indx,&d)) { free(indx);free(B);return -1;}
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) A[i+j*n]=0.0;
        A[j+j*n]=1.0;
        lubksb(B,n,indx,A+j*n);
    }
    free(indx);free(B);
    return 0;
}
/* solve linear equation -----------------------------------------------------*/
extern int solve(const char *tr, const double *A, const double *Y, int n,
    int m, double *X)
{
    double *B=mat(n,n);
    int info;

    matcpy(B,A,n,n);
    if (!(info=matinv(B,n))) matmul(tr[0]=='N'?"NN":"TN",n,m,n,1.0,B,Y,0.0,X);
    free(B);
    return info;
}
#endif
/* end of matrix routines ----------------------------------------------------*/

/* least square estimation -----------------------------------------------------
* least square estimation by solving normal equation (xp=(A*A')^-1*A*v)
* args  :double *A        I   transpose of (weighted) design matrix (n xp m)
*          double *v        I   (weighted) measurements (m xp 1)
*          int    n,m       I   number of parameters and measurements (n<=m)
*          double *xp        O   estmated parameters (n xp 1)
*          double *Q_        O   esimated parameters covariance matrix (n xp n)
* return:status (0:ok,0>:error)
* notes :for weighted least square,replace A and v by A*w and w*v (w=W^(1/2))
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern int lsq(const double *A,const double *y,int n,int m,double *x,
    double *Q)
{
    double *Ay;
    int info;

    if (m<n) return -1;
    Ay=mat(n,1);
    matmul("NN",n,1,m,1.0,A,y,0.0,Ay);/* Ay=A*v */
    matmul("NT",n,n,m,1.0,A,A,0.0,Q); /* Q_=A*A' */
    if (!(info=matinv(Q,n))) matmul("NN",n,1,n,1.0,Q,Ay,0.0,x);/* xp=Q_^-1*Ay */
    free(Ay);
    return info;
}
/* lsq filter in estProduct_d -----------------------------------------------------*/
extern int lsqfilter(double *x, double *P, const double *H, const double *v,
    const double *R, const int *index, int n, int m)
{
    double *dx,*Q_,*H_,*W,*F,*T;
    int i,j,k,info,*ix;

    ix=imat(n,1);
    for(i=k=0;i<n;i++) if(index[i]&1) ix[k++]=i;
    if (m<k) return -1;

    dx=mat(k,1);Q_=mat(k,k);H_=mat(k,m);
    W=mat(m,m);F=mat(k,m);T=mat(k,1);
    for (i=0;i<k;i++) for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
#if 0
    trace(2,"H_:\n"); tracemat(2,H_,k,m,6,3);
#endif

    matcpy(W,R,m,m);
    if (info=matinv(W,m)) return info;

    matmul("NN",k,m,m,1.0,H_,W,0.0,F);      /* F=H*W */
    matmul("NT",k,k,m,1.0,F,H_,0.0,Q_);     /* Q_=F*H' */
    if (!(info=matinv(Q_,k))) {
        matmul("NN",k,1,m,1.0,F,v,0.0,T);   /* T=F*v */
        matmul("NN",k,1,k,1.0,Q_,T,0.0,dx); /* dx=Q_*T */
    }
    for (i=0;i<k;i++) {
        x[ix[i]]+=dx[i];
        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Q_[i+j*k];
    }

    free(ix);free(dx);free(Q_);free(H_);free(W);free(F);free(T);
    return info;
}
/* lsq filter by block without constraint ------------------------------------*/
extern int lsqSingleBlock(double *x, double *P, const int *ix, const double *vl,
    const double *HL, const double *DL, int k, int nx, int nv)
{
    double *_HL,*WL,*FL,*PT,*FLT,*dx;
    int i,j,info;

    /* simplified designed matrix for measurement and constraint */
    _HL=mat(k,nv);
    for (i=0;i<k;i++) {
        for (j=0;j<nv;j++) _HL[i+j*k]=HL[ix[i]+j*nx];
    }
#if 0
    tracemat(2,_HL,k,nv,6,2);
#endif
    /* convert co-variance to weight */
    WL=mat(nv,nv);
    matcpy(WL,DL,nv,nv);
    if ((info=matinv(WL,nv))) {
        free(_HL);
        return info;
    }

    FL=mat(k,nv);PT=mat(k,k);
    matmul("NN",k,nv,nv,1.0,_HL, WL,0.0,FL);     /* FL=HL'*DL */
    matmul("NT",k, k,nv,1.0, FL,_HL,0.0,PT);     /* PT=HL'*DL*HL */

    FLT=mat(k,1);dx=mat(k,1);
    if (!(info=matinv(PT,k))) {
        matmul("NN",k,1,nv,1.0,FL, vl,0.0,FLT);   /* FLT=FL*v */
        matmul("NN",k,1, k,1.0,PT,FLT,0.0, dx);   /* dx=Q_*T */

        for (i=0;i<k;i++) {
            x[ix[i]]=dx[i];
            for (j=0;j<k;j++) P[ix[i]+ix[j]*nx]=PT[i+j*k];
        }
    }

    free(_HL);free(WL);
    free(FL);free(PT);free(FLT);free(dx);
    return info;
}
/* lsq filter by block -------------------------------------------------------*/
extern int lsqBlock(double *x, double *P, const int *index, const double *vl,
    const double *HL, const double *DL, const double *vp, const double *HP, 
    const double *DP, int nx, int nv, int nc)
{
    double *_HL,*_HP,*WL,*WP,*FL,*FP,*PT,*FLT,*dx;
    int i,j,k,info,*ix;

    ix=imat(nx,1);
    for(i=k=0;i<nx;i++) if(index[i]) ix[k++]=i;
    if ((nv+nc)<k) {
        free(ix);
        return -1;
    }

    if (nc==0) {
        info=lsqSingleBlock(x,P,ix,vl,HL,DL,k,nx,nv);
        free(ix);
        return info;
    }

    /* simplified designed matrix for measurement and constraint */
    _HL=mat(k,nv);_HP=mat(k,nc);
    for (i=0;i<k;i++) {
        for (j=0;j<nv;j++) _HL[i+j*k]=HL[ix[i]+j*nx];
        for (j=0;j<nc;j++) _HP[i+j*k]=HP[ix[i]+j*nx];
    }
#if 0
    trace(2,"HL:\n");
    tracemat(2,_HL,k,nv,6,2);
    trace(2,"HP:\n");
    tracemat(2,_HP,k,nc,6,2);
#endif
    /* convert co-variance to weight */
    WL=mat(nv,nv);WP=mat(nc,nc);
    matcpy(WL,DL,nv,nv);
    matcpy(WP,DP,nc,nc);

    if ((info=matinv(WL,nv))||(info=matinv(WP,nc))) {
        free(ix);free(_HL);free(_HP);free(WL);free(WP);
        return info;
    }

    FL=mat(k,nv);FP=mat(k,nc);PT=mat(k,k);
    matmul("NN",k,nv,nv,1.0,_HL, WL,0.0,FL);     /* FL=HL'*DL */
    matmul("NT",k, k,nv,1.0, FL,_HL,0.0,PT);     /* PT=HL'*DL*HL */
    matmul("NN",k,nc,nc,1.0,_HP, WP,0.0,FP);     /* FP=HP'*DP */
    matmul("NT",k, k,nc,1.0, FP,_HP,1.0,PT);     /* PT+=HL'*DL*HL */

    FLT=mat(k,1);dx=mat(k,1);
    if (!(info=matinv(PT,k))) {
        matmul("NN",k,1,nv,1.0,FL, vl,0.0,FLT);   /* FLT=FL*v */
        matmul("NN",k,1,nc,1.0,FP, vp,1.0,FLT);   /* FLT+=FP*vp */
        matmul("NN",k,1, k,1.0,PT,FLT,0.0, dx);   /* dx=Q_*T */

        for (i=0;i<k;i++) {
            x[ix[i]]=dx[i];
            for (j=0;j<k;j++) P[ix[i]+ix[j]*nx]=PT[i+j*k];
        }
    }

    free(ix);free(_HL);free(_HP);free(WL);free(WP);
    free(FL);free(FP);free(PT);free(FLT);free(dx);
    return info;
}
/* lsq filter in ppp ---------------------------------------------------------*/
extern int filterLSQ(double *x, double *P, const double *H, const double *v,
    const double *R, int n, int m)
{
    double *dx,*Q_,*H_,*W,*F,*T;
    int i,j,k,info,*ix;

    ix=imat(n,1);
    for (i=k=0;i<n;i++) if (x[i]!=0.0&&P[i+i*n]>0.0) ix[k++]=i;
    if (m<k) return -1;

    dx=mat(k,1);Q_=mat(k,k);H_=mat(k,m);
    W=mat(m,m);F=mat(k,m);T=mat(k,1);
    for (i=0;i<k;i++) for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
#if 0
    trace(2,"H_:\n"); tracemat(2,H_,k,m,6,3);
#endif

    matcpy(W,R,m,m);
    if (info=matinv(W,m)) return info;

    matmul("NN",k,m,m,1.0,H_,W,0.0,F);      /* F=H*W */
    matmul("NT",k,k,m,1.0,F,H_,0.0,Q_);     /* Q_=F*H' */
    if (!(info=matinv(Q_,k))) {
        matmul("NN",k,1,m,1.0,F,v,0.0,T);   /* T=F*v */
        matmul("NN",k,1,k,1.0,Q_,T,0.0,dx); /* dx=Q_*T */
    }
    for (i=0;i<k;i++) {
        x[ix[i]]+=dx[i];
        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Q_[i+j*k];
    }

    free(ix);free(dx);free(Q_);free(H_);free(W);free(F);free(T);
    return info;
}
/* kalman filter ---------------------------------------------------------------
* kalman filter state update as follows:
*
*   K=P*H*(H'*P*H+R)^-1,xp=xp+K*v,Q_=(I-K*H')*P
*
* args  :double *xp        I   states vector (n xp 1)
*          double *P        I   covariance matrix of states (n xp n)
*          double *H        I   transpose of design matrix (n xp m)
*          double *v        I   innovation (measurement-model) (m xp 1)
*          double *R        I   covariance matrix of measurement error (m xp m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n xp 1)
*          double *Q_       O   covariance matrix of states after update (n xp n)
* return:status (0:ok,<0:error)
* notes :matirix stored by column-major order (fortran convention)
*          if state xp[i]==0.0,not updates state xp[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
static int filter_(const double *x,const double *P,const double *H,
    const double *v,const double *R,int n,int m,
    double *xp,double *Pp)
{
    double *F=mat(n,m),*Q=mat(m,m),*K=mat(n,m),*I=eye(n);
    double *P1=mat(n,n),*T=mat(n,m);
    int info;

    matcpy(Q,R,m,m);
    matcpy(xp,x,n,1);
    matmul("NN",n,m,n,1.0,P,H,0.0,F);      /* F=P*H */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);      /* Q_ =H'*F+Q_ =H'*P*H+R */
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);  /* K=P*H*Q_^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp); /* xp=xp+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I); /* I=I-K*H' */
        matmul("NN",n,n,n,1.0,I,P,0.0,P1); /* P1=IP =(I-K*H')*P */
        matmul("NT",n,n,n,1.0,P1,I,0.0,Pp);/* Q_=P1*I=(I-KH')P(I-KH')' */
        matmul("NN",n,m,m,1.0,K,R,0.0,T);  /* T=K*R */
        matmul("NT",n,n,m,1.0,T,K,1.0,Pp); /* Q_=Q_+TK' =(I-KH')P(I-KH')'+KRK'*/
    }
    free(F);free(Q);free(K);free(I);
    free(P1);free(T);
    return info;
}
extern int filter(double *x, double *P, const double *H, const double *v,
    const double *R, int n, int m)
{
    double *x_,*xp_,*P_,*Pp_,*H_;
    int i,j,k,info,*ix;

    ix=imat(n,1); for (i=k=0;i<n;i++) if (x[i]!=0.0&&P[i+i*n]>0.0) ix[k++]=i;
    x_=mat(k,1);xp_=mat(k,1);P_=mat(k,k);Pp_=mat(k,k);H_=mat(k,m);
    for (i=0;i<k;i++) {
        x_[i]=x[ix[i]];
        for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
        for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
    }
    //tracemat(2,H_,k,m,6,3);
    info=filter_(x_,P_,H_,v,R,k,m,xp_,Pp_);
    for (i=0;i<k;i++) {
        x[ix[i]]=xp_[i];
        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp_[i+j*k];
    }
    free(ix);free(x_);free(xp_);free(P_);free(Pp_);free(H_);
    return info;
}
extern int filterKalman(double *x, double *P, const double *H, const double *v,
    const double *R, const int *index, int n, int m)
{
    double *x_,*xp_,*P_,*Pp_,*H_;
    int i,j,k,info,*ix;

    ix=imat(n,1); for(i=k=0;i<n;i++) if(index[i]==1) ix[k++]=i;
    //if (m<k) return -1;

    x_=mat(k,1);xp_=mat(k,1);P_=mat(k,k);Pp_=mat(k,k);H_=mat(k,m);
    for (i=0;i<k;i++) {
        x_[i]=x[ix[i]];
        for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
        for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
    }
    tracemat(2,H_,k,m,6,3);
    info=filter_(x_,P_,H_,v,R,k,m,xp_,Pp_);
    for (i=0;i<k;i++) {
        x[ix[i]]=xp_[i];
        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp_[i+j*k];
    }
    free(ix);free(x_);free(xp_);free(P_);free(Pp_);free(H_);
    return info;
}
/* smoother --------------------------------------------------------------------
* combine forward and backward filters by fixed-interval smoother as follows:
*
*   xs=Qs*(Qf^-1*xf+Qb^-1*xb),Qs=(Qf^-1+Qb^-1)^-1)
*
* args  :double *xf       I   forward solutions (n xp 1)
* args  :double *Qf       I   forward solutions covariance matrix (n xp n)
*          double *xb       I   backward solutions (n xp 1)
*          double *Qb       I   backward solutions covariance matrix (n xp n)
*          int    n         I   number of solutions
*          double *xs       O   smoothed solutions (n xp 1)
*          double *Qs       O   smoothed solutions covariance matrix (n xp n)
* return:status (0:ok,0>:error)
* notes :see reference [4] 5.2
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern int smoother(const double *xf,const double *Qf,const double *xb,
    const double *Qb,int n,double *xs,double *Qs)
{
    double *invQf=mat(n,n),*invQb=mat(n,n),*xx=mat(n,1);
    int i,info=-1;

    matcpy(invQf,Qf,n,n);
    matcpy(invQb,Qb,n,n);
    if (!matinv(invQf,n)&&!matinv(invQb,n)) {
        for (i=0;i<n*n;i++) Qs[i]=invQf[i]+invQb[i];
        if (!(info=matinv(Qs,n))) {
            matmul("NN",n,1,n,1.0,invQf,xf,0.0,xx);
            matmul("NN",n,1,n,1.0,invQb,xb,1.0,xx);
            matmul("NN",n,1,n,1.0,Qs,xx,0.0,xs);
        }
    }
    free(invQf);free(invQb);free(xx);
    return info;
}
/* print matrix ----------------------------------------------------------------
* print matrix to stdout
* args  :double *A        I   matrix A (n xp m)
*          int    n,m       I   number of rows and columns of A
*          int    p,q       I   total columns,columns under decimal point
*         (FILE  *fp        I   output file pointer)
* return:none
* notes :matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern void matfprint(const double A[],int n,int m,int p,int q,FILE *fp)
{
    int i,j;

    for (i=0;i<n;i++) {
        for (j=0;j<m;j++) fprintf(fp,"%*.*f",p,q,A[i+j*n]);
        fprintf(fp,"\n");
    }
}
extern void matprint(const double A[],int n,int m,int p,int q)
{
    matfprint(A,n,m,p,q,stdout);
}
/* print integer matrix ------------------------------------------------------*/
extern void imatfprint(const int A[],int n,int m,int p,FILE *fp)
{
    int i,j;

    for (i=0;i<n;i++) {
        for (j=0;j<m;j++) fprintf(fp," %*d",p,A[i+j*n]);
        fprintf(fp,"\n");
    }
}
extern void imatprint(const int A[],int n,int m,int p)
{
    imatfprint(A,n,m,p,stdout);
}

/* convert degree to deg-min-sec -----------------------------------------------
* convert degree to degree-minute-second
* args  :double deg       I   degree
*          double *dms      O   degree-minute-second {deg,min,sec}
*          int    ndec      I   number of decimals of second
* return:none
*-----------------------------------------------------------------------------*/
extern void deg2dms(double deg,double *dms,int ndec)
{
    double sign=deg<0.0?-1.0:1.0,a=fabs(deg);
    double unit=pow(0.1,ndec);
    dms[0]=floor(a);a=(a-dms[0])*60.0;
    dms[1]=floor(a);a=(a-dms[1])*60.0;
    dms[2]=floor(a/unit+0.5)*unit;
    if (dms[2]>=60.0) {
        dms[2]=0.0;
        dms[1]+=1.0;
        if (dms[1]>=60.0) {
            dms[1]=0.0;
            dms[0]+=1.0;
        }
    }
    dms[0]*=sign;
}
/* convert deg-min-sec to degree -----------------------------------------------
* convert degree-minute-second to degree
* args  :double *dms      I   degree-minute-second {deg,min,sec}
* return:degree
*-----------------------------------------------------------------------------*/
extern double dms2deg(const double *dms)
{
    double sign=dms[0]<0.0?-1.0:1.0;
    return sign*(fabs(dms[0])+dms[1]/60.0+dms[2]/3600.0);
}
/* get index of max double value in 1d array ---------------------------------*/
extern int d_max(const double *a, const int n)
{
    int i,j=-1;

    if (norm(a,n)==0.0) return -1;
    for (i=0;i<n;i++) {
        if (j<0||a[i]>a[j]) j=i;
    }
    return j;
}
/* get index of min double value in 1d array ---------------------------------*/
extern int d_min(const double *a, const int n)
{
    int i,j=-1;

    if (norm(a,n)==0.0) return -1;
    for (i=0;i<n;i++) {
        if (j<0||a[i]<a[j]) j=i;
    }
    return j;
}
/* int in int array or not --------------------------------------------------------*/
extern int intinta(const int a, const int *A, const int n)
{
    int i;
    for (i=0;i<n;i++) if(a==A[i]) return 1;
    return 0;
}
extern double normalCDF(double x)
{
    return erfc(-x/sqrt(2.0))/2.0;
}

#define LOG_PI          1.14472988584940017 /* log(pi) */
#define SQRT2           1.41421356237309510 /* sqrt(2) */
/* complementaty error function (ref [1] p.227-229) --------------------------*/
static double q_gamma(double a, double x, double log_gamma_a);
static double p_gamma(double a, double x, double log_gamma_a)
{
    double y,w;
    int i;

    if (x==0.0) return 0.0;
    if (x>=a+1.0) return 1.0-q_gamma(a,x,log_gamma_a);

    y=w=exp(a*log(x)-x-log_gamma_a)/a;

    for (i=1;i<100;i++) {
        w*=x/(a+i);
        y+=w;
        if (fabs(w)<1E-15) break;
    }
    return y;
}
static double q_gamma(double a, double x, double log_gamma_a)
{
    double y,w,la=1.0,lb=x+1.0-a,lc;
    int i;

    if (x<a+1.0) return 1.0-p_gamma(a,x,log_gamma_a);
    w=exp(-x+a*log(x)-log_gamma_a);
    y=w/lb;
    for (i=2;i<100;i++) {
        lc=((i-1-a)*(lb-la)+(i+x)*lb)/i;
        la=lb; lb=lc;
        w*=(i-1-a)/i;
        y+=w/la/lb;
        if (fabs(w/la/lb)<1E-15) break;
    }
    return y;
}
static double f_erfc(double x)
{
    return x>=0.0?q_gamma(0.5,x*x,LOG_PI/2.0):1.0+p_gamma(0.5,x*x,LOG_PI/2.0);
}
/* confidence function of integer ambiguity ----------------------------------*/
extern double conffunc(int N, double B, double sig)
{
    double x,p=1.0;
    int i;

    x=fabs(B-N);
    for (i=1;i<8;i++) {
        p-=f_erfc((i-x)/(SQRT2*sig))-f_erfc((i+x)/(SQRT2*sig));
    }
    return p;
}