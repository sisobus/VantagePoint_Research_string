
#define UNIXMACH

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#ifdef UNIXMACH
#include <strings.h>
#include <sys/time.h>
#else
#include <ctype.h>
#endif

#define true      1
#define false     0
#define LineSize  80
#define PerCent   0.01
#define MicroSec  0.001 

typedef enum {
    MaxCut=-1,
    EquCut= 0,
    UquCut= 2,
    DvdCut=-2,
    BoxCut=-3
  } prbset;

typedef enum {
    Gtype=0,
    Mtype=1
  } datset;
  
typedef enum {
    START=0,
    DATAIN,
    OPTIM,
    INTEG,
    ELAPS
  } tim_set;

typedef enum {
    CfcOk=0,
    CfcSpace,   /* fail to allocate required space */
    CfcIndef    /* indefinity is detected          */
  } cfc_sta;

typedef enum {
    OptFound=0,
    IterLimit,
    
    SysError=100,
    OutOfSpc,WriteFail,ReadFail,CholErr,
    MtxSzeErr,InvalCut,DataErr
  } xcode;
  
typedef struct {
    int    maxiter;
    int    prtlev;
    int    dynrho;
    int    varstep;
    double tolgap;
    double alph;
    double tolpiv;
    double bmfac;
    int    findvec;
    double rhofac;
    int    cachesize;
    int    cacheunit;
  } params;
  
typedef struct {
    int    nn0;
    int    *ja;
    double *an;
  } array;

typedef struct {
    int    maxnrow;
    int    nrow;
    int    non0;
    array  *rows;
  } smatx;

typedef struct {
    int    nrow;
    int    nnzo;
    double *diag;
    array  *roff;
 } symat;
 
typedef struct { 
    int    nrow; 
    int    nnzo; 
    array  *roff; 
 } syoff; 
  
typedef struct {
    int    mrow;     /* number of rows allocated                    */
    int    nrow;     /* number of rows used                         */
    
    int    snnz;     /* number of indices for nonzeros in S         */
    int    *shead;   /* position of first nonzero in row i of S     */
    int    *ssize;   /* number of non-zeros in row i of S below     */
                     /* the diagonal                                */
    int    *ssub;    /* column index buffer for non-zeros in S      */
    double *diag;    /* diagonal matrix D in the factorization      */
    
    int    unnz;     /* number of nonzeros in the upper factor      */
    int    ujnz;     /* number of column indices in the compressed  */
                     /* indices buffer ujsub                        */ 
    int    *ujbeg;   /* beginning position of indices in row i of U */
    int    *uhead;   /* position of first nonzero in row i of U     */
    int    *ujsze;   /* number of indices in row i of U             */ 
    int    *usub;    /* compressed column index buffer of U         */
    double *uval;    /* nonzero values in factor U                  */
    
    int    *perm;    /* permutation order                           */
    int    *invp;    /* inverse order of perm                       */
    
    int    nsnds;    /* number of supernodes                        */
    int    *subg;    /* index of the first column in supernode i    */
    int    ndens;    /* numer of dense rows                         */
    int    nsndn;    /* number supernodes in dense rows             */
    int    *dhead;   /* pointer first column in each dense row      */
    int    *dsub;    /* indices in dense rows                       */
    int    *dbeg;    /* beginning of column index                   */ 
    int    sdens;    /* separate dense row                          */
    int    upst;     /* specified index in uval for st-cut problem  */
    int    NegDiag;  /* Index of negative diagonal element. -1 if no*/
  } chfac;

typedef struct model {
    char   ss[60];  /* problem name                                     */
    char   spc[20]; /* specific file name                               */
    
    prbset ptyp;    /* problem type                                     */
    datset dtyp;    /* input data type                                  */
    
    int    iter;    /* iteration count                                  */
    int    nrow;    /* number of constraints                            */
    int    ncol;    /* number of variables                              */
    int    nnzo;    /* number of nonzeros in S                          */
    
    int    nblk;    /* number of diagonal blocks in S                   */
    int    maxblk;  /* maximum order of diagonal block in S             */
    int    minblk;  /* minimum order of diagonal block in S             */
       
    int    nbiv;    /* number of diagonal blocks in Sinv
    int    maxbiv;  /* maximum order of diagonal blocks in Sinv         */
    int    minbiv;  /* minimum order of diagonal blocks in Sinv         */
    
    int    dnfc;    /* dense decomposition count for linesearch         */
    int    dnck;    /* dense decomposition count for upper bound update */
    
    int    is;      /* specified index for st-cut problems              */
    int    it;      /* specified index for st-cut problems              */
    
    double rho;     /* penalty parameter in dual potential function     */
    double step;    /* the current step size                            */
    double pval;    /* the current P-value                              */
    
    double rgap;    /* the current primal dual gap                      */    
    double pobj;    /* the current primal objective value               */
    double dobj;    /* the current dual objective value                 */
    
    double bigM;    /* bigM for equal cut problems                      */
    double x0;      /* artificial variable for equal cut problems       */ 
    double s0;      /* dual slack for equal or unequal cut problems     */
    double lamda;   /* dual variable corresponding to x0                */
    double dl;      /* search direction component concerning lamda      */
    double kap;     /* cut value for unequal cut problems               */

    double *y;      /* dual variables                                   */
    double *dy;     /* search direction vector concerning y             */
    double *s;      /* dual slacks                                      */
    double *sinv;   /* inverse of S                                     */
    
    smatx *c;       /* symmetric matrix in full format                  */
    symat *cy;      /* symmetric matrix in triangular format            */
    syoff *st;      /* symmetric matrix in off-diagonal triangular form */
    
    chfac *sf;      /* cholesky factor for S                            */
    chfac *mf;      /* cholesky factor for M                            */
    
    double **u;     /* matrix for primal solution generation            */
    double **v;     /* factor of primal solution                        */
    
    int    *iw;     /* integer working buffer                           */
    double *rw;     /* double float working buffer                      */
    
    params par;     /* parameters                                       */
  } sdpdat;
  
extern sdpdat  *sdat;
extern FILE    *fout;
extern clock_t *stim;
extern clock_t lt,vt,qt,gt;

#include "SDPfun.h"
