#include "SDPdef.h"

static void McutInitVal(sdpdat *sdt)
{
  int i,n=sdt->ncol;
  
  sdt->rgap=0.0;
  sdt->dobj=0.0;
  
  for (i=0; i<n; i++) {
    sdt->rgap+=sdt->s[i];
    sdt->dobj+=sdt->y[i];
  }  
  sdt->pobj=sdt->dobj+sdt->rgap;
} /* McutInitVal */

static void BcutInitVal(sdpdat *sdt)
{
  int i,n=sdt->ncol;
  
  sdt->rgap=0.0;
  sdt->dobj=0.0;
  
  for (i=0; i<n; i++) {
    sdt->y[i] =min(sdt->y[i],-1.0);
    sdt->s[i] =sdt->cy->diag[i]-sdt->y[i];
    sdt->pobj+=sdt->cy->diag[i];
    sdt->dobj+=sdt->y[i];
  }
  sdt->pobj=sdt->pobj*0.9;
  sdt->rgap=sdt->pobj-sdt->dobj;
  
} /* BcutInitVal */

static void EcutInitVal(sdpdat *sdt)
{
  int    j,k,n=sdt->ncol;
  double r,rtmp;
  array *cj;
  
  sdt->lamda  =0.0;
  sdt->bigM   =(double)sdt->par.bmfac*n;
  sdt->x0     =1.0;
  sdt->s0     =sdt->bigM;
  sdt->kap    =0;
  
  r           =-1.0/(double)n;
  rtmp        =0.0;
  for (j=0; j<n; j++) {
    cj        =sdt->cy->roff+j;
    for (k=0; k<cj->nn0; k++)
      rtmp   +=2.0*cj->an[k];
  }
  
  rtmp       *=r;
  sdt->dobj   =0.0;
  for (j=0; j<n; j++) {
    rtmp     +=sdt->s[j];
    sdt->dobj+=sdt->y[j];
  }
  sdt->rgap   =rtmp+(double)sdt->bigM;
  sdt->pobj   =sdt->dobj+sdt->rgap;
} /* EcutInitVal */

static void UcutInitVal(sdpdat *sdt)
{
  int    i,j,k,n;
  double *x,*z;
  array  *cj;
  
  n=sdt->ncol;
  x=sdt->rw;
  z=x+n;
  k=(int)(n+sdt->kap)/2;
  
  for (j=0; j<k; j++) {
    x[j]=1.0;
    z[j]=sdt->cy->diag[j]*x[j];
  }
  for (j=k; j<n; j++) {
    x[j]=-1.0;
    z[j]=sdt->cy->diag[j]*x[j];
  }
  
  for (j=0; j<n; j++) {
    cj=sdt->cy->roff+j;
    for (k=0; k<cj->nn0; k++) {
      i    =cj->ja[k];
      z[j]+=cj->an[k]*x[i];
      z[i]+=cj->an[k]*x[j];
    }
    sdt->pobj+=z[j]*x[j];
    sdt->dobj+=sdt->y[j];
  }
  
  sdt->lamda =0.0;
  sdt->bigM  =0.0;
  sdt->x0    =0.0;
  sdt->s0    =sdt->bigM;
  sdt->kap  *=sdt->kap;
  sdt->rgap  =sdt->pobj-sdt->dobj;
  /*
  r           =(sdt->kap-n)/(double)(n*(n-1));
  rtmp        =0.0;
  for (j=0; j<n; j++) {
    cj        =sdt->cy->roff+j;
    for (k=0; k<cj->nn0; k++)
      rtmp   +=2.0*cj->an[k];
  }
  
  rtmp       *=r;
  sdt->dobj   =0.0;
  for (j=0; j<n; j++) {
    rtmp     +=sdt->s[j];
    sdt->dobj+=sdt->y[j];
  }
  
  sdt->rgap   =rtmp;
  sdt->pobj   =sdt->dobj+sdt->rgap;
  */
} /* EcutInitVal */

static void DcutInitVal(sdpdat *sdt)
{
  int    j,k,n=sdt->ncol;
  double rtmp;
  array *cj;
  
  sdt->lamda  =0.0;
  sdt->bigM   =(double)sdt->par.bmfac*n;
  sdt->x0     =1.0;
  sdt->s0     =sdt->bigM;
  sdt->kap    =0;
  rtmp        =0.0;
  
  cj          =sdt->cy->roff+sdt->is;
  for (k=0; k<cj->nn0; k++)
    if (cj->ja[k]==sdt->it)
      rtmp    =-2.0*cj->an[k];
  
  sdt->dobj   =0.0;
  for (j=0; j<n; j++) {
    rtmp     +=sdt->s[j];
    sdt->dobj+=sdt->y[j];
  }
  sdt->rgap   =rtmp+sdt->bigM;
  sdt->pobj   =sdt->dobj+sdt->rgap;
} /* DcutInitVal */

void InitSet(sdpdat *sdt)
{
  int    i,j,k,nnz,n=sdt->ncol;
  double spar,rtmp;
  array  *aj;
  symat  *cy=sdt->cy;
      
  for (j=0; j<n; j++) {
    sdt->y[j]=cy->diag[j]-1.0; 
    sdt->s[j]=cy->diag[j];
  } 

  for (j=0; j<n; j++) {
    aj =cy->roff+j;
    nnz=aj->nn0;
    for (k=0; k<nnz; k++) {
      i=aj->ja[k];
      sdt->y[j]-=fabs(aj->an[k]);
      sdt->y[i]-=fabs(aj->an[k]);
    }
  }
  
  for (j=0; j<n; j++)
    sdt->s[j]-=sdt->y[j];
    
  switch (sdt->ptyp) {
    case MaxCut:
      McutInitVal(sdt);
      sdt->rho =(double)5.0*n;
      break;
    case BoxCut:
      BcutInitVal(sdt);
      sdt->rho =(double)5.0*n;
      break;
    case EquCut:
      EcutInitVal(sdt);
      sdt->rho =(double)5.0*n;
      break;
    case UquCut:
      UcutInitVal(sdt);
      sdt->rho =(double)5.0*n;
      break;
    case DvdCut:
      DcutInitVal(sdt);
      sdt->rho =(double)5.0*n;
      break;
    default:
      ExitProc(SysError,"problem type error");
      break;
  }
  
  spar    =100.0*((2.0*sdt->cy->nnzo+n)/(double)(n))/(double)n;
  rtmp    =sdt->rgap/(1.0+fabs(sdt->dobj));
  
  ShowMsgTit(n,sdt,spar);
  ShowMsg(sdt,-1,rtmp,0.0);  
} /* InitSet */

