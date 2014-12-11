#include "SDPlib.h"

int InvMemo(sdpdat *sdt)
{
  symat *cy=sdt->cy;
  int   i,j,k,nnzo,nrow=cy->nrow;
  array *ci;
  
  nnzo       =nrow*(nrow-1)/2;  
  sdt->st    =SyoAlloc(nrow,nnzo,"st, InverAssign"); 
  sdt->sinv  =dAlloc(nrow,"sinv, InverAssign");
  
  nnzo=0;
  for (i=0; i<nrow; i++) {
    ci     =sdt->st->roff+i;
    ci->nn0=nrow-1-i;
    ci->ja=sdt->st->roff->ja+nnzo;
    ci->an=sdt->st->roff->an+nnzo;
    for (k=0,j=i+1; j<nrow; j++,k++)
      ci->ja[k]=j;
    nnzo+=k;
  }
  
  return true;
} /* InvMemo */

void FindSinv(chfac  *cf,
              syoff  *st,
              double *sinv,
              double *w,
              double *u,
              double lamda,
              prbset ptyp)
{
  int    i,j,k,nn0,n=cf->nrow;
  double *b,*x,sigma;
  array  *sj;
  
  nn0  =0;
  sigma=0.0;
  b    =w;
  x    =w+n;
  
  for (j=0; j<n; j++) {
    u[j]=0.0;
    x[j]=0.0;
    b[j]=0.0;
  }
    
  for (j=0; j<n; j++) {
    b[j]=1.0;
    
    ChlSolve(cf,b,x);
    
    sj     =st->roff+j;
    sj->ja=st->roff->ja+nn0;
    sj->an=st->roff->an+nn0;
        
    sinv[j]=x[j];
    for (k=0; k<sj->nn0; k++) {
      i=sj->ja[k];
      sj->an[k]=x[i];
    }
    
    nn0+=sj->nn0;
    
    if (ptyp==DvdCut) {
      u[j]=x[sdat->is]+x[sdat->it];
      for (i=0; i<n; i++) {
        x[i] =0.0;
        b[i] =0.0;
      }
    }
    else {
      for (i=0; i<n; i++) {
        u[j]+=x[i];
        x[i] =0.0;
        b[i] =0.0;
      }
      sigma+=u[j];
    }
  }
  
  if (ptyp==DvdCut) return;
  
  sigma=lamda/(1.0-lamda*sigma);
  
  for (j=0; j<n; j++) {
    sinv[j]+=sigma*u[j]*u[j];
    sj      =st->roff+j;
    for (k=0; k<sj->nn0; k++)
      sj->an[k]+=sigma*u[j]*u[sj->ja[k]];
  }
  
  for (j=0; j<n; j++)
    u[j]=0.0;
  
  for (j=0; j<n; j++) {
    u[j]+=sinv[j];
    sj   =st->roff+j;
    for (k=0; k<sj->nn0; k++) {
      i    =sj->ja[k];
      u[j]+=sj->an[k];
      u[i]+=sj->an[k];
    }
  }
} /* FindSinv */
