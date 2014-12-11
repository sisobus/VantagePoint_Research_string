#include "SDPdef.h"

double HouseTrans(int    n,
                  double *x,
                  double *v)
{
  int    i;
  double r,xmax,beta;
  
  xmax=0.0;
  for (i=0; i<n; i++) {
    r=fabs(x[i]);
    if (xmax<r)
      xmax=r;
  }
  
  r=0.0;
  for (i=0; i<n; i++) {
    v[i] =x[i]/xmax;
    r   +=v[i]*v[i];
  }
  
  r    =sqrt(r);
  beta =1.0/(r*(r+fabs(v[0])));
  v[0]+=sign(v[0])*r;
  
  return beta;
} /* HouseTrans */

int ChkPosEig(chfac  *sf,
              double lamda,
              double *rwk)
{
  int    i,n=sf->nrow;
  double *e,*v,*u,r,p;
  
  e=rwk;
  u=e+n;
  v=e;
  
  for (i=0; i<n; i++)
    e[i]=1.0;
    
  ForwSubst(sf,e,u);
  
  for (i=0; i<n; i++)
    e[i]=u[i]*sqrt(sf->diag[i]);
    
  for (i=0; i<n; i++)
    u[i]=e[sf->invp[i]];
    
  r=HouseTrans(n,u,v);
    
  r*=v[0];
  p =0.0;
  for (i=0; i<n; i++)
    p+=u[i]*v[i];
    
  p =u[0]-r*p;
  p*=p;
    
  if (1.0+lamda*p>1.0e-13)
    return true;
    
  else 
    return false;
} /* ChkPosEig */

static int isPosDef(chfac  *sf,
                    double rtmp,
                    double *rwk)
{
  int    i,j=0,n=sf->nrow;
  double *e,*u,r;
  
  e=rwk;
  u=e+n;
  
  for (i=0; i<n; i++)
    e[i]=1.0;
    
  ForwSubst(sf,e,u);
  
  for (i=0; i<n; i++) {
    r=sf->diag[i];
    if (r<1.0e-13) {
      j=i;
      r=fabs(r);
    }
    e[i]=u[i]*sqrt(r);
  }
  
  r=rtmp*(2*e[j]*e[j]-dDot(e,e,n))-1.0;

  if (r>0) return true;
  else return false;  
} /* isPosDef */

int ChkPosDef1(sdpdat *sdt,
               double *y,
               double lamda,
               double *dy,
               double dl,
               double *rwk)
{
  int    i,j,k,n=sdt->ncol;
  double rtmp,r,p;
  symat  *cy=sdt->cy;
  syoff  *st=sdt->st;
  array  *cj,*sj;
  chfac  *sf=sdt->sf,*mf=sdt->mf;
  
  /*
   * sparse check
   */
  rtmp=dl-lamda;
  
  if (sdt->ptyp==DvdCut) {
    for (j=0; j<n; j++)
      sf->diag[sf->invp[j]]=cy->diag[j]+(dy[j]-y[j]);
    dCopy(sf->unnz,rwk+3*n,sf->uval);
    
    sf->diag[sf->invp[sdt->is]]+=rtmp;
    sf->diag[sf->invp[sdt->it]]+=rtmp;
    sf->uval[sf->upst]         +=rtmp;
    
    if (CfcOk==ChlFact(sf,sdt->iw,rwk+n,true))
      return true;
    else
      return false; 
  }
  
  else if (sdt->ptyp==BoxCut) {
    for (j=0; j<n; j++) {
      rtmp=dy[j]-y[j];
      
      if (rtmp<1.0e-13)
        return false;
        
      sf->diag[sf->invp[j]]=cy->diag[j]+rtmp;
    }
    
    dCopy(sf->unnz,rwk+3*n,sf->uval);
  
    if (CfcOk==ChlFact(sf,sdt->iw,rwk+n,true))
      return true;
    else
      return false; 
  }
  
  for (j=0; j<n; j++) {
    p                    =cy->diag[j]+(dy[j]-y[j]);
    sf->diag[sf->invp[j]]=p;
    mf->diag[j]          =p+rtmp;
  }  
  dCopy(sf->unnz,rwk+3*n,sf->uval);
  
  if (sdt->ptyp==MaxCut) {
    if (CfcOk==ChlFact(sf,sdt->iw,rwk+n,true))
      return true;
    else
      return false; 
  }
  
  if (CfcOk==ChlFact(sf,sdt->iw,rwk+n,false)) {
    r=1.0e40;
    for (j=0; j<n; j++) {
      if (sf->diag[j]<1.0e-13&&r<1.0e-13)
        return false;
      
      if (r>sf->diag[j])
        r=sf->diag[j];
    }
    
    if (rtmp<1.0e-13&&r<1.0e-13)
      return false;
    
    else if (r>=1.0e-13)
      return (ChkPosEig(sf,rtmp,rwk));
    
    else
      return (isPosDef(sf,rtmp,rwk));
  }
  
  /*
   * dense check
   */
  sdt->dnck++;
  for (j=0; j<n; j++) {
    for (i=j+1; i<n; i++)
      rwk[i]=rtmp;
    
    cj=cy->roff+j;
    for (k=0; k<cj->nn0; k++)
      rwk[cj->ja[k]]+=cj->an[k];
    
    sj=st->roff+j;
    for (k=0; k<sj->nn0; k++)
      sj->an[k]=rwk[sj->ja[k]];
  }
  if (CfcOk!=ChlFact(mf,sdt->iw,rwk,true))
    return false;
  
  return true;
  
} /* ChkPosDef1 */

int ChkPosDef2(sdpdat *sdt,
               double *y,
               double lamda,
               double *rwk)
{
  int    i,j,k,n=sdt->ncol;
  double rtmp,r;
  symat  *cy=sdt->cy;
  syoff  *st=sdt->st;
  array  *cj,*sj;
  chfac  *sf=sdt->sf,*mf=sdt->mf;
  
  /*
   * sparse check
   */
  rtmp=-lamda;
  
  if (sdt->ptyp==DvdCut) {
    for (j=0; j<n; j++)
      sf->diag[sf->invp[j]]=cy->diag[j]-y[j];
    dCopy(sf->unnz,rwk+3*n,sf->uval);
    
    sf->diag[sf->invp[sdt->is]]+=rtmp;
    sf->diag[sf->invp[sdt->it]]+=rtmp;
    sf->uval[sf->upst]         +=rtmp;
    
    if (CfcOk==ChlFact(sf,sdt->iw,rwk+n,true))
      return true;
    else
      return false; 
  }
  
  for (j=0; j<n; j++) {
    r                    =cy->diag[j]-y[j];
    sf->diag[sf->invp[j]]=r;
    mf->diag[j]          =r+rtmp;
  }
  dCopy(sf->unnz,rwk+3*n,sf->uval);
  
  if (sdt->ptyp==MaxCut||sdt->ptyp==BoxCut) {
    if (CfcOk==ChlFact(sf,sdt->iw,rwk+n,true))
      return true;
    else
      return false;
  }
  
  if (CfcOk==ChlFact(sf,sdt->iw,rwk+n,false)) {
    r=1.0e40;
    for (j=0; j<n; j++) {
      if (sf->diag[j]<=1.0e-13&&r<=1.0e-13)
        return false;
      
      if (r>sf->diag[j])
        r=sf->diag[j];
    }
    
    if (rtmp<=1.0e-13&&r<=1.0e-13)
      return false;
    
    else if (r>1.0e-13)
      return (ChkPosEig(sf,rtmp,rwk));
    
    else
      return (isPosDef(sf,rtmp,rwk));
  }
  
  /*
   * dense check
   */
  sdt->dnfc++;
  
  for (j=0; j<n; j++) {
    
    for (i=j+1; i<n; i++)
      rwk[i]=rtmp;
    
    cj=cy->roff+j;
    for (k=0; k<cj->nn0; k++)
      rwk[cj->ja[k]]+=cj->an[k];
    
    sj=st->roff+j;
    for (k=0; k<sj->nn0; k++)
      sj->an[k]=rwk[sj->ja[k]];
  }
  
  if (CfcOk!=ChlFact(mf,sdt->iw,rwk+n,true))
    return false;
  
  return true;
  
} /* ChkPosDef */ 
