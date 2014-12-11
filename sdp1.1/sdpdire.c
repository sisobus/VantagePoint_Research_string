#include "SDPdef.h"

static void FormMcutDy(sdpdat *sdt, 
                       double *dy1,
                       double *dy2,
                       double *rw)
{
  int    i,n=sdt->ncol;
  double *b1,*b2;
  
  b1=rw;
  b2=rw+n;
  
  for (i=0; i<n; i++) {
    b1[i]=1.0;
    b2[i]=sdt->sinv[i];
  }
  
  ChlSolve(sdt->mf,b1,dy1);
  ChlSolve(sdt->mf,b2,dy2);  
} /* FormMcutDy */

static void FormBcutDy(sdpdat *sdt, 
                       double *dy1,
                       double *dy2,
                       double *rw)
{
  int    i,n=sdt->ncol;
  double *b1,*b2;
  
  b1=rw;
  b2=rw+n;
  
  for (i=0; i<n; i++) {
    b1[i]=1.0;
    b2[i]=sdt->sinv[i]-1.0/sdt->y[i];
  }
  
  ChlSolve(sdt->mf,b1,dy1);
  ChlSolve(sdt->mf,b2,dy2);  
} /* FormBcutDy */

static void FormEcutDy(sdpdat *sdt,
                       double *dy1,
                       double *dy2,
                       double *ese,
                       double *dl1,
                       double *dl2,
                       double *rwk)
{
  int    i,n=sdt->ncol;
  double *b,*p,*q,*u,beta,rtmp;
  chfac  *mf=sdt->mf;
  
  b=rwk;
  p=b+n;
  u=p+n;
  q=dy2;
  
  (*ese)   =0.0;
  for (i=0; i<n; i++) {
    (*ese)+=u[i];
    u[i]  *=u[i];
    b[i]   =1.0;
    dy1[i] =u[i];
  }
  
  ChlSolve(mf,dy1,q);
  ChlSolve(mf,b,p);
  
  rtmp=(*ese)*(*ese);
  rtmp+=1.0/(sdt->s0*sdt->s0);
    
  beta=sdt->kap;
  for (i=0; i<n; i++) {
    beta-=u[i]*p[i];
    rtmp-=u[i]*q[i];
  }
  
  (*dl1)=beta/rtmp;
  
  for (i=0; i<n; i++) {
    dy1[i]=p[i]-(*dl1)*q[i];
    b[i]  =sdt->sinv[i];
  }
  
  ChlSolve(mf,b,p);
  
  beta=(*ese)-1.0/sdt->s0;  
  for (i=0; i<n; i++)
    beta-=u[i]*p[i];
  
  (*dl2)=beta/rtmp;
  
  for (i=0; i<n; i++) 
    dy2[i]=p[i]-(*dl2)*q[i];
  
} /* FormEcutDy */

static void FormUcutDy(sdpdat *sdt,
                       double *dy1,
                       double *dy2,
                       double *ese,
                       double *dl1,
                       double *dl2,
                       double *rwk)
{
  int    i,n=sdt->ncol;
  double *b,*p,*q,*u,beta,rtmp;
  chfac  *mf=sdt->mf;
  
  b=rwk;
  p=b+n;
  u=p+n;
  q=dy2;
  
  (*ese)   =0.0;
  for (i=0; i<n; i++) {
    (*ese)+=u[i];
    u[i]  *=u[i];
    b[i]   =1.0;
    dy1[i] =u[i];
  }
  
  ChlSolve(mf,dy1,q);
  ChlSolve(mf,b,p);
  
  rtmp=(*ese)*(*ese);    
  beta=sdt->kap;
  for (i=0; i<n; i++) {
    beta-=u[i]*p[i];
    rtmp-=u[i]*q[i];
  }
  
  (*dl1)=beta/rtmp;
  
  for (i=0; i<n; i++) {
    dy1[i]=p[i]-(*dl1)*q[i];
    b[i]  =sdt->sinv[i];
  }
  
  ChlSolve(mf,b,p);
  
  beta=(*ese);  
  for (i=0; i<n; i++)
    beta-=u[i]*p[i];
  
  (*dl2)=beta/rtmp;
  
  for (i=0; i<n; i++) 
    dy2[i]=p[i]-(*dl2)*q[i];
  
} /* FormUcutDy */

static void FormDcutDy(sdpdat *sdt,
                       double *dy1,
                       double *dy2,
                       double *ese,
                       double *dl1,
                       double *dl2,
                       double *rwk)
{
  int    i,n=sdt->ncol;
  double *b,*p,*q,*u,beta,rtmp;
  chfac  *mf=sdt->mf;
  
  b=rwk;
  p=b+n;
  u=p+n;
  q=dy2;
  
  (*ese)   =u[sdt->is]+u[sdt->it];
  for (i=0; i<n; i++) {
    u[i]  *=u[i];
    b[i]   =1.0;
    dy1[i] =u[i];
  }
  
  ChlSolve(mf,dy1,q);
  ChlSolve(mf,b,p);
  
  rtmp=(*ese)*(*ese);
  rtmp+=1.0/(sdt->s0*sdt->s0);
    
  beta=sdt->kap;
  for (i=0; i<n; i++) {
    beta-=u[i]*p[i];
    rtmp-=u[i]*q[i];
  }
  
  (*dl1)=beta/rtmp;
  
  for (i=0; i<n; i++) {
    dy1[i]=p[i]-(*dl1)*q[i];
    b[i]  =sdt->sinv[i];
  }
  
  ChlSolve(mf,b,p);
  
  beta=(*ese)-1.0/sdt->s0;  
  for (i=0; i<n; i++)
    beta-=u[i]*p[i];
  
  (*dl2)=beta/rtmp;
  
  for (i=0; i<n; i++) 
    dy2[i]=p[i]-(*dl2)*q[i];
  
} /* FormDcutDy */

void FormDy(sdpdat *sdt,
            double *dy1,
            double *dy2,
            double *ese,
            double *dl1,
            double *dl2,
            double *rwk)
{
  switch (sdt->ptyp) {
    case MaxCut:
      FormMcutDy(sdt,dy1,dy2,rwk);
      break;
    case BoxCut:
      FormBcutDy(sdt,dy1,dy2,rwk);
      break;
    case EquCut:
      FormEcutDy(sdt,dy1,dy2,ese,dl1,dl2,rwk);
      break;
    case UquCut:
      FormUcutDy(sdt,dy1,dy2,ese,dl1,dl2,rwk);
      break;
    case DvdCut:
      FormDcutDy(sdt,dy1,dy2,ese,dl1,dl2,rwk);
      break;
    default:
      ExitProc(SysError,"problem type error");
      break;
  }
} /* FormDy */
