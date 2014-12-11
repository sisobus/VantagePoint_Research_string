#include "SDPdef.h"

void ModMcutUpp(sdpdat *sdt,
                double *dy,
                double *dy1,
                double *dy2,
                double *pnew,
                double gapnew)
{
  int    i,n=sdt->nrow;
  double rtmp;
  
  (*pnew)=0.0;
  rtmp   =sdt->rho/gapnew;
  for (i=0; i<n; i++) {
    dy[i]   =rtmp*dy1[i]-dy2[i];
    (*pnew)+=(rtmp-sdt->sinv[i])*dy[i];
  }
  (*pnew)=sqrt((*pnew));
  
} /* ModMcutUpp */

void ModBcutUpp(sdpdat *sdt,
                double *dy,
                double *dy1,
                double *dy2,
                double *pnew,
                double gapnew)
{
  int    i,n=sdt->nrow;
  double rtmp;
  
  (*pnew)=0.0;
  rtmp   =sdt->rho/gapnew;
  for (i=0; i<n; i++) {
    dy[i]   =rtmp*dy1[i]-dy2[i];
    (*pnew)+=(rtmp-sdt->sinv[i]+1.0/sdt->y[i])*dy[i];
  }
  (*pnew)=sqrt((*pnew));
  
} /* ModBcutUpp */

void ModEcutUpp(sdpdat *sdt,
                double *dy,
                double *dy1,
                double *dy2,
                double *dl,
                double dl1,
                double dl2,
                double *pnw,
                double gnw)
{
  int    i,n=sdt->ncol;
  double rtmp;
  
  sdt->x0  =sdt->rgap/sdt->rho*(sdt->s0-(*dl))/(sdt->s0*sdt->s0);
  rtmp     =sdt->rho/gnw;
  (*dl)    =rtmp*dl1-dl2;
  (*pnw)   =-(*dl)*((*pnw)-1.0/sdt->s0);
  for (i=0; i<n; i++) {
    dy[i]  =rtmp*dy1[i]-dy2[i];
    (*pnw)+=(rtmp-sdt->sinv[i])*dy[i];
  }
  (*pnw)   =sqrt((*pnw));
} /* ModEcutUpp */

void ModUcutUpp(sdpdat *sdt,
                double *dy,
                double *dy1,
                double *dy2,
                double *dl,
                double dl1,
                double dl2,
                double *pnw,
                double gnw)
{
  int    i,n=sdt->ncol;
  double rtmp;
  
  rtmp     =sdt->rho/gnw;
  (*dl)    =rtmp*dl1-dl2;
  (*pnw)   =(rtmp*sdt->kap-(*pnw))*(*dl);
  for (i=0; i<n; i++) {
    dy[i]  =rtmp*dy1[i]-dy2[i];
    (*pnw)+=(rtmp-sdt->sinv[i])*dy[i];
  }
  (*pnw)   =sqrt((*pnw));
} /* ModUcutUpp */

void ModDcutUpp(sdpdat *sdt,
                double *dy,
                double *dy1,
                double *dy2,
                double *dl,
                double dl1,
                double dl2,
                double *pnw,
                double gnw)
{
  int    i,n=sdt->ncol;
  double rtmp;
  
  sdt->x0  =sdt->rgap/sdt->rho*(sdt->s0-(*dl))/(sdt->s0*sdt->s0);
  rtmp     =sdt->rho/gnw;
  (*dl)    =rtmp*dl1-dl2;
  (*pnw)   =-(*dl)*((*pnw)-1.0/sdt->s0);
  for (i=0; i<n; i++) {
    dy[i]  =rtmp*dy1[i]-dy2[i];
    (*pnw)+=(rtmp-sdt->sinv[i])*dy[i];
  }
  (*pnw)   =sqrt((*pnw));
} /* ModDcutUpp */

void ModUppBnd(sdpdat *sdt,
               double *dy,
               double *dy1,
               double *dy2,
               double *dl,
               double dl1,
               double dl2,
               double *pnw,
               double gnw)
{
  switch (sdt->ptyp) {
    case MaxCut:
      ModMcutUpp(sdt,dy,dy1,dy2,pnw,gnw);
      break;
    case BoxCut:
      ModBcutUpp(sdt,dy,dy1,dy2,pnw,gnw);
      break;
    case EquCut:
      ModEcutUpp(sdt,dy,dy1,dy2,dl,dl1,dl2,pnw,gnw);
      break;
    case UquCut:
      ModUcutUpp(sdt,dy,dy1,dy2,dl,dl1,dl2,pnw,gnw);
      break;
    case DvdCut:
      ModDcutUpp(sdt,dy,dy1,dy2,dl,dl1,dl2,pnw,gnw);
      break;
    default:
      ExitProc(SysError, "problem type error");
      break;
  }
} /* ModUppBnd */

int GetStep(sdpdat *sdt,
            double *y,
            double *dy,
            double lamda,
            double dl,
            double alpha,
            double *lstp,
            int    *bstp,
            double *rwk)
{
  int    i,n=sdt->ncol,bpd,id;
  double *ynew,lnew,beta,rtmp;
  
  ynew=rwk;
  beta=alpha/sdt->pval;
  id  =false;
  
  if (sdt->par.varstep) {
    do {
      rtmp=(*lstp)*beta;
    
      lnew=lamda+rtmp*dl;
      for (i=0; i<n; i++)
        ynew[i]=y[i]+rtmp*dy[i];
            
      bpd=(sdt->bigM+lnew>1.0e-13)||
          (sdt->ptyp!=EquCut&&sdt->ptyp!=DvdCut);
      
      if (sdt->ptyp==BoxCut) {
        for (i=0; i<n; i++)
          bpd=bpd&&(ynew[i]<-1.0e-13);
      }
      
      if (bpd) {
        if (ChkPosDef2(sdt,ynew,lnew,rwk)) {
          sdt->step=rtmp;
          id       =true;
          (*bstp)  ++;
          if ((*bstp)>2) {
            (*lstp)=min(2.0*(*lstp),6.0);
            (*bstp)=0;
          }
          break;
        }
      }
    
      if ((*lstp)<=1.5) {
        lnew=lamda+beta*dl;
        
        bpd=(sdt->bigM+lnew>1.0e-13)||
            (sdt->ptyp!=EquCut&&sdt->ptyp!=DvdCut);
        
        if (bpd) sdt->step=beta;
        else     sdt->step=0.9999*(sdt->bigM+lamda)/dl;
        (*lstp)  =1.5;
        break;
      }
    
      else {
        (*lstp)  =max(0.5*(*lstp),1.5);
        (*bstp)  =0;
      }
    } while (!id);
  }
  
  else {
    rtmp=(*lstp)*beta;
    
    lnew=lamda+rtmp*dl;
    for (i=0; i<n; i++)
      ynew[i]=y[i]+rtmp*dy[i];
    
    if (sdt->ptyp==BoxCut) {
      for (i=0; i<n; i++) {
        if (ynew[i]>-1.0e-13) {
          sdt->step=beta;
          return id;
        }
      }
    }
    
    if (ChkPosDef2(sdt,ynew,lnew,rwk)) {
      id       =true;
      sdt->step=rtmp;
    }
    
    else sdt->step=beta;
  }
  
  return id;
} /* GetStep */

void UpdSol(sdpdat *sdt,
            double *y,
            double *s,
            double *dy,
            double dl)
{
  int    i,n=sdt->ncol;
  double rtmp,rdl,r=sdt->step;
  
  rdl         =r*dl;
  sdt->lamda +=rdl;
  sdt->s0    +=rdl;
  if (sdt->kap>1.0e-12) {
    rdl      *=sdt->kap;
    sdt->rgap-=rdl;
    sdt->dobj+=rdl;
  }
  for (i=0; i<n; i++) {
    rtmp      =r*dy[i];
    y[i]     +=rtmp;
	s[i]     -=rtmp;
    sdt->rgap-=rtmp;
    sdt->dobj+=rtmp;
  }
  sdt->pobj   =sdt->dobj+sdt->rgap;
} /* UpdSol */
