#include "SDPdef.h"

static void GetMcutVal(sdpdat *sdt,
                       double *dy,
                       double *dy1,
                       double *dy2,
                       double *gnw)
{
  int    i,n=sdt->ncol;
  double rtmp,pval;
  
  pval     =0.0;
  *gnw     =0.0;
  rtmp     =sdt->rho/sdt->rgap;
  
  for (i=0; i<n; i++) {
    dy[i]  =rtmp*dy1[i]-dy2[i];
    pval  +=(rtmp-sdt->sinv[i])*dy[i];
    (*gnw)+=dy[i]*sdt->sinv[i];
  }
  
  sdt->pval=sqrt(pval);
  (*gnw)   =((*gnw)+n)/rtmp;
} /* GetMcutVal */

static void GetBcutVal(sdpdat *sdt,
                       double *dy,
                       double *dy1,
                       double *dy2,
                       double *gnw)
{
  int    i,n=sdt->ncol;
  double rtmp,pval,invy;
  
  pval     =0.0;
  *gnw     =0.0;
  rtmp     =sdt->rho/sdt->rgap;
  
  for (i=0; i<n; i++) {
    invy   =1.0/sdt->y[i];
    dy[i]  =rtmp*dy1[i]-dy2[i];
    pval  +=(rtmp-sdt->sinv[i]+invy)*dy[i];
    (*gnw)+=dy[i]*(sdt->sinv[i]-invy);
  }
  
  sdt->pval=sqrt(pval);
  (*gnw)   =((*gnw)+2*n)/rtmp;
} /* GetBcutVal */

static void GetEcutVal(sdpdat *sdt,
                       double *dy,
                       double *dy1,
                       double *dy2,
                       double ese,
                       double *dl,
                       double dl1,
                       double dl2,
                       double *gapnew)
{
  int    i,n=sdt->ncol;
  double rtmp,pval,gpnw;
  
  rtmp      =sdt->rho/sdt->rgap;
  (*dl)     =rtmp*dl1-dl2;
  pval      =(rtmp*sdt->kap-(ese-1.0/sdt->s0))*(*dl);
  gpnw      =(ese-1.0/sdt->s0)*(*dl)+(double)n;
  
  for (i=0; i<n; i++) {
    dy[i]   =rtmp*dy1[i]-dy2[i];
    pval   +=(rtmp-sdt->sinv[i])*dy[i];
    gpnw   +=sdt->sinv[i]*dy[i];
  }
  sdt->pval =sqrt(pval);
  (*gapnew) =gpnw/rtmp+sdt->x0*sdt->s0;
} /* GetEcutVal */

static void GetUcutVal(sdpdat *sdt,
                       double *dy,
                       double *dy1,
                       double *dy2,
                       double ese,
                       double *dl,
                       double dl1,
                       double dl2,
                       double *gapnew)
{
  int    i,n=sdt->ncol;
  double rtmp,pval,gpnw;
  
  rtmp      =sdt->rho/sdt->rgap;
  (*dl)     =rtmp*dl1-dl2;
  pval      =(rtmp*sdt->kap-ese)*(*dl);
  gpnw      =ese*(*dl)+(double)n;
  
  for (i=0; i<n; i++) {
    dy[i]   =rtmp*dy1[i]-dy2[i];
    pval   +=(rtmp-sdt->sinv[i])*dy[i];
    gpnw   +=sdt->sinv[i]*dy[i];
  }
  sdt->pval =sqrt(pval);
  (*gapnew) =gpnw/rtmp;
} /* GetUcutVal */

static void GetDcutVal(sdpdat *sdt,
                       double *dy,
                       double *dy1,
                       double *dy2,
                       double ese,
                       double *dl,
                       double dl1,
                       double dl2,
                       double *gapnew)
{
  int    i,n=sdt->ncol;
  double rtmp,pval,gpnw;
  
  rtmp      =sdt->rho/sdt->rgap;
  (*dl)     =rtmp*dl1-dl2;
  pval      =-(ese-1.0/sdt->s0)*(*dl);
  gpnw      =(ese-1.0/sdt->s0)*(*dl)+(double)n;
  
  for (i=0; i<n; i++) {
    dy[i]   =rtmp*dy1[i]-dy2[i];
    pval   +=(rtmp-sdt->sinv[i])*dy[i];
    gpnw   +=sdt->sinv[i]*dy[i];
  }
  sdt->pval =sqrt(pval);
  (*gapnew) =gpnw/rtmp+sdt->x0*sdt->s0;
} /* GetDcutVal */

void GetValues(sdpdat *sdt,
               double *dy,
               double *dy1,
               double *dy2,
               double ese,
               double *dl,
               double dl1,
               double dl2,
               double *gapnew)
{
  switch (sdt->ptyp) {
    case MaxCut:
      GetMcutVal(sdt,dy,dy1,dy2,gapnew);
      break;
    case BoxCut:
      GetBcutVal(sdt,dy,dy1,dy2,gapnew);
      break;
    case EquCut:
      GetEcutVal(sdt,dy,dy1,dy2,ese,dl,dl1,dl2,gapnew);
      break;
    case UquCut:
      GetUcutVal(sdt,dy,dy1,dy2,ese,dl,dl1,dl2,gapnew);
      break;
    case DvdCut:
      GetDcutVal(sdt,dy,dy1,dy2,ese,dl,dl1,dl2,gapnew);
      break;
    default:
      ExitProc(SysError,"problem type error");
      break;
  }
} /* GetValues */
