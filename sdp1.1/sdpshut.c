#include "SDPdef.h"

void ShutDown(void)
{
  clock_t tim;
  double  rtmp;
  
  rtmp=fabs(sdat->pobj-sdat->dobj)/(1.0+fabs(sdat->dobj));
  if (rtmp>sdat->par.tolgap) {
    tim=GetTime()-stim[START];
  
#ifdef PCMACHINE
    printf(" time=%.2f\n\n",(double)tim/CLOCKS_PER_SEC);
#else
    printf(" time=%.2f\n\n",(double)0.01*tim);
#endif

    fclose(fout);
  }

  dFree(&sdat->y);
  dFree(&sdat->dy);
  dFree(&sdat->s);
  dFree(&sdat->sinv);

  SmtFree(&sdat->c);
  SymFree(&sdat->cy);
  SyoFree(&sdat->st);
  
  if (sdat->mf)
	sdat->mf->uval=NULL;
  CfcFree(&sdat->mf);
  CfcFree(&sdat->sf);
  
  dPtFree(&sdat->u);
  dPtFree(&sdat->v);
  
  iFree(&sdat->iw);
  dFree(&sdat->rw);
  
} /* ShutDown */
