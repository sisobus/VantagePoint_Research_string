#include "SDPdef.h"

static void copyChl(chfac *af,
                    chfac *bf)
{
  iFree(&af->shead);
  iFree(&af->ssize);
  iFree(&af->ssub);
  
  iFree(&bf->shead);
  iFree(&bf->ssize);
  iFree(&bf->ssub);
  
  bf->unnz =af->unnz;
  bf->ujnz =af->ujnz;
  bf->nsnds=af->nsnds;
  bf->ndens=af->ndens;
  bf->nsndn=af->nsndn;
  bf->sdens=af->sdens;
  bf->upst =af->upst;
  
  bf->usub =iAlloc(bf->ujnz,"usub, copyChl");
  bf->uval =dAlloc(bf->unnz,"uval, copyChl");
  
  iCopy(af->ujnz,af->usub,bf->usub);
  iCopy(af->nrow,af->ujbeg,bf->ujbeg);
  iCopy(af->nrow,af->uhead,bf->uhead);
  iCopy(af->nrow,af->ujsze,bf->ujsze);
  iCopy(af->nrow,af->perm,bf->perm);
  iCopy(af->nrow,af->invp,bf->invp);
  iCopy(af->nrow+1,af->subg,bf->subg);
  
  bf->dhead=iAlloc(af->ndens+1,"dhead, copyChl");
  iCopy(af->ndens+1,af->dhead,bf->dhead);
  
  if (af->nsnds) {
    bf->dbeg=iAlloc(af->nsnds,"dbeg, copyChl");
    bf->dsub=iAlloc(af->nsnds,"dsub, copyChl");
    iCopy(af->nsnds,af->dbeg,bf->dbeg);
    iCopy(af->nsnds,af->dsub,bf->dsub);
  }
} /* copyChl */

void GetLhat(sdpdat *sdt)
{
  int    i,j,k,n,itmp=0;
  double lamda,dl,rtmp,*q,*y,*dy,*uv;
  chfac  *cf,*mf;
  array  *cj,*sj;
  symat  *cy;
  syoff  *st;
  
  n=sdt->nrow;


  if ( sdt->ptyp!=UquCut) 
  {
    sdt->mf=NULL;
    CfcFree(&sdat->mf);
    SyoFree(&sdat->st);
    /* Reassign L factor for rank reduction */
    sdt->mf=sdt->mf=CfcAlloc(n,"mf, GetLhat");
    copyChl(sdt->sf,sdt->mf);
  }
  
  lamda=sdt->lamda;
  dl   =sdt->dl;
  cf   =sdt->sf;
  mf   =sdt->mf;
  y    =sdt->y;
  dy   =sdt->dy;
  q    =sdt->rw;
  uv   =sdt->rw+5*n;
  cy   =sdt->cy;
  st   =sdt->st;
  
  /*
   * find L such that LL^T=S+Diag(dy)+dl*e*e^T
   */
  if (sdt->ptyp==UquCut) 
  {
    rtmp=dl-lamda;
    
    for (j=0; j<n; j++)
      mf->diag[j]=cy->diag[j]+(dy[j]-y[j])+rtmp;
    
    for (j=0; j<n; j++) 
	{
      for (i=j+1; i<n; i++)
        q[i]=rtmp;
    
      cj=cy->roff+j;
      for (k=0; k<cj->nn0; k++)
        q[cj->ja[k]]+=cj->an[k];
    
      sj=st->roff+j;
      for (k=0; k<sj->nn0; k++)
        sj->an[k]=q[sj->ja[k]];
    }
  }  
  else
  {
	  for (j=0; j<n; j++)
	  {
         mf->diag[mf->invp[j]]=cy->diag[j]+(dy[j]-y[j]);
	  }

      dCopy(sdt->sf->unnz,uv,mf->uval);
   
      if (sdt->ptyp==DvdCut) 
	  {
        rtmp=dl-lamda;
        mf->diag[mf->invp[sdt->is]]+=rtmp;
        mf->diag[mf->invp[sdt->it]]+=rtmp;
        mf->uval[mf->upst]         +=rtmp;
	  }
  }


  /* If return is not 0 */
  if( sdt->ptyp != EquCut )
  {
	  if (CfcOk!=ChlFact(mf,sdt->iw,q,true))
	  { 
        ExitProc(SysError,NULL);    
	  }
  }
  else ChlFact(mf,sdt->iw,q,false);

  mf->NegDiag =-1;

  if( sdt->ptyp == EquCut )
  {
    for(j=0;j<n;j++) 
	{
	  
	  if( mf->diag[mf->invp[j]]<0 ) /* Why invp? */
	  {
		 
		  if( mf->NegDiag <0 ) 
		  {
			  mf->NegDiag = mf->invp[j];
			  itmp++;
		  }
          else mf->NegDiag =-2;
	  }
	}
  }

  if(mf->NegDiag<-1) 
  {
	  printf("\n Error in ecut rank reduction.");
      printf("\n More than one Negative number in diagonal.");
	  printf("\n There are %d negitives.",itmp);
      ShutDown();
  }

} /* GetLhat */

void GetVhat(double  *u,
             double  *q,
             sdpdat  *sdt)
{
  int    i,n;
  double ese,lamda,dl,rtmp,*r,*y,*dy,*uv;
  chfac  *cf,*mf;
  symat  *cy;
  syoff  *st;
  
  n    =sdt->nrow;
  lamda=sdt->lamda;
  dl   =sdt->dl;
  cf   =sdt->sf;
  mf   =sdt->mf;
  y    =sdt->y;
  dy   =sdt->dy;
  r    =sdt->rw+n;
  uv   =sdt->rw+5*n;
  cy   =sdt->cy;
  st   =sdt->st;
  
  /*
   * r=L*u;
   */
  GetUhat(mf,u,r);

  /*
   * r=S^-1*u
   */
  ChlSolve(cf,u,r);
  dCopy(n,r,u);
  
  if (sdt->ptyp==EquCut||sdt->ptyp==UquCut) 
  {
    rtmp=0.0;
    ese =0.0;
    for (i=0; i<n; i++) 
	{
      rtmp+=r[i];
      ese +=q[i];
    }
    ese=1.0-lamda*ese;
    
    rtmp=lamda*rtmp/ese;
    for (i=0; i<n; i++)
      u[i]=r[i]+rtmp*q[i];
  }  
} /* GetVhat */

