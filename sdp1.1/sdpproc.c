#include "SDPdef.h"

void SdpProc(sdpdat  *sdt,
             clock_t otim[])
{
  int    i,n=sdt->ncol,bprt,bupd,bfac,bstp,bs;
  double *y,*s,*sinv,*dy1,*dy2,*dy,dl1,dl2,dl,
         *rbuf,*u,*uv,*ynew,lnew,
         rgap0,gapnew,pnew,lstp,rtmp,beta,
         ese;
  char   sname[60];
  params *par=&sdt->par;
  chfac  *sf,*mf;
  symat  *cy;
  syoff  *st; 
  long   tm;
  FILE   *fp=NULL;
  
  WorkMem(sdt);
  InitSet(sdt);
  InvMemo(sdt);
  
  sdt->mf=CfcAlloc(n,"mf, PspSolver"); 
  mf     =sdt->mf; 
  InvSymbo(sdt->st,mf); 
  
  bprt =0;
  bstp =0;
  bs   =0;
  lstp =3.0;
  rgap0=sdt->rgap/(1.0+sdt->rgap);
  
  y    =sdt->y;
  s    =sdt->s;
  sinv =sdt->sinv;
  sf   =sdt->sf;
  st   =sdt->st;
  cy   =sdt->cy;
  
  dy   =sdt->dy;
  dy1  =sdt->rw;
  dy2  =dy1+n;
  rbuf =dy2+n;
  ynew =rbuf;
  u    =ynew+2*n;
  uv   =u+n;
  
  dl   =0.0;
  dl1  =0.0;
  dl2  =0.0;
  lnew =0.0;
    
  CfcInit(sf,cy,sdt->rw); 
  for (i=0; i<n; i++)
    sf->diag[sf->invp[i]]=s[i];
  dCopy(sf->unnz,sf->uval,uv);
  
  if (sdt->ptyp==DvdCut)
    FixIndex(sdt);
    
  if (CfcOk!=ChlFact(sf,sdt->iw,rbuf,true))
    ExitProc(CholErr,"EquCutSol");
  
  if (par->prtlev) {
    sprintf(sname,"coplsdp.swp");
    fp=fopen(sname,"w");
    if (!fp) ExitProc(WriteFail,"sname");
  }
  
  for (sdt->iter=0; sdt->iter<par->maxiter; sdt->iter++) {
    /*
     * form S inverse
     */
    FindSinv(sf,st,sinv,rbuf,u,sdt->lamda,sdt->ptyp);
    
    /*
     * form M and factorize it
     */
    DotProd(sinv,mf);
    if (sdt->ptyp==BoxCut)
      for (i=0; i<n; i++)
        mf->diag[i]+=1.0/(y[i]*y[i]);
        
    if (CfcOk!=ChlFact(mf,sdt->iw,rbuf,true)) {
      for (i=0; i<n; i++)
        if (mf->diag[i]<1.0e-13)
          printf(" %d %e\n",i+1,mf->diag[i]);
          
      ExitProc(CholErr,"PspSolver");
    }
     
    /*
     * compute dy1 and dy2
     */
    FormDy(sdt,dy1,dy2,&ese,&dl1,&dl2,rbuf);
    
    /*
     * compute values
     */
    GetValues(sdt,dy,dy1,dy2,ese,&dl,dl1,dl2,&gapnew);
    
    /*
     * check if primal upper bound needs to be updated
     */ 
    bupd=(dl<sdt->s0)||
         ((sdt->ptyp!=EquCut)&&(sdt->ptyp!=DvdCut));
    if (gapnew<sdt->rgap&&bupd) {
      if (ChkPosDef1(sdt,y,sdt->lamda,dy,dl,rbuf)) {
        if (par->prtlev) {
          rewind(fp);
          for (i=0; i<n; i++)
            fprintf(fp,"%+24.16e %+24.16e\n",y[i],dy[i]);
          fprintf(fp,"%+24.16e %+24.16e\n",sdt->lamda,dl);
          fprintf(fp,"%+24.16e %24d\n",
                     sdt->rgap/sdt->rho,sdt->ptyp);
          fprintf(fp,"%24d %24d\n",sdt->is+1,sdt->it+1);
        }
        
        pnew     =ese;
        ModUppBnd(sdt,dy,dy1,dy2,&dl,dl1,dl2,&pnew,gapnew);
        sdt->rgap=gapnew;
        sdt->pval=pnew;
      }
    }
    
    bfac=GetStep(sdt,y,dy,sdt->lamda,dl,
                 par->alph,&lstp,&bstp,rbuf);
    
    UpdSol(sdt,y,s,dy,dl);
    
    rtmp=sdt->rgap/(1.0+fabs(sdt->dobj));
    beta=par->alph/sdt->pval;
    
    ShowMsg(sdt,sdt->iter,rtmp, sdt->step/beta);
    
	 
    if (rtmp<par->tolgap) {
      tm=GetTime(); 
	  ese=TimeInSec(otim[START],tm);	  
      printf("\n Exit -- 0: relative gap %3.1e < %3.1e\n",
             rtmp,par->tolgap);      
      fprintf(fout,"\n Exit -- 0: relative gap %3.1e < %3.1e\n",
                   rtmp,par->tolgap);
      break;
    }
    
    if (!bfac) {
      if (sdt->ptyp==DvdCut) {
        for (i=0; i<n; i++)
          sf->diag[sf->invp[i]]=cy->diag[i]-y[i];
        dCopy(sf->unnz,uv,sf->uval);
        
        sf->diag[sf->invp[sdt->is]]-=sdt->lamda;
        sf->diag[sf->invp[sdt->it]]-=sdt->lamda;
        sf->uval[sf->upst]         -=sdt->lamda;
        
        if (CfcOk!=ChlFact(sf,sdt->iw,rbuf,false))
          ExitProc(CholErr,"EquCutSol");
      }
      else {
        for (i=0; i<n; i++)
          sf->diag[sf->invp[i]]=cy->diag[i]-y[i];
        dCopy(sf->unnz,uv,sf->uval);
        if (CfcOk!=ChlFact(sf,sdt->iw,rbuf,false)) {
          for (i=0; i<n; i++)
            if (fabs(mf->diag[i])<1.0e-13)
              printf(" %d %e\n",i+1,mf->diag[i]);
          
          ExitProc(CholErr,"EquCutSol");
        }
      }
    }
    
    if (par->dynrho) {
      if (sdt->ptyp==BoxCut) {
        if (rtmp>0.01) {
          if (rtmp>0.1) {
            if (!bs&&rgap0/rtmp>=1.05)
              sdt->rho=5.0*n*sqrt(rgap0/rtmp);
            else {
              sdt->rho=par->rhofac*n*sqrt(rgap0/rtmp);
              bs=true;
            }
          }
          else
            sdt->rho=2.5*n*sqrt(rgap0/rtmp);
          
          sdt->rho=min(10.0*n,sdt->rho);
          rgap0   =rtmp;
        }
        else
          sdt->rho=5.0*n;
      }
    
      else {
        if (rtmp>0.01) {
          if (rtmp<0.1)
            par->rhofac=max(par->rhofac,1.3);
          sdt->rho=par->rhofac*n*sqrt(rgap0/rtmp);
          rgap0   =rtmp;
        }
      
        else
          sdt->rho=3.0*n;
      }
    }
  } /* main loop */ 
  
  if (par->prtlev) {
    fclose(fp);
    
    dFree(&sdat->s);
    dFree(&sdat->sinv);
    
    SmtFree(&sdat->c);
    
    fp=fopen(sname,"r");
    if (!fp) ExitProc(ReadFail,"sname");
    
    for (i=0; i<n; i++)
      fscanf(fp,"%lf%lf",&y[i],&dy[i]);
    
    fscanf(fp,"%lf%lf",&sdt->lamda,&sdt->dl);
    fscanf(fp,"%lf%d",&sdt->rho,&i);
    fclose(fp);
      
    if (sdt->ptyp==DvdCut) {
      for (i=0; i<n; i++)
        sf->diag[sf->invp[i]]=cy->diag[i]-y[i];
      dCopy(sf->unnz,uv,sf->uval);
      
      sf->diag[sf->invp[sdt->is]]-=sdt->lamda;
      sf->diag[sf->invp[sdt->it]]-=sdt->lamda;
      sf->uval[sf->upst]         -=sdt->lamda;
        
      if (CfcOk!=ChlFact(sf,sdt->iw,rbuf,false))
        ExitProc(CholErr,"EquCutSol");
    }
    else {
      for (i=0; i<n; i++)
        sf->diag[sf->invp[i]]=cy->diag[i]-y[i];  
      dCopy(sf->unnz,uv,sf->uval);
      if (CfcOk!=ChlFact(sf,sdt->iw,rbuf,false))
        ExitProc(CholErr,"EquCutSol");
    }
  }
} /* SdpProc */
