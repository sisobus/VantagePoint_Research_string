#include "SDPlib.h"

int SmatxTrans(int   n,
               smatx *a,
               int   pid,
               int   *p,
               smatx **at)
{
  int   i,j,k,t,nn0,nrow;
  smatx *r;
  array *ai,*aj;
  
  nn0 =a->non0;
  nrow=a->nrow;
  r   =*at;
  
  if (r) {
    if (nn0>r->non0||n>r->maxnrow)
      ExitProc(MtxSzeErr,"a, at in SmatxTrans");
      
    r->nrow=n;
  }

  else
    r=SmtAlloc(n,nn0,"r, SmatxTrans");

  for (i=0; i<n; i++)
    r->rows[i].nn0=0;
  
  for (i=0; i<nrow; i++) {
    ai=a->rows+i;
    for (k=0; k<ai->nn0; k++)
      r->rows[ai->ja[k]].nn0++;
  }
  
  nn0=0;
  for (j=0; j<n; j++) {
    r->rows[j].ja=r->rows->ja+nn0;
    r->rows[j].an=r->rows->an+nn0;
    nn0+=r->rows[j].nn0;
    r->rows[j].nn0=0;
  }
     
  for (k=0; k<a->nrow; k++) {
    i=k;
    if (p) {
      ai=a->rows+p[k];
      if (!pid) i=p[k];
    }
    else ai=a->rows+k;
    
    for (t=0; t<ai->nn0; t++) {
      j  =ai->ja[t];
      aj =r->rows+j;
      nn0=aj->nn0;
      
      aj->ja[nn0]=i;
      aj->an[nn0]=ai->an[t];
      aj->nn0++;
    }
  }
  *at=r;
    
  return true;
} /* SmatxTrans */

double GetFnorm(int   n,
                smatx *a)
{
  int    i,k,nn0;
  array  *ai;
  double r=0;
  
  for (i=0; i<n; i++) {
    ai =a->rows+i;
    nn0=ai->nn0;
    for (k=0; k<nn0; k++)
      r+=ai->an[k]*ai->an[k];
  }
  r=sqrt(r);
  return (r);  
} /* GetFnorm */

void SmtCopy(int   n,
             smatx *a,
             smatx *b)
{
  int   i,j,k,nn0;
  array *arow,*brow;
  
  if (b->nrow<a->nrow||b->non0<a->non0)
    ExitProc(MtxSzeErr,"a, b in SmtCopy");
  
  b->nrow=a->nrow;  
  arow   =a->rows;
  brow   =b->rows;
  
  nn0        =arow[0].nn0;
  brow[0].nn0=nn0;
  for (k=0; k<nn0; k++) {
    brow[0].ja[k]=arow[0].ja[k];
    brow[0].an[k]=arow[0].an[k];
  }
  
  for (i=1; i<n; i++) {
    j          =i-1;
    nn0        =arow[i].nn0;
    k          =brow[j].nn0;
    brow[i].nn0=nn0;
    brow[i].ja=brow[j].ja+k;
    brow[i].an=brow[j].an+k;
    for (k=0; k<nn0; k++) {
      brow[i].ja[k]=arow[i].ja[k];
      brow[i].an[k]=arow[i].an[k];
    }
  }
} /* SmtCopy */

static void SolFwdSnode(chfac  *sf,
                         int    snde,
                         int    f,
                         int    l,
                         double x[])
{
  int    i,t,sze,*ls,*subg=sf->subg,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead,
         *usub=sf->usub;
  double xi,*l1,*diag=sf->diag,*uval=sf->uval;

  f += subg[snde];
  l += subg[snde];

  for(i=f; i<l; ++i)
  {
    x[i] /= diag[i];
    xi    = x[i];

    ls    = usub+ujbeg[i];
    l1    = uval+uhead[i];
    sze   = l-i-1;

    for(t=0; t<sze; ++t)
      x[ls[t]] -= l1[t]*xi;
  }
} /* SolFwdSnode */

static void SolBward(int    nrow,
                     double diag[],
                     double uval[],
                     int    fir[],
                     double x[])
{
  int    i,t,sze;
  double x1,x2,rtemp,
         *x0,*l1,*l2;

  for(i=nrow; i;) {
    for(; i>1; --i) {
          -- i;
      l1   = uval+fir[i-1]+1;
      l2   = uval+fir[i  ]+0;
      sze  = nrow-i-1;
      x1   = 0.0;
      x2   = 0.0;
      x0   = x+1+i;

      for(t=0; t<sze; ++t)
      {
        rtemp = x0[t];

        x1   += l1[t]*rtemp;
        x2   += l2[t]*rtemp;
      }

      x[i]   -= x2/diag[i];
      x[i-1] -= (uval[fir[i-1]]*x[i]+x1)/diag[i-1];
    }

    for(; i;) {
          -- i;
      l1   = uval+fir[i];
      sze  = nrow-i-1;
      x1   = 0.0;
      x0   = x+1+i;

      for(t=0; t<sze; ++t)
        x1 += l1[t]*x0[t];

      x[i] -= x1/diag[i];
    }
  }
} /* SolBward */

void ForwSubst(chfac  *sf,
               double b[],
               double x[])
{
  int    i,k,s,t,sze,f,l,itemp,*ls,
         *subg=sf->subg,*ujsze=sf->ujsze,*usub=sf->usub,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead;
  double rtemp1,rtemp2,rtemp3,rtemp4,
         rtemp5,rtemp6,rtemp7,rtemp8,
         *l1,*l3,*l2,*l4,*l5,*l6,*l7,*l8,
         *diag=sf->diag,*uval=sf->uval;
   
  for(i=0; i<sf->nrow; ++i)
    x[i] = b[sf->perm[i]];

  for(s=0; s<sf->nsnds; ++s) {
    f = subg[s];
    l = subg[s+1];

    SolFwdSnode(sf,s,0,l-f,x);

    itemp = l-f-1;
    ls    = usub+ujbeg[f]+itemp;
    sze   = ujsze[f]-itemp;
    k     = f;

    itemp = l-1;
    for(; k+7<l; k+=8) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);
      l5       = uval+uhead[k+4]+itemp-(k+4);
      l6       = uval+uhead[k+5]+itemp-(k+5);
      l7       = uval+uhead[k+6]+itemp-(k+6);
      l8       = uval+uhead[k+7]+itemp-(k+7);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];
      rtemp5   = x[k+4];
      rtemp6   = x[k+5];
      rtemp7   = x[k+6];
      rtemp8   = x[k+7];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t]
                    + rtemp5*l5[t]
                    + rtemp6*l6[t]
                    + rtemp7*l7[t]
                    + rtemp8*l8[t];
    }

    for(; k+3<l; k+=4) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t];
    }

    for(; k+1<l; k+=2) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t];
    }


    for(; k<l; ++k) {
      l1       = uval+uhead[k+0]+itemp-(k+0);

      rtemp1   = x[k+0];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t];
    }
  }
  
} /* ForwSubst */

void ChlSolve(chfac  *sf,
              double b[],
              double x[])
{
  int    i,k,s,t,sze,f,l,itemp,*ls,
         *subg=sf->subg,*ujsze=sf->ujsze,*usub=sf->usub,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead;
  double x1,x2,
         rtemp1,rtemp2,rtemp3,rtemp4,
         rtemp5,rtemp6,rtemp7,rtemp8,
         *l1,*l3,*l2,*l4,*l5,*l6,*l7,*l8,
         *diag=sf->diag,*uval=sf->uval;
   
  for(i=0; i<sf->nrow; ++i)
    x[i] = b[sf->perm[i]];

  for(s=0; s<sf->nsnds; ++s) {
    f = subg[s];
    l = subg[s+1];

    SolFwdSnode(sf,s,0,l-f,x);

    itemp = l-f-1;
    ls    = usub+ujbeg[f]+itemp;
    sze   = ujsze[f]-itemp;
    k     = f;

    itemp = l-1;
    for(; k+7<l; k+=8) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);
      l5       = uval+uhead[k+4]+itemp-(k+4);
      l6       = uval+uhead[k+5]+itemp-(k+5);
      l7       = uval+uhead[k+6]+itemp-(k+6);
      l8       = uval+uhead[k+7]+itemp-(k+7);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];
      rtemp5   = x[k+4];
      rtemp6   = x[k+5];
      rtemp7   = x[k+6];
      rtemp8   = x[k+7];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t]
                    + rtemp5*l5[t]
                    + rtemp6*l6[t]
                    + rtemp7*l7[t]
                    + rtemp8*l8[t];
    }

    for(; k+3<l; k+=4) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t];
    }

    for(; k+1<l; k+=2) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t];
    }


    for(; k<l; ++k) {
      l1       = uval+uhead[k+0]+itemp-(k+0);

      rtemp1   = x[k+0];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t];
    }
  }

  if (sf->nsnds) {
    s = sf->nsnds - 1;
    f = subg[s];
    l = subg[s+1];

    dCopy(l-f,x+f,b+f);

    SolBward(l-f,diag+f,uval,uhead+f,b+f);

    s = sf->nsnds-1;

    for(; s>=1; --s) {
      f = subg[s-1];
      l = subg[s];
      i = l;

      for(; i>1+f; --i) {
            -- i;
        ls   = usub+ujbeg[i];
        l1   = uval+uhead[i-1]+1;
        l2   = uval+uhead[i  ]+0;
        sze  = ujsze[i];
        x1   = 0.0;
        x2   = 0.0;

        for(t=0; t<sze; ++t) {
          rtemp1 = b[ls[t]];

          x1    += l1[t]*rtemp1;
          x2    += l2[t]*rtemp1;
        }

        b[i]   = x[i  ] -  x2  / diag[i];
        b[i-1] = x[i-1] - (x1 + uval[uhead[i-1]]*b[i]) / diag[i-1];
      }

      for(; i>f;) {
            -- i;
        l1   = uval+uhead[i];
        ls   = usub+ujbeg[i];
        sze  = ujsze[i];
        x1   = 0.0;

        for(t=0; t<sze; ++t)
          x1+= l1[t]*b[ls[t]];

        b[i] = x[i] - x1/diag[i];
      }
    }
  }

  for(i=0; i<sf->nrow; ++i)
    x[i] = b[sf->invp[i]];

} /* ChlSolve */

void FindInver(chfac  *cf,
               syoff  *xinv,
               double *rw,
               sdpdat *sdt)
{
  int    i,j,k,nn0,n=cf->nrow;
  double *b,*x;
  array  *rows;
  
  rows=xinv->roff;
  
  b   =rw;
  x   =rw+n; 
  
  for (j=0; j<n; j++) {
    
    for (i=0; i<n; i++) {
      x[i]=0.0;
      b[i]=0.0;
    }
    b[j]=1.0;
    
    ChlSolve(cf,b,x);
    
    if (j>0) {
      i          =j-1;
      rows[j].ja=rows[i].ja+rows[i].nn0;
      rows[j].an=rows[i].an+rows[i].nn0;
    }
    
    sdt->sinv[j]=x[j];
    nn0=rows[j].nn0;
    for (k=0; k<nn0; k++) {
      i=rows[j].ja[k];
      rows[j].an[k]=x[i];
    } 
  }
} /* FindInver */
 
void DotProd(double *sinv, 
             chfac  *mf) 
{ 
  int i,n=mf->nrow; 
 
  for (i=0; i<n; i++) 
	mf->diag[i]=sinv[i]*sinv[i]; 
 
  for (i=0; i<mf->unnz; i++) 
    mf->uval[i]*=mf->uval[i]; 
} /* DotProd */ 

void Smt2Sym(smatx *a,
             symat **ay)
{
  int   h,i,j,k,n,nn0,nnzo;
  array *ai,*aj,*ayj;
  symat *r=*ay;
  
  n  =a->nrow;
  if (!r) {
    nn0=0;
    for (i=0; i<n; i++) {
      ai=a->rows+i;
      for (k=0; k<ai->nn0; k++) {
        if (ai->ja[k]>i)
          nn0++;
      }
    }
    
    r=SymAlloc(n,nn0,"Smt2Sym");
  }
  
  if (r->nrow>n||r->nnzo>a->non0)
    ExitProc(MtxSzeErr,"Smt2Sym");
  
  nnzo=0;
  for (j=0; j<n; j++) {
  
    r->diag[j]=0.0;
    aj =a->rows+j;
    ayj=r->roff+j;
    
    ayj->ja=r->roff->ja+nnzo;
    ayj->an=r->roff->an+nnzo;
    
    h  =0;
    nn0=aj->nn0;
    for (k=0; k<nn0; k++) {
      i=aj->ja[k];
      if (i==j)
        r->diag[j]=aj->an[k];
      else if (i>j) {
        ayj->ja[h]=aj->ja[k];
        ayj->an[h]=aj->an[k];
        h++;
      }
    }
    
    ayj->nn0 =h;
    nnzo    +=h;
  }
  
  *ay=r;
} /* Smt2Sym */

void MtxTimesVct(int    nrow,
                 int    ncol,
                 smatx  *a,
                 double *x,
                 double *y,
                 int    btran)
{
  int   i,j,k,nn0;
  array *aj;
  
  if (nrow<0 || ncol<0)
    ExitProc(MtxSzeErr,"MtxTimesVct");
    
  if (btran) { /* x=x+A^T*y */    
    for (j=0; j<ncol; j++) {
      aj =a->rows+j;
      nn0=aj->nn0;
      for (k=0; k<nn0; k++) {
        i    =aj->ja[k];
        x[j]+=aj->an[k]*y[i];
      }
    }
  }
  
  else { /* y=y+A*x */
    for (j=0; j<ncol; j++) {
      aj =a->rows+j;
      nn0=aj->nn0;
      for (k=0; k<nn0; k++) {
        i    =aj->ja[k];
        y[i]+=aj->an[k]*x[j];
      }
    }  
  }
} /* MtxTimesVct */

double SymSum(symat  *a)
{
  int    i,k,n=a->nrow;
  double r=0.0,sum=0.0;
  array  *ai;
  
  for (i=0; i<n; i++) {
    sum+=a->diag[i];
    ai  =a->roff+i;
    for (k=0; k<ai->nn0; k++)
      r+=ai->an[k];
  }
  sum+=r+r;
  
  return sum;
} /* SymSum */

