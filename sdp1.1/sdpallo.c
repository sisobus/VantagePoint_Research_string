#include "SDPdef.h"

int *iAlloc(int  len,
            char *info)
{
  int *r=NULL;
  
  if (len) {
    r=(int*)calloc(len,sizeof(int));
    if (!r) ExitProc(OutOfSpc,info);
  }
  return r;
} /* iAlloc */

void iFree(int **x)
{
  int *r=*x;
  
  if (r) {
    free(r);
    *x=NULL;
  }
} /* iFree */

double *dAlloc(int  len,
               char *info)
{
  double *r=NULL;
  
  if (len) {
    r=(double*)calloc(len,sizeof(double));
    if (!r) ExitProc(OutOfSpc,info);
  }
  return r;
} /* dAlloc */

void dFree(double **x)
{
  double *r=*x;
  
  if (r) {
    free(r);
    *x=NULL;
  }
} /* dFree */

smatx *SmtAlloc(int  nrow,
                int  nnzo,
                char *info)
{
  smatx *r;
  
  r=(smatx*)calloc(1,sizeof(smatx));
  if (!r) ExitProc(OutOfSpc,info);
  
  if (nrow) {
    r->rows=(array*)calloc(nrow,sizeof(array));
    if (!r) ExitProc(OutOfSpc,info);
    
    if (nnzo) {
      r->rows->ja=iAlloc(nnzo,info);
      r->rows->an=dAlloc(nnzo,info);
    }
  }
  r->maxnrow=nrow;
  r->nrow=nrow;
  r->non0=nnzo;
  return r;
} /* SmatAlloc */

void SmtFree(smatx **a)
{
  smatx *r=*a;
  
  if (r) {
    if (r->rows) {
      iFree(&r->rows->ja);
      dFree(&r->rows->an);
      free(r->rows);
      r->rows=NULL;
    }
    
    r->maxnrow=0;
    r->nrow=0;
    free(r);
    *a=NULL;
  }
} /* SmatFree */

symat *SymAlloc(int  nrow,
                int  nnzo,
                char *info)
{
  symat *r;
  
  r=(symat*)calloc(1,sizeof(symat));
  if (!r) ExitProc(OutOfSpc,info);
  
  if (nrow) {
    r->diag=dAlloc(nrow,info);
    r->roff=(array*)calloc(nrow,sizeof(array));
    if (!r) ExitProc(OutOfSpc,info);
    
    if (nnzo) {
      r->roff->ja=iAlloc(nnzo,info);
      r->roff->an=dAlloc(nnzo,info);
    }
  }
  r->nrow=nrow;
  r->nnzo=nnzo;
  
  return r;
} /* SymAlloc */

syoff *SyoAlloc(int  nrow, 
                int  nnzo, 
                char *info) 
{ 
  syoff *r; 
   
  r=(syoff*)calloc(1,sizeof(syoff)); 
  if (!r) ExitProc(OutOfSpc,info); 
   
  if (nrow) { 
    r->roff=(array*)calloc(nrow,sizeof(array)); 
    if (!r) ExitProc(OutOfSpc,info); 
     
    if (nnzo) { 
      r->roff->ja=iAlloc(nnzo,info); 
      r->roff->an=dAlloc(nnzo,info); 
    } 
  } 
  r->nrow=nrow; 
  r->nnzo=nnzo; 
   
  return r; 
} /* SyoAlloc */ 
 
void SymFree(symat **a)
{
  symat *r=*a;
  
  if (r) {
    dFree(&r->diag);
    if (r->roff) {
      iFree(&r->roff->ja);
      dFree(&r->roff->an);
      free(r->roff);
      r->roff=NULL;
    }
    r->nrow=0;
    free(r);
    *a=NULL;
  }
} /* SymFree */

void SyoFree(syoff **a) 
{ 
  syoff *r=*a; 
   
  if (r) { 
    if (r->roff) { 
      iFree(&r->roff->ja); 
      dFree(&r->roff->an); 
      free(r->roff); 
      r->roff=NULL; 
    } 
    r->nrow=0; 
    free(r); 
    *a=NULL; 
  } 
} /* SyoFree */ 
 
int LvalAlloc(chfac *sf,
              char  *info)
{
  int nnz;
  
  nnz=iSum(sf->nrow,sf->ujsze);
  if ( nnz<=sf->unnz )
    return true;
  
  sf->unnz=0;
  if (sf->uval) dFree(&sf->uval);
  sf->uval=dAlloc(nnz,info);
  
  sf->unnz=nnz;
  return true;
} /* LvalAlloc */

chfac *CfcAlloc(int  maxrow,
                char *info)
{
  chfac *r=NULL;
  
  if (maxrow) {
    r=(chfac*)calloc(1,sizeof(chfac));
    if (!r) ExitProc(OutOfSpc,info);
    
    r->mrow =maxrow;
    r->nrow =maxrow;
      
    r->snnz =0;
    r->shead=iAlloc(maxrow,info);
    r->ssize=iAlloc(maxrow,info);
    r->ssub =NULL;
    r->diag =dAlloc(maxrow,info);
    r->unnz =0;
    r->ujnz =0;
    r->ujbeg=iAlloc(maxrow,info);
    r->uhead=iAlloc(maxrow,info);
    r->ujsze=iAlloc(maxrow,info);
    r->usub =NULL;
    r->uval =NULL;
    r->perm =iAlloc(maxrow,info);
    r->invp =iAlloc(maxrow,info);
    r->nsnds=0;
    r->subg =iAlloc(maxrow+1,info);      
  }
  return r;
} /* SchlAlloc */

void CfcFree(chfac **sf)
{
  chfac *r=*sf;
  
  if (*sf) {
    iFree(&r->shead);
    iFree(&r->ssize);
    iFree(&r->ssub);
    dFree(&r->diag);
    iFree(&r->ujbeg);
    iFree(&r->uhead);
    iFree(&r->ujsze);
    iFree(&r->usub);
    dFree(&r->uval);
    iFree(&r->perm);
    iFree(&r->invp);
    iFree(&r->subg);
    iFree(&r->dhead);
    iFree(&r->dbeg);
    iFree(&r->dsub);
    free(r);
  }
  *sf=NULL;
} /* CfcFree */

double **dPtAlloc(int  n,
                  char *info)
{
  int    i;
  double **r;
  
  r=NULL;
  if (!n) return r;
  
  r=(double **)calloc(n,sizeof(double*));
  if (!r) ExitProc(OutOfSpc,info);
  
  r[0]=dAlloc(n*n,info);
  for (i=1; i<n; i++)
    r[i]=r[i-1]+n;
  
  return r;
} /* dPtAlloc */

void dPtFree(double ***x)
{
  double **r=*x;
  
  if (r) {
    if (r[0])
      dFree(&r[0]);
    free(r);
    *x=NULL;
  }
} /* dPtFree */
