#include "SDPlib.h"

void PermSmatx(smatx *a,
               int   *p,
               int   *ip)
{
  int   h,i,j,k,nnz,
        n=a->nrow,
        n0=a->non0;
  smatx *r;
  array *aj,*rows;
  
  if (!p) return;
  
  if (!n0||!n) return;
    
  r  =SmtAlloc(n,n0,"PermSmatx");
    
  rows=r->rows;
  h   =0;
  for (i=0; i<n; i++) {
    j=p[i];
    aj=a->rows+j;
    nnz=aj->nn0;
    
    rows[i].nn0=nnz;
    rows[i].ja=rows->ja+h;
    rows[i].an=rows->an+h;
    
    for (k=0; k<nnz; k++) {
      rows->ja[h]=ip[aj->ja[k]];
      rows->an[h]=aj->an[k];
      h++;
    }
  }
  
  SmatxTrans(n,r,false,NULL,&a);
  
  for (i=0; i<n; i++) {
    p[i] =i;
    ip[i]=i;
  }
  
  SmtFree(&r);
} /* PermSmatx */
