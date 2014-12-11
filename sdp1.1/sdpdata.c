#include "SDPdef.h"

void ChkBlocks(sdpdat *sdt,
               smatx  *a)
{
  int   g,i,j,k,nn0,beblk,nblk,
        *last,*bvtx,nrow=a->nrow;
  array *cj;
  
  bvtx=iAlloc(2*nrow,"bvtx, InverAssign");
  last=bvtx+nrow;
  
  nblk=0;
  for (j=0; j<nrow; j++) {
  
    last[j]=0;
    bvtx[j]=-1;
    
    cj     =a->rows+j;
    nn0    =cj->nn0;
    for (k=0; k<nn0; k++)
      if (cj->ja[k]>last[j])
        last[j]=cj->ja[k];
        
    beblk=false;
    
    if (last[j]<=j) {
      for (k=0; k<j; k++)
        if (last[k]>j)
          break;
      if (k==j) beblk=true; 
    }
    
    if (beblk) {
      bvtx[nblk]=j+1;
      nblk++;
    }
  }
 
  sdt->nblk  =nblk;
  sdt->maxblk=0;
  sdt->minblk=nrow;
  
  g   =0;
  for (i=0; i<nblk; i++) {
    
    k=bvtx[i];
    if (i) g=bvtx[i-1];
    
    k-=g+1;
    
    if (sdt->maxblk<k) sdt->maxblk=k;
    if (sdt->minblk>k) sdt->minblk=k;
    
  }
  iFree(&bvtx);
  
} /* ChkBlocks */

void DataInput(sdpdat *sdt)
{
  FILE   *fp;
  int    h,i,j,k,n,nnzo,nrow,*fir,
         beg,end,*sub,*afir,*asub,
         *coli,badd;
  double r,*val,*aval;
  smatx  *c;
  char   ss[60];
  
  sprintf(ss,"%s",sdt->ss);  
  fp=fopen(ss,"r");
  if (!fp) ExitProc(ReadFail,ss);  
  
  fscanf(fp,"%d%d",&n,&nnzo);

  if (sdt->ptyp==EquCut&&n%2) {
    printf("\n Invalid data - cannot equally divide an odd number.\n");
    ExitProc(DataErr,NULL);
  }
  else if (sdt->ptyp==UquCut&&(
          (n-(int)sdt->kap)<0||
          (n+(int)sdt->kap)%2)) {
    printf("\n Invalid data - bad cut definition.\n");
    ExitProc(DataErr,NULL);
  }
  else if (sdt->ptyp==DvdCut&&(
           max(sdt->is,sdt->it)>=n||
           min(sdt->is,sdt->it)<0||
           sdt->is==sdt->it)) {
    printf("\n Invalid data - bad indices s, t.\n");
    ExitProc(DataErr,NULL);
  }
 
  if (n<=0||nnzo<=0)
    ExitProc(MtxSzeErr,NULL);
  
  fir=iAlloc(n+1,"fir, DataInput");
  sub=iAlloc(nnzo,"fir, DataInput");
  val=dAlloc(nnzo,"fir, DataInput");
  
  for (i=0; i<=n; i++) 
    fir[i]=0; 
  
  for (k=0; k<nnzo; k++) {
    fscanf(fp,"%d%d%lf",&i,&sub[k],&val[k]); 
    fir[i]=k+1;
  }
  
  for (i=1; i<=n; i++)
    if (fir[i]<fir[i-1]) 
      fir[i]=fir[i-1]; 

  fclose(fp); 

  asub=sub;
  aval=val;
  afir=fir;
  coli=iAlloc(n,"coli, sdpdata");
  
  for (i=0; i<n; i++)
    coli[i]=-1;
  
  h      =0;
  beg    =0;
  for (i=0; i<n; i++) {
    end=fir[i+1];
    for (k=beg; k<end; k++) {
      j=sub[k];
      if (coli[j-1]==i) {
        printf(" duplicate %d %d\n",i+1,j);
        continue;
      }
      coli[j-1]=i;
      asub[h]=j-1;
      aval[h]=val[k];
      h++;
    }
    beg=end;
    afir[i+1]=h;
  }
  
  iFree(&coli);
  nnzo=h;
  
  if (sdt->ptyp!=DvdCut)
    sdt->c=SmtAlloc(n,2*nnzo+n,"sdt->c, MemAssign");
  else {
    if (sdt->is>=n) sdt->is=n-1;
    if (sdt->it>=n) sdt->it=n-1;
    if (sdt->is<0)  sdt->is=0;
    if (sdt->it<0)  sdt->it=0;
    if (sdt->is>sdt->it) {
      i      =sdt->is;
      sdt->is=sdt->it;
      sdt->it=i;
    }

    if (sdt->is==sdt->it) {
      sdt->is=n/4;
      sdt->it=(3*n)/4;
    }
       
    sdt->c=SmtAlloc(n,2*nnzo+n+2,"sdt->c, MemAssign");
  }
  
  c=sdt->c;
  
  for (i=0; i<n; i++)
    c->rows[i].nn0=0;
  
  for (i=0; i<n; i++) {
    beg=fir[i];
    end=fir[i+1];
    for (k=beg; k<end; k++) {
      j=sub[k];
      c->rows[i].nn0++;
      if (j>i)
        c->rows[j].nn0++;
    }
  }  
  
  badd=true;
  if (sdt->ptyp==DvdCut) {
    beg=fir[sdt->is];
    end=fir[sdt->is+1];
    for (k=beg; k<end; k++)
      if (sub[k]==sdt->it)
        badd=false;
    if (badd) {
      c->rows[sdt->is].nn0++;
      c->rows[sdt->it].nn0++;
    }
  }
  
  h=0;
  for (i=0; i<n; i++) {
    c->rows[i].ja =c->rows->ja+h;
    c->rows[i].an =c->rows->an+h;
    h            +=c->rows[i].nn0;
    c->rows[i].nn0=0;
  }
  
  nnzo=0;
  nrow=0;
  
  for (i=0; i<n; i++) {
    beg=fir[i];
    end=fir[i+1];
      
    for (k=beg; k<end; k++) {
      j=sub[k];
      r=val[k];
      h=c->rows[i].nn0;
      c->rows[i].ja[h]=j;
      c->rows[i].an[h]=r;
      c->rows[i].nn0++;
      nnzo++;
      
      if (j>i) {
        h=c->rows[j].nn0;
        c->rows[j].ja[h]=i;
        c->rows[j].an[h]=r;
        c->rows[j].nn0++;
        nnzo++;
      }
      else nrow++;
    }
  }
      
  if (sdt->ptyp==DvdCut&&badd) {
    j=sdt->is;
    i=sdt->it;
    h=c->rows[j].nn0;
    c->rows[j].ja[h]=i;
    c->rows[j].an[h]=0.0;
    c->rows[j].nn0++;
    h=c->rows[i].nn0;
    c->rows[i].ja[h]=j;
    c->rows[i].an[h]=0.0;
    c->rows[i].nn0++;
    nnzo+=2;
  }
  
  iFree(&fir);
  iFree(&sub);
  dFree(&val);
  
  sdt->nrow=n;
  sdt->ncol=n;
  sdt->nnzo=nnzo+(n-nrow);
  
  ChkBlocks(sdt,sdt->c);
  
} /* DataInput */

