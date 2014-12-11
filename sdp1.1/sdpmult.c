#include "SDPdef.h"

static void mulSnod(chfac  *sf,
                    int    snde,
                    int    f,
                    int    l,
                    double *b,
                    double *x)
{
  int    i,t,sze,*ls,*subg,*ujbeg,*uhead,*usub;
  double xi,*l1,*diag,*uval;

  subg =sf->subg;
  ujbeg=sf->ujbeg;
  uhead=sf->uhead;
  usub =sf->usub;
  diag =sf->diag;
  uval =sf->uval;
  
  f += subg[snde];
  l += subg[snde];

  for(i=f; i<l; ++i) {
    xi   =b[i];
    ls   =usub+ujbeg[i];
    l1   =uval+uhead[i];
    sze  =l-i-1;
    x[i]+=xi*diag[i];
    for(t=0; t<sze; ++t)
      x[ls[t]]+=l1[t]*xi;
  }
} /* mulSnod */

void GetUhat(chfac  *sf,
             double *b,
             double *x)
{
  int    i,k,n,s,t,sze,f,l,itemp,*ls,
         *subg,*ujsze,*usub,*ujbeg,*uhead;
  double rtemp1,rtemp2,rtemp3,rtemp4,
         rtemp5,rtemp6,rtemp7,rtemp8,
         *l1,*l3,*l2,*l4,*l5,*l6,*l7,*l8,
         *diag,*uval;
  
  n    =sf->nrow; 
  subg =sf->subg;
  ujsze=sf->ujsze;
  usub =sf->usub;
  ujbeg=sf->ujbeg;
  uhead=sf->uhead;
  diag =sf->diag;
  uval =sf->uval;
  
  for (i=0; i<n; i++) {
    x[i]=b[i]/sqrt(diag[i]);
    b[i]=0.0;
  }
  
  for (s=0; s<sf->nsnds; s++) {
    f=subg[s];
    l=subg[s+1];
    
    mulSnod(sf,s,0,l-f,x,b);
    
    itemp=l-f-1;  
    ls   =usub+ujbeg[f]+itemp;
    sze  =ujsze[f]-itemp;
    k    =f;
    
    itemp=l-1;
    for(; k+7<l; k+=8) {
      l1    =uval+uhead[k+0]+itemp-(k+0);
      l2    =uval+uhead[k+1]+itemp-(k+1);
      l3    =uval+uhead[k+2]+itemp-(k+2);
      l4    =uval+uhead[k+3]+itemp-(k+3);
      l5    =uval+uhead[k+4]+itemp-(k+4);
      l6    =uval+uhead[k+5]+itemp-(k+5);
      l7    =uval+uhead[k+6]+itemp-(k+6);
      l8    =uval+uhead[k+7]+itemp-(k+7);
        
      rtemp1=x[k+0];
      rtemp2=x[k+1];
      rtemp3=x[k+2];
      rtemp4=x[k+3];
      rtemp5=x[k+4];
      rtemp6=x[k+5];
      rtemp7=x[k+6];
      rtemp8=x[k+7];
       
      for(t=0; t<sze; ++t)
        b[ls[t]]+= rtemp1*l1[t]
                  +rtemp2*l2[t]
                  +rtemp3*l3[t]
                  +rtemp4*l4[t]
                  +rtemp5*l5[t]
                  +rtemp6*l6[t]
                  +rtemp7*l7[t]
                  +rtemp8*l8[t];
    }
    

    for(; k+3<l; k+=4) {
      l1    =uval+uhead[k+0]+itemp-(k+0);
      l2    =uval+uhead[k+1]+itemp-(k+1);
      l3    =uval+uhead[k+2]+itemp-(k+2);
      l4    =uval+uhead[k+3]+itemp-(k+3);
       
      rtemp1=x[k+0];
      rtemp2=x[k+1];
      rtemp3=x[k+2];
      rtemp4=x[k+3];

      for(t=0; t<sze; ++t)
        b[ls[t]]+= rtemp1*l1[t]
                  +rtemp2*l2[t]
                  +rtemp3*l3[t]
                  +rtemp4*l4[t];
    }

    for(; k+1<l; k+=2) {
      l1    =uval+uhead[k+0]+itemp-(k+0);
      l2    =uval+uhead[k+1]+itemp-(k+1);

      rtemp1=x[k+0];
      rtemp2=x[k+1];

      for(t=0; t<sze; ++t)
        b[ls[t]]+= rtemp1*l1[t]
                  +rtemp2*l2[t];
    }


    for(; k<l; ++k) {
      l1    =uval+uhead[k+0]+itemp-(k+0);
      
      rtemp1=x[k+0];

      for(t=0; t<sze; ++t)
        b[ls[t]]+=rtemp1*l1[t];
    }
  }
  
  for (i=0; i<n; i++)
    x[sf->invp[i]]=b[i];
} /* GetUhat */
