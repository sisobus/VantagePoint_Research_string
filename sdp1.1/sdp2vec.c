
#include "SDP2Vec.h"
 
 
double EvalObjective( symat  objMatrix,  
                      double *nvec,  
                      int    n,
                      int    ptype)
/* Evaluate the quadratic function xQx where Q  is a symmetric
   sparse matrix and x is a vector.  This is done by calculating
   the inner product of Q with the outer product of the vector  */                       
{  double answer=0; 
   int i,j,k,nnzr; 
  
   for (i=0;i<n;i++)
   {
      answer += objMatrix.diag[i] * nvec[i] * nvec[i];
      nnzr = objMatrix.roff[i].nn0;
      for (j=0;j<nnzr;j++)
      {
         k = objMatrix.roff[i].ja[j];
         answer += objMatrix.roff[i].an[j] * nvec[i] * nvec[k] * 2;
      }
   } 
   if (ptype == -3) return (-answer);  /* box QP */
   else return (answer/-4);
   /* Recall that original problem asked for a maximum */
} 
 
static int cut_comp(const void *e1,const void *e2)
/* Used in qsort routine to sort the structure 'orderCut' by its value */ 
{ 
   double d1, d2;
   orderCut *a1,*a2;
   a1=(orderCut*)e1;
   a2=(orderCut*)e2;
   d1 = a1->val; 
   d2 = a2->val; 
 
   if (d1<d2) 
      return (-1); 
   else if (d1>d2) 
      return (1);
   return(0); 
} 
 
void DefineCut(double  *thisCut, 
               int      type, 
               orderCut *vhat, 
               int      n,
               int      sindex,
               int      tindex) 
/* This uses the 'vhat' to define a specific answer to the problem */
{ int i; 
  if (type == -1 || type == -2) /* Max Cut or s-t cut */ 
     for (i=0;i<n;i++) 
       thisCut[i] = sign(vhat[i].val); 
  else if (type == -3) 
     for (i=0;i<n;i++) 
        thisCut[i] = fabs(thisCut[i]) * sign(vhat[i].val); 
  else if (type == 0) 
  {  for (i=0;i<n;i++) vhat[i].index=i; 
     qsort( (void *)vhat, n, sizeof(orderCut), cut_comp); 
     for ( i=0;i<n/2;i++) 
       thisCut[vhat[i].index] = 1; 
     for ( i=n/2; i<n; i++) 
       thisCut[vhat[i].index] = -1;    
  } 
  else 
  {  for (i=0;i<n;i++)
        vhat[i].index=i; 
     qsort( (void *)vhat, n, sizeof(orderCut), cut_comp); 
     if ( vhat[n/2].val > 0) 
     {  for ( i=0;i<(n-type)/2;i++) 
          thisCut[vhat[i].index] = -1; 
        for ( i=(n-type)/2; i<n; i++) 
          thisCut[vhat[i].index] = 1; 
     } 
     else 
     {  for ( i=0;i<(n+type)/2;i++) 
          thisCut[vhat[i].index] = -1; 
        for ( i=(n+type)/2; i<n; i++) 
          thisCut[vhat[i].index] = 1; 
     } 
  } 
  if (type==-2 && sign(vhat[sindex].val)==sign(vhat[tindex].val) )
     vhat[sindex].val *= -1;
} 
 
void GetRandVec(double *randVec, int n) 
/* This routine obtains a vector from the n dimensional
   unit sphere.  Actually, if the vector were to be normalized
   to length one, it would be a randomly chosen vector from the
   unit sphere.                                                 */
   
{  int i; 
   double scal,x;
   scal=pow(2.0,15)-1.0;
   for (i=0;i<n;i++)
    { x = (( rand())/scal - .5);
     randVec[i] = tan(3.1415926*x); 
     }
}
 
void SaveIt(double *cut, int type, int n, char *outfile2)
/* Format and save the best answer.                      */
{
   int i;
   FILE *fp;
   char ss[30];
   
   if (type>0) 
     sprintf(ss,"%s.uct",outfile2);
   else if (type==0) 
     sprintf(ss,"%s.eqt",outfile2);
   else if (type==-1) 
     sprintf(ss,"%s.max",outfile2);
   else if (type==-2) 
     sprintf(ss,"%s.st",outfile2);
   else if (type==-3) 
     sprintf(ss,"%s.box",outfile2);
   fp=fopen(ss,"w");
   if (!fp)
     ExitProc(WriteFail,ss);
     
   if (type>-3)
      for (i=0; i<n; i++)
      {  if (cut[i] > 0) fprintf(fp," 1 ");
         else fprintf(fp,"-1 ");
         if ( i%10 == 9) fprintf(fp,"\n");
      }
   else
      for (i=0; i<n; i++)
      {  fprintf(fp,"%5.2f ",cut[i]);
         if ( i%10 == 9) fprintf(fp,"\n");
      }
  fclose(fp);
}      

orderCut *OrdAlloc(int  len,
                   char *info)
{
  orderCut *r=NULL;
  
  if (len) {
    r=(orderCut*)calloc(len,sizeof(orderCut));
    if (!r) ExitProc(OutOfSpc,info);
  }
  return r;
} /* OrdAlloc */

void OrdFree(orderCut **x)
{
  orderCut *r=*x;
  
  if (r) {
    free(r);
    *x=NULL;
  }
} /* ordFree */


void GetSinve(sdpdat *pd,double *sinve)
/* Solve the system of equations S*sinve = e where e is vec of
   all ones.   Also save the solution to (C-D(dy))*spSinve = e */
{
  int    i,n;
  double eSpSinve,rtmp,*spSinve,*e,lambda;
  
  n=pd->nrow;
  spSinve=pd->rw;
  e=spSinve+n;
  lambda=pd->lamda;
  
  for (i=0; i<n; i++)
    e[i]=1.0;
    
  ChlSolve(pd->sf,e,spSinve);
  if (pd->ptyp>=0) {
    eSpSinve =0.0;
    for (i=0; i<n; i++)
      eSpSinve +=spSinve[i];
    
    rtmp=1.0-lambda*eSpSinve;
    rtmp=lambda*eSpSinve/rtmp;
     
    for (i=0; i<n; i++)
      sinve[i]=(1.0+rtmp)*spSinve[i];
  }
  else {
    for (i=0; i<n; i++)
      sinve[i]=spSinve[i];
  }
} /* GetSinve */

void GetXrow(sdpdat *pd,double *sinve,double *xrowi,double *ei)
/* This routine finds one row of the primal matrix X */
{  int    i,j,n=pd->nrow;
   double rtmp,lambda,eSpSinve=0.0,*spSinve,*dySinvei,
	      *sDys,eTxrowi=0.0; 
   
   spSinve  = pd->rw;     /* inv(C-D(dy)) * e   */
   dySinvei = spSinve+n;  
   sDys     = dySinvei+n; /* Sinv*D(dy)*Sinv*ei */
   lambda   =pd->lamda;
    
    ChlSolve(pd->sf, ei, xrowi);     /* Solve (C-D(dy))*xrowi=ei */
    if (pd->ptyp==EquCut||pd->ptyp==UquCut) { /* Rank one update */  
    
      rtmp=0.0;
      eSpSinve =0.0;
      for (i=0; i<n; i++) {
        rtmp+=xrowi[i];
        eSpSinve +=spSinve[i];
      }
      eSpSinve=1.0-lambda*eSpSinve;
    
      rtmp=lambda*rtmp/eSpSinve;
      for (i=0; i<n; i++)
        xrowi[i]=xrowi[i]+rtmp*spSinve[i];
    }   /* Now  S*xrowi = ei */
    
    for (j=0;j<pd->nrow; j++)
       eTxrowi += xrowi[j];
    
    for (j=0;j<pd->nrow; j++)
        dySinvei[j] = xrowi[j]*pd->dy[j];
        
    ChlSolve(pd->sf, dySinvei, sDys);  
    if (pd->ptyp==EquCut||pd->ptyp==UquCut) {
      rtmp=0.0;
      for (i=0; i<n; i++) 
	  {
        rtmp+=sDys[i];
      }
    
      rtmp=lambda*rtmp/eSpSinve;
      for (i=0; i<n; i++)
        sDys[i]=sDys[i]+rtmp*spSinve[i];
    }
    
    rtmp=pd->lamda * eTxrowi;
    for (j=0;j<pd->nrow; j++)
       xrowi[j] +=   sDys[j] +  rtmp*sinve[j]; 
     
}  /* GetXrow */


void RankReduce(sdpdat *pd)
/* This routine generates a specific cut or solution to the original, 
   unrelaxed problem. */ 
{
   int    i,j,dim=pd->nrow, type = pd->ptyp;
   double bestObj=0.0,thisObj=0.0,bestRObj=0.0,
	      bestHObj=0.0,*thisCut,*bestCut,*randVec,
		  *dblPtr,*sinve,*ei,hth=0.0,beta=0.0,*u,rtmp,ese; 
   orderCut *vhat;

   thisCut = dAlloc(dim,"Insufficient memory to store cut"); 
   bestCut = dAlloc(dim,"Insufficient memory to store cut"); 
   randVec = dAlloc(dim,"Insufficient memory to store cut"); 
   sinve   = dAlloc(dim,"Insufficient memory to store cut");
   vhat    = OrdAlloc(dim,"Insufficient memory to store cut");
   ei      = dAlloc(dim,"Insufficient memory to store cut");  
   
   u = pd->rw+dim;
   
   srand(clock());   /*  Used to generate random vectors */

   if (type == -3)   /* For Box QP, use sqrt of the diagonal */ 
    for (i=0;i<dim;i++) 
     thisCut[i]=bestCut[i]=sqrt(1-pd->rho*(pd->dy[i]-pd->y[i])/(pd->y[i]*pd->y[i]));  
 
   GetSinve(pd,sinve);

   if (pd->par.findvec <= 0)
   {
     for (i=0; i<dim; i++) 
     { 
        memset(ei,0,sizeof(double)*dim);
        ei[i]=1.0; 
               
        GetXrow(pd, sinve, randVec, ei);
        for (j=0; j<dim; j++)
          vhat[j].val=randVec[j];
        DefineCut(thisCut,pd->ptyp,vhat,pd->nrow,pd->is,pd->it); 
        thisObj = EvalObjective(*(pd->cy),thisCut,dim,type); 
      
        if (i==0 || thisObj > bestHObj)
        { 
          bestHObj=thisObj; 
          dblPtr=bestCut; bestCut=thisCut; thisCut=dblPtr; 
        }
     }
   }

   /* We need to factor ( S + A(dy) ) when we want to
   do ecut rank reduction */

   if (pd->par.findvec >=0)
   {
	 /* In case of ecut, we have to save L^-1 * e to sinve */
     if(pd->ptyp == EquCut)
	 {
		 
	   for(i=0;i<dim;i++) { randVec[i]=2.0; u[i]=0.0; }

       FSubst( pd->mf,randVec,u , pd->mf->NegDiag);

       for(i=0;i<dim;i++){ hth += u[i]*u[i]; }
	
	   beta = hth*( pd->dl - pd->lamda );
	   beta +=1.0;
      
       if( beta < 0 ) 
	   {   
		   printf("\n { 1 + h'h*(dl-l) } is less than zero. Error in rank reduction."); 
           ShutDown();
	   }  
       beta = sqrt(beta)-1.0;  
       beta /= hth;
	 }

     bestRObj = -dim*dim;

     for (i=0; i<dim; i++) 
	 {
        GetRandVec(randVec,dim);
                                            
        if(pd->ptyp == EquCut)
		{
          if(pd->mf->NegDiag < 0) Uhat1(pd,u,randVec,ei,beta);
		  else Uhat2(pd,u,randVec,ei,beta,pd->mf->NegDiag,hth);

		  /* ChlSolve(cf,b,x); */
		  ChlSolve( pd->sf,randVec,ei);

          rtmp = 0.0;
		  ese  = 0.0;
		  
		  for(j=0;j<dim;j++)
		  {
			  rtmp += ei[j];
			  ese  += sinve[j];
		  }
		  ese = 1.0- pd->lamda*ese;

		  rtmp = pd->lamda*rtmp/ese;
          for(j=0;j<dim;j++) randVec[j] = ei[j] + rtmp*sinve[j];

		}
		else
		GetVhat(randVec,sinve,pd);
	    
		for (j=0; j<dim; j++) vhat[j].val=randVec[j];
                       
	    DefineCut(thisCut,pd->ptyp,vhat,pd->nrow,pd->is,pd->it); 
	    thisObj = EvalObjective(*(pd->cy),thisCut,dim,type); 
      
	    if (i==dim || thisObj > bestRObj)
		{ 
	      bestRObj=thisObj; 
	      dblPtr=bestCut; bestCut=thisCut; thisCut=dblPtr; 
		}
      
	 }
   }
   else bestRObj = bestHObj;

   bestObj = max(bestHObj, bestRObj);
   if (pd->par.prtlev) SaveIt(bestCut,pd->ptyp,dim,pd->ss);
   
   dFree(&thisCut); 
   dFree(&bestCut); 
   dFree(&randVec); 
   OrdFree(&vhat); 
   dFree(&sinve);
   dFree(&ei);
   
   if (pd->par.findvec == 0)
    printf("\n The best solution found is              : %8.6e \n",bestObj);
   if (pd->par.findvec <= 0) 
    printf(" The best solution (heuristic) found is  : %8.6e \n",bestHObj);
   if (pd->par.findvec >= 0)
    printf(" The best solution (randomized) found is : %8.6e \n",bestRObj);
   printf("\n"); 
   fprintf(fout," The best solution found is : %8.6e \n\n",bestObj); 

} 


void FSubst(chfac *sf,double *b,double *x, int negDiag )
{
  int    i,t,sze,*ls,*ujsze,*usub,
         *ujbeg,*uhead,*invp;
  double *uval,*diag,*li,xi;

  ujsze=sf->ujsze;
  usub =sf->usub;
  ujbeg=sf->ujbeg;
  uhead=sf->uhead;
  diag =sf->diag;
  uval =sf->uval;
  invp =sf->invp;

   for(i=0; i<sf->nrow; i++)
   {
     ls    = usub+ujbeg[i];
     li    = uval+uhead[i];
     sze   = ujsze[i];  
     x[i]  = b[i]/diag[i];
     xi    = x[i];

     for(t=0; t<sze; ++t)
	 {
        b[ls[t]] -= li[t]*xi;
	 }
   }
     for (i=0; i<sf->nrow; i++)
     {

       if( i== negDiag )
       {
          x[i] *= sqrt(-1.0*diag[i]);
       }
       else if (diag[i]<1.0E-20)
       { 
         printf("\n diagonal diag[%d] is %f.\n",i,diag[i]);     
         return;   
       }
       else   
       { x[i] *= sqrt(diag[i]);
      
       }
     }
    
} /* FSubst */



void L_vect(chfac *sf,double *b,double *x,int negDiag )

{
  int    i,t,sze,*ls,*ujsze,*usub,
         *ujbeg,*uhead,*invp;
  double *uval,*diag,*li,xi;

  ujsze=sf->ujsze;
  usub =sf->usub;
  ujbeg=sf->ujbeg;
  uhead=sf->uhead;
  diag =sf->diag;
  uval =sf->uval;
  invp =sf->invp;
  
  for (i=0; i<sf->nrow; i++)
  {
	if(i==negDiag){ x[i]=b[i]/sqrt(-1.0*diag[i]); }
    else if (diag[i]<1.0E-20)
    { 
      printf("\n Diagonal diag[%d] is %f.\n",i,diag[i]);    
      return;      
    }
    else x[i]=b[i]/sqrt(diag[i]);
    b[i]=0.0;
  }

  for(i=0; i<sf->nrow; i++)
  {
    ls    = usub+ujbeg[i];
    li    = uval+uhead[i];
    sze   = ujsze[i];
    xi    = x[i];
    b[i] += xi*diag[i];

    for(t=0; t<sze; ++t)
      b[ls[t]] += li[t]*xi;
  }
  
  for (i=0; i<sf->nrow; i++)
    x[i]=b[invp[i]];

} /* L_Vect */


void  Uhat1( sdpdat *sdt, double *h, double *randv, double *ans,
             double beta )
{
  int i;
  double temp=0.0;

  for(i=0; i< sdt->nrow ; i++) temp += h[i]*randv[i]; 
  temp *= beta;

  for(i=0; i< sdt->nrow ; i++) 
  ans[i] = randv[i] + temp*h[i];

  /* Then ans = D(0)*(I+Bhh')*U */

  L_vect( sdt->mf,ans,randv,-1);
 
}

void Uhat2( sdpdat *sdt, double *h, double *randv, double *ans,
            double beta ,int negDiag, double hth)
{

  int i;
  double temp=0.0,temp2=0.0;
  
  temp2 = beta/(hth*beta +1.0);
 
  for( i=0; i < sdt->nrow; i++)
  { 
	if(i==negDiag)
    {        
       ans[i] = 1 - temp2*h[negDiag]*h[i]; 
       temp += ans[i]*ans[i];
	}
    else
    {
       ans[i] =  -1.0*temp2*h[negDiag]*h[i];
	   temp += ans[i]*ans[i];
	}
  }

  if( temp > 0.5 ) 
  {   
     printf("\nError in rank reduction. 1-2ff' is less than zero.");
     printf("\nwhich is %f.",1.0-2.0*temp );
     return;
  } 
     
  temp = (-1.0 + sqrt( 1.0-2.0*temp))/temp; 
  temp2 = 0.0;       
  
  for( i=0; i < sdt->nrow; i++)
  {
     temp2 += randv[i]*ans[i];
  }
  temp2 *= temp;

  for( i=0; i < sdt->nrow; i++)
  {
     randv[i] += ans[i]*temp2;
  }
     
  temp2=0.0;
  for( i=0; i < sdt->nrow; i++)
  {
     temp2 += h[i]*randv[i];
  }
  temp2 *= beta;

  for( i=0; i < sdt->nrow; i++)
  {
     ans[i] = randv[i] + h[i]*temp2;
  }

  L_vect( sdt->mf,ans,randv,negDiag);

	 /* Solve S(s)*x = ans */  
}


