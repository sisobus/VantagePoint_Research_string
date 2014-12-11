#include "SDPdef.h"

void PrintHead(void)
{
  char Title[8][60]={
    " ==========================================\n",
    " Computational Optimization Program Library\n",
    "                                           \n",
    "     (  C  O  P  L  )     Version 1.0      \n",
    "                                           \n",
    " Computational Optimization Lab.           \n",
    " The University of Iowa (1997)             \n",
    " ==========================================\n"};
  int i;

  printf("\n\n\n");
  for (i=0; i<8; i++)
    printf("%s", Title[i]);
  printf("\n *************************** COPL STARTS ");
  printf("***************************\n\n");

  return;
} /* PrintHead */

int FprintHead(void)
{
  fprintf(fout,"========================\n");
  fprintf(fout,"C O P L  1.0  (Oct 1997)\n");
  fprintf(fout,"========================\n\n");
  return true;
} /* FprintHead */

void PrintEnd(clock_t tim[])
{
  double tmproc;
  char   ss[LineSize];
    
  fprintf(fout,"\nTime distribution");
  fprintf(fout,"\n-----------------\n");
  
  tmproc=TimeInSec(tim[START],tim[DATAIN]);
  sprintf(ss," %.2f sec",tmproc);
  LeftDots(ss,39);
  printf(" time for initial data input %s\n",ss);
  fprintf(fout," time for initial data input %s\n",ss);

  tmproc=TimeInSec(tim[DATAIN],tim[OPTIM]);
  sprintf(ss," %.2f sec",tmproc);
  LeftDots(ss,32);
  printf(" time to solve semidefinite program %s\n",ss);
  fprintf(fout," time to solve semidefinite program %s\n",ss);

#ifdef TEST
  fprintf(fres,"%8.2f\\\\\n",tmproc);
#endif

  tmproc=TimeInSec(tim[OPTIM],tim[INTEG]);
  sprintf(ss," %.2f sec",tmproc);
  LeftDots(ss,31);
  printf(" time for vector solution generation %s\n",ss);
  fprintf(fout," time for vector solution generation %s\n",ss);

  tmproc=TimeInSec(tim[START],tim[ELAPS]);
  sprintf(ss," %.2f sec",tmproc);
  LeftDots(ss,40);
  printf(" time for the whole process %s\n",ss);
  fprintf(fout," time for the whole process %s\n",ss);

  printf("\n **************************** COPL ENDS ");
  printf("****************************\n\n");  

  fprintf(fout,"\n------------------------------  EOF  ");
  fprintf(fout,"-------------------------------\n\n");
  
  return;
}/* PrintEnd */

void ShowMsgTit(int   n,
               sdpdat *sdt,
	       double dens)
{
  char ss[80];
  
  sprintf(ss," %d",n);
  LeftDots(ss,57);
  printf(" dimension %s\n",ss);
  fprintf(fout," dimension %s\n",ss);
  
  sprintf(ss," %.2f%s",dens,"%");
  LeftDots(ss,59);
  printf(" density %s\n\n",ss);
  fprintf(fout," density %s\n\n",ss);

  fprintf(fout,"\nIteration Information\n");
  fprintf(fout,"---------------------\n");
  
  switch (sdt->ptyp) {
    case MaxCut:
    case BoxCut:
      printf(" %-4s %6s  %14s %14s %7s %7s  %7s\n",
             "ITER","|P(z)|","POBJ    ","DOBJ     ",
             "RGAP ","MRHO","MSTEP");
      fprintf(fout," %-4s %6s  %14s %14s %7s %7s  %7s\n",
                   "ITER","|P(z)|","POBJ    ","DOBJ     ",
                   "RGAP ","MRHO","MSTEP");
      break;
      
    default:
      printf(" %-4s %6s %10s %10s %7s %7s %10s  %5s\n",
             "ITER","|P(z)|","POBJ  ","DOBJ  ",
             "RGAP ","X0  ","LAMBDA","MSTEP");
      fprintf(fout," %-4s %6s %10s %10s %7s %7s %10s  %5s\n",
              "ITER","|P(z)|","POBJ  ","DOBJ  ",
              "RGAP ","X0  ","LAMBDA","MSTEP");
      break;
  }
  
  printf(" ---------------------------------");
  printf("----------------------------------\n");
  fprintf(fout," ---------------------------------");
  fprintf(fout,"----------------------------------\n");
} /* ShowMsgTit */

void ShowMsg(sdpdat *sdt,
             int    iter,
             double rgap,
             double lstp)
{
  double pobj,dobj;
  
  if (sdt->ptyp==BoxCut) {
    pobj=-sdt->pobj;
    dobj=-sdt->dobj;
  }
  else {
    pobj=-sdt->pobj/4.0;
    dobj=-sdt->dobj/4.0;
  }
  
  switch (sdt->ptyp) {
    case MaxCut:
    case BoxCut:
#ifdef PCMACHINE
      printf(" %-4d %6.2f  %12.7e  %12.7e"
             " %7.1e %7.1e %5.2f\n",
             iter+1,sdt->pval,pobj,dobj,rgap,
             sdt->rho/(double)sdt->ncol,lstp);
      fprintf(fout," %-4d %6.2f  %13.7e  %13.7e"
                   " %7.1e %7.1e %5.2f\n",
                   iter+1,sdt->pval,pobj,dobj,rgap,
                   sdt->rho/(double)sdt->ncol,lstp);
#else
      printf(" %-4d %6.2f   %13.7e  %13.7e "
             " %7.1e  %7.1e %6.2f\n",
             iter+1,sdt->pval,pobj,dobj,rgap,
             sdt->rho/(double)sdt->ncol,lstp);
      fprintf(fout," %-4d %6.2f   %13.7e  %13.7e "
                   " %7.1e  %7.1e %6.2f\n",
                   iter+1,sdt->pval,pobj,dobj,rgap,
                   sdt->rho/(double)sdt->ncol,lstp);
#endif
      break;
    
    default:
#ifdef PCMACHINE
      printf(" %-4d %6.2f %10.3f %10.3f"
             " %7.1e %7.1e %+9.3f %5.2f\n",
             iter+1,sdt->pval,pobj,dobj,
             rgap,sdt->x0,sdt->lamda,lstp);
      fprintf(fout," %-4d %6.2f %10.3f %10.3f"
                   " %7.1e %7.1e %+9.3f %5.2f\n",
                   iter+1,sdt->pval,pobj,dobj,
                   rgap,sdt->x0,sdt->lamda,lstp);
#else
      printf(" %-4d %6.2f %10.3f %10.3f"
             "  %7.1e  %7.1e %+9.3f %5.2f\n",
             iter+1,sdt->pval,pobj,dobj,
             rgap,sdt->x0,sdt->lamda,lstp);
      fprintf(fout," %-4d %6.2f %10.3f %10.3f"
                   "  %7.1e  %7.1e %+9.3f %5.2f\n",
                   iter+1,sdt->pval,pobj,dobj,
                   rgap,sdt->x0,sdt->lamda,lstp);
#endif
      break;
  }
} /* ShowMsg */
