#include "SDPdef.h"

sdpdat  *sdat;
clock_t *stim;
FILE    *fout,*fres;

void InitParams(params *par)
{
  par->maxiter  =200;
  par->prtlev   =true;
  par->dynrho   =true;
  par->varstep  =true;
  par->tolgap   =1.0e-04;
  par->alph     =0.95;
  par->tolpiv   =1.0e-13;
  par->bmfac    =1.0;
  par->findvec  =0;
  par->cachesize=256;
  par->cacheunit=1000;
} /* SetParams */

int LocParam(char *ss)
{
  int i,k,id;
  
  k =strlen(ss);
  id=k;
  for (i=0; i<k; i++) {
    if (ss[i]==':') id=i;
    ss[i]=(char)tolower(ss[i]);
  }
  
  return id;
} /* LocParam */

void SetParams(sdpdat *sdt)
{
  int  i,k,h,ival;
  FILE *fp;
  char ss[80],sold[80]="",
       paramtitle[][30]={
         "problem_type",       "data_type",
         "maximum_iterations", "save",
         "dynamic_rho",        "variable_stepsize",
         "stopping_tolerance", "pivoting_tolerance",
         "trust_region_radius","cache_memory_size",
         "cache_memory_unit",  "cut_difference",
         "stcut_index",        "penalty_factor",
         "rho_factor",         "findvec",
         ""
       };
  double rval;
  params *par=&sdt->par;
  
  InitParams(par);
  
  if (*(sdt->spc)) {
    fp=fopen(sdt->spc,"r");
    if (!fp) ExitProc(ReadFail,sdt->spc);
  }
  
  else {
    fp=fopen("COPLSDP.spc","r");
    if (!fp) {
      if (sdt->ptyp==UquCut) par->rhofac=1.2;
      else if (sdt->ptyp==BoxCut) par->rhofac=2.1;
      else par->rhofac=1.5;
      return;
    }
  }
  
  fprintf(fout,"SPC file\n");
  fprintf(fout,"--------\n");
  fprintf(fout," BEGIN\n");
  
  while (!feof(fp)) {
    fgets(ss,80,fp);
    if (strcmp(ss,sold)) {
      fprintf(fout,"   ");
      fputs(ss,fout);
    }
    for (i=0; i<80; i++)
      sold[i]=ss[i];
      
    if (*ss=='*' ||*ss=='\n'||
        *ss=='\0'||*ss==' ') 
      continue;
      
    k=LocParam(ss)+1;
    for (i=0; i<17; i++)
      if (strstr(ss,paramtitle[i]))
        break;
        
    switch (i) {
      case 0:
        if (strstr(ss+k,"mcut"))
          sdt->ptyp=MaxCut;
        else if (strstr(ss+k,"ecut"))
          sdt->ptyp=EquCut;
        else if (strstr(ss+k,"ucut"))
          sdt->ptyp=UquCut;
        else if (strstr(ss+k,"dcut"))
          sdt->ptyp=DvdCut;
        else if (strstr(ss+k,"box"))
          sdt->ptyp=BoxCut;
        break;
      case 1:
        if (strstr(ss+k,"graph"))
          sdt->dtyp=Gtype;
        else if (strstr(ss+k,"matrix"))
          sdt->dtyp=Mtype;
        break;
      case 2:
        ival=atoi(ss+k);
        if (ival>0) par->maxiter=ival;
        break;
      case 3:
        if (strstr(ss+k,"yes"))
          par->prtlev=true;
        else if (strstr(ss+k,"no"))
          par->prtlev=false;
        break;
      case 4:
        if (strstr(ss+k,"on"))
          par->dynrho=true;
        else if (strstr(ss+k,"off"))
          par->dynrho=false;
        break;
      case 5:
        if (strstr(ss+k,"on"))
          par->varstep=true;
        else if (strstr(ss+k,"off"))
          par->varstep=false;
        break;
      case 6:
        rval=atof(ss+k);
        if (rval>0.0) par->tolgap=rval;
        break;
      case 7:
        rval=atof(ss+k);
        if (rval>0.0) par->tolpiv=rval;
        break;
      case 8:
        rval=atof(ss+k);
        if (rval>0.0) par->alph=rval;
        break;
      case 9:
        ival=atoi(ss+k);
        par->cachesize=ival;
        break;
      case 10:
        ival=atoi(ss+k);
        if (ival>0) par->cacheunit=ival;
        break;
      case 11:
        ival=atoi(ss+k);
        if (ival>0) sdt->kap=ival;
        break;
      case 12:
        ival=atoi(ss+k);
        if (ival>=0) sdt->is=ival-1;
        h=strlen(ss);
        
        for (i=k; i<h; i++)
          if (isdigit(ss[i]))
            break;
        for (k=i; k<h; k++)
          if (!isdigit(ss[k]))
            break;
        ival=atoi(ss+k)-1;
        sdt->it=max(ival,sdt->is);
        sdt->is=min(ival,sdt->is);
        break;
      case 13:
        rval=atof(ss+k);
        par->bmfac=max(rval,0.01);
        break;
      case 14:
        rval=atof(ss+k);
        par->rhofac=max(rval,1.0);
        break;
      case 15:
        if (strstr(ss+k,"quick"))
          par->findvec=-1;
        else if (strstr(ss+k,"both"))
          par->findvec=0;
        else if (strstr(ss+k,"random"))
          par->findvec=1;
        break;
      default:
        break;
    }
  }
  fprintf(fout," END\n");
  
  fclose(fp);
  
  if (!par->rhofac) {
    if (sdt->ptyp==UquCut) par->rhofac=1.2;
    else if (sdt->ptyp==BoxCut) par->rhofac=2.1;
    else par->rhofac=1.5;
  }
} /* SetParams */

int SdpSolver(sdpdat  *sdt,
              clock_t otim[])
{
  int     i,k,n,nnz;
  array   *ai;
  char    ss[80];
  
  SetParams(sdt);
  
  fprintf(fout,"\nProblem Statics\n");
  fprintf(fout,"---------------\n");
  
  sprintf(ss," %s",sdt->ss);
  LeftDots(ss,54);
  printf(" problem name %s\n",ss);
  fprintf(fout," problem name %s\n",ss);
  
  switch(sdt->ptyp) {
    case MaxCut:
      sprintf(ss," maximal cut");
      break;
    case EquCut:
      sprintf(ss," equal cut");
      break;
    case UquCut:
      sprintf(ss," unequal cut");
      break;
    case DvdCut:
      sprintf(ss," s-t cut");
      break;
    case BoxCut:
      sprintf(ss,"  box QP");
      break;
    default:
      ExitProc(SysError,"problem type error");
      break;
  }
  
  LeftDots(ss,54);
  printf(" problem type %s\n",ss);
  fprintf(fout," problem type %s\n",ss);

  DataInput(sdt);

#ifdef TEST
  fprintf(fres,"%4d&%5.2f%s&",sdt->ncol,
	      100.0*(((double)sdt->nnzo/
		  (double)sdt->ncol)/(double)sdt->ncol),
		  "\\%");
#endif

  if (sdt->ptyp==DvdCut) {
#ifdef TEST
    sdt->is=sdt->ncol/3+1;
	sdt->it=2*sdt->ncol/3+1;
#endif

    sprintf(ss," %d, %d",sdt->is+1,sdt->it+1);
    LeftDots(ss,55);
    printf(" s-t indices %s\n",ss);
    fprintf(fout," s-t indices %s\n",ss);
  }
  
  if (sdt->ptyp==UquCut) {
    if (sdt->kap<1.0)       sdt->kap=1;
    if (sdt->kap>sdt->ncol) sdt->kap=sdt->ncol;

#ifdef TEST
    sdt->kap=sdt->ncol/8;
#endif

    sprintf(ss," %d",(int)sdt->kap);
    LeftDots(ss,52);
    printf(" cut difference %s\n",ss);
    fprintf(fout," cut difference %s\n",ss);
  }
  
  SetTime(otim,DATAIN);
   
  Smt2Sym(sdt->c,&sdt->cy);
  n      =sdt->cy->nrow;
  sdt->sf=CfcAlloc(n,"sdt->sf, PspProc");
  PspSymbo(sdt->cy,sdt->sf);
  
  if (sdt->dtyp==Gtype) {
    for (i=0; i<n; i++) {
      sdt->cy->diag[i]=0.0;
      ai=sdt->c->rows+i;
      nnz=ai->nn0;
      for (k=0; k<nnz; k++) {
        if (ai->ja[k]==i) continue;
        sdt->cy->diag[i]-=ai->an[k];
      }
    }
  }
  
  SmtFree(&sdt->c);
  
  SdpProc(sdt,otim);
  
  return true;
} /* SdpProc */

int main(int argc,char *argv[])
{  
  clock_t otim[5];
  sdpdat  pd={0};
  int     i,beg,end;
  char    *cs,ss[60];

#ifdef TEST
  char   sn[][10]=
         {
	       "G1", "G2", "G3", "G4", "G5",
	       "G6", "G7", "G8", "G9", "G10",
	       "G11","G12","G13","G14","G15",
	       "G16","G17","G18","G19","G20",
	       "G21","G22","G23","G24","G25",
	       "G26","G27","G28","G29","G30",
	       "G31","G32","G33","G34","G35",
	       "G36","G37","G38","G39","G40",
	       "G41","G42","G43","G44","G45",
	       "G46","G47","G48","G49","G50",
		   "G51","G52","G53","G54","G55",
		   "G56","G57","G58","G59","G60",
		   "G61","G62","G63","G64"
         };
  int    iter;
#endif

  sdat=&pd;
  stim=otim;
  PrintHead();
  pd.ptyp=MaxCut;
  
  beg=CheckArgs(argc,argv,"-p");
  if (beg) {
    cs=argv[beg];
    end=strlen(cs);
    for (i=0; i<end; i++)
      cs[i]=(char)tolower(cs[i]);
      
    if (!strncmp(cs,"ecut",4))      pd.ptyp=EquCut;
    else if (!strncmp(cs,"ucut",4)) pd.ptyp=UquCut;
    else if (!strncmp(cs,"dcut",4)) pd.ptyp=DvdCut;
    else if (!strncmp(cs,"box",4))  pd.ptyp=BoxCut;
    else                            pd.ptyp=MaxCut;
  }
#ifdef TEST
  for (iter=54; iter<56; iter++)
  {
    sprintf(pd.ss,"%s",sn[iter]);
	fres=fopen("sdptest.res","a+");
	if (!fres)
	  ExitProc(WriteFail,"sdptest.res");
	fprintf(fres,"%-5s&",sn[iter]);

	pd.dobj=0.0;
	pd.pobj=0.0;
	pd.rgap=0.0;
    
	pd.nrow=0;
	pd.ncol=0;
	pd.nnzo=0;
#else
  beg=CheckArgs(argc,argv,"-f");
  if (beg) {
    cs=argv[beg];
    sprintf(pd.ss,"%s",cs);
  }

  if (!(*pd.ss)) {
    printf(" Please type the name of your data file: ");
    fflush(stdin);
    gets(pd.ss);
  }
#endif
  fout=fopen(pd.ss,"r");
  if (!fout)
    ExitProc(ReadFail,pd.ss);
    
  fclose(fout);
  
  beg=CheckArgs(argc,argv,"-s");
  if (beg) {
    cs=argv[beg];
    sprintf(pd.spc,"%s",cs);
  }

  sprintf(ss,"%s.out",pd.ss);
  fout=fopen(ss,"w");
  if (!fout) ExitProc(WriteFail,ss);
  
  FprintHead();
  SetTime(otim,START);
  SdpSolver(&pd,otim);
  SetTime(otim,OPTIM);
  
  if (pd.ptyp==BoxCut) {
    printf("\n Upper bound for solution              : %8.6e \n",-pd.dobj);
    fprintf(fout,"\n Upper bound for solution %8.6e \n",-pd.dobj);
  }
  else {
    printf("\n Upper bound for solution              : %8.6e \n",pd.dobj/-4.0);
    fprintf(fout,"\n Upper bound for solution %8.6e \n",pd.dobj/-4.0);
  }
  
  if (pd.par.findvec>=0&&pd.par.prtlev)
    GetLhat(&pd);

  if(pd.ptyp>0) pd.ptyp=(int)sqrt(pd.kap);
  if (pd.par.prtlev) RankReduce(&pd);

  SetTime(otim,INTEG);

#ifdef TEST
  if (pd.ptyp==BoxCut)
  fprintf(fres,"%12e&%12e&%8.1e&%3d&",
	      -pd.pobj,-pd.dobj,
		  pd.rgap/(1.0+fabs(pd.dobj)),
		  pd.iter);
  else
  fprintf(fres,"%12e&%12e&%8.1e&%3d&",
	      -pd.pobj/4.0,-pd.dobj/4.0,
		  pd.rgap/(1.0+fabs(pd.dobj)),
		  pd.iter);
#endif

  ShutDown();
  SetTime(otim,ELAPS);
  PrintEnd(otim);
  
  fclose(fout);

#ifdef TEST
  fclose(fres);
  }
#endif

  return CfcOk;
} /* main */
