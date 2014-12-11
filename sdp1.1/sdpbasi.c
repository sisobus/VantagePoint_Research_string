#include "SDPdef.h" 
 
clock_t GetTime(void)
{
  clock_t t;
#ifdef UNIXMACH
  struct tms {
    clock_t    tms_utime;
    clock_t    tms_stime;
    clock_t    tms_cutime;
    clock_t    tms_cstime;
  } tmrec;

  t = times(&tmrec);
  t = tmrec.tms_utime+tmrec.tms_stime;
#else
  t =clock();
#endif
  return t;
} /* GetTime */

int SetTime(clock_t tim[],
            int     phase)
{
  int i;
  
  tim[phase]=GetTime(); 
  for (i=phase+1; i<ELAPS; i++)
    tim[i]=tim[phase];
    
  return true;
}

int CheckArgs(int argc,char *argv[],char *ss)
{
  int i,id;
  
  id = 0;
  for (i=1; i<argc; i++) {
    if(strstr(argv[i],ss)) {
      argv[i]+=2;
      if (argv[i][0] == ':') {
        argv[i]++;
        id=i;
        break;
      }
      if (argv[i][0]) {
        id=i;
        break;
      }
      id = i+1;
      if (argv[id][0] == '-')
        id = 0;
      break;
    }
  }
  return id;
} /* Checkargs */

int LeftDots(char *ss,int len)
{
  char sz[60];
  int  i,j,k;
  
  k=strlen(ss);
  j=len-k;
  
  strcpy(sz,ss);
  for (i=0; i<j; i++)
    ss[i]='.';
  for (i=j; i<len; i++)
    ss[i]=sz[i-j];
  ss[len]='\0';
  return true;
}

double TimeInSec(clock_t head,
                 clock_t rear)
{
  double tscal;

#ifdef UNIXMACH
  tscal=0.01;
#else
  tscal=1.0/(double)CLOCKS_PER_SEC;
#endif

  return ((double)(rear-head)*tscal);
} /* TimeInSec */

int ExitProc(int  code,
             char *str)
{
  printf("\n Exit -- %d: ",code);
  
  switch (code) {
    case OptFound:
      printf("optimal solution found");
      return code;
      break;
    case OutOfSpc:
      printf("out of memory space");
      break;
    case IterLimit:
      printf("too many iterations");
      return code;
      break;
    case WriteFail:
      printf("can't create file");
      break;
    case ReadFail:
      printf("can't open file");
      break;
    default:
      break;
  }
  if (str) printf(", %s",str);
  ShutDown();
  exit(0);
} /* ExitProc */

void WorkMem(sdpdat *sdt)
{
  int n=sdt->ncol,isze,rsze;
  
  isze=2*n;
  rsze=5*n+sdt->sf->unnz;
  
  sdt->s =dAlloc(n,"s, WorkMem");
  sdt->y =dAlloc(n,"y, WorkMem");
  sdt->dy=dAlloc(n,"dy, WorkMem");
  sdt->iw=iAlloc(isze,"iw, WorkMem");
  sdt->rw=dAlloc(rsze,"rw, WorkMem");
} /* WorkMem */

