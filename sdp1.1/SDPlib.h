#include "SDPdef.h"

typedef struct {
    int idep;
    int last;
    int most;
    int cure;
    int loca;
    int lowp;
    int ntot;
    
    int *head;
    int *port;
    int *fwrd;
    int *bwrd;
  } xlist;
  
typedef struct {
    int nnod;
    int nn0;
    int raft;
    int head;
    int last;
    int ntot;
    
    int *adjn;
    int *rbeg;
    int *rexs;
    int *rlen;
    int *rend;
    int *pres;
    int *succ;
  } order;

int    InverAssign(sdpdat*);

order  *OdAlloc(int,int,char*);
void   OdFree(order**);
void   OdInit(order*,int*);
void   OdIndex(order*,int,int);
void   OdProc(order*,xlist*,int*,int*,int*,int*,int*,
              int*,int*,int*,int*,int*,int*,int*,int*);
int    GetOrder(order*,int*);

int    SmatxTrans(int,smatx*,int,int*,smatx**);
double GetFnorm(int,smatx*);
void   SmtCopy(int,smatx*,smatx*);
void   ChlSolve(chfac*,double*,double*);
void   FindInver(chfac*,syoff*,double*,sdpdat*);
void   DotProd(double*,chfac*);
void   Smt2Sym(smatx*,symat**);
void   MtxTimesVct(int,int,smatx*,double*,double*,int);
double SymSum(symat*);

void   CfcInit(chfac*,symat*,double*);
int    ChlFact(chfac*,int*,double*,int);

void   PermSmatx(smatx*,int*,int*);

xlist  *XtAlloc(int,int,char*);
void   XtFree(xlist**);
int    XtSucc(xlist*);
void   XtDel(xlist*,int);
void   XtPut(xlist*,int,int);
int    XtLeast(xlist*);
int    XtGet(xlist*,int*,int*);

void   IptAlloc(int,int,int**,char*);
void   IptFree(int,int**);
int    LocIntPos(int,int,int*);
void   PermTransSym(int,int*,int*,int*,int*,int,int*,int*,int*);
int    PspSymbo(symat*,chfac*);
int    InvSymbo(syoff*,chfac*); 

void   iSet(int,int,int*,int*);
void   iZero(int,int*,int*);
void   iFill(int,int,int*,int*,int*);
void   iSwap(int,int,int*);
void   iCopy(int,int*,int*);
int    iSum(int,int*);
void   dZero(int,double*,int*,int*);
void   dCopy(int,double*,double*);
void   dCat(int,int*,double*,double*);
double dDot(double*,double*,int);
double dSum(int,double*);
void   plusXs(int,int*,int*);
void   PlusByOne(int,int*,int*,int*);

