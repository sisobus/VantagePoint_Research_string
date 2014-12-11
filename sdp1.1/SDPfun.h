/*
 * common functions
 */
#if !defined (min)
#define min(a,b) ((a <= b)? (a) : (b))
#endif
#if !defined (max)
#define max(a,b) ((a >= b)? (a) : (b))
#endif
#if !defined (sign)
#define sign(a) ((a<0)? (-1) : (1))
#endif

/*
 * functions in sdpallo.c
 */
int     *iAlloc(int,char*);
void    iFree(int**);
double  *dAlloc(int,char*);
void    dFree(double**);
smatx   *SmtAlloc(int,int,char*);
void    SmtFree(smatx**);
symat   *SymAlloc(int,int,char*);
void    SymFree(symat**);
syoff   *SyoAlloc(int,int,char*); 
void    SyoFree(syoff**); 
int     LvalAlloc(chfac*,char*);
chfac   *CfcAlloc(int,char*);
void    CfcFree(chfac**);
double  **dPtAlloc(int,char*);
void    dPtFree(double***);

/*
 * functions in sdpbas.c
 */
clock_t GetTime(void);
int     SetTime(clock_t*,int);
double  TimeInSec(clock_t,clock_t);
int     CheckArgs(int,char**,char*);
int     LeftDots(char*,int);
int     ExitProc(int,char*);
void    WorkMem(sdpdat*);

/*
 * functions in sdpchec.c
 */
double  HouseTrans(int,double*,double*);
int     ChkPosEig(chfac*,double,double*);
int     ChkPosDef1(sdpdat*,double*,double,double*,double,double*);
int     ChkPosDef2(sdpdat*,double*,double,double*);

/*
 * functions in sdpdata.c
 */
void    ChkBlocks(sdpdat*,smatx*);
void    DataInput(sdpdat*);

/*
 * functions in sdpdire.c
 */
void    FormDy(sdpdat*,double*,double*,double*,double*,
               double*,double*);

/* 
 * functions in sdpinit.c 
 */ 
void    InitSet(sdpdat*);

/* 
 * functions in sdpmain.c 
 */ 
void    SetParams(sdpdat*);
int     SdpSolver(sdpdat*,clock_t*);
int     main(int argc,char *argv[]);
void    GetLhat(sdpdat*);
void    find_cut(sdpdat*);
void    GetVhat(double*,double*,sdpdat*);
void    GetUhat(chfac*,double*,double*);

/* 
 * functions in sdpprin.c 
 */ 
void    PrintHead(void);
int     FprintHead(void);
void    PrintEnd(clock_t*);
void    ShowMsgTit(int,sdpdat*,double);
void    ShowMsg(sdpdat*,int,double,double);

/*
 * functions in sdpproc.c
 */
void    SdpProc(sdpdat*,clock_t*);
void    PostProc(sdpdat*);

/*
 * functions in sdpshut.c
 */
void    ShutDown(void);

/*
 * functions in sdpupda.c
 */
void    GetValues(sdpdat*,double*,double*,double*,
                  double,double*,double,double,double*);
/*
 * functions in sdpupda.c
 */
void    ModMcutUpp(sdpdat*,double*,double*,double*,double*,double);
void    ModEcutUpp(sdpdat*,double*,double*,double*,double*,
                   double,double,double *,double);
void    ModUcutUpp(sdpdat*,double*,double*,double*,
                   double*,double,double,double*,double);
int     GetStep(sdpdat*,double*,double*,double,
                double,double,double*,int*,double*);
void    UpdSol(sdpdat*,double*,double*,double*,double);

/*
 * functions in sdplib.a
 */
int     InvMemo(sdpdat*);
void    FindSinv(chfac*,syoff*,double*,double*,double*,double,prbset);
void    FindInver(chfac*,syoff*,double*,sdpdat*);
void    DotProd(double*,chfac*);
void    Smt2Sym(smatx*,symat**);
void    ForwSubst(chfac*,double*,double*);
int     iSum(int,int*);
void    dCopy(int,double*,double*);
int     ChlFact(chfac*,int*,double*,int);
void    ChlSolve(chfac*,double*,double*);
int     PspSymbo(symat*,chfac*);
int     InvSymbo(syoff*,chfac*);
void    CfcInit(chfac*,symat*,double*); 
void    ModUppBnd(sdpdat*,double*,double*,double*,
                  double*,double,double,double*,double);
void    FixIndex(sdpdat*);
void    iCopy(int,int*,int*);
double  dDot(double*,double*,int);
void    RankReduce(sdpdat*);

/* Ecut modifications */
void FSubst(chfac*,double*,double*,int);
void GetEcut(double*,double*,sdpdat*);
void Uhat1( sdpdat*, double*, double*, double*, double);
void Uhat2( sdpdat*, double*, double*, double*, double, int, double);
