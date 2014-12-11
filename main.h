#include <cstdio>
#include <queue>
#include <ctime>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
using namespace std;
#define MAX_DIMENSION_POW 10
#define MAX_DIMENSION (1<<10)
#define INIT_SIZE 100
#define MAX_GENERATION 500
#define DOMINANCE_SIZE 0.1
#define MUTATION_CHANCE 0.01

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE !(FALSE)
#endif

#define MAX_CUT TRUE
#define EQU_CUT FALSE


double dabs(double a) {return a>0?a:-a;}
bool isRange(double a,double comp) {
    const double eps = 1e-9;
    if ( comp-eps <= a && a <= comp+eps ) return true;
    return false;
}
double f(const vector<string> &g,int n,int m) {
    double ret=0;
    for ( int i = 0 ; i < (int)g.size() ; i++ ) 
        for ( int j = 0 ; j < (int)g.size() ; j++ ) 
            if ( i != j ) {
                int now = 0;
                for ( int k = 0 ; k < (int)g[i].length() ; k++ ) 
                    now += (g[i][k] != g[j][k]);
                ret += dabs(sqrt((double)n/2.0f)-sqrt((double)now));
            }
    return ret;
}

struct Gene {
    int n,m;
    vector<string> g;
    Gene(int _n,int _m) {
        n = _n;
        m = _m;
        string now = "";
        for ( int i = 0 ; i < n ; i++ ) 
            now += "0";
        for ( int i = 0 ; i < m ; i++ ) 
            g.push_back(now);
    }
    inline bool operator < (const Gene& t) const {
        return f(g,n,m) < f(t.g,t.n,t.m);
    }
    inline bool operator > (const Gene& t) const {
        return f(g,n,m) > f(t.g,t.n,t.m);
    }
    inline bool operator == (const Gene& t) const {
        return f(g,n,m) == f(t.g,t.n,t.m);
    }
    void printGene(FILE *fp) {
        fprintf(fp,"f = %lf\n",f(g,n,m));
        for ( int i = 0 ; i < (int)g.size() ; i++ ) {
            for ( int j = 0 ; j < (int)g[i].length() ; j++ ) 
                fprintf(fp,"%c ",g[i][j]);
            fprintf(fp,"\n");
        }
        fprintf(fp,"\n");
    }
};

struct Graph {
    int u,v;
    double w;
    Graph(int _u,int _v,double _w) {
        u = _u;
        v = _v;
        w = _w;
    }
};

int noCut1,noCut2;

vector<Gene> twoPowOptimalPoints;
void makeTwoPowOptimalPoints();
int getDimensionPow(int now);

Gene initializingGene(int n,int m,int flag);
Gene makeGene0(int n,int m);
Gene makeGene1(int n,int m);
Gene makeGene2(int n,int m);
Gene makeGene3(int n,int m);
Gene makeGene4(int n,int m);
Gene makeGene5(int n,int m);
Gene makeGene6(int n,int m);
Gene makeGene7(int n,int m);
Gene makeGene8(int n,int m);
Gene makeRandomGene(int n,int m);

bool isSameVantagePointsInGene(Gene p,Gene other,int pos);
vector<Graph> makeGraph(Gene p,vector<int> v);
vector<Graph> makeGraph(Gene p,Gene other);
vector<Graph> makeGraph(Gene p);
vector<vector<int> > equalCut(int generation,int m,vector<Graph> &v,vector<int> &realCutData);
vector<vector<int> > maxCut(int generation,int m,vector<Graph> &v,vector<int> realCutData);
vector<Gene> crossOver(Gene p1,Gene p2,vector<vector<int> > &V1,vector<vector<int> > &V2);
void mutation(vector<Gene> &genes);
vector<int> changeVertex;
void maxCutForAdjustPoints(int generation,Gene p,vector<vector<int> > &now,int targetNumber,vector<int> realCutData);

double calculateNowGenesF(vector<Gene> &genes);
void deallocGenes(Gene gene);


vector<vector<int> > maxCut(int generation,int m,vector<Graph> &v);
