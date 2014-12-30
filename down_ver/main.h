#include <cstdio>
#include <cassert>
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

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE !(FALSE)
#endif

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

vector<Gene> twoPowOptimalPoints;
void makeTwoPowOptimalPoints();
double getDistanceVariance(Gene& ans, int pos);
int getDimensionPow(int now);

