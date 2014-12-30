#include "main.h"

int N;
//int N = 78;
//#define DEBUG

int main(int argc,char *argv[]) {
    assert(argc >= 2);

    sscanf(argv[1],"%d",&N);

    char filename[(1<<10)];
    sprintf(filename,"%d.dat",N);
    FILE *fp = fopen(filename,"w");

    makeTwoPowOptimalPoints();

#ifdef DEBUG
    printf("%d\n",(int)twoPowOptimalPoints.size());
    for ( int i = 0 ; i < 5 ; i++ ) {
        twoPowOptimalPoints[i].printGene(stdout);
    }
#endif

    // N : want 
    int iWantDeleteSomeColumn; // 2^n - N 
    int originalColumn; // N <=2^n
    int optimalPosition;
    for ( int i = 0 ; i < 10 ; i++ ) 
        if ( (1<<i) >= N ) {
            iWantDeleteSomeColumn = (1<<i)-N;
            originalColumn = (1<<i);
            optimalPosition = i;
            break;
        }
#ifdef DEBUG
    printf("%d\n",iWantDeleteSomeColumn);
    printf("%d\n",originalColumn);
    printf("%d\n",optimalPosition);
    printf("%d\n",twoPowOptimalPoints[optimalPosition].n);
#endif

    if ( iWantDeleteSomeColumn == 0 ) {
        twoPowOptimalPoints[optimalPosition].printGene(fp);
        return 0;
    }
    Gene now(N,originalColumn);
    for ( int i = 0 ; i < originalColumn ; i++ ) 
        for ( int j = 0 ; j < N ; j++ ) 
            now.g[i][j] = twoPowOptimalPoints[optimalPosition].g[i][j];

    Gene ans(N,N);
    bool b[MAX_DIMENSION];
    memset(b,false,sizeof b);
    int pos = 0;
    for ( int i = 0 ; i < ans.n ; i++ ) 
        ans.g[pos][i] = '0';
    pos++;
    b[0] = true;
    for ( int i = 0 ; i < ans.n/2 ; i++ ) 
        ans.g[pos][i] = '0';
    for ( int i = ans.n/2 ; i < ans.n ; i++ ) 
        ans.g[pos][i] = '1';
    pos++;

    for ( ; pos < N ; pos++ ) {
#ifdef DEBUG
        printf("current position : %d\n",pos);
#endif

        pair<double,int> pdi;
        pdi.first = 1234567.0f;
        pdi.second = 0;
        for ( int i = 0 ; i < now.m ; i++ ) {
            if ( b[i] ) continue;
            for ( int j = 0 ; j < now.n ; j++ ) 
                ans.g[pos][j] = now.g[i][j];
            double distVariance = getDistanceVariance(ans,pos);
            if ( distVariance < pdi.first ) {
                pdi.first = distVariance;
                pdi.second = i;
            }
        }
        for ( int j = 0 ; j < now.n ; j++ ) 
            ans.g[pos][j] = now.g[pdi.second][j];
        b[pdi.second] = true;
    }

    ans.printGene(fp);

    fclose(fp);
    return 0;
}

double getDistanceVariance(Gene& ans, int pos) {
    double ret = 0.0;
    for ( int i = 0 ; i <= pos ; i++ ) 
        for ( int j = 0 ; j <= pos ; j++ ) 
            if ( i != j ) {
                int now = 0;
                for ( int k = 0 ; k < ans.n ; k++ ) 
                    now += (ans.g[i][k] != ans.g[j][k]);
                ret += dabs(sqrt((double)ans.n/2.0f)-sqrt((double)now));
            }
    return ret;
}

void makeTwoPowOptimalPoints() {
    for ( int i = 0 ; i < MAX_DIMENSION_POW ; i++ ) {
        Gene now((1<<i),(1<<i));
        if ( i == 0 ) {
            now.g[0][0] = '0';
        } else {
            for ( int j = 0 ; j < twoPowOptimalPoints[i-1].m ; j++ ) {
                for ( int k = 0 ; k < twoPowOptimalPoints[i-1].n ; k++ ) 
                    now.g[j*2][k] = twoPowOptimalPoints[i-1].g[j][k];
                for ( int k = 0 ; k < twoPowOptimalPoints[i-1].n ; k++ ) 
                    now.g[j*2][k+twoPowOptimalPoints[i-1].m] = twoPowOptimalPoints[i-1].g[j][k];
                for ( int k = 0 ; k < twoPowOptimalPoints[i-1].n ; k++ ) 
                    now.g[j*2+1][k] = twoPowOptimalPoints[i-1].g[j][k];
                for ( int k = 0 ; k < twoPowOptimalPoints[i-1].n ; k++ ) 
                    now.g[j*2+1][k+twoPowOptimalPoints[i-1].m] = twoPowOptimalPoints[i-1].g[j][k]^1;
            }
        }
        twoPowOptimalPoints.push_back(now);
    }
}
