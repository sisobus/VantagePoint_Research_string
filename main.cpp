#include "main.h"

int main(int argc,char *argv[]) {
    srand((unsigned)time(NULL));
    int n,m;


    /*
     *  make 2^n optimal vantage points set
     *  vector<Gene> twoPowOptimalPoints
     */
    makeTwoPowOptimalPoints();
#ifdef DEBUG
    for ( int i = 0 ; i < MAX_DIMENSION_POW ; i++ ) {
        puts("");
        for ( int j = 0 ; j < twoPowOptimalPoints[i].m ; j++ ) {
            for ( int k = 0 ; k < twoPowOptimalPoints[i].n ; k++ ) 
                printf("%d ",twoPowOptimalPoints[i].g[j][k]);
            puts("");
        }
    }
#endif

    for ( int tot = 11 ; tot <= 35; tot+= 1 ) {
        /*
        system("rm graph/g*.max");
        system("rm graph/g*.out");
        system("rm graph/g*");
        */
        n = m = tot;
//        scanf("%d %d",&n,&m);
        /*
         *  make genes
         *  vector<Gene> genes;
         */
        vector<Gene> genes;
        for ( int i = 0 ; i < INIT_SIZE ; i++ ) 
            genes.push_back(initializingGene(n,m,i));
#ifdef DEBUG
        for ( int i = 0 ; i < INIT_SIZE ; i++ ) 
            genes[i].printGene(stdout);
#endif

        for ( int generation = 1 ; generation <= MAX_GENERATION ; generation++ ) {
            printf("current generation : %d\n",generation);
            /*
             *  print now generation information
             */
            char filename[222];
            FILE *fp;

            bool printOption = (!(generation%50));
//            bool printOption = true;
            //bool printOption = ( !(generation%1000) || (1 <= generation && generation <= 10) || generation == MAX_GENERATION);

            if ( printOption ) {
                sprintf(filename,"result/g_%d_%d.vp",n,generation);
                fp = fopen(filename,"w");
                fprintf(fp,"now generation = %d tot f = %lf\n",generation,calculateNowGenesF(genes));
                for ( int i = 0 ; i < genes.size() ; i++ ) 
                    genes[i].printGene(fp);
                fclose(fp);

                sprintf(filename,"result/g_%d_%d.result",n,generation);
                fp = fopen(filename,"w");
                fprintf(fp,"now generation = %d tot f = %lf\n",generation,calculateNowGenesF(genes));
            }

            /*
             *  min heap based on ff
             */
            priority_queue<Gene,vector<Gene>,greater<Gene> > pq;
            for ( int i = 0 ; i < genes.size() ; i++ ) {
                Gene now(genes[i].n,genes[i].m);
                for ( int j = 0 ; j < genes[i].m ; j++ ) 
                    for ( int k = 0 ; k < genes[i].n ; k++ ) 
                        now.g[j][k] = genes[i].g[j][k];
                pq.push(now);
            }
            genes.clear();

            if ( printOption ) {
                fprintf(fp,"min f = %lf\n",f(pq.top().g,pq.top().n,pq.top().m));

                fclose(fp);
            }

            for ( int i = 0 ; i < INIT_SIZE*DOMINANCE_SIZE ; i++ ) {
                Gene now = pq.top();pq.pop();

                Gene insertGene(now.n,now.m);
                for ( int j = 0 ; j < now.m ; j++ ) 
                    for ( int k = 0 ; k < now.n ; k++ ) 
                        insertGene.g[j][k] = now.g[j][k];
                genes.push_back(insertGene);
            }
            while ( !pq.empty() ) {
                if ( (int)pq.size() == 1 ) {
                    Gene tp1 = pq.top();pq.pop();
                    Gene p1(tp1.n,tp1.m);
                    for ( int i = 0 ; i < tp1.m; i++ ) 
                        for ( int j = 0 ; j < tp1.n ; j++ ) 
                            p1.g[i][j] = tp1.g[i][j];

                    genes.push_back(p1);
                    continue;
                }
                Gene tp1 = pq.top();pq.pop();
                Gene tp2 = pq.top();pq.pop();
                Gene p1(tp1.n,tp1.m);
                Gene p2(tp2.n,tp2.m);
                for ( int i = 0 ; i < tp1.m ; i++ ) 
                    for ( int j = 0 ; j < tp1.n ; j++ ) 
                        p1.g[i][j] = tp1.g[i][j];
                for ( int i = 0 ; i < tp2.m ; i++ ) 
                    for ( int j = 0 ; j < tp2.n ; j++ ) 
                        p2.g[i][j] = tp2.g[i][j];


                vector<Gene> nextGenes;
                nextGenes.push_back(Gene(p1.n,p1.m));
                nextGenes.push_back(Gene(p2.n,p2.m));
                int pos[2]={};

                Gene np1(p1.n,p1.m);
                Gene np2(p2.n,p2.m);
                int t_pos[2]={};

                for ( int i = 0 ; i < p1.m ; i++ ) {
                    if ( isSameVantagePointsInGene(p1,p2,i) ) {
                        for ( int j = 0 ; j < p1.n ; j++ ) 
                            nextGenes[0].g[pos[0]][j] = p1.g[i][j];
                        pos[0]++;
                    } else {
                        for ( int j = 0 ; j < p1.n; j++ ) 
                            np1.g[t_pos[0]][j] = p1.g[i][j];
                        t_pos[0]++;
                    }
                }
                np1.m = t_pos[0];
                for ( int i = 0 ; i < p2.m ; i++ ) {
                    if ( isSameVantagePointsInGene(p2,p1,i) ) {
                        for ( int j = 0 ; j < p2.n ; j++ ) 
                            nextGenes[1].g[pos[1]][j] = p2.g[i][j];
                        pos[1]++;
                    } else {
                        for ( int j = 0 ; j < p2.n ; j++ ) 
                            np2.g[t_pos[1]][j] = p2.g[i][j];
                        t_pos[1]++;
                    }
                }
                np2.m = t_pos[1];
                vector<Graph> g1,g2;
                g1 = makeGraph(np1);
                vector<vector<int> > V1,V2;

                V1 = maxCut(generation,np1.m,g1);
                for ( int i = 0 ; i < V1[0].size() ; i++ ) {
                    for ( int j = 0 ; j < np1.n ; j++ ) 
                        nextGenes[0].g[pos[0]][j] = np1.g[V1[0][i]][j];
                    pos[0]++;
                }
                for ( int i = 0 ; i < V1[1].size() ; i++ ) {
                    for ( int j = 0 ; j < np2.n ; j++ ) 
                        nextGenes[1].g[pos[1]][j] = np1.g[V1[1][i]][j];
                    pos[1]++;
                }

                g2 = makeGraph(np2);
                V2 = maxCut(generation,np2.m,g2);
                for ( int i = 0 ; i < V2[0].size() ; i++ ) {
                    for ( int j = 0 ; j < np2.n ; j++ ) 
                        nextGenes[0].g[pos[0]][j] = np2.g[V2[0][i]][j];
                    pos[0]++;
                }

                Gene now(np2.n,V2[1].size());
                for ( int i = 0 ; i < V2[1].size() ; i++ ) {
                    for ( int j = 0 ; j < np2.n; j++ ) 
                        now.g[i][j] = np2.g[V2[1][i]][j];
                }

                if ( pos[0] == m ) {
                    for ( int i = 0 ; i < now.m ; i++ ) {
                        for ( int j = 0 ; j < now.n ; j++ ) 
                            nextGenes[1].g[pos[1]][j] = now.g[i][j];
                        pos[1]++;
                    }
                }

                while ( pos[0] < m ) {
                    vector<Graph> ng = makeGraph(now);
                    vector<vector<int> > nV;
                    nV = maxCut(generation,now.m,ng);

                    for ( int i = 0 ; i < nV[0].size() ; i++ ) {
                        if ( pos[0] < m ) {
                            for ( int j = 0 ; j < now.n ; j++ ) 
                                nextGenes[0].g[pos[0]][j] = now.g[nV[0][i]][j];
                            pos[0]++;
                        } else {
                            nV[1].push_back(nV[0][i]);
                        }
                    }
                    if ( pos[0] == m ) {
                        for ( int i = 0 ; i < nV[1].size() ; i++ ) {
                            for ( int j = 0 ; j < now.n ; j++ ) 
                                nextGenes[1].g[pos[1]][j] = now.g[nV[1][i]][j];
                            pos[1]++;
                        }
                        break;
                    } else {
                        Gene next(now.n,nV[1].size());
                        int nPos = 0;
                        for ( int i = 0 ; i < nV[1].size() ; i++ ) {
                            for ( int j = 0 ; j < now.n ; j++ ) 
                                next.g[nPos][j] = now.g[nV[1][i]][j];
                            nPos++;
                        }
                        now.m = nV[1].size();
                        for ( int i = 0 ; i < next.m ; i++ ) 
                            for ( int j = 0 ; j < now.n ; j++ ) 
                                now.g[i][j] = next.g[i][j];
                    }
                }
                if ( f(nextGenes[0].g,nextGenes[0].n,nextGenes[0].m) <= 
                        f(nextGenes[1].g,nextGenes[1].n,nextGenes[1].m) ) {
                    genes.push_back(nextGenes[0]);
                    pq.push(nextGenes[1]);
                } else {
                    genes.push_back(nextGenes[1]);
                    pq.push(nextGenes[0]);
                }

                g1.clear();
                g2.clear();
                for ( int i = 0 ; i < V1.size() ; i++ ) 
                    V1[i].clear();
                V1.clear();
                for ( int i = 0 ; i < V2.size() ; i++ ) 
                    V2[i].clear();
                V2.clear();
            }
            mutation(genes);
        }
    }

    return 0;
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
int getDimensionPow(int now) {
    int ret = 0 ;
    for ( int i = 1 ;  (i<<1) <= now ; i <<= 1 )
        ret ++;
    return ret;
}

Gene initializingGene(int n,int m,int flag) {
    Gene (*makeGeneFunction[10])(int n,int m) = {
        makeGene0,makeGene1,makeGene2,makeGene3,makeGene4,makeGene5,
        makeGene6,makeGene7,makeGene8,makeRandomGene
    };

    /*
       if ( flag >= 10 ) return makeGeneFunction[9](n,m);
       else return makeGeneFunction[flag](n,m);
       */
    //if ( flag == 0 ) return makeGeneFunction[0](n,m);
    if ( flag < 10 ) return makeGeneFunction[0](n,m);
    return makeGeneFunction[9](n,m);
}
Gene makeGene0(int n,int m) {
    Gene ret(n,m);

    int pos = 0;
    int nowDimension = n;
    int twoPow = (1<<getDimensionPow(n));
    int totDimensionPow = getDimensionPow(twoPow);

    while ( nowDimension ) {
        if ( nowDimension < twoPow ) {
            twoPow /= 2;
            continue;
        }
        int nowDimensionPow = getDimensionPow(twoPow);
        Gene nowGene = twoPowOptimalPoints[nowDimensionPow];


        for ( int i = 0 ; i < nowGene.m ; i++ ) 
            for ( int j = 0 ; j < nowGene.n ; j++ ) 
                ret.g[i][pos+j] = nowGene.g[i][j];

        int nowN = nowGene.n;
        int posN = nowGene.n;
        for ( int i = 0 ; i < totDimensionPow - nowDimensionPow ; i++ ) {
            for ( int j = 0 ; j < nowN ; j++ ) 
                for ( int k = 0 ; k < nowGene.m ; k++ ) 
                    ret.g[posN+j][pos+k] = ret.g[j][pos+k]^1;
            posN += nowN;
            nowN *= 2;
        }

        pos += nowGene.m;
        nowDimension -= twoPow;
        twoPow /= 2;
    }
    int nn = (1<<getDimensionPow(m));
    for ( int i = 0 ; i < m-nn ; i++ ) 
        for ( int j = 0 ; j < n ; j++ ) 
            ret.g[nn+i][j] = ret.g[i][j]^1;
    return ret;
}
Gene makeGene1(int n,int m) {

}
Gene makeGene2(int n,int m) {

}
Gene makeGene3(int n,int m) {

}
Gene makeGene4(int n,int m) {

}
Gene makeGene5(int n,int m) {

}
Gene makeGene6(int n,int m) {

}
Gene makeGene7(int n,int m) {

}
Gene makeGene8(int n,int m) {

}
Gene makeRandomGene(int n,int m) {
    Gene ret(n,m);

    for ( int i = 0 ; i < m ; i++ ) {
hell:;
     for ( int j = 0 ; j < n ; j++ ) 
         ret.g[i][j] = rand()%2+'0';
     for ( int j = 0 ; j < i ; j++ ) {
         bool ok = false;
         for ( int k = 0 ; k < n ; k++ ) 
             if ( ret.g[i][k] != ret.g[j][k] ) 
                 ok = true;
         if ( !ok ) goto hell;
     }
    }

    return ret;
}

bool isSameVantagePointsInGene(Gene p,Gene other,int pos) {
    for ( int i = 0 ; i < other.m ; i++ ) {
        bool isSame = true;
        for ( int j = 0 ; j < other.n ; j++ ) 
            if ( p.g[pos][j] != other.g[i][j] ) 
                isSame = false;
        if ( isSame ) return true;
    }
    return false;
}
vector<Graph> makeGraph(Gene p,vector<int> v) {
    vector<Graph> ret;

    for ( int i = 0 ; i < (int)v.size() ; i++ )
        for ( int j = i+1 ; j < (int)v.size() ; j++ ) {
            int cnt = 0;
            for ( int k = 0 ; k < p.n ; k++ )
                cnt += (p.g[v[i]][k] != p.g[v[j]][k]);
            double weight = dabs(sqrt((double)p.n/2.0f)-sqrt((double)cnt));
            if ( isRange(weight,0.0f) ) weight = 0.1f;
            ret.push_back(Graph(i+1,j+1,weight));
        }
    return ret;
}
vector<Graph> makeGraph(Gene p,Gene other) {
    vector<Graph> ret;
    vector<int> v;
    for ( int i = 0 ; i < p.m ; i++ ) 
        if ( !isSameVantagePointsInGene(p,other,i) ) 
            v.push_back(i);
    for ( int i = 0 ; i < (int)v.size() ; i++ ) 
        for ( int j = i+1 ; j < (int)v.size() ; j++ ) {
            int cnt = 0;
            for ( int k = 0 ; k < p.n ; k++ ) 
                cnt += (p.g[v[i]][k] != p.g[v[j]][k]);
            double weight = dabs(sqrt((double)p.n/2.0f)-sqrt((double)cnt));
            if ( isRange(weight,0.0f) ) weight = 0.1f;
            ret.push_back(Graph(i+1,j+1,weight));
        }
    return ret;
}
vector<Graph> makeGraph(Gene p) {
    vector<Graph> ret;
    for ( int i = 0 ; i < p.m ; i++ ) 
        for ( int j = i+1 ; j < p.m ; j++ ) {
            int cnt = 0;
            for ( int k = 0 ; k < p.n ; k++ ) 
                cnt += (p.g[i][k] != p.g[j][k]);
            double weight = dabs(sqrt((double)p.n/2.0f)-sqrt((double)cnt));
            if ( isRange(weight,0.0f) ) weight = 0.1f;
            ret.push_back(Graph(i+1,j+1,weight));
        }
    return ret;
}

vector<vector<int> > equalCut(int generation,int m,vector<Graph> &v,vector<int> &realCutData) {
    vector<vector<int> > ret;
    ret.resize(2);

    char *filename = (char *)malloc(sizeof(char)*222);
    sprintf(filename,"graph/g%d",generation);
    FILE *fp = fopen(filename,"w");
    fprintf(fp,"%d %d\n",m,(int)v.size());
    for ( int i = 0 ; i < v.size() ; i++ ) 
        fprintf(fp,"%d %d %d\n",v[i].u,v[i].v,(int)(v[i].w+1.00001));
    fclose(fp);

    char *command = (char *)malloc(sizeof(char)*222);
    sprintf(command,"sdp1.1/sdp -p ecut -f %s",filename);
    system(command);

    sprintf(filename,"graph/g%d.eqt",generation);
    fp = fopen(filename,"r");

    for ( int i = 0 ; i < m ; i++ ) {
        int t;
        fscanf(fp,"%d",&t);
        if ( ~t ) ret[0].push_back(realCutData[i]);
        else ret[1].push_back(realCutData[i]);
    }
    fclose(fp);

    free(filename);
    free(command);

    return ret;
}
vector<vector<int> > maxCut(int generation,int m,vector<Graph> &v,vector<int> realCutData) {
    vector<vector<int> > ret;
    ret.resize(2);

    char *filename = (char *)malloc(sizeof(char)*222);
    sprintf(filename,"graph/g%d",generation);
    FILE *fp = fopen(filename,"w");
    fprintf(fp,"%d %d\n",m,(int)v.size());
    for ( int i = 0 ; i < v.size() ; i++ ) 
        fprintf(fp,"%d %d %d\n",v[i].u,v[i].v,(int)(v[i].w+1.00001));
    fclose(fp);

    char *command = (char *)malloc(sizeof(char)*222);
    sprintf(command,"sdp1.1/sdp -f %s 1>/dev/null",filename);
    system(command);

    sprintf(filename,"graph/g%d.max",generation);
    fp = fopen(filename,"r");

    vector<int> now;
    int oneCnt = 0;
    for ( int i = 0 ; i < m ; i++ ) {
        int t;
        fscanf(fp,"%d",&t);
        now.push_back(t);
        oneCnt += (t == 1);
    }
    fclose(fp);

    int minusOneCnt = m-oneCnt;
    if ( minusOneCnt < oneCnt ) {
        for ( int i = 0 ; i < m ; i++ ) {
            if ( now[i] == -1 ) ret[0].push_back(i);
            else ret[1].push_back(i);
        }
    } else {
        for ( int i = 0 ; i < m ; i++ ) {
            if ( now[i] == 1 ) ret[0].push_back(i);
            else ret[1].push_back(i);
        }
    }

    now.clear();
    free(filename);
    free(command);

    return ret;
}

vector<Gene> crossOver(Gene p1,Gene p2,vector<vector<int> > &V1,vector<vector<int> > &V2) {
    vector<Gene> ret;
    ret.push_back(Gene(p1.n,p1.m));
    ret.push_back(Gene(p2.n,p2.m));

    vector<int> samePointsPosition[2];
    for ( int i = 0 ; i < p1.m ; i++ ) 
        if ( isSameVantagePointsInGene(p1,p2,i) ) 
            samePointsPosition[0].push_back(i);
    if ( ~noCut1 ) 
        samePointsPosition[0].push_back(noCut1);
    for ( int i = 0 ; i < p2.m ; i++ ) 
        if ( isSameVantagePointsInGene(p2,p1,i) ) 
            samePointsPosition[1].push_back(i);
    if ( ~noCut2 ) 
        samePointsPosition[1].push_back(noCut2);

    int pos = 0;
    for ( int i = 0 ; i < V1[0].size() ; i++,pos++ ) 
        for ( int j = 0 ; j < p1.n ; j++ ) 
            ret[0].g[pos][j] = p1.g[V1[0][i]][j];
    for ( int i = 0 ; i < V2[1].size() ; i++,pos++ ) 
        for ( int j = 0 ; j < p2.n ; j++ ) 
            ret[0].g[pos][j] = p2.g[V2[1][i]][j];
    for ( int i = 0 ; i < samePointsPosition[0].size() ; i++,pos++ ) 
        for ( int j = 0 ; j < p1.n ; j++ ) 
            ret[0].g[pos][j] = p1.g[samePointsPosition[0][i]][j];

    pos = 0;
    for ( int i = 0 ; i < V1[1].size() ; i++,pos++ ) 
        for ( int j = 0 ; j < p1.n ; j++ ) 
            ret[1].g[pos][j] = p1.g[V1[1][i]][j];
    for ( int i = 0 ; i < V2[0].size() ; i++,pos++ ) 
        for ( int j = 0 ; j < p2.n ; j++ ) 
            ret[1].g[pos][j] = p2.g[V2[0][i]][j];
    for ( int i = 0 ; i < samePointsPosition[1].size() ; i++,pos++ ) 
        for ( int j = 0 ; j < p2.n ; j++ ) 
            ret[1].g[pos][j] = p2.g[samePointsPosition[1][i]][j];

    for ( int i = 0 ; i < 2 ; i++ ) 
        samePointsPosition[i].clear();
    return ret;
}
void mutation(vector<Gene> &genes) {
    for ( int i = 1 ; i < (int)genes.size() ; i++ )
        for ( int j = 0 ; j < genes[i].m ; j++ ) {
            double r = (double)rand()/RAND_MAX;
            if ( r >= MUTATION_ROW_CHANCE )
                continue;
            for ( int k = 0 ; k < genes[i].n ; k++ ) {
                r = (double)rand()/RAND_MAX;
                if ( r >= MUTATION_BIT_CHANCE )
                    continue;
                genes[i].g[j][k] ^= 1;
            }
        }
}
/*
void mutation(vector<Gene> &genes) {
    for ( int i = 1 ; i < genes.size() ; i++ ) 
        for ( int j = 0 ; j < genes[i].m ; j++ ) 
            for ( int k = 0 ; k < genes[i].n ; k++ ) {
                double r = (double)rand()/RAND_MAX;
                if ( r < MUTATION_CHANCE ) 
                    genes[i].g[j][k]^=1;
                bool eq = false;
                for ( int jj = 0 ; jj < genes[i].m ; jj++ ) {
                    if ( jj == j ) continue;
                    bool now = true;
                    for ( int kk = 0 ; kk < genes[i].n ; kk++ ) 
                        if ( genes[i].g[j][k] != genes[i].g[jj][kk] ) 
                            now = false;
                    if ( now ) {
                        genes[i].g[j][k] ^= 1;
                        break;
                    }
                }
            }
}
*/
void maxCutForAdjustPoints(int generation,Gene p,vector<vector<int> > &now,int targetNumber,vector<int> realCutData) {
    if ( (int)changeVertex.size() == targetNumber ) return;
    if ( (int)now[0].size()+(int)changeVertex.size() <= targetNumber ) {
        for ( int i = 0 ; i < (int)now[0].size() ; i++ ) 
            changeVertex.push_back(realCutData[now[0][i]]);
    }
    if ( (int)now[1].size()+(int)changeVertex.size() <= targetNumber) {
        for ( int i = 0 ; i < (int)now[1].size() ; i++ ) 
            changeVertex.push_back(realCutData[now[1][i]]);
    }
    if ( (int)changeVertex.size() == targetNumber ) return;
    Gene np1(p.n,now[0].size());
    Gene np2(p.n,now[1].size());

    int pos = 0;
    for ( int i = 0 ; i < now[0].size() ; i++,pos++ ) 
        for ( int j = 0 ; j < np1.n ; j++ ) 
            np1.g[pos][j] = p.g[realCutData[now[0][i]]][j];
    pos = 0;
    for ( int i = 0 ; i < now[1].size() ; i++,pos++ ) 
        for ( int j = 0 ; j < np2.n ; j++ ) 
            np2.g[pos][j] = p.g[realCutData[now[1][i]]][j];
    vector<Graph> ng1,ng2;
    ng1 = makeGraph(np1);
    ng2 = makeGraph(np2);
    vector<vector<int> > next[2];
    next[0] = maxCut(generation,np1.m,ng1,realCutData);
    next[1] = maxCut(generation,np2.m,ng2,realCutData);

    maxCutForAdjustPoints(generation,np1,next[0],targetNumber,realCutData);
    maxCutForAdjustPoints(generation,np2,next[1],targetNumber,realCutData);
}

double calculateNowGenesF(vector<Gene> &genes) {
    double ret = 0;
    for ( int i = 0 ; i < genes.size() ; i++ ) 
        ret += f(genes[i].g,genes[i].n,genes[i].m);
    return ret;
}
/*
void deallocGenes(Gene gene) {
    if ( gene.g == NULL ) return ;
    for ( int i = 0 ; i < gene.m ; i++ ) 
        free(gene.g[i]);
    free(gene.g);
}
*/

vector<vector<int> > maxCut(int generation,int m,vector<Graph> &v) {
    vector<vector<int> > ret;
    ret.resize(2);

    char *filename = (char *)malloc(sizeof(char)*222);
    sprintf(filename,"graph/g%d",generation);
    FILE *fp = fopen(filename,"w");
    fprintf(fp,"%d %d\n",m,(int)v.size());
    for ( int i = 0 ; i < v.size() ; i++ ) 
        fprintf(fp,"%d %d %d\n",v[i].u,v[i].v,(int)(v[i].w+1.00001));
    fclose(fp);

    char *command = (char *)malloc(sizeof(char)*222);
    sprintf(command,"sdp1.1/sdp -f %s > dev",filename);
    system(command);

    sprintf(filename,"graph/g%d.max",generation);
    fp = fopen(filename,"r");

    vector<int> now;
    int oneCnt = 0;
    for ( int i = 0 ; i < m ; i++ ) {
        int t;
        fscanf(fp,"%d",&t);
        now.push_back(t);
        oneCnt += (t == 1);
    }
    fclose(fp);

    int minusOneCnt = m-oneCnt;
    if ( minusOneCnt < oneCnt ) {
        for ( int i = 0 ; i < m ; i++ ) {
            if ( now[i] == -1 ) ret[0].push_back(i);
            else ret[1].push_back(i);
        }
    } else {
        for ( int i = 0 ; i < m ; i++ ) {
            if ( now[i] == 1 ) ret[0].push_back(i);
            else ret[1].push_back(i);
        }
    }

    now.clear();
    free(filename);
    free(command);

    return ret;
}
