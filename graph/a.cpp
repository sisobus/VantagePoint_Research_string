#include <cstdio>
int main() {
    FILE *fp = fopen("g685","r");
    FILE *fout = fopen("gg685","w");
    int n,m;
    fscanf(fp,"%d %d",&n,&m);
    fprintf(fout,"%d %d\n",n,m);
    for ( int i = 0 ; i < m ; i++ ) {
        int a,b;
        double c;
        fscanf(fp,"%d %d %lf",&a,&b,&c);
        fprintf(fout,"%d %d %d\n",a,b,(int)c);
    }
    fclose(fp);
    fclose(fout);
    return 0;
}
