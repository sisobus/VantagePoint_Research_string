#include <cstdio>
#include <cstring>
#include <string>
#include <cstdlib>
using namespace std;

int main() {
    for ( int i = 50 ; i <= 50 ; i++ ) {
        char filename[1111];
        sprintf(filename,"g_%d_1.vp",i);
        FILE *fp = fopen(filename,"r");
        char s[1111];
        fscanf(fp,"%[^\n]\n",s);
        fscanf(fp,"%[^\n]\n",s);

        sprintf(filename,"d%d.vp",i);
        FILE *fout = fopen(filename,"w");
        for ( int j = 0 ; j < i ; j++ ) {
            fscanf(fp,"%[^\n]\n",s);
            fprintf(fout,"%s\n",s);
        }
        fclose(fout);
        fclose(fp);
    }
    return 0;
}
