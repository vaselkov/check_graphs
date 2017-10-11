#include <iostream>
#include <cstdio>
#include <time.h>
#include <cstdlib>
#include <algorithm>
#include <climits>
#include <cmath>
#include <string.h>

#include <vector>
#include <set>
#include <map>

#include "Matrix.h"
#include "isomorphi.h"

using namespace std;

const char* USAGE="\n\
    check [-g[R]|-s|-l|q]\n\
    [--count C] [--zprob P] [--left L] [--right R]\n\
    N infile";
const char* HELPTEXT="\
    Check invariants of semigroups or graphs.\n\n\
    Infile is the input file with graph or semigroups matrices.\n\n\
    Options:\n\
       -g: check invariants of graps\n\
       -s: check invariants of semigroups\n\
       -l: read all matrices from list and calculate invariants\n\
       -q: suppress auxiliary output\n\
       -R: generate C paris random graphs \n\
       --count C: set the number of generated graphs\n\
       --zprob P: set the probability(in percents) of zero in matrix\n\
       --left L and --right R: set the left and right border\n\
            of the number of edges\n\
       If -R is not set --count,--zprob,--left,--right will be ignored\n\
    \n\
       -R is incompatible with -s and -l";
const char* SEPARATOR="==============================\n";

const int MAX_ITERATIONS=100000;
const int ZPROB=50;
const int LEFT=1;
const int RIGHT=7;

int det(const Matrix&);
int Diameter(const Matrix&);
int Wiener(const Matrix&);
int edges(const Matrix&);

void getIODeg(const Matrix &X, vector<int>& odeg);

void nCount(const Matrix &a, vector<int>& vv);
void inv2(const Matrix &a, vector<int>& vv);


struct svInv{
    int (*Inv)(const Matrix&);
    void (*vInv)(const Matrix&, vector<int>&);
    const char* name;
    const char* format;
};

bool checkVInv(const Matrix&, const Matrix&, svInv, bool);
bool checkInv(const Matrix&, const Matrix&, svInv, bool);


svInv GraphInvs[]={
    //scalar        vector     name                  format
    {Wiener,        NULL,      "Wiener index    ",   "w(%c) = "      },
    {Diameter,      NULL,      "diameter        ",   "d(%c) = "      },
    {edges,         NULL,      "number of edges ",   "m(%c) = "      },
    {det,           NULL,      "determinant     ",   "det(%c) = "    },
    {NULL,          getIODeg,  "vertices degree ",   "deg(%c) = "    },
};

svInv SGrpInvs[]={
    //scalar        vector     name                  format
    {NULL,          nCount,    "inv1 ",              "inv1 of %c = " },
    {NULL,          inv2,      "inv2 ",              "inv2 of %c = " },
};


int gauss(vector<vector<double> >& in){
    int swaps=1;
    for(uint i=0;i<in.size()-1;i++){
        if(in[i][i]==0){
            uint ii=i+1;
            for(;ii<in.size() && in[ii][i]==0;ii++);
            if(ii<in.size()){
                swap(in[i], in[ii]);
                swaps*=-1;
            }
            else continue;
        }
        for(uint ii=i+1;ii<in.size();ii++){
            double c=-in[ii][i]/in[i][i];
            for(uint j=i;j<in.size();j++)
                in[ii][j]+=c*in[i][j];
        }
    }
    return swaps;
}

int det(const Matrix& in){
    double retval=1;
    vector<vector<double> > iin(in.size(), vector<double>(in.size()));
    for(uint i=0;i<in.size();i++)
        for(uint j=0;j<in.size();j++)
            iin[i][j]=(double)in[i][j];
    retval*=gauss(iin);
    for(uint i=0;i<in.size();i++)
        retval*=iin[i][i];
    return (int) retval;
}

void vReinit(vector<int>& vv, unsigned size){
    fill(vv.begin(), vv.end(), 0);
    if(vv.size()!=size)
        vv.resize(size, 0);
}

//vector<int> nCount(const Matrix &a){
void nCount(const Matrix &a, vector<int>& vv){
    //vector<int> vv(a.size(), 0);
    vReinit(vv, a.size());
    for(unsigned i=0;i<a.size();i++)
        for(unsigned j=0;j<a.size();j++)
            vv[a[i][j]-1]++;
    //return vv;
}

//vector<int> inv2(const Matrix &a){
void inv2(const Matrix &a, vector<int>& vv){
    vReinit(vv, a.size());
    //vector<int> vv(a.size(), 0);
    for(unsigned i=0;i<a.size();i++){
        set<int> s;
        unsigned its=0; unsigned right=i;
        do{
            its++;
            s.insert(a[right][i]-1);
            right=a[right][i]-1;
        }while(s.size()==its);
        vv[i]=s.size();
    }
    //return vv;
}

int npow(int base, int exp){
    return (int) pow((double) base, (double) exp);
}

//vector<int> getIODeg(const Matrix &X){
void getIODeg(const Matrix &X, vector<int>& odeg){
    //vector<int> odeg(X.size(), 0);
    vReinit(odeg, X.size());
    vector<int> ideg(X.size(), 0);
    int max=0;
    for(unsigned i=0;i<X.size();i++)
        for(unsigned j=0;j<X.size();j++){
            odeg[i]+=X[i][j];
            ideg[j]+=X[i][j];
            if(ideg[j]>max) max=ideg[j];
        }

    int n=(int) log10(max)+1;
    for(unsigned i=0;i<X.size();i++){
        odeg[i]=odeg[i]*npow(10, n)+ideg[i];
    }
    //return odeg;
}


bool CompareGraphs(const Matrix &a, const Matrix &b,
                   const vector<int>& f){
    for(unsigned i=0;i<a.size();i++)
        for(unsigned j=0;j<a.size();j++)
            if(a[i][j]!=b[f[i]][f[j]])
                return false;
    return true;
}


Matrix matrix_permutation(const Matrix &a, const vector<int>& f){
    Matrix b(a.size());
    for(unsigned i=0;i<a.size();i++)
        for(unsigned j=0;j<a.size();j++)
            b[f[i]][f[j]]=f[a[i][j]-1]+1;
    return b;
}

bool CompareSemigroups(const Matrix &a, const Matrix &b, const vector<int>& f){
    for(unsigned i=0;i<a.size();i++)
        for(unsigned j=0;j<a.size();j++)
            if(f[a[i][j]-1]!=b[f[i]][f[j]]-1)
                return false;
    return true;
}

int Min(int A, int B) {
    int Result = (A < B) ? A : B;
    if((A < 0) && (B >= 0)) Result = B;
    if((B < 0) && (A >= 0)) Result = A;
    if((A < 0) && (B < 0)) Result = -1;
    return Result;
}

void floyd_alg(const Matrix &a, Matrix &d){
    d.create(a.size());
    for(unsigned i=0;i<a.size();i++){
            for(unsigned j=0;j<a.size();j++)
                if(i!=j) d[i][j]=a[i][j] ? 1 : -1;
                else d[i][j]=0;
    }

    for(unsigned k=0;k<a.size();k++)
        for(unsigned i=0;i<a.size();i++)
            for(unsigned j=0;j<a.size();j++)
                if(d[i][k]!=-1 && d[k][j]!=-1){
                    d[i][j]=Min(d[i][j], d[i][k]+d[k][j]);
                }
}

int Diameter(const Matrix &a){
    int rtv=-1;
    for(unsigned i=0;i<a.size();i++)
        for(unsigned j=0;j<a.size();j++)
            if(a[i][j]>rtv) rtv=a[i][j];
    return rtv;
}

int Wiener(const Matrix &a){
    int rtv=0;
    for(unsigned i=0;i<a.size();i++)
        for(unsigned j=0;j<a.size();j++)
            if(a[i][j]!=-1) rtv+=a[i][j];
    return rtv;
}

int edges(const Matrix &a){
    return Wiener(a);
}

void printVector(const vector<int>& v){
    printf("(");
    for(unsigned i=0;i<v.size()-1;i++)
        printf("%d, ", v[i]);
    printf("%d)", v[v.size()-1]);
}

bool checkVInv(const Matrix& a, const Matrix& b,
               svInv inv, bool print){
    if(inv.Inv!=NULL)
        return checkInv(a, b, inv, print);
    if(inv.Inv==NULL && inv.vInv==NULL)
        return false;
    vector<int> i1,i2;
    inv.vInv(a, i1);
    inv.vInv(b, i2);
    sort(i1.begin(), i1.end());
    sort(i2.begin(), i2.end());
    if(print){
        printf(inv.format, 'A');
        printVector(i1);
        printf("\n");
        printf(inv.format, 'B');
        printVector(i2);
        printf("\n");
    }

    if(i1!=i2) return false;
    else return true;
}

bool checkInv(const Matrix& a, const Matrix& b,
               svInv inv, bool print){
    if(print) printf("\n");
    if(inv.vInv!=NULL)
        return checkVInv(a, b, inv, print);
    if(inv.Inv==NULL && inv.vInv==NULL)
        return false;
    int i1=inv.Inv(a);
    int i2=inv.Inv(b);

    if(print){
        printf(inv.format, 'A');
        printf("%d\n", i1);
        printf(inv.format, 'B');
        printf("%d\n", i2);
    }

    if(i1!=i2) return false;
    else return true;
}



void checkGraphs(const Matrix& A, const Matrix& B, bool minout){
    const char* sep = minout? "": SEPARATOR;
    Matrix dA,dB;
    floyd_alg(A, dA);
    floyd_alg(B, dB);

    for(unsigned i=0; i<sizeof(GraphInvs)/sizeof(GraphInvs[0]); i++){
        printf("%sChecking %s: ", sep, GraphInvs[i].name);
        bool a;
        if(i==0 || i==1)
            a=checkInv(dA, dB, GraphInvs[i], !minout);
        else a=checkInv(A, B, GraphInvs[i], !minout);
        if(a) printf("OK\n");
        else  printf("NO\n");
    }

    printf("%sChecking graphs isomorphism: ", sep);
    if(isomorphi(A, B, CompareGraphs, getIODeg))
        printf("YES\n");
    else printf("NO\n");
}

void toGraph(const Matrix& a, Matrix& b){
    b.create(a.size(), 0);
    for(unsigned i=0;i<a.size();i++)
        for(unsigned j=0;j<a.size();j++)
            b[i][a[i][j]-1]++;
}


void checkSGrps(const Matrix A, const Matrix B, bool minout){
    const char* sep = minout? "": SEPARATOR;
    for(unsigned i=0; i<sizeof(SGrpInvs)/sizeof(SGrpInvs[0]); i++){
        printf("%sChecking %s: ", sep, SGrpInvs[i].name);
        bool a;
        a=checkInv(A, B, SGrpInvs[i], !minout);
        if(a) printf("OK\n");
        else printf("NO");
    }

    printf("%sChecking semigroups as graphs:\n", sep);
    Matrix gA,gB;
    toGraph(A, gA);
    toGraph(B, gB);
    if(!minout){
        printf("gA=\n");
        gA.show();
        printf("gB=\n");
        gB.show();
        printf("\n");
    }

    checkGraphs(gA, gB, minout);

    printf("%sChecking groups isomorphism: ", sep);
    //if(isomorph(A, B, CompareSemigroups)) cout<<"- YES"<<endl;
    //else cout<<"- NO"<<endl;
    if(isomorphi(A, B, CompareSemigroups, inv2)) printf("YES\n");
    else printf("NO\n");
}

bool nextMatrix(FILE* in, Matrix& a){
    //a.create(n);
    unsigned i=0,j=0;
    char c;
    while(fscanf(in, "%c", &c)!=EOF){
        if(c!=' ' && c!='\n'){
            a[i][j]=(int) c-'0';
            if(++j==a.size()){
                j=0; i++;
            }
            if(i==a.size()){
                i=0;
                return true;
            }
        }
    }
    return false;
}

void randGraph(Matrix& a, const unsigned size,
               const int zeroPercent,
               int left, int right){
    a.create(size);

    if(right<left) swap(left, right);
    for(unsigned i=0;i<size;i++)
        for(unsigned j=0;j<size;j++)
            if((rand()%100) < zeroPercent)
                a[i][j]=0;
            else a[i][j]=rand() % (right-left+1) + left;
}

void findEqualGraphs(FILE* in, const unsigned N){
    vector<Matrix> sgrps;
    vector<Matrix> graphs;

    Matrix A(N);
    int count=0;
    for(;nextMatrix(in, A);count++){
        Matrix gA(A.size());
        sgrps.push_back(A);
        toGraph(A, gA);
        graphs.push_back(gA);
    }

    printf("finding graphs\n");
    for(int i=0;i<count;i++)
        for(int j=i+1;j<count;j++){
            if(isomorphi(graphs[i], graphs[j], CompareGraphs, getIODeg)){
                printf("A=\n");
                sgrps[i].show();
                printf("\ngA=\n");
                graphs[i].show();
                printf("\nB=\n");
                sgrps[j].show();
                printf("\ngB=\n");
                graphs[j].show();
                printf("\n");
                //return;
            }
        }
}


void checkRandGraphs(const unsigned size, const int iterations,
                     const int zeroPercent,
                     const int left, const int right){
    const unsigned N=sizeof(GraphInvs)/sizeof(GraphInvs[0]);

    vector<int> sameInv(N);
    int same=0;
    Matrix A(size), B(size);

    srand(time(NULL));
    for(int I=0;I<iterations;I++){
        randGraph(A, size, zeroPercent, left, right);
        randGraph(B, size, zeroPercent, left, right);
        //A.show();
        //printf("\n");
        //B.show();
        //printf("\n");
        Matrix dA(size), dB(size);
        floyd_alg(A, dA);
        floyd_alg(B, dB);
        bool eq=true;
        bool curEq[N];
        for(unsigned i=0; i<N; i++){
            bool a;
            if(i==0 || i==1)
                a=checkInv(dA, dB, GraphInvs[i], false);
            else a=checkInv(A, B, GraphInvs[i], false);
            if(a) curEq[i]=true;
            else  eq=curEq[i]=false;
        }
        if(eq && isomorphi(A, B, CompareGraphs, getIODeg))
            same++;
        for(unsigned i=0;i<N;i++)
            if(curEq[i]) sameInv[i]++;
    }

    for(unsigned i=0;i<N;i++)
        printf("Same %s: %5d, probability: %.5f\n",
                GraphInvs[i].name, sameInv[i],
                (double)sameInv[i]/(double)(iterations-same));
    printf("Isomorphs graphs     : %5d\n", same);
}

template<typename SV>
double getProb(map<SV, int>& inm, int N){
    if(N<2) return 0;
    double q=((double)N)*(N-1);
    double p=0;
    for(auto it: inm)
        p+=((double)it.second)*(it.second-1);
    return p/q;
}

void graphsFromList(FILE* in, unsigned N){
    unsigned CN=sizeof(GraphInvs)/sizeof(GraphInvs[0]);
    map<int, int> sInvs[CN];
    map<vector<int>, int> vInvs[CN];

    int count=0;
    Matrix A(N), d(N);
    vector<int> v(N);
    while(nextMatrix(in, A)){
        count++;
        floyd_alg(A, d);
        for(unsigned i=0;i<CN;i++){
            if(GraphInvs[i].vInv==NULL && GraphInvs[i].Inv!=NULL){
                int k;
                if(i==0 || i==1) k=GraphInvs[i].Inv(d);
                else k=GraphInvs[i].Inv(A);
                sInvs[i][k]++;
            }
            if(GraphInvs[i].vInv!=NULL && GraphInvs[i].Inv==NULL){
                GraphInvs[i].vInv(A, v);
                sort(v.begin(), v.end());
                vInvs[i][v]++;
            }
        }
    }


    printf("Total: %d\n", count);
    for(unsigned i=0;i<CN;i++){
        printf("Different %s: ", GraphInvs[i].name);
        if(GraphInvs[i].vInv==NULL && GraphInvs[i].Inv!=NULL)
            printf("%5lu, probability: %f\n", sInvs[i].size(),
                   getProb<int>(sInvs[i], count));
        if(GraphInvs[i].vInv!=NULL && GraphInvs[i].Inv==NULL)
            printf("%5lu, probability: %f\n", vInvs[i].size(),
                   getProb<vector<int> >(vInvs[i], count));
    }
}

void sgroupsFromList(FILE* in, unsigned N){
    unsigned gCN=sizeof(GraphInvs)/sizeof(GraphInvs[0]);
    map<int, int> gsInvs[gCN];
    map<vector<int>, int> gvInvs[gCN];

    unsigned sCN=sizeof(SGrpInvs)/sizeof(SGrpInvs[0]);
    map<int, int> SsInvs[sCN];
    map<vector<int>, int> SvInvs[sCN];

    int count=0;
    Matrix A(N);
    Matrix gA(N), d(N);
    vector<int> v(N);
    while(nextMatrix(in, A)){
        count++;
        for(unsigned i=0;i<sCN;i++){
            if(SGrpInvs[i].vInv==NULL && SGrpInvs[i].Inv!=NULL)
                SsInvs[i][SGrpInvs[i].Inv(A)]++;
            if(SGrpInvs[i].vInv!=NULL && SGrpInvs[i].Inv==NULL){
                SGrpInvs[i].vInv(A, v);
                sort(v.begin(), v.end());
                //SvInvs[i].insert(v);
                SvInvs[i][v]++;
            }
        }

        toGraph(A, gA);
        floyd_alg(gA, d);
        for(unsigned i=0;i<gCN;i++){
            if(GraphInvs[i].vInv==NULL && GraphInvs[i].Inv!=NULL){
                int k;
                if(i==0 || i==1) k=GraphInvs[i].Inv(d);
                else k=GraphInvs[i].Inv(gA);
                gsInvs[i][k]++;
            }
            if(GraphInvs[i].vInv!=NULL && GraphInvs[i].Inv==NULL){
                GraphInvs[i].vInv(gA, v);
                sort(v.begin(), v.end());
                gvInvs[i][v]++;
            }
        }

    }

    printf("total: %d\n", count);
    for(unsigned i=0;i<sCN;i++){
        printf("Different %s: ", SGrpInvs[i].name);
        //if(SGrpInvs[i].vInv==NULL && SGrpInvs[i].Inv!=NULL)
        //    printf("%5lu, probability: %f\n", SsInvs[i].size(),
        //            getProb<int>(SsInvs[i], count));
        if(SGrpInvs[i].vInv!=NULL && SGrpInvs[i].Inv==NULL){
            printf("%5lu, probability: %f\n", SvInvs[i].size(),
                    getProb<vector<int> >(SvInvs[i], count));
        }
    }
    printf("Graph invs:\n");
    for(unsigned i=0;i<gCN;i++){
        printf("Different %s: ", GraphInvs[i].name);
        if(GraphInvs[i].vInv==NULL && GraphInvs[i].Inv!=NULL)
            printf("%5lu, probability: %f\n", gsInvs[i].size(),
                    getProb<int>(gsInvs[i], count));
        if(GraphInvs[i].vInv!=NULL && GraphInvs[i].Inv==NULL)
            printf("%5lu, probability: %f\n", gvInvs[i].size(),
                    getProb<vector<int> >(gvInvs[i], count));
    }

    /*for(auto i: SvInvs[1]){
        printVector(i.first);
        printf(" - %d\n", i.second);
    }*/
}

int main(int argc, char *argv[]){
    bool gsw, ssw, lsw, Rsw, qsw;
    gsw=ssw=lsw=Rsw=qsw=false;

    bool cfl,zfl,lfl,rfl;
    cfl=zfl=lfl=rfl=false;

    int N=0;
    char* filename;
    int right=RIGHT,left=LEFT;
    int zprob=ZPROB;
    int count=MAX_ITERATIONS;

    if(argc>1 && (strcmp(argv[1], "-help")==0 ||
                strcmp(argv[1], "--help")==0)){
        printf("Usage: %s\n\n%s\n", USAGE, HELPTEXT);
        return 0;
    }

    bool badargs=false;
    bool czlr=false;
    int argnum=0;
    for(int i=1;i<argc && !badargs;i++){
        char* arg=argv[i];

        if(arg[0]!='-' &&
                (cfl || zfl || lfl || rfl)){
            if(cfl) count=atoi(arg);
            if(zfl) zprob=atoi(arg);
            if(lfl) left=atoi(arg);
            if(rfl) right=atoi(arg);
            czlr=true;
            cfl=zfl=lfl=rfl=false;
        }

        if(strcmp(arg, "--count")==0) cfl=true;
        if(strcmp(arg, "--zprob")==0) zfl=true;
        if(strcmp(arg, "--left")==0) lfl=true;
        if(strcmp(arg, "--right")==0) rfl=true;

        if(!(cfl || zfl || lfl || rfl) && !czlr){
            if(arg[0]=='-'){
                int j=0;
                while(arg[++j]!='\0'){
                    switch(arg[j]){
                        case('g'): gsw=true; break;
                        case('s'): ssw=true; break;
                        case('l'): lsw=true; break;
                        case('R'): Rsw=true; break;
                        case('q'): qsw=true; break;
                        default: badargs=true;
                    }
                }
            }
            else{
                if(argnum==0) N=atoi(arg);
                else if(argnum==1) filename=arg;
                else badargs=true;
                argnum++;
            }
        }
        czlr=false;
    }

    if(badargs){
        printf("Usage: %s\n", USAGE);
        printf("Use check -help to see a list of the options\n");
        return 0;
    }

    if(!(gsw != ssw)){
        printf("-g or -s are incompatible but one of them MUST be set\n");
        return 0;
    }

    if(Rsw && (ssw || lsw)){
        printf("-R is incompatible with -a and -l\n");
        return 0;
    }

    if(N<=0){
        printf("N is not correct\n");
        return 0;
    }

    if(!Rsw){
        FILE* input=fopen(filename, "r");
        if(input==NULL){
            printf("Unable to open file: %s\n", filename);
            return 0;
        }

        if(!lsw){
            Matrix A(N), B(N);
            nextMatrix(input, A);
            nextMatrix(input, B);
            if(gsw) checkGraphs(A, B, qsw);
            else if(ssw) checkSGrps(A, B, qsw);
        }
        else{
            if(gsw) graphsFromList(input, N);
            else if(ssw) sgroupsFromList(input, N);
        }
    }
    else{
        if(right<=left || right<0 || left<0){
            printf("R or L is not correct\n");
            return 0;
        }
        if(count<=0){
            printf("C is not correct\n");
            return 0;
        }
        if(zprob>100 || zprob<0){
            printf("P is not corect\n");
            return 0;
        }

        checkRandGraphs(N, count, zprob, left, right);
    }

    //findEqualGraphs(input, n);

    return 0;
}
