#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <map>

#include "Matrix.h"
#include "isomorphi.h"

using namespace std;

template<class BidirIt>
void Reverse(BidirIt first, BidirIt last){
    while ((first != last) && (first != --last)) {
        std::swap(**first++, **last);
    }
}

template<class BidirIt>
bool Next_permutation(BidirIt first, BidirIt last){
    if (first == last) return false;
    BidirIt i = last;
    if (first == --i) return false;

    while (true) {
        BidirIt i1, i2;
        i1 = i;
        if (**--i < **i1) {
            i2 = last;
            while (!(**i < **--i2));
            std::iter_swap(*i, *i2);
            Reverse(i1, last);
            return true;
        }
        if (i == first) {
            Reverse(first, last);
            return false;
        }
    }
}

bool isomorphi(const Matrix &X, const Matrix &Y,
               bool (*Compare)(const Matrix&, const Matrix&,
                               const vector<int>&),
               void (*getDeg)(const Matrix&, vector<int>&)){
    if(X.size()!=Y.size())
        return false;

    vector<int> xdeg,ydeg;
    getDeg(X, xdeg);
    getDeg(Y, ydeg);
    //vector<int> xdeg=getDeg(X);
    //vector<int> ydeg=getDeg(Y);

    vector<int> f(X.size());
    map<int, vector<int*> > xsets;
    map<int, vector<int> > ysets;
    for(unsigned i=0;i<X.size();i++){
        if(xsets.find(xdeg[i])==xsets.end())
            xsets[xdeg[i]]=vector<int*>(1, &f[i]);
        else xsets.find(xdeg[i])->second.push_back(&f[i]);

        if(ysets.find(ydeg[i])==ysets.end())
            ysets[ydeg[i]]=vector<int>(1, i);
        else ysets.find(ydeg[i])->second.push_back(i);
    }

    if(xsets.size()!=ysets.size())
        return false;
    auto it1=xsets.begin(); auto it2=ysets.begin();
    for(;it1!=xsets.end() && it2!=ysets.end(); it1++, it2++){
        if(it1->first!=it2->first)
            return false;
        if(it1->second.size()!=it2->second.size())
            return false;
        for(unsigned i=0;i<it1->second.size();i++){
            *it1->second[i]=it2->second[i];
        }
    }

    auto it=xsets.begin();
    while(true){
        if(Compare(X, Y, f))
            return true;
        if(!Next_permutation(it->second.begin(), it->second.end())){
            auto ei=++xsets.begin();
            for(;ei!=xsets.end() &&
                 !Next_permutation(ei->second.begin(),
                                   ei->second.end());ei++);
            if(ei==xsets.end())
                return false;
        }

    }
}

bool isomorph(const Matrix &X, const Matrix &Y, 
              bool (*Compare)(const Matrix&, const Matrix&, 
                              const vector<int>&)){
    if(X.size()!=Y.size())
        return false;
    
    vector<int> f(X.size());
    for(unsigned i=0;i<X.size();i++)
        f[i]=i;
    
    do{
        if(Compare(X, Y, f))
            return true;
    }while(next_permutation(f.begin(), f.end()));
    return false;
}

/*#define N 8
int main(){
    //ifstream inp("input.txt");
    ifstream inp("graphs.txt");
    Matrix A(N, inp),B(N, inp);
    //A.show(); cout<<endl; B.show(); cout<<endl;
    if(isomorphi(A, B, CompareGraphs, getIODeg)) cout<<"YES"<<endl;
    else cout<<"NO"<<endl;

    return 0;
}


        for(auto i:xsets){
            cout<<i.first<<": ";
            for(auto j:i.second)
                cout<<"["<<j<<", "<<*j<<"] ";
            cout<<endl;
        }
        for(auto i:f)
            cout<<i<<' ';
        cout<<endl;

    for(auto i:xdeg) cout<<i<<' ';
    cout<<endl;
    for(auto i:ydeg) cout<<i<<' ';
    cout<<endl;*/
