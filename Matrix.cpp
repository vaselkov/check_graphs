#include <iostream>
#include <fstream>

#include "Matrix.h"

using namespace std;

Matrix::Matrix(){
    msize=0; data=NULL;
}

Matrix::Matrix(unsigned n){
    msize=n;
    data=new int[msize*msize];
}

Matrix::Matrix(const Matrix& a){
    msize=a.size();
    data=new int[msize*msize];
    for(unsigned i=0;i<msize;i++)
        for(unsigned j=0;j<msize;j++)
            data[i*msize+j]=a[i][j];
}

Matrix::~Matrix(){
    delete[] data;
}

int* Matrix::operator[](unsigned row) const{
    return data+row*msize;
}

unsigned Matrix::size() const{
    return msize;
}

void Matrix::show(){
    if(!msize) std::cout<<"empty"<<endl;
    for(unsigned i=0;i<msize;i++){
        for(unsigned j=0;j<msize;j++)
            std::cout<<data[i*msize+j]<<' ';
        std::cout<<std::endl;
    }
}

void Matrix::show(ofstream& out){
    if(!msize) out<<"empty"<<std::endl;
    for(unsigned i=0;i<msize;i++){
        for(unsigned j=0;j<msize;j++)
            out<<data[i*msize+j]<<' ';
        out<<std::endl;
    }
}

void Matrix::clear(){
    msize=0;
    if(data) delete[] data; data=NULL;
}

void Matrix::create(unsigned n){
    if(n!=msize){
        clear();
        msize=n;
        data=new int[msize*msize];
    }
}

void Matrix::create(unsigned n, int vl){
    create(n);
    for(unsigned i=0;i<msize;i++)
        for(unsigned j=0;j<msize;j++)
            data[i*msize+j]=vl;
}
