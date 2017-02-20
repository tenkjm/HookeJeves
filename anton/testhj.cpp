/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testhj.cpp
 * Author: anton
 *
 * Created on 31 января 2017 г., 0:42
 */

/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <iostream>
#include <oneobj/contboxconstr/dejong.hpp>
#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
#include <methods/lins/quadls/quadls.hpp>

#include "hookejeeves.h"

 
 
  bool stopnow(double xdiff, double fdiff, double fval, int n) {
        //        std::cout << "n = " << n << "fdiff = " << fdiff << "\n";
        mCnt++;
        return false;
    }

    int mCnt = 0;

/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 10;
    OPTITEST::DejongProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);

   
    LOCSEARCH::HookeJeevesSearchAnton<double> desc(*mpp, stopnow);
 
    
    desc.getOptions().mInc = 1;
    desc.getOptions().mLambdaLB = desc.getOptions().mLambda;

    double x[n];

    for (int i = 0; i < n; i++)
        x[i] = i * 100 + 1;
    double v;
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "In " << mCnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);
    return 0;
}


