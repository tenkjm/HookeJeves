/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   hjser.cpp
 * Author: anton
 *
 * Created on 30 января 2017 г., 14:49
 */

#include <methods/hookejeeves/coorhjexplorer.hpp>
 
#include <methods/hookejeeves/hookjeeves.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include "crystproblemfact.hpp"
#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
#include <methods/lins/quadls/quadls.hpp>
#include <methods/gfsdesc/gfsdesc.hpp>
#include <methods/coordesc/coordesc.hpp>
#include <methods/varcoordesc/varcoordesc.hpp>
#include <anton/hj.h>
#include <methods/hookejeeves/coorhjexplorer.hpp>
#include <methods/hookejeeves/rndhjexplorer.hpp>
#include <methods/hookejeeves/hookjeeves.hpp>





using namespace std;


class HJStopper : public LOCSEARCH::MHookeJeeves<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double fval, int n) {
        //        std::cout << "n = " << n << "fdiff = " << fdiff << "\n";
        mCnt++;
        return false;
    }

    int mCnt = 0;
};



/*
 * 
 */
int main(int argc, char** argv) {
    HJStopper hjstp;
    
    CrystallProblemFactory cpf(argv[1]);
    COMPI::MPProblem<double>& mpp = *cpf.get();
    int cnt = 0;
    auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n) {
        cnt++;
        //std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
        if (cnt > 7000)
            return true;
        else
            return false;
    };
    
    auto stopper1 = [&](double xdiff, double fdiff,  double fval, int n) {
        cnt++;
        //std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
        if (cnt > 7000)
            return true;
        else
            return false;
    };
    LOCSEARCH::CoorDesc<double> desc(mpp, stopper);
    const int n = 12;
    double *x = new double[n];
    snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
    double v;
    v = mpp.mObjectives[0]->func(x);
    std::cout << "Initial v = " << v << "\n";
    std::cout << "Initial x = " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "In " << cnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
    
    //////////////////////////////////////////////////////////////////////////////////////
    cnt = 0;
    
    
    snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
    
    v = mpp.mObjectives[0]->func(x);
    
    
    LOCSEARCH::HookeJeevesSearchAnton<double>* desc1 = new LOCSEARCH::HookeJeevesSearchAnton< double>(mpp, stopper1);
 
    
    //desc1->firstPhase->getMOptions().mInc = 1;
    desc1->firstPhase->getMOptions().mHInit = 0.01;
    desc1->firstPhase->getMOptions().mHLB = 1e-6;
    desc1->firstPhase->getMOptions().mResetEveryTime = true;
    
    desc1->getOptions().mLambda = 1;
     
    rv = desc1->search(x, v);
    std::cout << desc1->about() << "\n";
    std::cout << "In " << cnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);
    
    LOCSEARCH::CoorHJExplorer<double> explr(mpp);
    explr.getOptions().mHInit = 0.01;
    explr.getOptions().mHLB = 1e-6;
    explr.getOptions().mResetEveryTime = true;
    
    snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
    
    v = mpp.mObjectives[0]->func(x);
    
    LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, hjstp, explr);
    hjdesc.getOptions().mLambda = 1;
    
    hjstp.mCnt = 0;
    hjdesc.search(x, v);
    
    std::cout << hjdesc.about() << "\n";
    std::cout << "In " << cnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);
    
    
    
    return 0;
}

