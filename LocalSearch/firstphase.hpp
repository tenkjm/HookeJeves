/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   firstphase.h
 * Author: anton
 *
 * Created on 31 января 2017 г., 2:16
 */

#ifndef FIRSTPHASE_H
#define FIRSTPHASE_H

#include <thread>
#include <common/locsearch.hpp>
#include <common/lineseach.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>

namespace LOCSEARCH {
template <typename FT> class HookeJeevesFirstPhase {
    public:        
        
    struct Options {
           /**
            * Initial value of granularity
            */
           FT mHInit = 1;

           /**
            * Increase in the case of success
            */
           FT mInc = 1;
           /**
            * Decrease in the case of failure
            */
           FT mDec = 0.5;
           /**
            * Lower bound for granularity
            */
           FT mHLB = 1e-08;
           /**
            * Upper bound on granularity
            */
           FT mHUB = 1e+02;
           /**
            * Reset parameter in each invocation
            */
           bool mResetEveryTime = false;
       };
    HookeJeevesFirstPhase() {
        mH=this->mOptions.mHInit;
    }
    
    
    
    FT explore(FT* x, const COMPI::MPProblem<FT>& mProb)
    {
        COMPI::Functor<FT>* obj = mProb.mObjectives.at(0);
        int n = mProb.mVarTypes.size();
        const snowgoose::Box<double>& box = *(mProb.mBox);
        FT fcur = obj->func(x);
        FT fold = fcur;       
        FT* mincoordinate = new FT[n];
        FT maxcoordinate[n];
          
        if (mOptions.mResetEveryTime)
                reset();
        for (;;) {                       
                for (int i = 0; i < n; i++) {                    
                    snowgoose::VecUtils::vecCopy(n, x, mincoordinate);
                    snowgoose::VecUtils::vecCopy(n, x, maxcoordinate);                   
                    mincoordinate[i] = mincoordinate[i] - mH;
                    maxcoordinate[i] = maxcoordinate[i] + mH;
                                        
                    if (mincoordinate[i] < box.mA[i]) {
                        mincoordinate[i] = box.mA[i];
                    }                   
                    if (maxcoordinate[i] > box.mB[i]) {
                        maxcoordinate[i] = box.mB[i];
                    }                   
                    FT fnmin;                       
                    std::thread thr(HookeJeevesFirstPhase<FT>::exploreOneSide, mincoordinate, obj,  std::ref(fnmin));                    
                    FT fnmax = obj->func(maxcoordinate);                      
                    thr.join();

                    if(fnmin<fnmax){
                        if (fnmin <= fcur) {
                            x[i] = mincoordinate[i];
                            fcur = fnmin;
                            continue;
                        }                  
                    }else{
                         if (fnmax <= fcur) {
                            x[i] = maxcoordinate[i]; 
                            fcur = fnmax;
                            continue;
                        }                           
                    }                                                          
                }
                
                if (fcur < fold) {
                    mH *= mOptions.mInc;
                    mH = SGMIN(mOptions.mHUB, mH);
                    break;
                } else {
                    mH *= mOptions.mDec;
                    if (mH <= mOptions.mHLB)
                        break;
                }
            }  
       
        return fcur;               
    }
   
    static void exploreOneSide( FT* x , COMPI::Functor<FT>* obj, FT& fn)
    {       
        fn = obj->func(x);       
    }
    Options& getMOptions()  {
        return mOptions;
    }

    void setMOptions(Options mOptions) {
        this->mOptions = mOptions;
    }

     /**
         * Reset the current vicinity size
         */
    void reset() {
        mH = mOptions.mHInit;
    }
    
    private:
        Options mOptions;
        FT mH;
    
};
}


#endif /* FIRSTPHASE_H */

