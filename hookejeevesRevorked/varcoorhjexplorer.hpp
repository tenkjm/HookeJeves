/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   varcoorhjexplores.h
 * Author: anton
 *
 * Created on 26 марта 2017 г., 21:03
 */

#ifndef VARCOORHJEXPLORES_H
#define VARCOORHJEXPLORES_H

#include <sstream>
#include <stdlib.h>
#include <mpproblem.hpp>
#include "hjexplorer.hpp"
#include <common/utilmacro.hpp>
#include <iostream>
#include <thread>
#include <mputils.hpp>

namespace LOCSEARCH {

    /**
     * Traditional coordinate-based exploration used in Hooke and Jeeves Method
     */
    template <class FT> class VarCoorHJExplorer : public HJExplorer <FT> {
    public:

        struct Options {
            /**
             * Initial value of granularity
             */
            FT mHInit = 0.01;

            /**
             * Increase in the case of success
             */
            FT mInc = 1.75;
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
            
            std::vector<FT> mShifts;
        };

        VarCoorHJExplorer(const COMPI::MPProblem<FT>& prob) : mProb(prob) {
            this->mH = mOptions.mHInit;
            mOptions.mShifts.assign(prob.mVarTypes.size(), 1);
        }

        /**
         * Explore the vicinity of the 'x' point to find better value
         * @param x start vector on entry, resulting vector on exit 
         * @return new obtained value 
         */


        //         FT explore(FT* x)
        //    {
        //        COMPI::Functor<FT>* obj = mProb.mObjectives.at(0);
        //        int n = mProb.mVarTypes.size();
        //        const snowgoose::Box<double>& box = *(mProb.mBox);
        //        FT fcur = obj->func(x);
        //        FT fold = fcur;       
        //        FT* mincoordinate = new FT[n];
        //        FT maxcoordinate[n];
        //          
        //        if (mOptions.mResetEveryTime)
        //                reset();
        //        for (;;) {                       
        //                for (int i = 0; i < n; i++) {                    
        //                    snowgoose::VecUtils::vecCopy(n, x, mincoordinate);
        //                    snowgoose::VecUtils::vecCopy(n, x, maxcoordinate);  
        //                 
        //                    mincoordinate[i] = mincoordinate[i] - this->mH;
        //                    maxcoordinate[i] = maxcoordinate[i] + this->mH;
        //                                        
        //                    if (mincoordinate[i] < box.mA[i]) {
        //                        mincoordinate[i] = box.mA[i];
        //                    }                   
        //                    if (maxcoordinate[i] > box.mB[i]) {
        //                        maxcoordinate[i] = box.mB[i];
        //                    }                   
        //                    FT fnmin;                       
        //                    std::thread thr(exploreOneSide, mincoordinate, obj,  std::ref(fnmin));                    
        //                    FT fnmax = obj->func(maxcoordinate);                      
        //                    thr.join();
        //
        //                    if(fnmin<fnmax){
        //                        if (fnmin <= fcur) {
        //                            x[i] = mincoordinate[i];
        //                            fcur = fnmin;
        //                            continue;
        //                        }                  
        //                    }else{
        //                         if (fnmax <= fcur) {
        //                            x[i] = maxcoordinate[i]; 
        //                            fcur = fnmax;
        //                            continue;
        //                        }                           
        //                    }                                                          
        //                }
        //                
        //                if (fcur < fold) {
        //                    this->mH *= mOptions.mInc;
        //                    this->mH = SGMIN(mOptions.mHUB, this->mH);
        //                    break;
        //                } else {
        //                    this->mH *= mOptions.mDec;
        //                    if (this->mH <= mOptions.mHLB)
        //                        break;
        //                }
        //            }  
        //       
        //        return fcur;               
        //}

        FT explore(FT* x, FT fcur) {
            COMPI::Functor<FT>* obj = mProb.mObjectives.at(0);
            int n = mProb.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProb.mBox);
            FT fold = fcur;
            std::vector<FT> sft;
            sft.assign(n, mOptions.mHInit);
            if (mOptions.mResetEveryTime)
                reset();
            for (;;) {
                for (int i = 0; i < n; i++) {
                    FT h = sft[i];
                    FT dh = h *  mOptions.mShifts[i];
                    FT y = x[i] - dh;
                    if (y < box.mA[i]) {
                        y = box.mA[i];
                    }
                    FT tmp = x[i];
                    x[i] = y;
                    FT fn = obj->func(x);
                    if (fn >= fcur) {
                         FT t = h * mOptions.mDec;
                        
                        sft[i] = t;
                        x[i] = tmp;
                    } else {
                        FT t = h * mOptions.mInc;
                       
                        sft[i] = t;
                        fcur = fn;
                        continue;
                    }
                    y = x[i] + dh;
                    if (y > box.mB[i]) {
                        y = box.mB[i];
                    }
                    tmp = x[i];
                    x[i] = y;
                    fn = obj->func(x);
                    if (fn >= fcur) {
                        FT t = h * mOptions.mDec;
                       
                        sft[i] = t;
                        x[i] = tmp;
                    } else {
                        FT t = h * mOptions.mInc;
                       
                        sft[i] = t;
                        fcur = fn;
                        continue;
                    }
                }
               
                     
                    FT H = snowgoose::VecUtils::maxAbs(n, sft.data(), nullptr);
                    if (H <= mOptions.mHLB)
                        break;
                
            }
            return fcur;
        }

        void decMH() {
            this->mH *= mOptions.mDec;
        }

        static void exploreOneSide(FT* x, COMPI::Functor<FT>* obj, FT& fn) {
            fn = obj->func(x);
        }

        std::string about() const {
            std::ostringstream os;
            os << "Basic coordinate based exporer modifed for two side\n";
            os << "Initial granularity: " << mOptions.mHInit << "\n";
            os << "Increment coefficient: " << mOptions.mInc << "\n";
            os << "Decrement coefficient: " << mOptions.mDec << "\n";
            os << "Lower bound on granularity: " << mOptions.mHLB << "\n";
            os << "Upper bound on granularity: " << mOptions.mHUB << "\n";
            if (mOptions.mResetEveryTime)
                os << "Reset on every iteration\n";
            return os.str();
        }

        /**
         * Retrieve options
         * @return reference to options
         */
        Options& getOptions() {
            return mOptions;
        }

        /**
         * Reset the current vicinity size
         */
        void reset() {
            this->mH = mOptions.mHInit;
        }



    private:

        const COMPI::MPProblem<FT>& mProb;
        Options mOptions;


    };

};


#endif /* VARCOORHJEXPLORES_H */

