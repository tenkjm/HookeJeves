/* 
 * File:   hjexplorer.hpp
 * Author: medved
 *
 * Created on March 5, 2016, 3:23 PM
 */

#ifndef RNDHJEXPLORER_HPP
#define	RNDHJEXPLORER_HPP

#include <stdlib.h>
#include <mpproblem.hpp>
#include "hjexplorer.hpp"
#include <thread>
#include <mputils.hpp>
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
     * Random  exploration to be used in Hooke and Jeeves Method
     */
    template <class FT> class RndHJExplorer : public HJExplorer <FT> {
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
            FT mHLB = 1e-01;
            /**
             * Upper bound on granularity
             */
            FT mHUB = 1e+02;
            /**
             * maximal number of tries before improvement
             */
            unsigned int mMaxTries = 16;
            /**
             * Reset parameter in each invocation
             */
            bool mResetEveryTime = false;

        };

        RndHJExplorer(const COMPI::MPProblem<FT>& prob) : mProb(prob) {
            this->mH = mOptions.mHInit;
        }


	 static void exploreOneSide( FT* x , COMPI::Functor<FT>* obj, FT& fn)
	    {       
		fn = obj->func(x);       
	    }

        /**
         * Explore the vicinity of the 'x' point to find better value
         * @param x start vector on entry, resulting vector on exit 
         * @return new obtained value 
         */

//	 FT explore(FT* x)
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
//		FT r = 0.5 - (FT) rand() / (FT) RAND_MAX;                 
//                    snowgoose::VecUtils::vecCopy(n, x, mincoordinate);
//                    snowgoose::VecUtils::vecCopy(n, x, maxcoordinate);  
//                 
//                    mincoordinate[i] = mincoordinate[i] - this->mH*r;
//                    maxcoordinate[i] = maxcoordinate[i] + this->mH*r;
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
//    }

/**
         * Explore the vicinity of the 'x' point to find better value
         * @param x start vector on entry, resulting vector on exit 
         * @return new obtained value 
         */
        FT explore(FT* x) {
            COMPI::Functor<FT>* obj = mProb.mObjectives.at(0);
            int n = mProb.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProb.mBox);
            FT fcur = obj->func(x);
            FT fold = fcur;
            if (mOptions.mResetEveryTime)
                reset();
            for (;;) {
                for (int i = 0; i < n; i++) {
                    FT y = x[i] - this->mH;
                    if (y < box.mA[i]) {
                        y = box.mA[i];
                    }
                    FT tmp = x[i];
                    x[i] = y;
                    FT fn = obj->func(x);
                    if (fn >= fcur) {
                        x[i] = tmp;
                    } else {
                        fcur = fn;
                        continue;
                    }

                    y = x[i] + this->mH;
                    if (y > box.mB[i]) {
                        y = box.mB[i];
                    }
                    tmp = x[i];
                    x[i] = y;
                    fn = obj->func(x);
                    if (fn >= fcur) {
                        x[i] = tmp;
                    } else {
                        fcur = fn;
                        continue;
                    }
                }
                if (fcur < fold) {
                    this->mH *= mOptions.mInc;
                    this->mH = SGMIN(mOptions.mHUB, this->mH);
                    break;
                } else {
                    this->mH *= mOptions.mDec;
                    if (this->mH <= mOptions.mHLB)
                        break;
                }
            }
            return fcur;
        }



       
        std::string about() const {
            std::ostringstream os;
            os << "Randomized exporer\n";
            os << "Initial granularity: " << mOptions.mHInit << "\n";
            os << "Increment coefficient: " << mOptions.mInc << "\n";
            os << "Decrement coefficient: " << mOptions.mDec << "\n";
            os << "Lower bound on granularity: " << mOptions.mHLB << "\n";
            os << "Upper bound on granularity: " << mOptions.mHUB << "\n";
            os << "Max Tries: " << mOptions.mMaxTries << "\n";
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
         void decMH()
     {
         
         
         this->mH *= mOptions.mDec;
         
         
         
     }
    
    private:

        const COMPI::MPProblem<FT>& mProb;
        Options mOptions;
        

    };

};


#endif	/* HJEXPLORER_HPP */

