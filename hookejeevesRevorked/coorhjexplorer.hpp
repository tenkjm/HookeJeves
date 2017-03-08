/* 
 * File:   hjexplorer.hpp
 * Author: medved
 *
 * Created on March 5, 2016, 3:23 PM
 */

#ifndef COORHJEXPLORER_HPP
#define	COORHJEXPLORER_HPP

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
    template <class FT> class CoorHJExplorer : public HJExplorer <FT> {
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

        CoorHJExplorer(const COMPI::MPProblem<FT>& prob) : mProb(prob) {
            mH = mOptions.mHInit;
        }

        /**
         * Explore the vicinity of the 'x' point to find better value
         * @param x start vector on entry, resulting vector on exit 
         * @return new obtained value 
         */
          FT explore(FT* x)
    {
        COMPI::Functor<FT>* obj = mProb.mObjectives.at(0);
        int n = mProb.mVarTypes.size();
        const snowgoose::Box<double>& box = *(mProb.mBox);
        FT fcur = obj->func(x);
        FT fold = fcur;       
        FT* mincoordinate = new FT[n];
        FT maxcoordinate[n];
          
       // if (mOptions.mResetEveryTime)
       //         reset();
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
                    std::thread thr(exploreOneSide, mincoordinate, obj,  std::ref(fnmin));                    
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
     void decMH()
     {
         
         
         mH *= mOptions.mDec;
         
         
         
     }
     
     
    static void exploreOneSide( FT* x , COMPI::Functor<FT>* obj, FT& fn)
    {       
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
            mH = mOptions.mHInit;
        }



    private:

        const COMPI::MPProblem<FT>& mProb;
        Options mOptions;
        FT mH;

    };

};


#endif	/* HJEXPLORER_HPP */

