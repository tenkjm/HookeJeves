/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   hj.h
 * Author: anton
 *
 * Created on 31 января 2017 г., 2:17
 */

#ifndef HJ_H
#define HJ_H
#include "firstphase.hpp"
#include <common/locsearch.hpp>
#include <sstream>
#include <functional>
#include <common/locsearch.hpp>
#include <common/lineseach.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>
// basic file operations
#include <iostream>
#include <fstream>







 

namespace LOCSEARCH {
    template <typename FT> class HookeJeevesSearchAnton : public LocalSearch <FT> {
    

        public:
           
            
          struct Options {
            /**
             * Speedup multiplier
             */
            FT mLambda = 1;
            /**
             * Increment multiplier
             */
            FT mInc = 1;
            /**
             * Decrement multiplier
             */
            FT mDec = 0.5;
            /**
             * Lower bound on lambda
             */
            FT mLambdaLB = 0.01;
            
             
        };
        
         /**
         * Determines stopping conditions
         * @param xdiff - distance between next and previous x
         * @param fdiff - difference between next and previous f value
         * @param gram - current granularity
         * @param n - current step number
         */
        typedef std::function<bool(FT xdiff, FT fdiff, FT fval, int n) > Stopper;    
        HookeJeevesFirstPhase<FT>* firstPhase;    
         
        HookeJeevesSearchAnton(const COMPI::MPProblem<FT>& prob, const Stopper& stopper) :
        mProblem(prob), stopper(stopper) {
            
            
                        unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
            
            firstPhase = new HookeJeevesFirstPhase<FT>(); 
            
        }
        std::ofstream myfile;
        void openFile(char* filename)
        {
            
            
            myfile.open (filename);
            
        }
        
        void closeFile()
        {     
            myfile.close();
        }
        
        void appendToFile(const char* message)
        {
            myfile<<message<<"\n";
        }
        
        
         
        bool search(FT* x, FT& v){
            openFile("optimisation_process.tx");
            bool rv = false;
            
            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            snowgoose::BoxUtils::project(x, *(mProblem.mBox));
            FT fcur = obj->func(x);
            int n = mProblem.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            int sn = 0;
            FT lam = mOptions.mLambda;
            FT dir[n];
            FT xold[n];
            FT y[n];
            
            auto step = [&] (FT* x1, FT* x2, FT * x3) {
                
                    for (int i = 0; i < n; i++) {
                        x3[i] = x2[i] + lam * (x2[i] - x1[i]);
                    }
                
            };


            snowgoose::VecUtils::vecCopy(n, x, y);

            for (;;) {
                sn++;
                FT fnew = firstPhase->explore(y, mProblem);
                if (fnew < fcur) {
                    rv = true;

                    FT fdiff = fcur - fnew;
                    fcur = fnew;
                    snowgoose::VecUtils::vecCopy(n, x, xold);
                    snowgoose::VecUtils::vecCopy(n, y, x);
                    appendToFile(snowgoose::VecUtils::vecPrint(n, x).c_str());
                    FT xdiff = snowgoose::VecUtils::vecDist(n, xold, x);
                    step(xold, x, y);
                    lam *= mOptions.mInc;
                    if (stopper(xdiff, fdiff, fcur, sn))
                        break;

                } else {
                    if (lam > mOptions.mLambdaLB) {
                        lam *= mOptions.mDec;
                        snowgoose::VecUtils::vecCopy(n, x, y);
                    } else
                        break;
                }
            }
            v = fcur;
            closeFile();
            return rv;
            
        }

         
        std::string about() const {
            return "Hooke -Jeeves method for testing";
        }
     Options& getOptions() {
            return mOptions;
        }    

    private:         
        const COMPI::MPProblem<FT>& mProblem;
        Stopper stopper;
        Options mOptions;
        
    };
};






#endif /* HJ_H */

