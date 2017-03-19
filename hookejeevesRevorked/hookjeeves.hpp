/* 
 * File:   bbboxdesc.hpp
 * Author: medved
 *
 * Created on November 3, 2015, 5:05 PM
 */

#ifndef HOOKEJEEVES_HPP
#define  HOOKEJEEVES_HPP

#include <sstream>
#include <functional>
#include <solver.hpp>
#include <common/lineseach.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>

#include "hjexplorer.hpp"


namespace LOCSEARCH {

    /**
     * Modified Hooke-Jeeves Method
     */
    template <typename FT> class MHookeJeeves : public COMPI::Solver <FT> {
    public:

        struct Options {
            /**
             * Speedup multiplier
             */
            FT mLambda = 1;
            /**
             * Increment multiplier
             */
            FT mInc = 1.2;
            /**
             * Decrement multiplier
             */
            FT mDec = 0.5;
            /**
             * Lower bound on lambda
             */
            FT mLambdaLB = 0.1;
        };

        
        typedef std::function<bool(FT v, const FT* x)> Stopper;
        
        

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param explorer - reference to the explorer
         * @param ls - pointer to the line search
         */
        MHookeJeeves(const COMPI::MPProblem<FT>& prob, HJExplorer<FT>& explorer, LineSearch<FT>* ls = nullptr) :
        mProblem(prob),
        mExplorer(explorer),
        mLS(ls),
        mStopper(defaultStopper)
        {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
        }

        /**
         * Perform search
         * @param x start point and result
         * @param v  the resulting value
         * @return true if search converged and false otherwise
         */
        bool search(FT* x, FT& v) {
            bool rv = false;

            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            snowgoose::BoxUtils::project(x, *(mProblem.mBox));
            FT fcur = obj->func(x);
            int n = mProblem.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProblem.mBox);

            FT lam = mOptions.mLambda;
            FT dir[n];
            FT xold[n];
            FT y[n];

            auto step = [&] (FT* x1, FT* x2, FT * x3) {
                if (mLS == nullptr) {
                    for (int i = 0; i < n; i++) {
                        x3[i] = x2[i] + lam * (x2[i] - x1[i]);
                    }
                } else {
                    FT vv;
                    snowgoose::VecUtils::vecSaxpy(n, x2, x1, -1., dir);
                    snowgoose::VecUtils::vecCopy(n, x2, x3);
                    mLS->search(dir, x3, vv);
                }
            };


            snowgoose::VecUtils::vecCopy(n, x, y);

            for (;;) {
                if(mStopper(fcur, x))
                    break;
                FT fnew = mExplorer.explore(y, fcur);
                if (fnew < fcur) {
                    fcur = fnew;
                    rv = true;
                    if (mOptions.mLambda > 0) {
                        snowgoose::VecUtils::vecCopy(n, x, xold);
                        snowgoose::VecUtils::vecCopy(n, y, x);
                        step(xold, x, y);
                        snowgoose::BoxUtils::project(y, box);
                        FT fstep = obj->func(y);
                        if (fstep > fcur) {
                            if (lam >= mOptions.mLambdaLB) {
                                lam *= mOptions.mDec;
                            }
                            snowgoose::VecUtils::vecCopy(n, x, y);
                        } else {
                            lam *= mOptions.mInc;
                            fcur = fstep;
                        }
                    }
                } else {
                    break;
                }
            }


            v = fcur;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Modified Hooke-Jeeves method\n";
            os << "Initial speedup coefficient: " << mOptions.mLambda << "\n";
            os << "Lower bound on speedup coefficient: " << mOptions.mLambdaLB << "\n";
            os << "Decrement multiplier: " << mOptions.mDec << "\n";
            os << "Increment multiplier: " << mOptions.mInc << "\n";
            os << "Explorer: \n" << mExplorer.about() << "\n";
            if (mLS != nullptr)
                os << "Line search: " << mLS->about() << "\n";
            return os.str();
        }

        Options& getOptions() {
            return mOptions;
        }

        /**
         * Setup stopper
         * @param stopper 
         */
        void setStopper(Stopper&& stopper) {
            mStopper = stopper;
            
        }
    private:

        static bool defaultStopper(FT v, const FT* x) {
            return false;
        }

        Stopper mStopper;
        HJExplorer<FT> &mExplorer;
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
        LineSearch<FT>* mLS;

    };
}

#endif 

