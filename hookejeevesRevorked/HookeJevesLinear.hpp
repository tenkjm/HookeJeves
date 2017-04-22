/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HookeJevesLinear.hpp
 * Author: anton
 *
 * Created on 26 марта 2017 г., 22:18
 */

#ifndef HOOKEJEVESLINEAR_HPP
#define HOOKEJEVESLINEAR_HPP
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
#include <fstream>
#include "hjexplorer.hpp"


namespace LOCSEARCH {

    /**
     * Modified Hooke-Jeeves Method
     */
    template <typename FT> class MHookeJeevesLinear : public COMPI::Solver <FT> {
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
        std::ofstream myfile;
        char * filename;
        void openFile(char* filename) {


            myfile.open(this->filename);

        }

        void closeFile() {
            myfile.close();
        }

        void appendToFile(const char* message) {
            myfile << message << "\n";
        }

        typedef std::function<bool(FT v, const FT* x) > Stopper;

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param explorer - reference to the explorer
         * @param ls - pointer to the line search
         */
        MHookeJeevesLinear(const COMPI::MPProblem<FT>& prob, HJExplorer<FT>& explorer, LineSearch<FT>* ls = nullptr) :
        mProblem(prob),
        mExplorer(explorer),
        mLS(ls),
        mStopper(defaultStopper) {
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
               openFile("optimisation_process_linear_Hooke_jeves.tx");
            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            snowgoose::BoxUtils::project(x, *(mProblem.mBox));
            FT fcur = obj->func(x);
            int n = mProblem.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProblem.mBox);

            FT lam = mOptions.mLambda;
            FT dir[n];
            FT xold[n];
            FT y[n];
            
            FT xtemp[n];
            
            auto step = [&] (FT* x1, FT* x2, FT * x3) {
                 
                if (mLS == nullptr) {
                    FT last_value = obj->func(x2);

                    for (int i = 0; i < n; i++) {
                        x3[i] = x2[i] + lam * (x2[i] - x1[i]);
                    }
                    if (lam > 0) {
                    //  std::cout<<"before loop\n";
                       while (last_value > obj->func(x3)) {
                            {
                                snowgoose::VecUtils::vecCopy(n, x3, xtemp);
                                lam *= mOptions.mInc;
                               // std::cout<<"inc2\n";
                                if(snowgoose::BoxUtils::isIn(x3,box)){
                                    last_value = obj->func(x3); 
                            //        std::cout<<"ok!\n";
                                }
                                else
                                {   
                              ///      std::cout<<"break2\n";
                                    break;
                                }
                                 for (int i = 0; i < n; i++) {
                                            x3[i] = x2[i] + lam * (x2[i] - x1[i]);
                                 }
                                snowgoose::BoxUtils::project(x3, box);
                            }
                            snowgoose::VecUtils::vecCopy(n, xtemp, x3);
                        }                                           
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
                    if (mStopper(fcur, x))
                        break;
                    FT fnew = mExplorer.explore(y, fcur);
                     //std::cout<<"after explore\n";
                    if (fnew < fcur) {
                        fcur = fnew;
                        
                        rv = true;
                        if (mOptions.mLambda > 0) {
                            snowgoose::VecUtils::vecCopy(n, x, xold);
                            snowgoose::VecUtils::vecCopy(n, y, x);
                           // std::cout<<"in loop\n";
			    appendToFile(snowgoose::VecUtils::vecPrint(n, x).c_str());
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

                closeFile();
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

            Options & getOptions() {
                return mOptions;
            }

            /**
             * Setup stopper
             * @param stopper 
             */
            void setStopper(Stopper && stopper) {
                mStopper = stopper;

            }
            private:

            static bool defaultStopper(FT v, const FT * x) {
                return false;
            }

            Stopper mStopper;
            HJExplorer<FT> &mExplorer;
            const COMPI::MPProblem<FT>& mProblem;
            Options mOptions;
            LineSearch<FT>* mLS;

        };
    
}



#endif /* HOOKEJEVESLINEAR_HPP */

