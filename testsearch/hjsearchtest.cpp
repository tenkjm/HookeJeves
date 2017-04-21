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

//#include <methods/hookejeeves/coorhjexplorer.hpp>

//#include <methods/hookejeeves/hookjeeves.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include "crystproblemfact.hpp"

#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
//#include <methods/lins/quadls/quadls.hpp>
//#include <methods/gfsdesc/gfsdesc.hpp>
#include <methods/coordesc/coordesc.hpp>
#include <methods/varcoordesc/varcoordesc.hpp>
//#include <methods/varcoordesc/varcoordesc.hpp>


#include <hookejeevesRevorked/coorhjexplorer.hpp>
#include <hookejeevesRevorked/rndhjexplorer.hpp>
#include <hookejeevesRevorked/hookjeeves.hpp>
#include <hookejeevesRevorked/HookeJevesLinear.hpp>
#include <hookejeevesRevorked/hjexplorer.hpp>
#include <hookejeevesRevorked/varcoorhjexplorer.hpp>
#include <methods/lins/quadls/quadls.hpp>
/*#include <methods/hookejeeves/coorhjexplorer.hpp>
#include <methods/hookejeeves/rndhjexplorer.hpp>
#include <methods/hookejeeves/hookjeeves.hpp>
 */
#include <stdlib.h>




using namespace std;

class RosenbrockFunction : public COMPI::Functor <double> {
public:

    double func(const double* x) {

        double x1 = x[0];
        double y1 = x[1];
        double a = 1.0 - x1;
        double b = y1 - x1* x1;
        return a * a + b * b * 100;
    }

    /**
     * Retrieve energy 
     * @return reference to energy
     */


};

class HJTester {
public:
    COMPI::MPProblem<double>& mpp;
    int n = 0;

    class CoorStopper {
    public:

        bool operator()(double xdiff, double fdiff, const std::vector<double>& gran, double fval, int n) {
            mCnt++;
            return false;
        }

        int mCnt = 0;
    };

    class MyStopper : public LOCSEARCH::QuadLS<double>::Stopper {
    public:

        bool stopnow(double s, int k, double vo, double vn) {
            mCnt++;
            if (s < 1e-3)
                return true;
            else if (k > 16)
                return true;
            else
                return false;
        }

        int mCnt = 0;
    };

    HJTester(COMPI::MPProblem<double>& _mpp) : mpp(_mpp) {

        n = mpp.mBox->mDim;

    }

  /*  void TestLS(char** argv) {


        int cnt = 0;


        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
        double v;
        LOCSEARCH::CoorHJExplorer<double> explr(mpp);

        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);

        v = mpp.mObjectives[0]->func(x);

        MyStopper stp;
        LOCSEARCH::QuadLS<double>* ls = new LOCSEARCH::QuadLS<double>(mpp, stp);
        LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, explr, ls);

        hjdesc.getOptions().mLambda = 0.2;

        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> (mpp.mObjectives[0]);
        obj->reset();
        hjdesc.search(x, v);
        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        SG_ASSERT(v <= 0.01);


    }
*/
    void TestHJ(char** argv, double lambda, double inc, double dec) {
        int cnt = 0;

        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
        double v;
        LOCSEARCH::CoorHJExplorer<double> explr(mpp);
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);

        v = mpp.mObjectives[0]->func(x);

        LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, explr);
        hjdesc.getOptions().mLambda = lambda;
        hjdesc.getOptions().mInc = inc;
        hjdesc.getOptions().mDec = dec;
	hjdesc.filename = "TestHJstandart.csv";
        hjdesc.setStopper([&](double v, const double* x) {
            return false;
        });

        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> (mpp.mObjectives[0]);
        obj->reset();
        hjdesc.search(x, v);

        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
        SG_ASSERT(v <= 0.01);

    }
    
    void TestHJLinear(char** argv, double lambda, double inc, double dec) {
        int cnt = 0;

        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
        double v;
        LOCSEARCH::CoorHJExplorer<double> explr(mpp);
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);

        v = mpp.mObjectives[0]->func(x);

        LOCSEARCH::MHookeJeevesLinear<double> hjdesc(mpp, explr);
        hjdesc.getOptions().mLambda = lambda;
        hjdesc.getOptions().mInc = inc;
        hjdesc.getOptions().mDec = dec;
	hjdesc.filename = "TestHJlinear.csv";
        hjdesc.setStopper([&](double v, const double* x) {
            return false;
        });

        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> (mpp.mObjectives[0]);
        obj->reset();
        hjdesc.search(x, v);

        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
        SG_ASSERT(v <= 0.01);

    }

    void TestVarHJ(char** argv, double lambda, double inc, double dec) {
        int cnt = 0;

        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
        double v;
        LOCSEARCH::VarCoorHJExplorer<double> explr(mpp);
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);

        v = mpp.mObjectives[0]->func(x);

        LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, explr);
        hjdesc.getOptions().mLambda = lambda;
        hjdesc.getOptions().mInc = inc;
        hjdesc.getOptions().mDec = dec;
        hjdesc.setStopper([&](double v, const double* x) {
            return false;
        });
	hjdesc.filename = "TestHJVar.csv";
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> (mpp.mObjectives[0]);
        obj->reset();
        std::cout<<"search\n";
        hjdesc.search(x, v);

        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
        SG_ASSERT(v <= 0.01);

    }

    void TestRnd(char** argv, double lambda, double inc, double dec) {

        int cnt = 0;


        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
        double v;
        LOCSEARCH::RndHJExplorer<double> explr(mpp);
        explr.getOptions().mHInit = 0.01;
        explr.getOptions().mHLB = 1e-6;
        explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);

        v = mpp.mObjectives[0]->func(x);

        LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, explr);
        hjdesc.getOptions().mLambda = lambda;
        hjdesc.getOptions().mInc = inc;
        hjdesc.getOptions().mDec = dec;
        hjdesc.filename = "TestHJRnd.csv";
        hjdesc.search(x, v);

        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << " iterations found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
        SG_ASSERT(v <= 0.01);



    }

    void TestCoorDesk() {
        int cnt = 0;
        auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n) {
            cnt++;
            //std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
            if (cnt > 7000)
                return true;
            else
                return false;
        };


        //LOCSEARCH::CoorDesc<double> desc(mpp, stopper);
        LOCSEARCH::CoorDesc<double> desc(mpp);

        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
        double v;
        v = mpp.mObjectives[0]->func(x);
        std::cout << "Initial v = " << v << "\n";
        std::cout << "Initial x = " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> (mpp.mObjectives[0]);
        obj->reset();
        bool rv = desc.search(x, v);
        std::cout << desc.about() << "\n";
        std::cout << "In " << cnt << " iterations and " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
    }

    void TestVarCoorDesk() {
        int cnt = 0;



        CoorStopper stp;
        LOCSEARCH::VarCoorDesc<double> desc(mpp);

        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
        double v;
        v = mpp.mObjectives[0]->func(x);
        std::cout << "Initial v = " << v << "\n";
        std::cout << "Initial x = " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> (mpp.mObjectives[0]);
        obj->reset();
        bool rv = desc.search(x, v);
        std::cout << desc.about() << "\n";
        std::cout << "In " << cnt << " iterations and " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
    }



};

COMPI::MPProblem<double>* getRosenbrock() {

    COMPI::MPProblem<double>* mpp = new COMPI::MPProblem<double>();


    mpp->mObjectives.push_back(new RosenbrockFunction());

    const int n = 2;
    mpp->mVarTypes.assign(n, COMPI::MPProblem<double>::VariableTypes::GENERIC);
    mpp->mBox = new snowgoose::Box<double>(n);

    for (int i = 0; i < n; i++) {

        // height
        mpp->mBox->mA[i] = -5;
        mpp->mBox->mB[i] = 5;

    }

    return mpp;
}

/*
 * 
 */
int main(int argc, char** argv) {

    COMPI::MPProblem<double>* mpp;
    double lambda = 1, inc = 1, dec = 1;
    if (argc == 1) {
        mpp = getRosenbrock();
    } else {
        CrystallProblemFactory cpf12(argv[1]);
        mpp = cpf12.get();

	lambda = std::stod(argv[2]);
	inc = std::stod(argv[3]);
	dec = std::stod(argv[4]);

        std::cout<<"lambda = "<< lambda << "\n";
    }
    mpp->mObjectives[0] = new COMPI::FuncCnt<double>(*(mpp->mObjectives.at(0)));

    HJTester hjtester(*mpp);

//    std::cout << "================================================================================CoorDesk Test" << endl;
//    hjtester.TestCoorDesk();
    
//    std::cout << "================================================================================VarCoorDesk Test" << endl;
//    hjtester.TestVarCoorDesk();

    std::cout << "================================================================================VarCoorDesk HJ Test" << endl;
    hjtester.TestVarHJ(argv, lambda, inc, dec);

 std::cout << "================================================================================HJ std" << endl;
 hjtester.TestHJ(argv, lambda, inc, dec);
    
    std::cout << "================================================================================HJ linear cofficient" << endl;
    hjtester.TestHJLinear(argv, lambda, inc, dec);
    
     
    
//#if 0    
    std::cout << "===============================================================HJ rnd" << endl;
    hjtester.TestRnd(argv , lambda, inc, dec);
//#endif

   // std::cout << "===============================================================HJ linear" << endl;
   // hjtester.TestLS(argv);


    return 0;
}

