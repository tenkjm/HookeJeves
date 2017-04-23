/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   utestglobopt.hpp
 * Author: alusov
 *
 * Created on February 27, 2017, 11:00 AM
 */

#ifndef UTESTGLOBOPT_HPP
#define UTESTGLOBOPT_HPP


#include <limits>
#include "gtest/gtest.h"
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include "crystproblemfact.hpp"

#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
#include <methods/lins/quadls/quadls.hpp>
//#include <methods/gfsdesc/gfsdesc.hpp>
//#include <methods/coordesc/coordesc.hpp>
//#include <methods/varcoordesc/varcoordesc.hpp>
//#include <methods/varcoordesc/varcoordesc.hpp>


#include <hookejeevesRevorked/coorhjexplorer.hpp>
#include <hookejeevesRevorked/rndhjexplorer.hpp>
#include <hookejeevesRevorked/hookjeeves.hpp>
#include <hookejeevesRevorked/HookeJevesLinear.hpp>
#include <hookejeevesRevorked/hjexplorer.hpp>
#include <hookejeevesRevorked/varcoorhjexplorer.hpp>
#include <methods/lins/quadls/quadls.hpp>
#include "testfuncs/testfuncs.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"
#include "descfunc/descfunc.hpp"
#include "descfunc/keys.hpp"
#include "oneobj/contboxconstr/exprfunc.hpp"

#include <box/boxutils.hpp>

#define EPSILON 0.001
#define MAX_COUNT 10000

const char* JSONPATH;

using namespace snowgoose::expression;
using namespace OPTITEST;

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

class GlobOptTest : public ::testing::Test {
protected:

    GlobOptTest() : dfr(JSONPATH) {
    }

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

    void testglobopt(const std::string& key, const Expr<double>& expr, const Expr<Interval<double>> &exprInterval, int customDim = 0, double epsilon = EPSILON) {

        
        double lambda = 0.6;
        double inc = 1.1;
        double dec = 0.9;
        testgloboptCoorHJExplorer(key, expr, exprInterval, customDim, lambda, inc, dec);
        testgloboptTestVarHJ(key, expr, exprInterval, customDim, lambda, inc, dec);
        testgloboptTestHJLinear(key, expr, exprInterval, customDim, lambda, inc, dec);
        testgloboptTestHJRND(key, expr, exprInterval, customDim, lambda, inc, dec);
        ///testgloboptTestCoorDesk(key, expr, exprInterval, customDim);
       /// testgloboptTestVarCoorDesk(key, expr, exprInterval, customDim);
    }
    
    void testgloboptCoorHJExplorer(const std::string& key, const Expr<double>& expr, const Expr<Interval<double>> &exprInterval, 
     int customDim = 0,  double lambda=0, double inc=0, double dec=0, double epsilon = EPSILON) {

        cout << "HkeJeeves Coord explorer\n";
        auto desc = dfr.getdesr(key);
        int dim = desc.anyDim ? customDim : desc.dim;
        Bounds bounds = desc.anyDim ? Bounds(dim, desc.bounds[0]) : desc.bounds;

        // Setup problem
        OPTITEST::ExprProblemFactory fact(expr, bounds);
        COMPI::MPProblem<double> *mpp = fact.getProblem();

        double expected = desc.globMinY;
      
        int cnt = 0;
        int n = mpp->mBox->mDim;
        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        double v;
        LOCSEARCH::CoorHJExplorer<double> explr(*(mpp));
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;
        std::cout <<"search0\n";
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        mpp->mObjectives[0] = new COMPI::FuncCnt<double>(*(mpp->mObjectives.at(0)));
        std::cout <<" before v\n";
        v = mpp->mObjectives[0]->func(x);
        std::cout <<" before problem\n";
        COMPI::MPProblem<double> pr = *mpp;
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> ((pr).mObjectives[0]);
        
        obj->reset();
        std::cout <<" before search0\n";
        LOCSEARCH::MHookeJeeves<double> hjdesc(*(mpp), explr);
          std::cout <<"search1\n";
        hjdesc.getOptions().mLambda = lambda;
        hjdesc.getOptions().mInc = inc;
        hjdesc.getOptions().mDec = dec;
        hjdesc.filename = "TestHJstandart.csv";
        hjdesc.setStopper([&](double v, const double* x) {
            return false;
        });
        std::cout <<"search\n";
        hjdesc.search(x, v);
        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        

        ASSERT_NEAR(expected, v, epsilon);
    }
    
    
    
    
    
    void testgloboptTestVarHJ(const std::string& key, const Expr<double>& expr, const Expr<Interval<double>> &exprInterval, 
     int customDim = 0,double lambda=0, double inc=0, double dec=0, double epsilon = EPSILON) {

        cout << "HkeJeeves VarCoord explorer\n";
        auto desc = dfr.getdesr(key);
        int dim = desc.anyDim ? customDim : desc.dim;
        Bounds bounds = desc.anyDim ? Bounds(dim, desc.bounds[0]) : desc.bounds;

        // Setup problem
        OPTITEST::ExprProblemFactory fact(expr, bounds);
        COMPI::MPProblem<double> *mpp = fact.getProblem();

        double expected = desc.globMinY;
        
        int cnt = 0;
        int n = mpp->mBox->mDim;
        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        double v;
        LOCSEARCH::VarCoorHJExplorer<double> explr(*(mpp));
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        mpp->mObjectives[0] = new COMPI::FuncCnt<double>(*(mpp->mObjectives.at(0)));
        v = mpp->mObjectives[0]->func(x);
        
        COMPI::MPProblem<double> pr = *mpp;
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> ((pr).mObjectives[0]);
       
        obj->reset();
        
        LOCSEARCH::MHookeJeeves<double> hjdesc(*(mpp), explr);
        hjdesc.getOptions().mLambda = lambda;
        hjdesc.getOptions().mInc = inc;
        hjdesc.getOptions().mDec = dec;
        hjdesc.setStopper([&](double v, const double* x) {
            return false;
        });
        hjdesc.filename = "TestHJVar.csv";
        hjdesc.search(x, v);
        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        

        ASSERT_NEAR(expected, v, epsilon);
    }

    
    void testgloboptTestHJLinear(const std::string& key, const Expr<double>& expr, const Expr<Interval<double>> &exprInterval, 
   int customDim = 0, double lambda=0, double inc=0, double dec=0, double epsilon = EPSILON) {

        cout << "HkeJeeves Linear explorer\n";
        auto desc = dfr.getdesr(key);
        int dim = desc.anyDim ? customDim : desc.dim;
        Bounds bounds = desc.anyDim ? Bounds(dim, desc.bounds[0]) : desc.bounds;

        // Setup problem
        OPTITEST::ExprProblemFactory fact(expr, bounds);
        COMPI::MPProblem<double> *mpp = fact.getProblem();

        double expected = desc.globMinY;
        
        int cnt = 0;
        int n = mpp->mBox->mDim;
        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        double v;
        LOCSEARCH::CoorHJExplorer<double>  explr(*(mpp));
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        mpp->mObjectives[0] = new COMPI::FuncCnt<double>(*(mpp->mObjectives.at(0)));
        v = mpp->mObjectives[0]->func(x);
        
        COMPI::MPProblem<double> pr = *mpp;
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> ((pr).mObjectives[0]);
        
        obj->reset();
        
        LOCSEARCH::MHookeJeevesLinear<double> hjdesc(*(mpp), explr);
        hjdesc.getOptions().mLambda = lambda;
        hjdesc.getOptions().mInc = inc;
        hjdesc.getOptions().mDec = dec;
        hjdesc.filename = "TestHJlinear.csv";
        hjdesc.setStopper([&](double v, const double* x) {
            return false;
        });
        
        hjdesc.search(x, v);
        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        cout << "inside5\n";

        ASSERT_NEAR(expected, v, epsilon);
    }
    
    void testgloboptTestHJRND(const std::string& key, const Expr<double>& expr, const Expr<Interval<double>> &exprInterval, 
    int customDim = 0,double lambda=0, double inc=0, double dec=0, double epsilon = EPSILON) {

        cout << "HkeJeeves Rnd explorer\n";
        auto desc = dfr.getdesr(key);
        int dim = desc.anyDim ? customDim : desc.dim;
        Bounds bounds = desc.anyDim ? Bounds(dim, desc.bounds[0]) : desc.bounds;

        // Setup problem
        OPTITEST::ExprProblemFactory fact(expr, bounds);
        COMPI::MPProblem<double> *mpp = fact.getProblem();

        double expected = desc.globMinY;
        
        int cnt = 0;
        int n = mpp->mBox->mDim;
        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        double v;
        LOCSEARCH::RndHJExplorer<double>  explr(*(mpp));
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        mpp->mObjectives[0] = new COMPI::FuncCnt<double>(*(mpp->mObjectives.at(0)));
        v = mpp->mObjectives[0]->func(x);
        
        COMPI::MPProblem<double> pr = *mpp;
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> ((pr).mObjectives[0]);
        
        obj->reset();
        
        LOCSEARCH::MHookeJeevesLinear<double> hjdesc(*(mpp), explr);
        hjdesc.getOptions().mLambda = lambda;
        hjdesc.getOptions().mInc = inc;
        hjdesc.getOptions().mDec = dec;
         hjdesc.filename = "TestHJRnd.csv";
        hjdesc.setStopper([&](double v, const double* x) {
            return false;
        });
        
        hjdesc.search(x, v);
        std::cout << hjdesc.about() << "\n";
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        cout << "inside5\n";

        ASSERT_NEAR(expected, v, epsilon);
    }

/*
    void testgloboptTestCoorDesk(const std::string& key, const Expr<double>& expr, const Expr<Interval<double>> &exprInterval, int customDim = 0, double epsilon = EPSILON) {

        
        
        
        int cnt = 0;
        auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n) {
            cnt++;
            //std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
            if (cnt > 7000)
                return true;
            else
                return false;
        };
        cout << "TestCoorDesk explorer\n";
        auto desc = dfr.getdesr(key);
        int dim = desc.anyDim ? customDim : desc.dim;
        Bounds bounds = desc.anyDim ? Bounds(dim, desc.bounds[0]) : desc.bounds;

        // Setup problem
        OPTITEST::ExprProblemFactory fact(expr, bounds);
        COMPI::MPProblem<double> *mpp = fact.getProblem();

        double expected = desc.globMinY;
        cout << "inside2\n";
         
        int n = mpp->mBox->mDim;
        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        double v;
        LOCSEARCH::CoorDesc<double> desk(*mpp);
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        mpp->mObjectives[0] = new COMPI::FuncCnt<double>(*(mpp->mObjectives.at(0)));
        v = mpp->mObjectives[0]->func(x);
        cout << "inside3\n";
        
        cout << "inside33\n";
        COMPI::MPProblem<double> pr = *mpp;
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> ((pr).mObjectives[0]);
        cout << "inside331\n";
        obj->reset();
        cout << "inside333\n";
        
        
        desk.search(x, v);
         
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        cout << "inside5\n";

        ASSERT_NEAR(expected, v, epsilon);
    }
    
     void testgloboptTestVarCoorDesk(const std::string& key, const Expr<double>& expr, const Expr<Interval<double>> &exprInterval, int customDim = 0, double epsilon = EPSILON) {

        
        
        
        int cnt = 0;
        auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n) {
            cnt++;
            //std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
            if (cnt > 7000)
                return true;
            else
                return false;
        };
        cout << "TestVarCoorDesk explorer\n";
        auto desc = dfr.getdesr(key);
        int dim = desc.anyDim ? customDim : desc.dim;
        Bounds bounds = desc.anyDim ? Bounds(dim, desc.bounds[0]) : desc.bounds;

        // Setup problem
        OPTITEST::ExprProblemFactory fact(expr, bounds);
        COMPI::MPProblem<double> *mpp = fact.getProblem();

        double expected = desc.globMinY;
        cout << "inside2\n";
         
        int n = mpp->mBox->mDim;
        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        double v;
        LOCSEARCH::VarCoorDesc<double> desk(*mpp);
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        mpp->mObjectives[0] = new COMPI::FuncCnt<double>(*(mpp->mObjectives.at(0)));
        v = mpp->mObjectives[0]->func(x);
        cout << "inside3\n";
        
        cout << "inside33\n";
        COMPI::MPProblem<double> pr = *mpp;
        COMPI::FuncCnt<double> *obj = dynamic_cast<COMPI::FuncCnt<double>*> ((pr).mObjectives[0]);
        cout << "inside331\n";
        obj->reset();
        cout << "inside333\n";
        
        
        desk.search(x, v);
         
        std::cout << "In " << obj->mCounters.mFuncCalls << " function calls found v = " << v << "\n";
        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
        cout << "inside5\n";

        ASSERT_NEAR(expected, v, epsilon);
    }
    */
    DescFuncReader dfr;
};

//TEST_F(GlobOptTest, TestGlobOptAckley1) {
//    int N = 3;
//    auto expr = Ackley1<double>(N);
//    auto exprInterval = Ackley1<Interval<double>>(N);
//    testglobopt(K.Ackley1, expr, exprInterval, N);
//}

//TEST_F(GlobOptTest, TestGlobOptAckley2) {
//    int N = 4;
//    auto expr = Ackley2<double>(N);
//    auto exprInterval = Ackley2<Interval<double>>(N);
//    testglobopt(K.Ackley2, expr, exprInterval, N);
//}
//
TEST_F(GlobOptTest, TestGlobOptAckley3)
{
	auto expr = Ackley3<double>();
    auto exprInterval = Ackley3<Interval<double>>();
	testglobopt(K.Ackley3, expr, exprInterval);
}


//TEST_F(GlobOptTest, TestGlobOptAdjiman)
//{
//    auto expr = Adjiman<double>();
//    auto exprInterval = Adjiman<Interval<double>>();
//	testglobopt(K.Adjiman, expr, exprInterval);
//}
//
//TEST_F(GlobOptTest, TestGlobOptAlpine1)
//{  
//    int N = 3;
//    auto expr = Alpine1<double>(N);
//    auto exprInterval = Alpine1<Interval<double>>(N);
//    testglobopt(K.Alpine1, expr, exprInterval, N);
//}

#endif /* UTESTGLOBOPT_HPP */

