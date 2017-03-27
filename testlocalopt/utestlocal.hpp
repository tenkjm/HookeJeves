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

  GlobOptTest() : dfr(JSONPATH)
  {
  }
  
  virtual void SetUp() 
  {
		
  }
  
  virtual void TearDown() 
  {
      
  }

  void testglobopt(const std::string& key,const Expr<double>& expr, const Expr<Interval<double>> & exprInterval, int customDim = 0, double epsilon = EPSILON) 
  {

    cout<<"inside\n";
    auto desc = dfr.getdesr(key);
    int dim = desc.anyDim ? customDim : desc.dim;
    Bounds bounds = desc.anyDim ? Bounds(dim, desc.bounds[0]) : desc.bounds;

    // Setup problem
    OPTITEST::ExprProblemFactory fact(expr, bounds);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
		
    double expected = desc.globMinY;
    cout<<"inside2\n";
     int cnt = 0;
     int n = mpp->mBox->mDim;
        double *x = new double[n];
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
        double v;
        LOCSEARCH::CoorHJExplorer<double> explr(*(mpp));
        //explr.getOptions().mHInit = 0.01;
        //explr.getOptions().mHLB = 1e-6;
        //explr.getOptions().mResetEveryTime = false;

        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);

        v = mpp->mObjectives[0]->func(x);
cout<<"inside3\n";
        LOCSEARCH::MHookeJeeves<double> hjdesc(*(mpp), explr);
        hjdesc.getOptions().mLambda = 0;
        hjdesc.getOptions().mInc = 1;
        hjdesc.getOptions().mDec = 1;
        hjdesc.setStopper([&](double v, const double* x) {
            return false;
        });

        hjdesc.search(x, v);
        
        cout<<"inside5\n";
        
    ASSERT_NEAR(expected, v, epsilon);
}

DescFuncReader dfr;
};


	 
TEST_F(GlobOptTest, TestGlobOptAckley1)
{
    int N = 3;
	auto expr = Ackley1<double>(N);
    auto exprInterval = Ackley1<Interval<double>>(N);
    testglobopt(K.Ackley1, expr, exprInterval, N);
}

TEST_F(GlobOptTest, TestGlobOptAckley2)
{
	int N = 4;
	auto expr = Ackley2<double>(N);
    auto exprInterval = Ackley2<Interval<double>>(N);
	testglobopt(K.Ackley2, expr, exprInterval, N);
}



#endif /* UTESTGLOBOPT_HPP */

