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
#include <methods/lins/quadls/quadls.hpp>
//#include <methods/gfsdesc/gfsdesc.hpp>
#include <methods/coordesc/coordesc.hpp>
//#include <methods/varcoordesc/varcoordesc.hpp>


#include <hookejeevesRevorked/coorhjexplorer.hpp>
#include <hookejeevesRevorked/rndhjexplorer.hpp>
#include <hookejeevesRevorked/hookjeeves.hpp>
#include <hookejeevesRevorked/hjexplorer.hpp>
#include <methods/lins/quadls/quadls.hpp>
/*#include <methods/hookejeeves/coorhjexplorer.hpp>
#include <methods/hookejeeves/rndhjexplorer.hpp>
#include <methods/hookejeeves/hookjeeves.hpp>
*/




using namespace std;


class HJStopper : public LOCSEARCH::MHookeJeeves<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double fval, int n) {
        //        std::cout << "n = " << n << "fdiff = " << fdiff << "\n";
        mCnt++;
       // std::cout<<n<<endl;
        
        if(xdiff<0.0000001f){
            
            return true;
        
        }
        if(n>7000)
        {
            
            return true;
        }
        
        return false;
    }

    int mCnt = 0;
};


 class RosenbrockFunction : public COMPI::Functor <double> {
    public:
        
      
        double func(const double* x) {
            
            double x1 = x[0];
            double y1 = x[1];
            double a = 1.0 - x1;
            double b = y1 - x1* x1;
            return a*a + b*b*100;
        }

        /**
         * Retrieve energy 
         * @return reference to energy
         */
       

    };


 class HJTester   {
    public:
       COMPI::MPProblem<double>& mpp;
       
      
    class MyStopper : public LOCSEARCH::QuadLS<double>::Stopper {
public:

    bool stopnow(double s, int k, double vo, double vn){
        mCnt++;
        if (s < 1e-3)
            return true;
        else if(k > 16)
            return true;
        else
            return false;
    }

    int mCnt = 0;
};

         HJTester(COMPI::MPProblem<double>& _mpp):mpp(_mpp){
           
            
           
       }


	void TestHJ(char** argv)
	{
           
                
            
		 
		    
		    int cnt = 0;
		  
                        HJStopper hjstp;
                        int n = 9;
                        double *x = new double[n];
                        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
                        double v;
			LOCSEARCH::CoorHJExplorer<double> explr(mpp);
                        explr.getOptions().mHInit = 0.01;
                        explr.getOptions().mHLB = 1e-6;
                        explr.getOptions().mResetEveryTime = false;

                        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);

                        v = mpp.mObjectives[0]->func(x);

                        LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, hjstp, explr);
                        hjdesc.getOptions().mLambda = 1;

                        hjstp.mCnt = 0;
                        hjdesc.search(x, v);

                        std::cout << hjdesc.about() << "\n";
                        std::cout << "In " << hjstp.mCnt << " iterations found v = " << v << "\n";
                        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
                        //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
                        SG_ASSERT(v <= 0.01);

	}


	void TestLS(char** argv)
	{

                    
		    int cnt = 0;
		     
                        HJStopper hjstp;
                        int n = 9;
                        double *x = new double[n];
                        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
                        double v;
			LOCSEARCH::CoorHJExplorer<double> explr(mpp);
                        explr.getOptions().mHInit = 0.01;
                        explr.getOptions().mHLB = 1e-6;
                        explr.getOptions().mResetEveryTime = false;

                        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);

                        v = mpp.mObjectives[0]->func(x);
                        
                        MyStopper stp;
                        LOCSEARCH::QuadLS<double>* ls = new LOCSEARCH::QuadLS<double>(mpp, stp);
                        LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, hjstp, explr, ls);
                        
                         hjdesc.getOptions().mLambda = 1;

                        hjstp.mCnt = 0;
                        hjdesc.search(x, v);

                        std::cout << hjdesc.about() << "\n";
                        std::cout << "In " << hjstp.mCnt << " iterations found v = " << v << "\n";
                        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
                        //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
                        SG_ASSERT(v <= 0.01);


	}

	void TestRnd(char** argv)
	{
                 
		    int cnt = 0;
		     
                        HJStopper hjstp;
                        int n = 9;
                        double *x = new double[n];
                        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
                        double v;
			LOCSEARCH::RndHJExplorer<double> explr(mpp);
                        explr.getOptions().mHInit = 0.01;
                        explr.getOptions().mHLB = 1e-6;
                        explr.getOptions().mResetEveryTime = false;

                        snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);

                        v = mpp.mObjectives[0]->func(x);

                        LOCSEARCH::MHookeJeeves<double> hjdesc(mpp, hjstp, explr);
                        hjdesc.getOptions().mLambda = 1;

                        hjstp.mCnt = 0;
                        hjdesc.search(x, v);

                        std::cout << hjdesc.about() << "\n";
                        std::cout << "In " << hjstp.mCnt << " iterations found v = " << v << "\n";
                        std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
                        //std::cout << "Number of objective calls is " << mpp.mObjectives[0]->mCounters.mFuncCalls << "\n";
                        SG_ASSERT(v <= 0.01);



        }
        void TestCoorDesk()
        {
            
            
            
            
             
            int cnt = 0;
            auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n) {
                cnt++;
                //std::cout << "cnt = " << cnt << ", fval =" << fval << "\n";
                if (cnt > 7000)
                    return true;
                else
                    return false;
            };


            LOCSEARCH::CoorDesc<double> desc(mpp, stopper);
            const int n = 12;
            double *x = new double[n];
            snowgoose::BoxUtils::getCenter(*(mpp.mBox), x);
            double v;
            v = mpp.mObjectives[0]->func(x);
            std::cout << "Initial v = " << v << "\n";
            std::cout << "Initial x = " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";
            bool rv = desc.search(x, v);
            std::cout << desc.about() << "\n";
            std::cout << "In " << cnt << " iterations found v = " << v << "\n";
            std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x, 10) << "\n";


            
            
        }


    };





COMPI::MPProblem<double>* getRosenbrock() {

        COMPI::MPProblem<double>* mpp = new COMPI::MPProblem<double>();

        
        mpp->mObjectives.push_back(new  RosenbrockFunction() );

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
    
    
    

    
    
    CrystallProblemFactory cpf12(argv[1]);
    COMPI::MPProblem<double>& mpp1 = *cpf12.get();
    
    //HJTester hjtester(mpp1);
    HJTester hjtester(*getRosenbrock());
    
    std::cout<<"================================================================================CoorDesk Test"<<endl;
    hjtester.TestCoorDesk();
    
    std::cout<<"================================================================================HJ std"<<endl;
    hjtester.TestHJ(argv);
    std::cout<<"===============================================================HJ rnd"<<endl;
    hjtester.TestRnd(argv);
    std::cout<<"===============================================================HJ linear"<<endl;
    hjtester.TestLS(argv);
    return 0;
}

