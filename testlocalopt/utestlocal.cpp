/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   utestglobopt.cpp
 * Author: alusov
 *
 * Created on February 27, 2017, 10:31 AM
 */

#include <limits.h>
#include "gtest/gtest.h"
#include "utestlocal.hpp"

 
int main(int argc, char **argv) 
{
 if(argc==1){
   std::cout << "usage: " << argv[0] << ' ' << "'path to json file'" << '\n';
   return 1;
 }
 JSONPATH = argv[1];
 ::testing::InitGoogleTest(&argc, argv);
 return RUN_ALL_TESTS();
}

