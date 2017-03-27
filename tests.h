/* 
 * File:   tests.h
 * Author: mdesana
 *
 * Created on May 9, 2014, 11:48 AM
 */

#pragma once

#include "Common.h"
#include "SPGM.h"

void TestEM1();
void TestEM2();
void TestEM3();
void TestEM4();
void TestChowLiu();
void TestMixMaxTreesSPGM(Params p);
void TestMixMaxTreesSPGM_vsChowLiu(Params p);
void TestKnapsack();
void TestRuntime(bool useNltcs, int NmaxPaths, int nTimes);
void TestRuntime2(int nVars, int nSamples);
void TestKMeans(int K);
void TestSPGM_mixture(Params p);
void TestSPGM_mixture2(Params p);
void TestSPGM_mixture3(Params p);
void TestSPGM_mixture3_rand(Params p);
void TestSPGM_mixture4(Params p);
void TestSPGM_mixture5(Params p);
void TestEigen(int nSamples, int nVals);
void TestCopyConstructor(Params p);