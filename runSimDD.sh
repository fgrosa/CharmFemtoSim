#!/bin/bash
parallel -j20 "root -b -q -l 'SimulateDDstarCorrelation.cc(10000, kCRMode2, kSoftQCD, 14000, "{}", \"AnalysisResults_"{}".root\")'" ::: {101..120}
