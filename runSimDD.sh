#!/bin/bash
NEVENTSPERJOB=625000 #1250000
parallel -j96 "root -b -q -l 'SimulateDDstarCorrelation.cc($NEVENTSPERJOB, kCRMode2, kHardQCD, 14000, "{}", \"hardqcd/AnalysisResults_vsY_"{}".root\")'" ::: {1..96}
