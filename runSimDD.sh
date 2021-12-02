#!/bin/bash
parallel -j80 "root -b -q -l 'SimulateDDstarCorrelation.cc(2500000, kCRMode2, kSoftQCD, 14000, "{}", \"AnalysisResults_vsY_"{}".root\")'" ::: {1..80}
