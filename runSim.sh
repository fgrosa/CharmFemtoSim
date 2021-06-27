#!/bin/bash
parallel -j20 "root -b -q -l 'SimulateCharmLightCorrelation.cc(1000000, kCRMode2, kSoftQCD, 13000, "{}", \"AnalysisResults_"{}".root\")'" ::: {1..20}
