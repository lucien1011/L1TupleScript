{

  gROOT->ProcessLine(".x ../../L1Ntuples/macros/initL1Analysis.C++");
  std::cout << "L1Ntuple library loading ..." <<std::endl;
  gROOT->ProcessLine(".L trigTiming.C++");
  std::cout << "L1Ntuple Running ..." <<std::endl;

  RunL1("246908");

  std::cout << " Done ..." <<std::endl;
}
