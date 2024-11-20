#include <stdio.h>
#include <iostream>
#include "Riostream.h"
#include "TBranch.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include "TProfile.h"
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <TROOT.h>
#include <TMath.h>

// Macro to make a BH flat ntuple. John Hakala 12/1/2015

void BHflatTuplizer_origin(std::string inFilename, std::string outFilename, std::string metListFilename);
float dR(float eta1, float phi1, float eta2, float phi2);
std::map<unsigned, std::set<unsigned> > readEventList(char const* _fileName);

void BHflatTuplizer_origin(std::string inFilename, std::string outFilename, std::string metListFilename) {
  std::map<unsigned, std::set<unsigned> > list = readEventList(metListFilename.c_str());
  bool isData        = true  ;
  bool debugFlag     = false  ;
  int  eventsToDump  = 25    ;  // if debugFlag is true, then stop once the number of dumped events reaches eventsToDump
  bool dumpBigEvents = true  ;
  bool dumpIsoInfo   = false ;
  int  nDumpedEvents = 0     ;
  bool useMETcut     = false ;
  int  nBin	     = 400   ; // 100 for 100 GeV bin, 1000 for 10 GeV bin in ST histograms
  int  STlow	     = 0   ;
  //int  STlow	     = 2500   ;    // Lower bound of ST histogram: 500 GeV or 0 GeV
  bool is2015D       = true;    // switches for updated trigger name

  // define output textfile
  ofstream outTextFile;
  std::string outTextFilename  = outFilename+"_log.tx";
  outTextFile.open(outTextFilename.c_str());

  // define output histograms
  TH2F METvsMHT                            = TH2F("METvsMHT"                            ,  "METvsMHT"                        ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2                        = TH2F("METvsMHTinc2"                        ,  "METvsMHTinc2"                    ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasMuon                 = TH2F("METvsMHTinc2hasMuon"                 ,  "METvsMHTinc2hasMuon"             ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasPhoton               = TH2F("METvsMHTinc2hasPhoton"               ,  "METvsMHTinc2hasPhoton"           ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasElectron             = TH2F("METvsMHTinc2hasElectron"             ,  "METvsMHTinc2hasElectron"         ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2onlyJets                = TH2F("METvsMHTinc2onlyJets"                ,  "METvsMHTinc2onlyJets"            ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHT_tight                      = TH2F("METvsMHT_tight"                      ,  "METvsMHT_tight"                  ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2_tight                  = TH2F("METvsMHTinc2_tight"                  ,  "METvsMHTinc2_tight"              ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasMuon_tight           = TH2F("METvsMHTinc2hasMuon_tight"           ,  "METvsMHTinc2hasMuon_tight"       ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasPhoton_tight         = TH2F("METvsMHTinc2hasPhoton_tight"         ,  "METvsMHTinc2hasPhoton_tight"     ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasElectron_tight       = TH2F("METvsMHTinc2hasElectron_tight"       ,  "METvsMHTinc2hasElectron_tight"   ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2onlyJets_tight          = TH2F("METvsMHTinc2onlyJets_tight"          ,  "METvsMHTinc2onlyJets_tight"      ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH1F METoverSumET                      = TH1F("METoverSumET"                      ,  "METoverSumET"                      ,  300,  0.,  3.);
  TH1F METoverSumETinc2                  = TH1F("METoverSumETinc2"                  ,  "METoverSumETinc2"                  ,  300,  0.,  3.);
  TH1F METoverSumETinc2hasMuon           = TH1F("METoverSumETinc2hasMuon"           ,  "METoverSumETinc2hasMuon"           ,  300,  0.,  3.);
  TH1F METoverSumETinc2hasPhoton         = TH1F("METoverSumETinc2hasPhoton"         ,  "METoverSumETinc2hasPhoton"         ,  300,  0.,  3.);
  TH1F METoverSumETinc2hasElectron       = TH1F("METoverSumETinc2hasElectron"       ,  "METoverSumETinc2hasElectron"       ,  300,  0.,  3.);
  TH1F METoverSumETinc2onlyJets          = TH1F("METoverSumETinc2onlyJets"          ,  "METoverSumETinc2onlyJets"          ,  300,  0.,  3.);
  TH1F METoverSumET_tight                = TH1F("METoverSumET_tight"                ,  "METoverSumET_tight"                ,  300,  0.,  3.);
  TH1F METoverSumETinc2_tight            = TH1F("METoverSumETinc2_tight"            ,  "METoverSumETinc2_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc3_tight            = TH1F("METoverSumETinc3_tight"            ,  "METoverSumETinc3_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc4_tight            = TH1F("METoverSumETinc4_tight"            ,  "METoverSumETinc4_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc5_tight            = TH1F("METoverSumETinc5_tight"            ,  "METoverSumETinc5_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc6_tight            = TH1F("METoverSumETinc6_tight"            ,  "METoverSumETinc6_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc7_tight            = TH1F("METoverSumETinc7_tight"            ,  "METoverSumETinc7_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc8_tight            = TH1F("METoverSumETinc8_tight"            ,  "METoverSumETinc8_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc9_tight            = TH1F("METoverSumETinc9_tight"            ,  "METoverSumETinc9_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc10_tight            = TH1F("METoverSumETinc10_tight"            ,  "METoverSumETinc10_tight"            ,  300,  0.,  3.);
  TH1F METoverSumETinc2hasMuon_tight     = TH1F("METoverSumETinc2hasMuon_tight"     ,  "METoverSumETinc2hasMuon_tight"     ,  300,  0.,  3.);
  TH1F METoverSumETinc2hasPhoton_tight   = TH1F("METoverSumETinc2hasPhoton_tight"   ,  "METoverSumETinc2hasPhoton_tight"   ,  300,  0.,  3.);
  TH1F METoverSumETinc2hasElectron_tight = TH1F("METoverSumETinc2hasElectron_tight" ,  "METoverSumETinc2hasElectron_tight" ,  300,  0.,  3.);
  TH1F METoverSumETinc2onlyJets_tight    = TH1F("METoverSumETinc2onlyJets_tight"    ,  "METoverSumETinc2onlyJets_tight"    ,  300,  0.,  3.);

  TH1F MuonJetIso1             = TH1F("MuonJetIso1"             ,  "MuonJetIso1"             ,  300,   0.,  3);
  TH1F MuonJetIso2             = TH1F("MuonJetIso2"             ,  "MuonJetIso2"             ,  300,   0.,  3);
  TH1F MuonJetIso3             = TH1F("MuonJetIso3"             ,  "MuonJetIso3"             ,  300,   0.,  3);
  TH1F MuonJetIso4             = TH1F("MuonJetIso4"             ,  "MuonJetIso4"             ,  300,   0.,  3);
  TH1F MuonJetoverlapdR1       = TH1F("MuonJetIso4overlapdR1"   ,  "MuonJetIso4overlapdR1"   ,  300,   0., .3);
  TH1F MuonJetoverlapdR2       = TH1F("MuonJetIso4overlapdR2"   ,  "MuonJetIso4overlapdR2"   ,  300,   0., .3);
  TH1F MuonJetoverlapdR3       = TH1F("MuonJetIso4overlapdR3"   ,  "MuonJetIso4overlapdR3"   ,  300,   0., .3);
  TH1F MuonJetoverlapdR4       = TH1F("MuonJetIso4overlapdR4"   ,  "MuonJetIso4overlapdR4"   ,  300,   0., .3);
  TH1F ElectronJetIso1         = TH1F("ElectronJetIso1"         ,  "ElectronJetIso1"         ,  300,   0.,  3);
  TH1F ElectronJetIso2         = TH1F("ElectronJetIso2"         ,  "ElectronJetIso2"         ,  300,   0.,  3);
  TH1F ElectronJetIso3         = TH1F("ElectronJetIso3"         ,  "ElectronJetIso3"         ,  300,   0.,  3);
  TH1F ElectronJetIso4         = TH1F("ElectronJetIso4"         ,  "ElectronJetIso4"         ,  300,   0.,  3);
  TH1F ElectronJetoverlapdR1   = TH1F("ElectronJetoverlapdR1"   ,  "ElectronJetoverlapdR1"   ,  300,   0., .3);
  TH1F ElectronJetoverlapdR2   = TH1F("ElectronJetoverlapdR2"   ,  "ElectronJetoverlapdR2"   ,  300,   0., .3);
  TH1F ElectronJetoverlapdR3   = TH1F("ElectronJetoverlapdR3"   ,  "ElectronJetoverlapdR3"   ,  300,   0., .3);
  TH1F ElectronJetoverlapdR4   = TH1F("ElectronJetoverlapdR4"   ,  "ElectronJetoverlapdR4"   ,  300,   0., .3);
  TH1F PhotonJetIso1           = TH1F("PhotonJetIso1"           ,  "PhotonJetIso1"           ,  300,   0.,  3);
  TH1F PhotonJetIso2           = TH1F("PhotonJetIso2"           ,  "PhotonJetIso2"           ,  300,   0.,  3);
  TH1F PhotonJetIso3           = TH1F("PhotonJetIso3"           ,  "PhotonJetIso3"           ,  300,   0.,  3);
  TH1F PhotonJetIso4           = TH1F("PhotonJetIso4"           ,  "PhotonJetIso4"           ,  300,   0.,  3);
  TH1F PhotonJetoverlapdR1     = TH1F("PhotonJetoverlapdR1"     ,  "PhotonJetoverlapdR1"     ,  300,   0., .3);
  TH1F PhotonJetoverlapdR2     = TH1F("PhotonJetoverlapdR2"     ,  "PhotonJetoverlapdR2"     ,  300,   0., .3);
  TH1F PhotonJetoverlapdR3     = TH1F("PhotonJetoverlapdR3"     ,  "PhotonJetoverlapdR3"     ,  300,   0., .3);
  TH1F PhotonJetoverlapdR4     = TH1F("PhotonJetoverlapdR4"     ,  "PhotonJetoverlapdR4"     ,  300,   0., .3);

  // loop to create ST histograms for inclusive and exclusive multiplicities from 2 up to multMax
  TH1F stHist = TH1F("stHist", "ST", nBin, STlow, 20000);
  TH1F stHist_tight = TH1F("stHist_tight", "ST_tight", nBin, STlow, 20000);
  int mult=2;
  int multMax = 12;
  TH1F *stIncHist[multMax-2];
  TH1F *stExcHist[multMax-2];
  TH1F *stExcHist_2leadjets[multMax-2];
  TH1F *stExcHist_noISR[multMax-2];
  TH1F *Gen71stExcHist[multMax-2];
  TH1F *GenInteractstExcHist[multMax-2];
  TH1F *GenUnderlystExcHist[multMax-2];
  TH1F *GenJetstExcHist[multMax-2];
  TH1F *GenJetNoISRstExcHist[multMax-2];
  TH1F *stIncHist_tight[multMax-2];
  TH1F *stExcHist_tight[multMax-2];
  TH1F stHistMHT = TH1F("stHistMHT", "ST using MHT", nBin, STlow, 20000);
  TH1F stHistMHT_tight = TH1F("stHistMHT_tight", "ST_tight using MHT_tight", nBin, STlow, 20000);
  TH1F *stIncHistMHT[multMax-2];
  TH1F *stExcHistMHT[multMax-2];
  TH1F *stIncHistMHT_tight[multMax-2];
  TH1F *stExcHistMHT_tight[multMax-2];
  TH1F JetEtaExc08_3500Up = TH1F("JetEtaExc08_3500Up","JetEtaExc08_3500Up", 100, -5, 5);
  TH1F JetPtExc08_3500Up = TH1F("JetPtExc08_3500Up","JetPtExc08_3500Up", 130, 0, 1300);
  TH1F JetEtaExc08_2300to2400 = TH1F("JetEtaExc08_2300to2400","JetEtaExc08_2300to2400", 100, -5, 5);
  TH1F JetPtExc08_2300to2400 = TH1F("JetPtExc08_2300to2400","JetPtExc08_2300to2400", 130, 0, 1300);
  TH1F JetEtaExc03_1800to1900 = TH1F("JetEtaExc03_1800to1900","JetEtaExc03_1800to1900", 100, -5, 5);
  TH1F JetPtExc03_1800to1900 = TH1F("JetPtExc03_1800to1900","JetPtExc03_1800to1900", 100, 0, 2000);
  TH1F FinalGenPDGIDExc08_3500Up = TH1F("FinalGenPDGIDExc08_3500Up","FinalGenPDGIDExc08_3500Up", 20000, -10000, -10000);
  TH1F FinalGenPDGIDExc08_2300to2400 = TH1F("FinalGenPDGIDExc08_2300to2400","FinalGenPDGIDExc08_2300to2400", 20000, -10000, -10000);
  TH1F FinalGenPDGIDExc03_1800to1900 = TH1F("FinalGenPDGIDExc03_1800to1900","FinalGenPDGIDExc03_1800to1900", 20000000, -10000000, -10000000);
  TH1F FinalGenEtExc08_3500Up = TH1F("FinalGenEtExc08_3500Up","FinalGenEtExc08_3500Up", 100, 0, 500);
  TH1F FinalGenEtExc08_2300to2400 = TH1F("FinalGenEtExc08_2300to2400","FinalGenEtExc08_2300to2400", 100, 0, 500);
  TH1F FinalGenEtExc03_1800to1900 = TH1F("FinalGenEtExc03_1800to1900","FinalGenEtExc03_1800to1900", 100, 0, 2000);
  TH1F FinalGenEtaExc08_3500Up = TH1F("FinalGenEtaExc08_3500Up","FinalGenEtaExc08_3500Up", 100, -5, 5);
  TH1F FinalGenEtaExc08_2300to2400 = TH1F("FinalGenEtaExc08_2300to2400","FinalGenEtaExc08_2300to2400", 100, -5, 5);
  TH1F FinalGenEtaExc03_1800to1900 = TH1F("FinalGenEtaExc03_1800to1900","FinalGenEtaExc03_1800to1900", 100, -5, 5);
  TH1F FinalGenNExc08_3500Up = TH1F("FinalGenNExc08_3500Up","FinalGenNExc08_3500Up", 100, 0, 100);
  TH1F FinalGenNExc08_2300to2400 = TH1F("FinalGenNExc08_2300to2400","FinalGenNExc08_2300to2400", 100, 0, 100);
  TH1F FinalGenNExc03_1800to1900 = TH1F("FinalGenNExc03_1800to1900","FinalGenNExc03_1800to1900", 100, 0, 100);
  TH1F FinalGenSTExc03 = TH1F("FinalGenSTExc03","FinalGenSTExc03", 100, 0, 2000);
  TH1F FinalGenSTExc08 = TH1F("FinalGenSTExc08","FinalGenSTExc08", 100, 0, 2000);
  TH1F FinalGenST71Exc03 = TH1F("FinalGenST71Exc03","FinalGenST71Exc03", 40, 0, 4000);
  TH1F FinalGenST71Exc02 = TH1F("FinalGenST71Exc02","FinalGenST71Exc02", 40, 0, 4000);
  TH1F FinalGenST71Exc08 = TH1F("FinalGenST71Exc08","FinalGenST71Exc08", 40, 0, 4000);
  //TH1F from1500to3000_ST_N2 = TH1F("from1500to3000_ST_N2","from1500to3000_ST_N2", 40, 0, 4000);
  //TH1F from1500to3000_ST_N8 = TH1F("from1500to3000_ST_N8","from1500to3000_ST_N8", 40, 0, 4000);
  TH1F FinalGenN71Exc02 = TH1F("FinalGenN71Exc02","FinalGenN71Exc02", 200, 0, 200);
  TH1F FinalGenN71Exc08 = TH1F("FinalGenN71Exc08","FinalGenN71Exc08", 200, 0, 200);
  TH1F FinalGenPDGID71Exc02 = TH1F("FinalGenPDGID71Exc02","FinalGenPDGID71Exc02", 80, -40, 40);
  TH1F FinalGenPDGID71Exc08 = TH1F("FinalGenPDGID71Exc08","FinalGenPDGID71Exc08", 80, -40, 40);
  TH1F FinalGenEta71Exc02 = TH1F("FinalGenEta71Exc02","FinalGenEta71Exc02", 100, -5, 5);
  TH1F FinalGenEta71Exc08 = TH1F("FinalGenEta71Exc08","FinalGenEta71Exc08", 100, -5, 5);
  TH1F FinalGenEt71Exc02 = TH1F("FinalGenEt71Exc02","FinalGenEt71Exc02", 40, 0, 400);
  TH1F FinalGenEt71Exc08 = TH1F("FinalGenEt71Exc08","FinalGenEt71Exc08", 40, 0, 400);
  TH1F FinalGenPDGID45Exc02 = TH1F("FinalGenPDGID45Exc02","FinalGenPDGID45Exc02", 80, -40, 40);
  TH1F FinalGenPDGID45Exc08 = TH1F("FinalGenPDGID45Exc08","FinalGenPDGID45Exc08", 80, -40, 40);
  TH1F FinalGenEta45Exc02 = TH1F("FinalGenEta45Exc02","FinalGenEta45Exc02", 100, -5, 5);
  TH1F FinalGenEta45Exc08 = TH1F("FinalGenEta45Exc08","FinalGenEta45Exc08", 100, -5, 5);
  TH1F FinalGenEt45Exc02 = TH1F("FinalGenEt45Exc02","FinalGenEt45Exc02", 100, 500, 1500);
  TH1F FinalGenEt45Exc08 = TH1F("FinalGenEt45Exc08","FinalGenEt45Exc08", 100, 500, 1500);
  TH1F FinalGenSTExc03_interaction = TH1F("FinalGenSTExc03_interaction","FinalGenSTExc03_interaction", 40, 0, 4000);
  TH1F FinalGenSTExc02_interaction = TH1F("FinalGenSTExc02_interaction","FinalGenSTExc02_interaction", 40, 0, 4000);
  TH1F FinalGenSTExc08_interaction = TH1F("FinalGenSTExc08_interaction","FinalGenSTExc08_interaction", 40, 0, 4000);
  TH1F FinalGenPDGIDExc02_interaction = TH1F("FinalGenPDGIDExc02_interaction","FinalGenPDGIDExc02_interaction", 80, -40, 40);
  TH1F FinalGenPDGIDExc08_interaction = TH1F("FinalGenPDGIDExc08_interaction","FinalGenPDGIDExc08_interaction", 80, -40, 40);
  TH1F FinalGenEtaExc02_interaction = TH1F("FinalGenEtaExc02_interaction","FinalGenEtaExc02_interaction", 100, -5, 5);
  TH1F FinalGenEtaExc08_interaction = TH1F("FinalGenEtaExc08_interaction","FinalGenEtaExc08_interaction", 100, -5, 5);
  TH1F FinalGenEtExc02_interaction = TH1F("FinalGenEtExc02_interaction","FinalGenEtExc02_interaction", 40, 0, 400);
  TH1F FinalGenEtExc08_interaction = TH1F("FinalGenEtExc08_interaction","FinalGenEtExc08_interaction", 40, 0, 400);
  TH1F FinalGenSTExc03_radiation = TH1F("FinalGenSTExc03_radiation","FinalGenSTExc03_radiation", 40, 0, 4000);
  TH1F FinalGenSTExc02_radiation = TH1F("FinalGenSTExc02_radiation","FinalGenSTExc02_radiation", 40, 0, 4000);
  TH1F FinalGenSTExc08_radiation = TH1F("FinalGenSTExc08_radiation","FinalGenSTExc08_radiation", 40, 0, 4000);
  TH1F FinalGenPDGIDExc02_radiation = TH1F("FinalGenPDGIDExc02_radiation","FinalGenPDGIDExc02_radiation", 80, -40, 40);
  TH1F FinalGenPDGIDExc08_radiation = TH1F("FinalGenPDGIDExc08_radiation","FinalGenPDGIDExc08_radiation", 80, -40, 40);
  TH1F FinalGenEtaExc02_radiation = TH1F("FinalGenEtaExc02_radiation","FinalGenEtaExc02_radiation", 100, -5, 5);
  TH1F FinalGenEtaExc08_radiation = TH1F("FinalGenEtaExc08_radiation","FinalGenEtaExc08_radiation", 100, -5, 5);
  TH1F FinalGenEtExc02_radiation = TH1F("FinalGenEtExc02_radiation","FinalGenEtExc02_radiation", 40, 0, 400);
  TH1F FinalGenEtExc08_radiation = TH1F("FinalGenEtExc08_radiation","FinalGenEtExc08_radiation", 40, 0, 400);
  TH1F Energy_Eta_radiation_N2 = TH1F("Energy_Eta_radiation_N2","Energy_Eta_radiation_N2", 100, -5, 5);
  TH1F Energy_Eta_radiation_N8 = TH1F("Energy_Eta_radiation_N8","Energy_Eta_radiation_N8", 100, -5, 5);
  TH1F Energy_Eta_interaction_N2 = TH1F("Energy_Eta_interaction_N2","Energy_Eta_interaction_N2", 100, -5, 5);
  TH1F Energy_Eta_interaction_N8 = TH1F("Energy_Eta_interaction_N8","Energy_Eta_interaction_N8", 100, -5, 5);
  TH1F ST_Sphericity_p05top2_in = TH1F("ST_Sphericity_p05top2_in","ST_Sphericity_p05top2_in", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out = TH1F("ST_Sphericity_p05top2_out","ST_Sphericity_p05top2_out", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N3 = TH1F("ST_Sphericity_p05top2_in_N3","ST_Sphericity_p05top2_in_N3", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N3 = TH1F("ST_Sphericity_p05top2_out_N3","ST_Sphericity_p05top2_out_N3", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N4 = TH1F("ST_Sphericity_p05top2_in_N4","ST_Sphericity_p05top2_in_N4", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N4 = TH1F("ST_Sphericity_p05top2_out_N4","ST_Sphericity_p05top2_out_N4", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N5 = TH1F("ST_Sphericity_p05top2_in_N5","ST_Sphericity_p05top2_in_N5", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N5 = TH1F("ST_Sphericity_p05top2_out_N5","ST_Sphericity_p05top2_out_N5", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N6 = TH1F("ST_Sphericity_p05top2_in_N6","ST_Sphericity_p05top2_in_N6", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N6 = TH1F("ST_Sphericity_p05top2_out_N6","ST_Sphericity_p05top2_out_N6", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N7 = TH1F("ST_Sphericity_p05top2_in_N7","ST_Sphericity_p05top2_in_N7", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N7 = TH1F("ST_Sphericity_p05top2_out_N7","ST_Sphericity_p05top2_out_N7", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N8 = TH1F("ST_Sphericity_p05top2_in_N8","ST_Sphericity_p05top2_in_N8", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N8 = TH1F("ST_Sphericity_p05top2_out_N8","ST_Sphericity_p05top2_out_N8", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N9 = TH1F("ST_Sphericity_p05top2_in_N9","ST_Sphericity_p05top2_in_N9", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N9 = TH1F("ST_Sphericity_p05top2_out_N9","ST_Sphericity_p05top2_out_N9", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N10 = TH1F("ST_Sphericity_p05top2_in_N10","ST_Sphericity_p05top2_in_N10", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N10 = TH1F("ST_Sphericity_p05top2_out_N10","ST_Sphericity_p05top2_out_N10", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_in_N11 = TH1F("ST_Sphericity_p05top2_in_N11","ST_Sphericity_p05top2_in_N11", 80, 0, 4000);
  TH1F ST_Sphericity_p05top2_out_N11 = TH1F("ST_Sphericity_p05top2_out_N11","ST_Sphericity_p05top2_out_N11", 80, 0, 4000);
  TH1F Hist_sphericity = TH1F("Hist_sphericity","Hist_sphericity", 100, 0, 1);
  TH1F Hist_sphericity_ST4p5cut = TH1F("Hist_sphericity_ST4p5cut","Hist_sphericity_ST4p5cut", 100, 0, 1);
  TH1F Hist_sphericity_ST4p5cut_incN3 = TH1F("Hist_sphericity_ST4p5cut_incN3","Hist_sphericity_ST4p5cut_incN3", 100, 0, 1);
  TH1F Hist_sphericity_ST4p5cut_incN4 = TH1F("Hist_sphericity_ST4p5cut_incN4","Hist_sphericity_ST4p5cut_incN4", 100, 0, 1);
  TH1F Hist_sphericity_ST4p5cut_incN5 = TH1F("Hist_sphericity_ST4p5cut_incN5","Hist_sphericity_ST4p5cut_incN5", 100, 0, 1);
  TH1F Hist_sphericity_ST4p5cut_incN6 = TH1F("Hist_sphericity_ST4p5cut_incN6","Hist_sphericity_ST4p5cut_incN6", 100, 0, 1);
  TH1F Hist_sphericity_ST4p5cut_incN7 = TH1F("Hist_sphericity_ST4p5cut_incN7","Hist_sphericity_ST4p5cut_incN7", 100, 0, 1);
  TH1F Hist_sphericity_ST4p5cut_incN8 = TH1F("Hist_sphericity_ST4p5cut_incN8","Hist_sphericity_ST4p5cut_incN8", 100, 0, 1);
  TH1F ST_Sphericity_p1less = TH1F("ST_Sphericity_p1less","ST_Sphericity_p1less", 80, 2000, 10000);
  TH1F ST_Sphericity_p1less_N3 = TH1F("ST_Sphericity_p1less_N3","ST_Sphericity_p1less_N3", 80, 2000, 10000);
  TH1F ST_Sphericity_p1less_N4 = TH1F("ST_Sphericity_p1less_N4","ST_Sphericity_p1less_N4", 80, 2000, 10000);
  TH1F ST_Sphericity_p1less_N5 = TH1F("ST_Sphericity_p1less_N5","ST_Sphericity_p1less_N5", 80, 2000, 10000);
  TH1F ST_Sphericity_p1less_N6 = TH1F("ST_Sphericity_p1less_N6","ST_Sphericity_p1less_N6", 80, 2000, 10000);
  TH1F ST_Sphericity_p1less_N7 = TH1F("ST_Sphericity_p1less_N7","ST_Sphericity_p1less_N7", 80, 2000, 10000);
  TH1F ST_Sphericity_p1less_N8 = TH1F("ST_Sphericity_p1less_N8","ST_Sphericity_p1less_N8", 80, 2000, 10000);
  TH1F ST_Sphericity_p1less_N9 = TH1F("ST_Sphericity_p1less_N9","ST_Sphericity_p1less_N9", 80, 2000, 10000);
  TH1F ST_Sphericity_p1more = TH1F("ST_Sphericity_p1more","ST_Sphericity_p1more", 80, 2000, 10000);
  TH1F ST_Sphericity_p1more_N3 = TH1F("ST_Sphericity_p1more_N3","ST_Sphericity_p1more_N3", 80, 2000, 10000);
  TH1F ST_Sphericity_p1more_N4 = TH1F("ST_Sphericity_p1more_N4","ST_Sphericity_p1more_N4", 80, 2000, 10000);
  TH1F ST_Sphericity_p1more_N5 = TH1F("ST_Sphericity_p1more_N5","ST_Sphericity_p1more_N5", 80, 2000, 10000);
  TH1F ST_Sphericity_p1more_N6 = TH1F("ST_Sphericity_p1more_N6","ST_Sphericity_p1more_N6", 80, 2000, 10000);
  TH1F ST_Sphericity_p1more_N7 = TH1F("ST_Sphericity_p1more_N7","ST_Sphericity_p1more_N7", 80, 2000, 10000);
  TH1F ST_Sphericity_p1more_N8 = TH1F("ST_Sphericity_p1more_N8","ST_Sphericity_p1more_N8", 80, 2000, 10000);
  TH1F ST_Sphericity_p1more_N9 = TH1F("ST_Sphericity_p1more_N9","ST_Sphericity_p1more_N9", 80, 2000, 10000);

  char *histTitle = new char[20];
  // These use pat::slimmedMETs
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    sprintf(histTitle, "stInc%02dHist", mult);
    stIncHist[iHist] = new TH1F(histTitle, "Inclusive ST", nBin, STlow, 20000);
    sprintf(histTitle, "stExc%02dHist", mult);
    stExcHist[iHist] = new TH1F(histTitle, "Exclusive ST", nBin, STlow, 20000);
    sprintf(histTitle, "stExc%02dHist_2leadjets", mult);
    stExcHist_2leadjets[iHist] = new TH1F(histTitle, "Exclusive ST with 2 leading jets", nBin, STlow, 20000);
    sprintf(histTitle, "stExc%02dHist_noISR", mult);
    stExcHist_noISR[iHist] = new TH1F(histTitle, "Exclusive ST with no ISR events", nBin, STlow, 20000);
    sprintf(histTitle, "Gen71stExc%02dHist", mult);
    Gen71stExcHist[iHist] = new TH1F(histTitle, "Exclusive Gen71 ST", nBin, STlow, 20000);
    sprintf(histTitle, "GenInteractstExc%02dHist", mult);
    GenInteractstExcHist[iHist] = new TH1F(histTitle, "Exclusive Gen Interact ST", nBin, STlow, 20000);
    sprintf(histTitle, "GenUnderlystExc%02dHist", mult);
    GenUnderlystExcHist[iHist] = new TH1F(histTitle, "Exclusive Gen Underly ST", nBin, STlow, 20000);
    sprintf(histTitle, "GenJetstExc%02dHist", mult);
    GenJetstExcHist[iHist] = new TH1F(histTitle, "Exclusive Gen Jet ST", nBin, STlow, 20000);
    sprintf(histTitle, "GenJetNoISRstExc%02dHist", mult);
    GenJetNoISRstExcHist[iHist] = new TH1F(histTitle, "Exclusive Gen Jet No ISR ST", nBin, STlow, 20000);
    sprintf(histTitle, "stInc%02dHist_tight", mult);
    stIncHist_tight[iHist] = new TH1F(histTitle, "Inclusive ST_tight", nBin, STlow, 20000);
    sprintf(histTitle, "stExc%02dHist_tight", mult);
    stExcHist_tight[iHist] = new TH1F(histTitle, "Exclusive ST_tight", nBin, STlow, 20000);
    sprintf(histTitle, "stInc%02dHistMHT", mult);
    stIncHistMHT[iHist] = new TH1F(histTitle, "Inclusive ST using MHT", nBin, STlow, 20000);
    sprintf(histTitle, "stExc%02dHistMHT", mult);
    stExcHistMHT[iHist] = new TH1F(histTitle, "Exclusive ST using MHT", nBin, STlow, 20000);
    sprintf(histTitle, "stInc%02dHistMHT_tight", mult);
    stIncHistMHT_tight[iHist] = new TH1F(histTitle, "Inclusive ST_tight using MHT_tight", nBin, STlow, 20000);
    sprintf(histTitle, "stExc%02dHistMHT_tight", mult);
    stExcHistMHT_tight[iHist] = new TH1F(histTitle, "Exclusive ST_tight using MHT_tight", nBin, STlow, 20000);
    ++mult;
  }

  TH2D *Hist2D_N_ST = new TH2D("Hist2D_N_ST","ST to N", 25, 0, 25, 35, 3000, 10000);

  //TH1F sameGenST_GenPDGIDN2 = TH1F("SameGenST_GenPDGIDN2","SameGenST_GenPDGIDN2", 20000, -10000, 10000);
  //TH1F SameGenST_GenETN2 = TH1F("SameGenST_GenETN2","SameGenST_GenETN2", 100, 0, 5000);
  //TH1F SameGenST_GenEtaN2 = TH1F("SameGenST_GenEtaN2","SameGenST_GenEtaN2", 100, -5, 5);
  //TH1F SameGenST_GenPhiN2 = TH1F("SameGenST_GenPhiN2","SameGenST_GenPhiN2", 100, -5, 5);
  //TH1F SameGenST_GenPDGIDN8 = TH1F("SameGenST_GenPDGIDN8","SameGenST_GenPDGIDN8", 20000, -10000, 10000);
  //TH1F SameGenST_GenETN8 = TH1F("SameGenST_GenETN8","SameGenST_GenETN8", 100, 0, 5000);
  //TH1F SameGenST_GenEtaN8 = TH1F("SameGenST_GenEtaN8","SameGenST_GenEtaN8", 100, -5, 5);
  //TH1F SameGenST_GenPhiN8 = TH1F("SameGenST_GenPhiN8","SameGenST_GenPhiN8", 100, -5, 5);


  // variables calculated in the loop
  float OurMet           = 0.            ;
  float Px               = 0.            ;
  float Py               = 0.            ;
  float ST               = 0.            ;
  float ST_2leadjets     = 0.            ;
  float STMHTnoMET       = 0.            ;
  int multiplicity       = 0             ;
  bool passIso           = true          ;
  float OurMet_tight     = 0.            ;
  float Px_tight         = 0.            ;
  float Py_tight         = 0.            ;
  float ST_tight         = 0.            ;
  float STMHTnoMET_tight = 0.            ;
  int multiplicity_tight = 0             ;
  bool passIso_tight     = true          ;
  bool passMetCut        = true          ;
  bool passMetCut_tight  = true          ;
  char *messageBuffer    = new char[400] ;
  bool eventHasMuon      = false         ;
  bool eventHasPhoton    = false         ;
  bool eventHasElectron  = false         ;
  bool  TightJets[25]                    ;
  bool  isTightJet       = false         ;
  float JetMuonEt        = 0.            ;
  float JetElectronEt    = 0.            ;
  float JetPhotonEt      = 0.            ;

  // variables accessed from the tree
  //TODO
  Bool_t     firedHLT_PFHT800       ;
  Bool_t     firedHLT_PFHT900       ;
  //Bool_t     firedHLT_PFHT800_v1       ;
  //Bool_t     passed_CSCTightHaloFilter ;
  Bool_t     firedHLT_CaloJet500_NoJetID          ;
  Bool_t     firedHLT_AK8PFJet450      ;
  Bool_t     firedHLT_PFJet450         ;
  Bool_t     passed_globalTightHalo2016Filter ;
  Bool_t     passed_goodVertices       ;
  Bool_t     passed_eeBadScFilter      ;
  Bool_t     passed_HBHENoiseFilter      ;
  Bool_t     passed_HBHENoiseIsoFilter      ;
  Bool_t     passed_filterbadPFMuon      ;
  //Bool_t     passed_BadPFMuonFilter     ; //for UL samples
  Bool_t     passed_filterbadChCandidate      ;
  //Bool_t     passed_BadChargedCandidateFilter    ;//for UL samples
  Bool_t     passed_EcalDeadCellTriggerPrimitiveFilter      ;
  //Bool_t     passed_duplicateMuons      ;
  //Bool_t     passed_globalSuperTightHalo2016Filter      ;
  int        runno                     ;
  long long  evtno                     ;
  int        lumiblock                 ;
  //int        Gen_Multiplicity          ;
  float      JetEt [25]                ;
  float      JetPx [25]                ;
  float      JetPy [25]                ;
  float      JetEta[25]                ;
  float      JetPhi[25]                ;
  float      JetNeutHadFrac[25]        ;
  float      JetNeutEMFrac[25]         ;
  float      JetChgHadFrac[25]         ;
  float      JetMuFrac[25]             ;
  float      JetChgEMFrac[25]          ;
  int        JetNConstituents[25]      ;
  int        JetNNeutConstituents[25]  ;
  int        JetNChgConstituents[25]   ;
  float      EleEt[25]                 ;
  float      ElePx[25]                 ;
  float      ElePy[25]                 ;
  float      EleEta[25]                ;
  float      ElePhi[25]                ;
  float      PhEt[25]                  ;
  float      PhPx[25]                  ;
  float      PhPy[25]                  ;
  float      PhEta[25]                 ;
  float      PhPhi[25]                 ;
  float      MuEt[25]                  ;
  float      MuPx[25]                  ;
  float      MuPy[25]                  ;
  float      MuEta[25]                 ;
  float      MuPhi[25]                 ;
  float      MuPFdBiso[25]             ;
  //int        FinalGenPDGID[100]        ;
  //float      FinalGenEt[100]           ;
  //float      FinalGenPt[100]           ;
  //float      FinalGenEta[100]          ;
  //float      FinalGenPhi[100]          ;
  float      Met                       ;
  float      Sphericity                ;
  //float      FinalGenST;
  //float      GenST71;
  //float      HT                       ;

  int PDGID_Gen4;
  int PDGID_Gen5;
  int N_Gen71;
  int PDGID_Gen71[150];
  float Eta_Gen4;
  float Eta_Gen5;
  float Et_Gen4;
  float Et_Gen5;
  float GenST71;
  float Eta_Gen71[150];
  float Et_Gen71[150];

  float GenST_interaction;
  int PDGID_interaction[100];
  float Eta_interaction[100];
  float Et_interaction[100];

  float GenST_radiation;
  int PDGID_radiation[100];
  float Eta_radiation[100];
  float Et_radiation[100];

  float ST_GenJet;
  float ST_GenJet_noISR;
  int count2jets;


  // tree branches
  //TODO
  TBranch  *b_firedHLT_PFHT800       ;
  TBranch  *b_firedHLT_PFHT900       ;
  //TBranch  *b_firedHLT_PFHT800       ;
  //TBranch  *b_passed_CSCTightHaloFilter ;
  TBranch  *b_firedHLT_CaloJet500_NoJetID       ;
  TBranch  *b_firedHLT_AK8PFJet450   ;
  TBranch  *b_firedHLT_PFJet450      ;
  TBranch  *b_passed_globalTightHalo2016Filter ;
  TBranch  *b_passed_goodVertices       ;
  TBranch  *b_passed_eeBadScFilter      ;
  TBranch  *b_passed_HBHENoiseFilter      ;
  TBranch  *b_passed_HBHENoiseIsoFilter      ;
  TBranch  *b_passed_filterbadPFMuon      ;
  //TBranch  *b_passed_BadPFMuonFilter     ;
  TBranch  *b_passed_filterbadChCandidate      ;
  //TBranch  *b_passed_BadChargedCandidateFilter    ;
  TBranch  *b_passed_EcalDeadCellTriggerPrimitiveFilter      ;
  //TBranch  *b_passed_duplicateMuons      ;
  //TBranch  *b_passed_globalSuperTightHalo2016Filter      ;
  //TBranch  *b_HT                     ;
  TBranch  *b_JetEt                     ;
  TBranch  *b_JetPx                     ;
  TBranch  *b_JetPy                     ;
  TBranch  *b_JetEta                    ;
  TBranch  *b_JetPhi                    ;
  TBranch  *b_JetNeutHadFrac            ;
  TBranch  *b_JetNeutEMFrac             ;
  TBranch  *b_JetChgHadFrac             ;
  TBranch  *b_JetMuFrac                 ;
  TBranch  *b_JetChgEMFrac              ;
  TBranch  *b_JetNConstituents          ;
  TBranch  *b_JetNNeutConstituents      ;
  TBranch  *b_JetNChgConstituents       ;
  TBranch  *b_EleEt                     ;
  TBranch  *b_ElePx                     ;
  TBranch  *b_ElePy                     ;
  TBranch  *b_EleEta                    ;
  TBranch  *b_ElePhi                    ;
  TBranch  *b_PhEt                      ;
  TBranch  *b_PhPx                      ;
  TBranch  *b_PhPy                      ;
  TBranch  *b_PhEta                     ;
  TBranch  *b_PhPhi                     ;
  TBranch  *b_MuEt                      ;
  TBranch  *b_MuPx                      ;
  TBranch  *b_MuPy                      ;
  TBranch  *b_MuEta                     ;
  TBranch  *b_MuPhi                     ;
  TBranch  *b_MuPFdBiso                 ;
  TBranch  *b_Met                       ;
  TBranch  *b_runno                     ;
  TBranch  *b_evtno                     ;
  TBranch  *b_lumiblock                 ;
  TBranch  *b_Sphericity                ;
  //TBranch  *b_Gen_Multiplicity          ;
  //TBranch  *b_FinalGenST                ;
  //TBranch  *b_GenST71              ;
  //TBranch  *b_FinalGenPDGID             ;
  //TBranch  *b_FinalGenEt                ;
  //TBranch  *b_FinalGenPt                ;
  //TBranch  *b_FinalGenEta               ;
  //TBranch  *b_FinalGenPhi               ;

  TBranch  *b_PDGID_Gen4;
  TBranch  *b_PDGID_Gen5;
  TBranch  *b_N_Gen71;
  TBranch  *b_PDGID_Gen71;
  TBranch  *b_Eta_Gen4;
  TBranch  *b_Eta_Gen5;
  TBranch  *b_Et_Gen4;
  TBranch  *b_Et_Gen5;
  TBranch  *b_GenST71;
  TBranch  *b_Eta_Gen71;
  TBranch  *b_Et_Gen71;
  TBranch  *b_PDGID_interaction;
  TBranch  *b_GenST_interaction;
  TBranch  *b_Eta_interaction;
  TBranch  *b_Et_interaction;
  TBranch  *b_PDGID_radiation;
  TBranch  *b_GenST_radiation;
  TBranch  *b_Eta_radiation;
  TBranch  *b_Et_radiation;
  TBranch  *b_ST_GenJet;
  TBranch  *b_ST_GenJet_noISR;


  //create a chain by looping over the input filename
  //TChain chain("demo/t");
  TChain chain("bhana/t");
  //ifstream infile;
  //infile.open(inFilename.c_str());
  //std::string buffer;
  //const char *eosURL = "root://eoscms.cern.ch/";
  //chain.SetMakeClass(1);
  //while (std::getline(infile, buffer)) {
  //  std::string ntupleURL = eosURL + buffer;
  chain.Add(inFilename.c_str());
  //}

  cout << "Opened chain: " << chain.GetName() << endl;

  // set all branch addresses
  // TODO
  if(isData && is2015D){
      chain.SetBranchAddress( "firedHLT_PFHT800"       ,  &firedHLT_PFHT800       ,  &b_firedHLT_PFHT800       );
  } else{
      chain.SetBranchAddress( "firedHLT_PFHT800"       ,  &firedHLT_PFHT800       ,  &b_firedHLT_PFHT800       );
  }
  chain.SetBranchAddress( "firedHLT_PFHT900"       ,  &firedHLT_PFHT900       ,  &b_firedHLT_PFHT900       );
  //chain.SetBranchAddress( "passed_duplicateMuons"       ,  &passed_duplicateMuons       ,  &b_passed_duplicateMuons       );
  //chain.SetBranchAddress( "passed_globalSuperTightHalo2016Filter"       ,  &passed_globalSuperTightHalo2016Filter       ,  &b_passed_globalSuperTightHalo2016Filter       );
  //chain.SetBranchAddress( "passed_CSCTightHaloFilter" ,  &passed_CSCTightHaloFilter ,  &b_passed_CSCTightHaloFilter );
  chain.SetBranchAddress( "passed_globalTightHalo2016Filter" ,  &passed_globalTightHalo2016Filter ,  &b_passed_globalTightHalo2016Filter );
  chain.SetBranchAddress( "passed_goodVertices"       ,  &passed_goodVertices       ,  &b_passed_goodVertices       );
  chain.SetBranchAddress( "passed_eeBadScFilter"      ,  &passed_eeBadScFilter      ,  &b_passed_eeBadScFilter      );
  chain.SetBranchAddress( "passed_HBHENoiseFilter"      ,  &passed_HBHENoiseFilter      ,  &b_passed_HBHENoiseFilter      );
  chain.SetBranchAddress( "passed_HBHENoiseIsoFilter"      ,  &passed_HBHENoiseIsoFilter      ,  &b_passed_HBHENoiseIsoFilter      );
  chain.SetBranchAddress( "passed_filterbadPFMuon"      ,  &passed_filterbadPFMuon      ,  &b_passed_filterbadPFMuon      );
  //chain.SetBranchAddress( "passed_BadPFMuonFilter"      ,  &passed_BadPFMuonFilter      ,  &b_passed_BadPFMuonFilter      );
  chain.SetBranchAddress( "passed_filterbadChCandidate"      ,  &passed_filterbadChCandidate      ,  &b_passed_filterbadChCandidate      );
  //chain.SetBranchAddress( "passed_BadChargedCandidateFilter"      ,  &passed_BadChargedCandidateFilter      ,  &b_passed_BadChargedCandidateFilter      );
  chain.SetBranchAddress( "passed_EcalDeadCellTriggerPrimitiveFilter"      ,  &passed_EcalDeadCellTriggerPrimitiveFilter      ,  &b_passed_EcalDeadCellTriggerPrimitiveFilter      );
  chain.SetBranchAddress( "firedHLT_CaloJet500_NoJetID"       ,  &firedHLT_CaloJet500_NoJetID       ,  &b_firedHLT_CaloJet500_NoJetID       );
  chain.SetBranchAddress( "firedHLT_AK8PFJet450"       ,  &firedHLT_AK8PFJet450       ,  &b_firedHLT_AK8PFJet450       );
  chain.SetBranchAddress( "firedHLT_PFJet450"       ,  &firedHLT_PFJet450       ,  &b_firedHLT_PFJet450       );
  chain.SetBranchAddress( "runno"                     ,  &runno                     ,  &b_runno                     );
  chain.SetBranchAddress( "lumiblock"                 ,  &lumiblock                 ,  &b_lumiblock                 );
  chain.SetBranchAddress( "evtno"                     ,  &evtno                     ,  &b_evtno                     );
  chain.SetBranchAddress( "JetEt"                     ,  JetEt                      ,  &b_JetEt                     );
  chain.SetBranchAddress( "JetPx"                     ,  JetPx                      ,  &b_JetPx                     );
  chain.SetBranchAddress( "JetPy"                     ,  JetPy                      ,  &b_JetPy                     );
  chain.SetBranchAddress( "JetEta"                    ,  JetEta                     ,  &b_JetEta                    );
  chain.SetBranchAddress( "JetPhi"                    ,  JetPhi                     ,  &b_JetPhi                    );
  chain.SetBranchAddress( "JetNeutHadFrac"            ,  JetNeutHadFrac             ,  &b_JetNeutHadFrac            );
  chain.SetBranchAddress( "JetNeutEMFrac"             ,  JetNeutEMFrac              ,  &b_JetNeutEMFrac             );
  chain.SetBranchAddress( "JetChgHadFrac"             ,  JetChgHadFrac              ,  &b_JetChgHadFrac             );
  chain.SetBranchAddress( "JetMuFrac"                 ,  JetMuFrac                  ,  &b_JetMuFrac                 );
  chain.SetBranchAddress( "JetChgEMFrac"              ,  JetChgEMFrac               ,  &b_JetChgEMFrac              );
  chain.SetBranchAddress( "JetNConstituents"          ,  JetNConstituents           ,  &b_JetNConstituents          );
  chain.SetBranchAddress( "JetNNeutConstituents"      ,  JetNNeutConstituents       ,  &b_JetNNeutConstituents      );
  chain.SetBranchAddress( "JetNChgConstituents"       ,  JetNChgConstituents        ,  &b_JetNChgConstituents       );
  chain.SetBranchAddress( "EleEt"                     ,  EleEt                      ,  &b_EleEt                     );
  chain.SetBranchAddress( "ElePx"                     ,  ElePx                      ,  &b_ElePx                     );
  chain.SetBranchAddress( "ElePy"                     ,  ElePy                      ,  &b_ElePy                     );
  chain.SetBranchAddress( "EleEta"                    ,  EleEta                     ,  &b_EleEta                    );
  chain.SetBranchAddress( "ElePhi"                    ,  ElePhi                     ,  &b_ElePhi                    );
  chain.SetBranchAddress( "PhEt"                      ,  PhEt                       ,  &b_PhEt                      );
  chain.SetBranchAddress( "PhPx"                      ,  PhPx                       ,  &b_PhPx                      );
  chain.SetBranchAddress( "PhPy"                      ,  PhPy                       ,  &b_PhPy                      );
  chain.SetBranchAddress( "PhEta"                     ,  PhEta                      ,  &b_PhEta                     );
  chain.SetBranchAddress( "PhPhi"                     ,  PhPhi                      ,  &b_PhPhi                     );
  chain.SetBranchAddress( "MuEt"                      ,  MuEt                       ,  &b_MuEt                      );
  chain.SetBranchAddress( "MuPx"                      ,  MuPx                       ,  &b_MuPx                      );
  chain.SetBranchAddress( "MuPy"                      ,  MuPy                       ,  &b_MuPy                      );
  chain.SetBranchAddress( "MuEta"                     ,  MuEta                      ,  &b_MuEta                     );
  chain.SetBranchAddress( "MuPhi"                     ,  MuPhi                      ,  &b_MuPhi                     );
  chain.SetBranchAddress( "MuPFdBiso"                 , MuPFdBiso                   ,  &b_MuPFdBiso                 );
  chain.SetBranchAddress( "Met"                       ,  &Met                       ,  &b_Met                       );
  chain.SetBranchAddress( "Sphericity"                ,  &Sphericity                ,  &b_Sphericity                );
  //chain.SetBranchAddress( "Gen_Multiplicity"          ,  &Gen_Multiplicity          ,  &b_Gen_Multiplicity          );
  //chain.SetBranchAddress( "FinalGenST"                ,  &FinalGenST                ,  &b_FinalGenST                );
  //chain.SetBranchAddress( "GenST71"                   ,  &GenST71                   ,  &b_GenST71                   );
  //chain.SetBranchAddress( "FinalGenPDGID"             ,  FinalGenPDGID              ,  &b_FinalGenPDGID             );
  //chain.SetBranchAddress( "FinalGenEt"                ,  FinalGenEt                 ,  &b_FinalGenEt                );
  //chain.SetBranchAddress( "FinalGenPt"                ,  FinalGenPt                 ,  &b_FinalGenPt                );
  //chain.SetBranchAddress( "FinalGenEta"               ,  FinalGenEta                ,  &b_FinalGenEta               );
  //chain.SetBranchAddress( "FinalGenPhi"               ,  FinalGenPhi                ,  &b_FinalGenPhi               );
  //chain.SetBranchAddress( "HT"                       ,  &HT                       ,  &b_HT                       );

  chain.SetBranchAddress( "PDGID_Gen4"                       ,  &PDGID_Gen4                       ,  &b_PDGID_Gen4                       );
  chain.SetBranchAddress( "PDGID_Gen5"                       ,  &PDGID_Gen5                       ,  &b_PDGID_Gen5                       );
  chain.SetBranchAddress( "N_Gen71"                       ,  &N_Gen71                       ,  &b_N_Gen71                       );
  chain.SetBranchAddress( "PDGID_Gen71"                       ,  PDGID_Gen71                       ,  &b_PDGID_Gen71                       );
  chain.SetBranchAddress( "Eta_Gen4"                       ,  &Eta_Gen4                       ,  &b_Eta_Gen4                       );
  chain.SetBranchAddress( "Eta_Gen5"                       ,  &Eta_Gen5                       ,  &b_Eta_Gen5                       );
  chain.SetBranchAddress( "Et_Gen4"                       ,  &Et_Gen4                       ,  &b_Et_Gen4                       );
  chain.SetBranchAddress( "Et_Gen5"                       ,  &Et_Gen5                       ,  &b_Et_Gen5                       );
  chain.SetBranchAddress( "GenST71"                       ,  &GenST71                       ,  &b_GenST71                       );
  chain.SetBranchAddress( "Eta_Gen71"                       ,  Eta_Gen71                       ,  &b_Eta_Gen71                       );
  chain.SetBranchAddress( "Et_Gen71"                       ,  Et_Gen71                       ,  &b_Et_Gen71                       );
  chain.SetBranchAddress( "GenST_interaction"          ,  &GenST_interaction       ,  &b_GenST_interaction        );
  chain.SetBranchAddress( "PDGID_interaction"           ,  PDGID_interaction        ,  &b_PDGID_interaction            );
  chain.SetBranchAddress( "Eta_interaction"             ,  Eta_interaction               ,  &b_Eta_interaction            );
  chain.SetBranchAddress( "Et_interaction"                       ,  Et_interaction               ,  &b_Et_interaction             );
  chain.SetBranchAddress( "GenST_radiation"          ,  &GenST_radiation       ,  &b_GenST_radiation        );
  chain.SetBranchAddress( "PDGID_radiation"           ,  PDGID_radiation       ,  &b_PDGID_radiation            );
  chain.SetBranchAddress( "Eta_radiation"             ,  Eta_radiation              ,  &b_Eta_radiation         );
  chain.SetBranchAddress( "Et_radiation"               ,  Et_radiation             ,  &b_Et_radiation           );
  chain.SetBranchAddress( "ST_GenJet"          ,  &ST_GenJet       ,  &b_ST_GenJet        );
  chain.SetBranchAddress( "ST_GenJet_noISR"          ,  &ST_GenJet_noISR       ,  &b_ST_GenJet_noISR        );


  const int nEvents = chain.GetEntries();
  int   N_failpassBadChCand   = 0;
  cout << "Number of events in chain is: " << nEvents << endl;
  bool passMETfilterList = true;
  // loop over all events
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (iEvent%50000==0) {
      cout << std::fixed << std::setw(3) << std::setprecision(1) << (float(iEvent)/float(nEvents))*100 << "% done: Scanned " << iEvent << " events." << endl;
    }

    // reset variables
    isTightJet          = false ;
    OurMet              = 0.    ;
    Px                  = 0.    ;
    Py                  = 0.    ;
    ST                  = 0.    ;
    ST_2leadjets        = 0.    ;
    multiplicity        = 0     ;
    OurMet_tight        = 0.    ;
    Px_tight            = 0.    ;
    Py_tight            = 0.    ;
    ST_tight            = 0.    ;
    multiplicity_tight  = 0     ;
    eventHasMuon        = false ;
    eventHasPhoton      = false ;
    eventHasElectron    = false ;
    count2jets          = 0;
    std::fill(std::begin( TightJets ), std::end( TightJets ), false );

    chain.GetEntry(iEvent);
    //std::cout << " event number: " << evtno <<std::endl;
    // apply trigger and filter requirements
    //TODO
    //if ( isData && is2015D ){
          //if (    !firedHLT_PFHT800 || !passed_globalTightHalo2016Filter || !passed_goodVertices || !passed_eeBadScFilter     ) {
          //if (!firedHLT_PFHT900 || (!firedHLT_CaloJet500_NoJetID && !firedHLT_AK8PFJet450 && !firedHLT_PFJet450)) {
	//	cout <<"firedHLT_PFHT800_v3 is"<<firedHLT_PFHT800_v3<<"  halofilterResult is "<<passed_globalTightHalo2016Filter<< "  vertices = "<< passed_goodVertices<<"  ee="<<passed_eeBadScFilter<<endl;
		//continue;
	  //}
    //}
    //if (isData && !is2015D){
          //if (    !firedHLT_PFHT800 || !passed_globalTightHalo2016Filter || !passed_goodVertices || !passed_eeBadScFilter     ) continue;
          //if (!firedHLT_PFHT900 || (!firedHLT_CaloJet500_NoJetID && !firedHLT_AK8PFJet450 && !firedHLT_PFJet450)) continue;
    //}
    //if (HT<1050) continue;
    if (!passed_filterbadChCandidate){N_failpassBadChCand++;}
    if (!firedHLT_PFHT800) continue;
    if (!passed_globalTightHalo2016Filter || !passed_goodVertices || !passed_HBHENoiseFilter || !passed_HBHENoiseIsoFilter || !passed_filterbadPFMuon || !passed_filterbadChCandidate || !passed_EcalDeadCellTriggerPrimitiveFilter) continue; //for 2016 samples
    //if (!passed_globalSuperTightHalo2016Filter || !passed_goodVertices || !passed_HBHENoiseFilter || !passed_HBHENoiseIsoFilter || !passed_BadPFMuonFilter || !passed_EcalDeadCellTriggerPrimitiveFilter || !passed_eeBadScFilter) continue; //for UL samples
    //if (!passed_globalSuperTightHalo2016Filter || !passed_goodVertices || !passed_HBHENoiseFilter || !passed_HBHENoiseIsoFilter || !passed_filterbadPFMuon || !passed_filterbadChCandidate || !passed_EcalDeadCellTriggerPrimitiveFilter) continue;
        // use Yutaro's method for applying the event filter
        passMETfilterList=true;
        auto rItr(list.find(runno));
        if (rItr != list.end()) {
          if (rItr->second.find(evtno) != rItr->second.end()){
            if (dumpBigEvents && debugFlag) {
              sprintf(messageBuffer, "Event in MET list skipped: run number %d lumi section %d event number %lld\n", runno, lumiblock, evtno);
              outTextFile << messageBuffer;
            }
            passMETfilterList = false;
            continue;
          }
        }
        if (!passMETfilterList) cout << "ERROR! This event should be filtered!" << endl;
        if ( runno == 254790 && (lumiblock==211 || lumiblock==395) ) {
          sprintf(messageBuffer, "Event in lumiblock that could not be filtered skipped: run number %d lumi section %d event number %%lld\n", runno, lumiblock, evtno);
          outTextFile << messageBuffer;
          continue;
        }

        Hist_sphericity.Fill(Sphericity);

        // apply isolation requirement and calculate ST and MHT.
        //Jets
        for (int iJet = 0; iJet < 25; ++iJet) {
          passIso=true;
          passIso_tight=true;
          isTightJet=false;
          JetMuonEt     =0;
          JetElectronEt =0;
          JetPhotonEt   =0;
          if (fabs(JetEta[iJet])<=3 && JetNeutHadFrac[iJet]<0.9 && JetNeutEMFrac[iJet]<0.9 && JetNConstituents[iJet]>1 && JetMuFrac[iJet]<0.8) {
            isTightJet=true;
            if (fabs(JetEta[iJet])<=2.4) {
              if ( JetNChgConstituents[iJet] > 0 && JetChgHadFrac[iJet] > 0 && JetChgEMFrac[iJet]<0.9) isTightJet=true;
              else isTightJet=false;
            }
          }
          if (fabs(JetEta[iJet])>3 && JetNeutEMFrac[iJet] < 0.9 && JetNNeutConstituents[iJet] > 10 ) isTightJet=true;
          TightJets[iJet]=isTightJet;
          if (JetEt[iJet]>70.) {
            //std::cout << " Jet Et: "<< JetEt[iJet]<< " Eta: " << JetEta[iJet] << " Phi: " << JetPhi[iJet] <<std::endl;
            for (int iMuon = 0; iMuon < 25; ++iMuon ) {
              if (MuEt[iMuon]>70 && MuPFdBiso[iMuon]<0.15) {
                eventHasMuon = true;
                if (JetEt[iJet] && dR(JetEta[iJet],JetPhi[iJet], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
                  JetMuonEt+=MuEt[iMuon];
                  if (MuEt[iMuon]<150) {
                    MuonJetIso1.Fill(MuEt[iMuon]/JetEt[iJet]);
                    MuonJetoverlapdR1.Fill(dR(JetEta[iJet],JetPhi[iJet],MuEta[iMuon],MuPhi[iMuon]));
                  }
                  if (150<=MuEt[iMuon] && MuEt[iMuon]<250) {
                    MuonJetIso2.Fill(MuEt[iMuon]/JetEt[iJet]);
                    MuonJetoverlapdR2.Fill(dR(JetEta[iJet],JetPhi[iJet],MuEta[iMuon],MuPhi[iMuon]));
                  }
                  if (250<=MuEt[iMuon] && MuEt[iMuon]<400) {
                    MuonJetIso3.Fill(MuEt[iMuon]/JetEt[iJet]);
                    MuonJetoverlapdR3.Fill(dR(JetEta[iJet],JetPhi[iJet],MuEta[iMuon],MuPhi[iMuon]));
                  }
                  if (400<=MuEt[iMuon]) {
                    MuonJetIso4.Fill(MuEt[iMuon]/JetEt[iJet]);
                    MuonJetoverlapdR4.Fill(dR(JetEta[iJet],JetPhi[iJet],MuEta[iMuon],MuPhi[iMuon]));
                  }
                  if (JetMuonEt>0.8*JetEt[iJet]) {
                    passIso = false;
                    if (isTightJet) {
                      passIso_tight=false;
                      if (dumpIsoInfo) {
                        sprintf(messageBuffer, "Jet number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %lld\n", iJet, iMuon, runno, lumiblock, evtno);
                        outTextFile << messageBuffer;
                      }
                    }
                    break;
                  }
                }
              }
            }
            for (int iElectron = 0; iElectron < 25; ++iElectron ) {
              if (EleEt[iElectron]>70) {
                eventHasElectron = true;
                if (dR(JetEta[iJet],JetPhi[iJet], EleEta[iElectron], ElePhi[iElectron]) < 0.3) {
                  JetElectronEt+=EleEt[iElectron];
                  if (EleEt[iElectron]<150) {
                    ElectronJetIso1.Fill(EleEt[iElectron]/JetEt[iJet]);
                    ElectronJetoverlapdR1.Fill(dR(JetEta[iJet],JetPhi[iJet],EleEta[iElectron],ElePhi[iElectron]));
                  }
                  if (150<=EleEt[iElectron] && EleEt[iElectron]<250) {
                    ElectronJetIso2.Fill(EleEt[iElectron]/JetEt[iJet]);
                    ElectronJetoverlapdR2.Fill(dR(JetEta[iJet],JetPhi[iJet],EleEta[iElectron],ElePhi[iElectron]));
                  }
                  if (250<=EleEt[iElectron] && EleEt[iElectron]<400) {
                    ElectronJetIso3.Fill(EleEt[iElectron]/JetEt[iJet]);
                    ElectronJetoverlapdR3.Fill(dR(JetEta[iJet],JetPhi[iJet],EleEta[iElectron],ElePhi[iElectron]));
                  }
                  if (400<=EleEt[iElectron]) {
                    ElectronJetIso4.Fill(EleEt[iElectron]/JetEt[iJet]);
                    ElectronJetoverlapdR4.Fill(dR(JetEta[iJet],JetPhi[iJet],EleEta[iElectron],ElePhi[iElectron]));
                  }
                  if (JetElectronEt > 0.7*JetEt[iJet] ) {
                    passIso = false;
                    if (isTightJet) {
                      passIso_tight=false;
                      if (dumpIsoInfo) {
                        sprintf(messageBuffer, "Jet number %d failed isolation with Electron number %d  in run number %d lumi section %d event number %lld\n", iJet, iElectron, runno, lumiblock, evtno);
                        outTextFile << messageBuffer;
                      }
                    }
                    break;
                  }
                }
              }
            }
            for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
              if (PhEt[iPhoton]>70) {
                eventHasPhoton = true;
                if (dR(JetEta[iJet],JetPhi[iJet], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
                  JetPhotonEt+=PhEt[iPhoton];
                  if (PhEt[iPhoton]<150) {
                    PhotonJetIso1.Fill(PhEt[iPhoton]/JetEt[iJet]);
                    PhotonJetoverlapdR1.Fill(dR(JetEta[iJet],JetPhi[iJet],PhEta[iPhoton],PhPhi[iPhoton]));
                  }
                  if (150<=PhEt[iPhoton] && PhEt[iPhoton]<250) {
                    PhotonJetIso2.Fill(PhEt[iPhoton]/JetEt[iJet]);
                    PhotonJetoverlapdR2.Fill(dR(JetEta[iJet],JetPhi[iJet],PhEta[iPhoton],PhPhi[iPhoton]));
                  }
                  if (250<=PhEt[iPhoton] && PhEt[iPhoton]<400) {
                    PhotonJetIso3.Fill(PhEt[iPhoton]/JetEt[iJet]);
                    PhotonJetoverlapdR3.Fill(dR(JetEta[iJet],JetPhi[iJet],PhEta[iPhoton],PhPhi[iPhoton]));
                  }
                  if (400<=PhEt[iPhoton]) {
                    PhotonJetIso4.Fill(PhEt[iPhoton]/JetEt[iJet]);
                    PhotonJetoverlapdR4.Fill(dR(JetEta[iJet],JetPhi[iJet],PhEta[iPhoton],PhPhi[iPhoton]));
                  }
                  if (JetPhotonEt>0.5*JetEt[iJet] ) {
                    passIso = false;
                    if (isTightJet) {
                      passIso_tight=false;
                      if (dumpIsoInfo) {
                        sprintf(messageBuffer, "Jet number %d failed isolation with Photon number %d  in run number %d lumi section %d event number %lld\n", iJet, iPhoton, runno, lumiblock, evtno);
                        outTextFile << messageBuffer;
                      }
                    }
                    break;
                  }
                }
              }
            }
            if (!passIso) continue;

            if (debugFlag) outTextFile << "    JetEt for jet number " << iJet << " is: " << JetEt[iJet] << endl;

            if (iJet!=0 || isTightJet){
              ST += JetEt[iJet];
              multiplicity+=1;
              Px += JetPx[iJet];
              Py += JetPy[iJet];
              if(count2jets < 2){ST_2leadjets += JetEt[iJet];}
              ++count2jets;
            }
            /*
            ST += JetEt[iJet];
            multiplicity+=1;
            Px += JetPx[iJet];
            Py += JetPy[iJet];
            */
            /*
            if (iJet==0 && isTightJet){
              ST += JetEt[iJet];
              multiplicity+=1;
              Px += JetPx[iJet];
              Py += JetPy[iJet];
            }
            if(ijet!=0){
              ST += JetEt[iJet];
              multiplicity+=1;
              Px += JetPx[iJet];
              Py += JetPy[iJet];
            }
            */
            if(isTightJet) {
              ST_tight += JetEt[iJet];
              multiplicity_tight+=1;
              Px_tight += JetPx[iJet];
              Py_tight += JetPy[iJet];
            }
            if (debugFlag && dumpIsoInfo) {
              sprintf(messageBuffer, "Jet number %d passed isolation in run number %d lumi section %d event number %lld.\n       It had Px=%f and Py=%f\n", iJet, runno, lumiblock, evtno, JetPx[iJet], JetPy[iJet]);
              outTextFile << messageBuffer;
              sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
              outTextFile << messageBuffer;
            }
          }
          else break;
        }
        //std::cout << " Jet Multiplicity: "<< multiplicity << " ST: "<< ST <<std::endl;

        //Electrons
        if (eventHasElectron) {
          for (int iElectron = 0; iElectron < 25; ++iElectron) {
            passIso=true;
            passIso_tight=true;
            if (EleEt[iElectron]>150.) {
              for (int iJet = 0; iJet < 25; ++iJet ) {
                if (JetEt[iJet]>70 && dR(EleEta[iElectron],ElePhi[iElectron], JetEta[iJet], JetPhi[iJet]) < 0.3) {
                  if (EleEt[iElectron]<0.7*JetEt[iJet]) {
                    passIso = false;
                    if(TightJets[iJet]) {
                      passIso_tight=false;
                      if (dumpIsoInfo) {
                        sprintf(messageBuffer, "Electron number %d failed isolation with Jet number %d  in run number %d lumi section %d event number %lld\n", iElectron, iJet, runno, lumiblock, evtno);
                        outTextFile << messageBuffer;
                      }
                      break;
                    }
                  }
                }
              }
              if (!passIso_tight) continue;

              // Throw away electron if there's an electron/muon overlap.
              for (int iMuon = 0; iMuon < 25; ++iMuon ) {
                if (MuEt[iMuon]>150 && MuPFdBiso[iMuon]<0.15 && dR(EleEta[iElectron],ElePhi[iElectron], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
                  passIso = false;
                  passIso_tight = false;
                  if (dumpIsoInfo) {
                    sprintf(messageBuffer, "Electron number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %lld\n", iElectron, iMuon, runno, lumiblock, evtno);
                    outTextFile << messageBuffer;
                  }
                  break;
                }
              }
              if (!passIso_tight) continue;

              if (debugFlag) cout << "    EleEt for electron number " << iElectron << " is: " << EleEt[iElectron] << endl;
              ST_tight += EleEt[iElectron];
              multiplicity_tight+=1;
              Px_tight += ElePx[iElectron];
              Py_tight += ElePy[iElectron];
              if (passIso) {
                ST += EleEt[iElectron];
                multiplicity+=1;
                Px += ElePx[iElectron];
                Py += ElePy[iElectron];
              }
              if (debugFlag && dumpIsoInfo) {
                sprintf(messageBuffer, "Ele number %d passed isolation in run number %d lumi section %d event number %lld.      \n It had Px=%f and Py=%f\n", iElectron, runno, lumiblock, evtno, ElePx[iElectron], ElePy[iElectron]);
                outTextFile << messageBuffer;
                sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
                outTextFile << messageBuffer;
              }
            }
            else break;
          }
        }
        //std::cout << " Jet+Electron Multiplicity: "<< multiplicity << " ST: "<< ST <<std::endl;

        //Photons
        if (eventHasPhoton) {
          for (int iPhoton = 0; iPhoton < 25; ++iPhoton) {
            passIso=true;
            if (PhEt[iPhoton]>150.) {
              for (int iJet = 0; iJet < 25; ++iJet ) {
                if (JetEt[iJet]>70 && dR(PhEta[iPhoton],PhPhi[iPhoton], JetEta[iJet], JetPhi[iJet]) < 0.3) {
                  if (PhEt[iPhoton]<0.5*JetEt[iJet]) {
                    passIso = false;
                    if (TightJets[iJet]) {
                      passIso_tight=false;
                      if (dumpIsoInfo) {
                        sprintf(messageBuffer, "Photon number %d failed isolation with Jet number %d  in run number %d lumi section %d event number %lld\n", iPhoton, iJet, runno, lumiblock, evtno);
                        outTextFile << messageBuffer;
                      }
                      break;
                    }
                  }
                }
              }
              if (!passIso_tight) continue;

              // Throw out photon if there's a photon/muon overlap
              for (int iMuon = 0; iMuon < 25; ++iMuon ) {
                if (MuEt[iMuon]>70 && MuPFdBiso[iMuon]<0.15 && dR(PhEta[iPhoton], PhPhi[iPhoton], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
                  if (dumpIsoInfo) {
                    sprintf(messageBuffer, "Photon number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %lld\n", iPhoton, iMuon, runno, lumiblock, evtno);
                    outTextFile << messageBuffer;
                  }
                  passIso = false;
                  passIso_tight = false;
                  break;
                }
              }
              if (!passIso_tight) continue;

              // Throw out photon if there's a photon/electron overlap
              for (int iElectron = 0; iElectron < 25; ++iElectron ) {
                if (EleEt[iElectron]>70 && dR(PhEta[iPhoton], PhPhi[iPhoton], EleEta[iElectron], ElePhi[iElectron]) < 0.3) {
                  if (dumpIsoInfo) {
                    sprintf(messageBuffer, "Photon number %d failed isolation with Electron number %d  in run number %d lumi section %d event number %lld\n", iPhoton, iElectron, runno, lumiblock, evtno);
                    outTextFile << messageBuffer;
                  }
                  passIso = false;
                  passIso_tight = false;
                  break;
                }
              }
              if (!passIso_tight) continue;

              if (debugFlag) cout << "    PhEt for photon number " << iPhoton << " is: " << PhEt[iPhoton] << endl;

              ST_tight += PhEt[iPhoton];
              multiplicity_tight+=1;
              Px_tight += PhPx[iPhoton];
              Py_tight += PhPy[iPhoton];
              if (passIso) {
                ST += PhEt[iPhoton];
                multiplicity+=1;
                Px += PhPx[iPhoton];
                Py += PhPy[iPhoton];
              }
              if (debugFlag && dumpIsoInfo) {
                sprintf(messageBuffer, "Photon number %d passed isolation in run number %d lumi section %d event number %lld.\n      It had Px=%f and Py=%f\n", iPhoton, runno, lumiblock, evtno, PhPx[iPhoton], PhPy[iPhoton]);
                outTextFile << messageBuffer;
                sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
                outTextFile << messageBuffer;
              }
            }
            else break;
          }
        }
        //std::cout << " Jet+Electron+Photon Multiplicity: "<< multiplicity << " ST: "<< ST <<std::endl;

        //Muons
        if (eventHasMuon) {
          for (int iMuon = 0; iMuon < 25; ++iMuon) {
            passIso=true;
            passIso_tight=true;
            if (MuEt[iMuon]>150. && MuPFdBiso[iMuon]<0.15) {
              if (debugFlag) cout << "    MuEt for muon number " << iMuon << " is: " << MuEt[iMuon] << endl;
              ST += MuEt[iMuon];
              multiplicity+=1;
              Px += MuPx[iMuon];
              Py += MuPy[iMuon];
              ST_tight += MuEt[iMuon];
              multiplicity_tight+=1;
              Px_tight += MuPx[iMuon];
              Py_tight += MuPy[iMuon];
              if (debugFlag && dumpIsoInfo) {
                sprintf(messageBuffer, "Muon number %d passed isolation in run number %d lumi section %d event number %lld.\n       It had Px=%f and Py=%f\n", iMuon, runno, lumiblock, evtno, MuPx[iMuon], MuPy[iMuon]);
                outTextFile << messageBuffer;
                sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
                outTextFile << messageBuffer;
              }
            }
            else break;
          }
        }
        //std::cout << " Jet+Electron+Photon+Muon Multiplicity: "<< multiplicity << " ST: "<< ST <<std::endl;


        //debug info and big ST printing
        if (debugFlag) cout << "    Met from PAT collection is: " << Met << endl;
        OurMet = std::sqrt(Px*Px + Py*Py);
        OurMet_tight = std::sqrt(Px_tight*Px_tight + Py_tight*Py_tight);
        if (debugFlag) cout << "    Met calculated according to my recipe is: " << OurMet << endl;
        if (0.5*ST>Met || !useMETcut) passMetCut=true;
        else passMetCut=false;
        if (0.5*ST_tight>Met || !useMETcut) {
          passMetCut_tight=true;
        }
        else {
          sprintf(messageBuffer, "Event number %lld failed the MET cut in run number %d lumi section %d.\n       It had MET/HT=%f \n", evtno, runno, lumiblock, (Met/ST_tight));
          outTextFile << messageBuffer;
          passMetCut_tight=false;
        }
        if (ST>1500) {
          if (ST<3000) METoverSumET.Fill(Met/ST);
          if (ST_tight>1500 && ST_tight<3000) METoverSumET_tight.Fill(Met/ST_tight);
          if (multiplicity>=2){
            METoverSumETinc2.Fill(Met/ST);
            if (multiplicity_tight >=2 && ST_tight>1500 && ST_tight<3000)METoverSumETinc2_tight.Fill(Met/ST_tight);
            if (multiplicity_tight >=3 && ST_tight>1500 && ST_tight<3000)METoverSumETinc3_tight.Fill(Met/ST_tight);
            if (multiplicity_tight >=4 && ST_tight>1500 && ST_tight<3000)METoverSumETinc4_tight.Fill(Met/ST_tight);
            if (multiplicity_tight >=5 && ST_tight>1500 && ST_tight<3000)METoverSumETinc5_tight.Fill(Met/ST_tight);
            if (multiplicity_tight >=6 && ST_tight>1500 && ST_tight<3000)METoverSumETinc6_tight.Fill(Met/ST_tight);
            if (multiplicity_tight >=7 && ST_tight>1500 && ST_tight<3000)METoverSumETinc7_tight.Fill(Met/ST_tight);
            if (multiplicity_tight >=8 && ST_tight>1500 && ST_tight<3000)METoverSumETinc8_tight.Fill(Met/ST_tight);
            if (multiplicity_tight >=9 && ST_tight>1500 && ST_tight<3000)METoverSumETinc9_tight.Fill(Met/ST_tight);
            if (multiplicity_tight >=10 && ST_tight>1500 && ST_tight<3000)METoverSumETinc10_tight.Fill(Met/ST_tight);
            if (eventHasMuon)                                           METoverSumETinc2hasMuon.Fill(Met/ST);
            if (eventHasPhoton)                                         METoverSumETinc2hasPhoton.Fill(Met/ST);
            if (eventHasElectron)                                       METoverSumETinc2hasElectron.Fill(Met/ST);
            if (!eventHasMuon && !eventHasPhoton && !eventHasElectron)  METoverSumETinc2onlyJets.Fill(Met/ST);
            if (multiplicity_tight>=2 && ST_tight>1500 && ST_tight<3000) {
              if (eventHasMuon)                                           METoverSumETinc2hasMuon_tight.Fill(Met/ST_tight);
              if (eventHasPhoton)                                         METoverSumETinc2hasPhoton_tight.Fill(Met/ST_tight);
              if (eventHasElectron)                                       METoverSumETinc2hasElectron_tight.Fill(Met/ST_tight);
              if (!eventHasMuon && !eventHasPhoton && !eventHasElectron)  METoverSumETinc2onlyJets_tight.Fill(Met/ST_tight);
            }
          }
        }
        STMHTnoMET = ST + OurMet;
        STMHTnoMET_tight = ST_tight + OurMet_tight;
        ST += Met;
        ST_tight += Met;
        if (passMetCut){
          stHist.Fill(ST);
          stHistMHT.Fill(STMHTnoMET);
        }
        if (passMetCut_tight){
          stHist_tight.Fill(ST_tight);
          stHistMHT_tight.Fill(STMHTnoMET_tight);
        }

        for (int iHist = 0; iHist<multMax-2; ++iHist) {
          if (multiplicity == iHist+2 && passMetCut) {
            stExcHist[iHist]->Fill(ST);
            Gen71stExcHist[iHist]->Fill(GenST71);
            GenInteractstExcHist[iHist]->Fill(GenST_interaction);
            GenUnderlystExcHist[iHist]->Fill(GenST_radiation);
            GenJetstExcHist[iHist]->Fill(ST_GenJet);
            GenJetNoISRstExcHist[iHist]->Fill(ST_GenJet_noISR);
            stExcHist_2leadjets[iHist]->Fill(ST_2leadjets);
            if (ST_GenJet - ST_GenJet_noISR < 70){stExcHist_noISR[iHist]->Fill(ST);}
          }
          if (multiplicity >= iHist+2 && passMetCut) stIncHist[iHist]->Fill(ST);
          if (multiplicity_tight == iHist+2 && passMetCut_tight) stExcHist_tight[iHist]->Fill(ST_tight);
          if (multiplicity_tight >= iHist+2 && passMetCut_tight) stIncHist_tight[iHist]->Fill(ST_tight);
        }


        for (int iHist = 0; iHist<multMax-2; ++iHist) {
          if (multiplicity == iHist+2 && passMetCut) stExcHistMHT[iHist]->Fill(STMHTnoMET);
          if (multiplicity >= iHist+2 && passMetCut) stIncHistMHT[iHist]->Fill(STMHTnoMET);
          if (multiplicity_tight == iHist+2 && passMetCut_tight) stExcHistMHT_tight[iHist]->Fill(STMHTnoMET_tight);
          if (multiplicity_tight >= iHist+2 && passMetCut_tight) stIncHistMHT_tight[iHist]->Fill(STMHTnoMET_tight);
        }
        METvsMHT.Fill(OurMet,Met);
        METvsMHT_tight.Fill(OurMet_tight,Met);
        if (multiplicity>=2){
          METvsMHTinc2.Fill(OurMet,Met);
          METvsMHTinc2_tight.Fill(OurMet_tight,Met);
          if (eventHasMuon)                                           METvsMHTinc2hasMuon.Fill(OurMet, Met);
          if (eventHasPhoton)                                         METvsMHTinc2hasPhoton.Fill(OurMet, Met);
          if (eventHasElectron)                                       METvsMHTinc2hasElectron.Fill(OurMet, Met);
          if (!eventHasMuon && !eventHasPhoton && !eventHasElectron)  METvsMHTinc2onlyJets.Fill(OurMet, Met);
          if (multiplicity_tight>=2) {
            if (eventHasMuon)                                           METvsMHTinc2hasMuon_tight.Fill(OurMet_tight, Met);
            if (eventHasPhoton)                                         METvsMHTinc2hasPhoton_tight.Fill(OurMet_tight, Met);
            if (eventHasElectron)                                       METvsMHTinc2hasElectron_tight.Fill(OurMet_tight, Met);
            if (!eventHasMuon && !eventHasPhoton && !eventHasElectron)  METvsMHTinc2onlyJets_tight.Fill(OurMet_tight, Met);
          }
        }
        if (dumpIsoInfo && fabs(OurMet-Met)>300) {
          sprintf(messageBuffer, "MET-MHT is %f in run number %d lumi section %d event number %lld. ST is %f and multiplicity is %d\n", Met-OurMet, runno, lumiblock, evtno, ST, multiplicity);
          outTextFile << messageBuffer;
          if (debugFlag) cout << messageBuffer;
        }

        Hist2D_N_ST->Fill(multiplicity,ST);


        // dump info on events with very big ST
        /*
        if (multiplicity==8){
          FinalGenSTExc08.Fill(FinalGenST);
          if (ST>3500 && dumpBigEvents) {
            //sprintf(messageBuffer, "In run number %d lumi section %d event number %lld: ST is %f, ST_tight is %f, and multiplicity is %d\n", runno, lumiblock, evtno, ST, ST_tight, multiplicity);
            sprintf(messageBuffer, "In run number %d lumi section %d event number %llu: ST is %f, ST_tight is %f, and multiplicity is %d\n", runno, lumiblock, evtno, ST, ST_tight, multiplicity);
            std::cout << "N=8 3.5 TeV up: In run number "<< runno<< " lumi section "<< lumiblock<< " event number "<< evtno<<" :ST is " << ST<< " and multiplicity is " <<multiplicity<<"\n"<<std::endl;
            outTextFile << messageBuffer;
            for (int j=0; j<25; ++j) {
              if(JetEt[j]>0.000) {
                sprintf(messageBuffer, "    Jet %d has TightJet=%d Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, TightJets[j],  JetEt[j], JetPx[j], JetPy[j], JetEta[j], JetPhi[j]);
                std::cout << " Jet: "<< j<< " has TightJet= "<< TightJets[j]<< " Et= "<< JetEt[j]<<", Px= " << JetPx[j]<< ", Py=" <<JetPy[j]<< ", Eta=" <<JetEta[j]<< ", Phi=" <<JetPhi[j]<<std::endl;
                outTextFile << messageBuffer;
                JetPtExc08_3500Up.Fill(JetEt[j]);
              }
              if(EleEt[j]>0.000) {std::cout << " Electron: "<< j<<" Et= "<< EleEt[j]<<", Px= " << ElePx[j]<< ", Py=" <<ElePy[j]<< ", Eta=" <<EleEta[j]<< ", Phi=" <<ElePhi[j]<<std::endl;}
              if(MuEt[j]>0.000) {std::cout << " Muon: "<< j<<" Et= "<< MuEt[j]<<", Px= " << MuPx[j]<< ", Py=" <<MuPy[j]<< ", Eta=" <<MuEta[j]<< ", Phi=" <<MuPhi[j]<<std::endl;}
              if(PhEt[j]>0.000) {std::cout << " Photon: "<< j<<" Et= "<< PhEt[j]<<", Px= " << PhPx[j]<< ", Py=" <<PhPy[j]<< ", Eta=" <<PhEta[j]<< ", Phi=" <<PhPhi[j]<<std::endl;}

              JetEtaExc08_3500Up.Fill(JetEta[j]);
              if (debugFlag) cout  << messageBuffer;
            }
            for (int j=0; j<100; ++j) {
              if(FinalGenPDGID[j]!=0){FinalGenPDGIDExc08_3500Up.Fill(FinalGenPDGID[j]);}
              if(FinalGenEt[j]>0.0){FinalGenEtExc08_3500Up.Fill(FinalGenEt[j]);}
              if(FinalGenEta[j]>-9000.0){FinalGenEtaExc08_3500Up.Fill(FinalGenEta[j]);}
            }
            FinalGenNExc08_3500Up.Fill(Gen_Multiplicity);
          }
            if (ST>2300 && ST<2400 && dumpBigEvents) {
              //sprintf(messageBuffer, "In run number %d lumi section %d event number %lld: ST is %f, ST_tight is %f, and multiplicity is %d\n", runno, lumiblock, evtno, ST, ST_tight, multiplicity);
              sprintf(messageBuffer, "In run number %d lumi section %d event number %llu: ST is %f, ST_tight is %f, and multiplicity is %d\n", runno, lumiblock, evtno, ST, ST_tight, multiplicity);
              std::cout << " N=8 2.3 to 2.4 TeV: In run number "<< runno<< " lumi section "<< lumiblock<< " event number "<< evtno<<" :ST is " << ST<< " and multiplicity is " <<multiplicity<<"\n"<<std::endl;
              outTextFile << messageBuffer;
              for (int j=0; j<25; ++j) {
                if(JetEt[j]>0.000) {
                  sprintf(messageBuffer, " Jet %d has TightJet=%d Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, TightJets[j],  JetEt[j], JetPx[j], JetPy[j], JetEta[j], JetPhi[j]);
                  std::cout << " Jet "<< j<< " has TightJet= "<< TightJets[j]<< " Et= "<< JetEt[j]<<", Px= " << JetPx[j]<< ", Py=" <<JetPy[j]<< ", Eta=" <<JetEta[j]<< ", Phi=" <<JetPhi[j]<<std::endl;
                  outTextFile << messageBuffer;
                  JetPtExc08_2300to2400.Fill(JetEt[j]);
                }
                if(EleEt[j]>0.000) {std::cout << " Electron: "<< j<<" Et= "<< EleEt[j]<<", Px= " << ElePx[j]<< ", Py=" <<ElePy[j]<< ", Eta=" <<EleEta[j]<< ", Phi=" <<ElePhi[j]<<std::endl;}
                if(MuEt[j]>0.000) {std::cout << " Muon: "<< j<<" Et= "<< MuEt[j]<<", Px= " << MuPx[j]<< ", Py=" <<MuPy[j]<< ", Eta=" <<MuEta[j]<< ", Phi=" <<MuPhi[j]<<std::endl;}
                if(PhEt[j]>0.000) {std::cout << " Photon: "<< j<<" Et= "<< PhEt[j]<<", Px= " << PhPx[j]<< ", Py=" <<PhPy[j]<< ", Eta=" <<PhEta[j]<< ", Phi=" <<PhPhi[j]<<std::endl;}
                if(JetEta[j]> -9000) {JetEtaExc08_2300to2400.Fill(JetEta[j]);}
                if (debugFlag) cout  << messageBuffer;
              }
              for (int j=0; j<100; ++j) {
                if(FinalGenPDGID[j]!=0){FinalGenPDGIDExc08_2300to2400.Fill(FinalGenPDGID[j]);}
                if(FinalGenEt[j]>0.0){FinalGenEtExc08_2300to2400.Fill(FinalGenEt[j]);}
                if(FinalGenEta[j]>-9000.0){FinalGenEtaExc08_2300to2400.Fill(FinalGenEta[j]);}
              }
              FinalGenNExc08_2300to2400.Fill(Gen_Multiplicity);
            }
        }


        if (multiplicity==2){
          FinalGenSTExc03.Fill(FinalGenST);
          if ( ST>1800 && ST<1900){
            std::cout << " N=3 1.8 to 1.9 TeV: In run number "<< runno<< " lumi section "<< lumiblock<< " event number "<< evtno<<" :ST is " << ST<< " and multiplicity is " <<multiplicity<<"\n"<<std::endl;
            for (int j=0; j<25; ++j) {
              if(JetEt[j]>0.000) {
                sprintf(messageBuffer, "    Jet %d has TightJet=%d Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, TightJets[j],  JetEt[j], JetPx[j], JetPy[j], JetEta[j], JetPhi[j]);
                std::cout << " Jet "<< j<< " has TightJet= "<< TightJets[j]<< " Et= "<< JetEt[j]<<", Px= " << JetPx[j]<< ", Py=" <<JetPy[j]<< ", Eta=" <<JetEta[j]<< ", Phi=" <<JetPhi[j]<<std::endl;
                outTextFile << messageBuffer;
                JetPtExc03_1800to1900.Fill(JetEt[j]);
              }
              if(EleEt[j]>0.000) {std::cout << " Electron: "<< j<<" Et= "<< EleEt[j]<<", Px= " << ElePx[j]<< ", Py=" <<ElePy[j]<< ", Eta=" <<EleEta[j]<< ", Phi=" <<ElePhi[j]<<std::endl;}
              if(MuEt[j]>0.000) {std::cout << " Muon: "<< j<<" Et= "<< MuEt[j]<<", Px= " << MuPx[j]<< ", Py=" <<MuPy[j]<< ", Eta=" <<MuEta[j]<< ", Phi=" <<MuPhi[j]<<std::endl;}
              if(PhEt[j]>0.000) {std::cout << " Photon: "<< j<<" Et= "<< PhEt[j]<<", Px= " << PhPx[j]<< ", Py=" <<PhPy[j]<< ", Eta=" <<PhEta[j]<< ", Phi=" <<PhPhi[j]<<std::endl;}
              if(JetEta[j]> -9000) {JetEtaExc03_1800to1900.Fill(JetEta[j]);}
              if (debugFlag) cout  << messageBuffer;
            }
            for (int j=0; j<100; ++j) {
              if(FinalGenPDGID[j]!=0){FinalGenPDGIDExc03_1800to1900.Fill(FinalGenPDGID[j]);}
              if(FinalGenEt[j]>0.0){FinalGenEtExc03_1800to1900.Fill(FinalGenEt[j]);}
              if(FinalGenEta[j]>-9000.0){FinalGenEtaExc03_1800to1900.Fill(FinalGenEta[j]);}
            }
            FinalGenNExc03_1800to1900.Fill(Gen_Multiplicity);
          }
        }
        */
        /*
        if (multiplicity==2){
          FinalGenST71Exc02.Fill(GenST71);
        }
        if (multiplicity==3){
          FinalGenST71Exc03.Fill(GenST71);
        }
        if (multiplicity==8){
          FinalGenST71Exc08.Fill(GenST71);
        }
        */

        if (  Sphericity > 0.05 && Sphericity < 0.3){
          //std::cout << " Multiplicity: "<< multiplicity << " ST: "<< ST << " MET: " << Met <<std::endl;
          //std::cout << " Sphericity: "<< Sphericity  <<std::endl;
          ST_Sphericity_p05top2_in.Fill(ST);
          if (multiplicity == 3){ST_Sphericity_p05top2_in_N3.Fill(ST);}
          if (multiplicity == 4){ST_Sphericity_p05top2_in_N4.Fill(ST);}
          if (multiplicity == 5){ST_Sphericity_p05top2_in_N5.Fill(ST);}
          if (multiplicity == 6){ST_Sphericity_p05top2_in_N6.Fill(ST);}
          if (multiplicity == 7){ST_Sphericity_p05top2_in_N7.Fill(ST);}
          if (multiplicity == 8){ST_Sphericity_p05top2_in_N8.Fill(ST);}
          if (multiplicity == 9){ST_Sphericity_p05top2_in_N9.Fill(ST);}
          if (multiplicity == 10){ST_Sphericity_p05top2_in_N10.Fill(ST);}
          if (multiplicity == 11){ST_Sphericity_p05top2_in_N11.Fill(ST);}
        }
        if (  Sphericity < 0.05 || Sphericity > 0.3){
          //std::cout << " Multiplicity: "<< multiplicity << " ST: "<< ST << " MET: " << Met <<std::endl;
          //std::cout << " Sphericity: "<< Sphericity  <<std::endl;
          ST_Sphericity_p05top2_out.Fill(ST);
          if (multiplicity == 3){ST_Sphericity_p05top2_out_N3.Fill(ST);}
          if (multiplicity == 4){ST_Sphericity_p05top2_out_N4.Fill(ST);}
          if (multiplicity == 5){ST_Sphericity_p05top2_out_N5.Fill(ST);}
          if (multiplicity == 6){ST_Sphericity_p05top2_out_N6.Fill(ST);}
          if (multiplicity == 7){ST_Sphericity_p05top2_out_N7.Fill(ST);}
          if (multiplicity == 8){ST_Sphericity_p05top2_out_N8.Fill(ST);}
          if (multiplicity == 9){ST_Sphericity_p05top2_out_N9.Fill(ST);}
          if (multiplicity == 10){ST_Sphericity_p05top2_out_N10.Fill(ST);}
          if (multiplicity == 11){ST_Sphericity_p05top2_out_N11.Fill(ST);}
        }

        if ( Sphericity < 0.1 ){
          ST_Sphericity_p1less.Fill(ST);
          if (multiplicity == 3){ST_Sphericity_p1less_N3.Fill(ST);}
          if (multiplicity == 4){ST_Sphericity_p1less_N4.Fill(ST);}
          if (multiplicity == 5){ST_Sphericity_p1less_N5.Fill(ST);}
          if (multiplicity == 6){ST_Sphericity_p1less_N6.Fill(ST);}
          if (multiplicity == 7){ST_Sphericity_p1less_N7.Fill(ST);}
          if (multiplicity == 8){ST_Sphericity_p1less_N8.Fill(ST);}
          if (multiplicity == 9){ST_Sphericity_p1less_N9.Fill(ST);}
        }

        if ( Sphericity > 0.1 ){
          ST_Sphericity_p1more.Fill(ST);
          if (multiplicity == 3){ST_Sphericity_p1more_N3.Fill(ST);}
          if (multiplicity == 4){ST_Sphericity_p1more_N4.Fill(ST);}
          if (multiplicity == 5){ST_Sphericity_p1more_N5.Fill(ST);}
          if (multiplicity == 6){ST_Sphericity_p1more_N6.Fill(ST);}
          if (multiplicity == 7){ST_Sphericity_p1more_N7.Fill(ST);}
          if (multiplicity == 8){ST_Sphericity_p1more_N8.Fill(ST);}
          if (multiplicity == 9){ST_Sphericity_p1more_N9.Fill(ST);}
        }

        /*
        if (FinalGenST < 600. && FinalGenST > 500.){
          //if (multiplicity==2){
          //  std::cout << " Same Gen ST, reco N=2: "<<std::endl;
          //  std::cout << " Gen N: "<< Gen_Multiplicity << " Gen ST: "<< FinalGenST <<std::endl;
          //  for (int j=0; j<100; ++j) {
          //    if(FinalGenEta[j]>-9000.0){
          //      std::cout << " Gen Particles "<< j<< " has PDGID = "<< FinalGenPDGID[j] << ", Et= "<< FinalGenEt[j] << ", Eta=" <<FinalGenEta[j]<< ", Phi=" <<FinalGenPhi[j]<<std::endl;
          //      FinalGenPDGIDExc08_3500Up.Fill(FinalGenPDGID[j]);
          //      FinalGenEtExc08_3500Up.Fill(FinalGenEt[j]);
          //      FinalGenEtaExc08_3500Up.Fill(FinalGenEta[j]);
          //      FinalGenNExc08_3500Up.Fill(Gen_Multiplicity);
          //    }
          //  }
          //  for (int j=0; j<25; ++j) {
          //    if(JetEt[j]>0.000){
          //      std::cout << " Reco Jets "<< j<< " has Et= "<< JetEt[j] << ", Eta=" <<JetEta[j]<< ", Phi=" <<JetPhi[j]<<std::endl;
          //    }
          //  }
          //}
          if (multiplicity==8){
            //std::cout << " Same Gen ST, reco N=8: "<<std::endl;
            //std::cout << " Gen N: "<< Gen_Multiplicity << " Gen ST: "<< FinalGenST <<std::endl;
            if (ST>2000){
              for (int j=0; j<100; ++j) {
                if(FinalGenEta[j]>-9000.0){
                  std::cout << " Gen Particles "<< j<< " has PDGID = "<< FinalGenPDGID[j] << ", Et= "<< FinalGenEt[j] << ", Eta=" <<FinalGenEta[j]<< ", Phi=" <<FinalGenPhi[j]<<std::endl;
                  FinalGenPDGIDExc08_3500Up.Fill(FinalGenPDGID[j]);
                  FinalGenEtExc08_3500Up.Fill(FinalGenEt[j]);
                  FinalGenEtaExc08_3500Up.Fill(FinalGenEta[j]);
                  FinalGenNExc08_3500Up.Fill(Gen_Multiplicity);
                }
              }
              for (int j=0; j<25; ++j) {
                if(JetEt[j]>0.000){
                  std::cout << " Reco Jets "<< j<< " has Et= "<< JetEt[j] << ", Eta=" <<JetEta[j]<< ", Phi=" <<JetPhi[j]<<std::endl;
                  JetPtExc08_3500Up.Fill(JetEt[j]);
                }
              }
            }
            if (ST<2000){
              for (int j=0; j<100; ++j) {
                if(FinalGenEta[j]>-9000.0){
                  //std::cout << " Gen Particles "<< j<< " has PDGID = "<< FinalGenPDGID[j] << ", Et= "<< FinalGenEt[j] << ", Eta=" <<FinalGenEta[j]<< ", Phi=" <<FinalGenPhi[j]<<std::endl;
                  FinalGenPDGIDExc08_2300to2400.Fill(FinalGenPDGID[j]);
                  FinalGenEtExc08_2300to2400.Fill(FinalGenEt[j]);
                  FinalGenEtaExc08_2300to2400.Fill(FinalGenEta[j]);
                  FinalGenNExc08_2300to2400.Fill(Gen_Multiplicity);
                }
              }
              for (int j=0; j<25; ++j) {
                if(JetEt[j]>0.000){
                  std::cout << " Reco Jets "<< j<< " has Et= "<< JetEt[j] << ", Eta=" <<JetEta[j]<< ", Phi=" <<JetPhi[j]<<std::endl;
                  JetPtExc08_2300to2400.Fill(JetEt[j]);
                }
              }
            }
          }
        }
        */
        /*
        if (multiplicity==2){
          FinalGenST71Exc02.Fill(GenST71);
          FinalGenN71Exc02.Fill(N_Gen71);
          FinalGenPDGID45Exc02.Fill(PDGID_Gen4);
          FinalGenPDGID45Exc02.Fill(PDGID_Gen5);
          FinalGenEta45Exc02.Fill(Eta_Gen4);
          FinalGenEta45Exc02.Fill(Eta_Gen5);
          FinalGenEt45Exc02.Fill(Et_Gen4);
          FinalGenEt45Exc02.Fill(Et_Gen5);
          FinalGenSTExc02_interaction.Fill(GenST_interaction);
          FinalGenSTExc02_radiation.Fill(GenST_radiation);
          for (int j=0; j<150; ++j) {
            if(Et_Gen71[j]>0.000){
              FinalGenPDGID71Exc02.Fill(PDGID_Gen71[j]);
              FinalGenEta71Exc02.Fill(Eta_Gen71[j]);
              FinalGenEt71Exc02.Fill(Et_Gen71[j]);
            }
          }
          for (int j=0; j<100; ++j) {
            if(Et_interaction[j]>0.000){
              FinalGenPDGIDExc02_interaction.Fill(PDGID_interaction[j]);
              FinalGenEtaExc02_interaction.Fill(Eta_interaction[j]);
              FinalGenEtExc02_interaction.Fill(Et_interaction[j]);
              float weight = Et_interaction[j];
              Energy_Eta_interaction_N2.Fill(Eta_interaction[j], weight);
            }
            if(Et_radiation[j]>0.000){
              FinalGenPDGIDExc02_radiation.Fill(PDGID_radiation[j]);
              FinalGenEtaExc02_radiation.Fill(Eta_radiation[j]);
              FinalGenEtExc02_radiation.Fill(Et_radiation[j]);
              float weight = Et_radiation[j];
              Energy_Eta_radiation_N2.Fill(Eta_radiation[j], weight);
            }
          }
        }
        if (multiplicity==3){
          FinalGenSTExc03_interaction.Fill(GenST_interaction);
          FinalGenSTExc03_radiation.Fill(GenST_radiation);
          FinalGenST71Exc03.Fill(GenST71);
        }
        if (multiplicity==8){
          FinalGenST71Exc08.Fill(GenST71);
          FinalGenN71Exc08.Fill(N_Gen71);
          FinalGenPDGID45Exc08.Fill(PDGID_Gen4);
          FinalGenPDGID45Exc08.Fill(PDGID_Gen5);
          FinalGenEta45Exc08.Fill(Eta_Gen4);
          FinalGenEta45Exc08.Fill(Eta_Gen5);
          FinalGenEt45Exc08.Fill(Et_Gen4);
          FinalGenEt45Exc08.Fill(Et_Gen5);
          FinalGenSTExc08_interaction.Fill(GenST_interaction);
          FinalGenSTExc08_radiation.Fill(GenST_radiation);
          for (int j=0; j<150; ++j) {
            if(Et_Gen71[j]>0.000){
              FinalGenPDGID71Exc08.Fill(PDGID_Gen71[j]);
              FinalGenEta71Exc08.Fill(Eta_Gen71[j]);
              FinalGenEt71Exc08.Fill(Et_Gen71[j]);
            }
          }
          for (int j=0; j<100; ++j) {
            if(Et_interaction[j]>0.000){
              FinalGenPDGIDExc08_interaction.Fill(PDGID_interaction[j]);
              FinalGenEtaExc08_interaction.Fill(Eta_interaction[j]);
              FinalGenEtExc08_interaction.Fill(Et_interaction[j]);
              float weight = Et_interaction[j];
              Energy_Eta_interaction_N8.Fill(Eta_interaction[j], weight);
            }
            if(Et_radiation[j]>0.000){
              FinalGenPDGIDExc08_radiation.Fill(PDGID_radiation[j]);
              FinalGenEtaExc08_radiation.Fill(Eta_radiation[j]);
              FinalGenEtExc08_radiation.Fill(Et_radiation[j]);
              float weight = Et_radiation[j];
              Energy_Eta_radiation_N8.Fill(Eta_radiation[j], weight);
            }
          }
        }
        */

        /*
        if (ST<5000){
          if (multiplicity==2){
            from1500to3000_ST_N2.Fill(ST);
          }
          if (multiplicity==8){
            from1500to3000_ST_N8.Fill(ST);
          }
        }
        */

        if (ST>4500){
          Hist_sphericity_ST4p5cut.Fill(Sphericity);
          if (multiplicity>2){
            Hist_sphericity_ST4p5cut_incN3.Fill(Sphericity);
          }
          if (multiplicity>3){
            Hist_sphericity_ST4p5cut_incN4.Fill(Sphericity);
          }
          if (multiplicity>4){
            Hist_sphericity_ST4p5cut_incN5.Fill(Sphericity);
          }
          if (multiplicity>5){
            Hist_sphericity_ST4p5cut_incN6.Fill(Sphericity);
          }
          if (multiplicity>6){
            Hist_sphericity_ST4p5cut_incN7.Fill(Sphericity);
          }
          if (multiplicity>7){
            Hist_sphericity_ST4p5cut_incN8.Fill(Sphericity);
          }
        }

        if (debugFlag && ( STMHTnoMET>4000 || Met > 2000 || OurMet > 2000 || fabs(OurMet-Met)>100) && multiplicity>=2) {
          if (debugFlag) cout << "In run number " << runno << " lumi section " << lumiblock << " event number " << evtno << " sT is:" << ST << endl;
          if (dumpIsoInfo) {
            sprintf(messageBuffer, "In run number %d lumi section %d event number %lld ST is %f and multiplicity is %d\n", runno, lumiblock, evtno, ST, multiplicity);
            outTextFile << messageBuffer;
          }
          if (debugFlag) cout << messageBuffer;

          // dump all object info
          for (int j=0; j<25; ++j) {
            if(debugFlag && dumpIsoInfo && JetEt[j]>0.000) {
              sprintf(messageBuffer, "    Jet %d has TightJet=%d Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, TightJets[j],  JetEt[j], JetPx[j], JetPy[j], JetEta[j], JetPhi[j]);
              outTextFile << messageBuffer;
            }
            if (debugFlag) cout  << messageBuffer;
          }
          for (int j=0; j<25; ++j) {
            if(debugFlag && dumpIsoInfo && EleEt[j]>0.000) {
              sprintf(messageBuffer, "    Ele %d has Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, EleEt[j], ElePx[j], ElePy[j], EleEta[j], ElePhi[j]);
              outTextFile << messageBuffer;
            }
            if (debugFlag) cout  << messageBuffer;
          }
          for (int j=0; j<25; ++j) {
            if(debugFlag && dumpIsoInfo && PhEt[j]>0.000) {
              sprintf(messageBuffer, "    Ph %d has Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, PhEt[j], PhPx[j], PhPy[j], PhEta[j], PhPhi[j]);
              outTextFile << messageBuffer;
            }
            if (debugFlag) cout  << messageBuffer;
          }
          for (int j=0; j<25; ++j) {
            if(debugFlag && dumpIsoInfo && MuEt[j]>0.000) {
              sprintf(messageBuffer, "    Mu %d has Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, MuEt[j], MuPx[j], MuPy[j], MuEta[j], MuPhi[j]);
              outTextFile << messageBuffer;
            }
            if (debugFlag) cout  << messageBuffer;
          }
          if (debugFlag && dumpIsoInfo) {
            sprintf(messageBuffer, "    our Px is=%f\n", Px);
            outTextFile << messageBuffer;
            sprintf(messageBuffer, "    our Py is=%f\n", Py);
            outTextFile << messageBuffer;
            sprintf(messageBuffer, "    our MHT is=%f\n", OurMet);
            outTextFile << messageBuffer;
            sprintf(messageBuffer, "    MET is=%f\n", Met);
            outTextFile << messageBuffer;
            if (debugFlag) cout  << messageBuffer;
            sprintf(messageBuffer, "\n\n\n\n");
            outTextFile << messageBuffer;
          }
        }
        nDumpedEvents+=1;
        if (debugFlag && nDumpedEvents==eventsToDump) break;
  }
  cout << "Number of events in chain is: " << N_failpassBadChCand << endl;
  // write output textfile
  outTextFile.close();
  // write output root file
  TFile* outRootFile = new TFile(outFilename.c_str(), "RECREATE");
  outRootFile->cd();
  outRootFile->mkdir("ST");
  outRootFile->mkdir("ST_tight");
  outRootFile->mkdir("MET-MHT");
  outRootFile->mkdir("Isolation");

  outRootFile->cd("ST");
  stHist.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHist[iHist]->Write();
    Gen71stExcHist[iHist]->Write();
    GenInteractstExcHist[iHist]->Write();
    GenUnderlystExcHist[iHist]->Write();
    GenJetstExcHist[iHist]->Write();
    GenJetNoISRstExcHist[iHist]->Write();
    stIncHist[iHist]->Write();
    stExcHist_2leadjets[iHist]->Write();
    stExcHist_noISR[iHist]->Write();
  }
  Hist2D_N_ST->Write();
  stHistMHT.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHistMHT[iHist]->Write();
    stIncHistMHT[iHist]->Write();
  }
  //JetEtaExc08_3500Up.Write();
  JetPtExc08_3500Up.Write();
  //JetEtaExc08_2300to2400.Write();
  JetPtExc08_2300to2400.Write();
  //JetEtaExc03_1800to1900.Write();
  //JetPtExc03_1800to1900.Write();
  FinalGenPDGIDExc08_3500Up.Write();
  FinalGenPDGIDExc08_2300to2400.Write();
  //FinalGenPDGIDExc03_1800to1900.Write();
  FinalGenEtExc08_3500Up.Write();
  FinalGenEtExc08_2300to2400.Write();
  //FinalGenEtExc03_1800to1900.Write();
  FinalGenEtaExc08_3500Up.Write();
  FinalGenEtaExc08_2300to2400.Write();
  //FinalGenEtaExc03_1800to1900.Write();
  FinalGenNExc08_3500Up.Write();
  FinalGenNExc08_2300to2400.Write();
  //FinalGenNExc03_1800to1900.Write();
  //FinalGenSTExc03.Write();
  //FinalGenSTExc08.Write();
  //SameGenST_Gen_PDGID_N2.Write();
  //SameGenST_Gen_ET_N2.Write();
  //SameGenST_Gen_Eta_N2.Write();
  //SameGenST_Gen_Phi_N2.Write();
  //SameGenST_Gen_PDGID_N8.Write();
  //SameGenST_Gen_ET_N8.Write();
  //SameGenST_Gen_Eta_N8.Write();
  //SameGenST_Gen_Phi_N8.Write();
  FinalGenST71Exc02.Write();
  FinalGenST71Exc03.Write();
  FinalGenST71Exc08.Write();
  //from1500to3000_ST_N2.Write();
  //from1500to3000_ST_N8.Write();
  FinalGenN71Exc02.Write();
  FinalGenPDGID71Exc02.Write();
  FinalGenEta71Exc02.Write();
  FinalGenEt71Exc02.Write();
  FinalGenPDGID45Exc02.Write();
  FinalGenEta45Exc02.Write();
  FinalGenEt45Exc02.Write();
  FinalGenN71Exc08.Write();
  FinalGenPDGID71Exc08.Write();
  FinalGenEta71Exc08.Write();
  FinalGenEt71Exc08.Write();
  FinalGenPDGID45Exc08.Write();
  FinalGenEta45Exc08.Write();
  FinalGenEt45Exc08.Write();
  FinalGenSTExc02_interaction.Write();
  FinalGenPDGIDExc02_interaction.Write();
  FinalGenEtaExc02_interaction.Write();
  FinalGenEtExc02_interaction.Write();
  FinalGenSTExc02_radiation.Write();
  FinalGenPDGIDExc02_radiation.Write();
  FinalGenEtaExc02_radiation.Write();
  FinalGenEtExc02_radiation.Write();
  FinalGenSTExc08_interaction.Write();
  FinalGenPDGIDExc08_interaction.Write();
  FinalGenEtaExc08_interaction.Write();
  FinalGenEtExc08_interaction.Write();
  FinalGenSTExc08_radiation.Write();
  FinalGenPDGIDExc08_radiation.Write();
  FinalGenEtaExc08_radiation.Write();
  FinalGenEtExc08_radiation.Write();
  FinalGenSTExc03_interaction.Write();
  FinalGenSTExc03_radiation.Write();
  Energy_Eta_radiation_N2.Write();
  Energy_Eta_radiation_N8.Write();
  Energy_Eta_interaction_N2.Write();
  Energy_Eta_interaction_N8.Write();
  ST_Sphericity_p05top2_in.Write();
  ST_Sphericity_p05top2_out.Write();
  ST_Sphericity_p05top2_in_N3.Write();
  ST_Sphericity_p05top2_in_N4.Write();
  ST_Sphericity_p05top2_in_N5.Write();
  ST_Sphericity_p05top2_in_N6.Write();
  ST_Sphericity_p05top2_in_N7.Write();
  ST_Sphericity_p05top2_in_N8.Write();
  ST_Sphericity_p05top2_in_N9.Write();
  ST_Sphericity_p05top2_in_N10.Write();
  ST_Sphericity_p05top2_in_N11.Write();
  ST_Sphericity_p05top2_out_N3.Write();
  ST_Sphericity_p05top2_out_N4.Write();
  ST_Sphericity_p05top2_out_N5.Write();
  ST_Sphericity_p05top2_out_N6.Write();
  ST_Sphericity_p05top2_out_N7.Write();
  ST_Sphericity_p05top2_out_N8.Write();
  ST_Sphericity_p05top2_out_N9.Write();
  ST_Sphericity_p05top2_out_N10.Write();
  ST_Sphericity_p05top2_out_N11.Write();
  Hist_sphericity.Write();
  Hist_sphericity_ST4p5cut.Write();
  Hist_sphericity_ST4p5cut_incN3.Write();
  Hist_sphericity_ST4p5cut_incN4.Write();
  Hist_sphericity_ST4p5cut_incN5.Write();
  Hist_sphericity_ST4p5cut_incN6.Write();
  Hist_sphericity_ST4p5cut_incN7.Write();
  Hist_sphericity_ST4p5cut_incN8.Write();
  ST_Sphericity_p1less.Write();
  ST_Sphericity_p1less_N3.Write();
  ST_Sphericity_p1less_N4.Write();
  ST_Sphericity_p1less_N5.Write();
  ST_Sphericity_p1less_N6.Write();
  ST_Sphericity_p1less_N7.Write();
  ST_Sphericity_p1less_N8.Write();
  ST_Sphericity_p1less_N9.Write();
  ST_Sphericity_p1more.Write();
  ST_Sphericity_p1more_N3.Write();
  ST_Sphericity_p1more_N4.Write();
  ST_Sphericity_p1more_N5.Write();
  ST_Sphericity_p1more_N6.Write();
  ST_Sphericity_p1more_N7.Write();
  ST_Sphericity_p1more_N8.Write();
  ST_Sphericity_p1more_N9.Write();


  outRootFile->cd("ST_tight");
  stHist_tight.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHist_tight[iHist]->Write();
    stIncHist_tight[iHist]->Write();
  }
  stHistMHT_tight.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHistMHT_tight[iHist]->Write();
    stIncHistMHT_tight[iHist]->Write();
  }
  FinalGenST71Exc02.Write();
  FinalGenST71Exc03.Write();
  FinalGenST71Exc08.Write();

  outRootFile->cd("MET-MHT");
  METvsMHT.Write();
  METvsMHTinc2.Write();
  METvsMHTinc2hasMuon.Write();
  METvsMHTinc2hasPhoton.Write();
  METvsMHTinc2hasElectron.Write();
  METvsMHTinc2onlyJets.Write();
  METvsMHT_tight.Write();
  METvsMHTinc2_tight.Write();
  METvsMHTinc2hasMuon_tight.Write();
  METvsMHTinc2hasPhoton_tight.Write();
  METvsMHTinc2hasElectron_tight.Write();
  METvsMHTinc2onlyJets_tight.Write();
  METoverSumET.Write();
  METoverSumETinc2.Write();
  METoverSumETinc2hasMuon.Write();
  METoverSumETinc2hasPhoton.Write();
  METoverSumETinc2hasElectron.Write();
  METoverSumETinc2onlyJets.Write();
  METoverSumET_tight.Write();
  METoverSumETinc2_tight.Write();
  METoverSumETinc3_tight.Write();
  METoverSumETinc4_tight.Write();
  METoverSumETinc5_tight.Write();
  METoverSumETinc6_tight.Write();
  METoverSumETinc7_tight.Write();
  METoverSumETinc8_tight.Write();
  METoverSumETinc9_tight.Write();
  METoverSumETinc10_tight.Write();
  METoverSumETinc2hasMuon_tight.Write();
  METoverSumETinc2hasPhoton_tight.Write();
  METoverSumETinc2hasElectron_tight.Write();
  METoverSumETinc2onlyJets_tight.Write();

  outRootFile->cd("Isolation");
  MuonJetIso1.Write();
  MuonJetIso2.Write();
  MuonJetIso3.Write();
  MuonJetIso4.Write();
  MuonJetoverlapdR1.Write();
  MuonJetoverlapdR2.Write();
  MuonJetoverlapdR3.Write();
  MuonJetoverlapdR4.Write();
  ElectronJetIso1.Write();
  ElectronJetIso2.Write();
  ElectronJetIso3.Write();
  ElectronJetIso4.Write();
  ElectronJetoverlapdR1.Write();
  ElectronJetoverlapdR2.Write();
  ElectronJetoverlapdR3.Write();
  ElectronJetoverlapdR4.Write();
  PhotonJetIso1.Write();
  PhotonJetIso2.Write();
  PhotonJetIso3.Write();
  PhotonJetIso4.Write();
  PhotonJetoverlapdR1.Write();
  PhotonJetoverlapdR2.Write();
  PhotonJetoverlapdR3.Write();
  PhotonJetoverlapdR4.Write();

  std::cout << " End of event "<<std::endl;

  outRootFile->Close();
}



// function to calculate dR between two objects
float dR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt( ( eta1 - eta2 )*( eta1 - eta2 ) + std::pow(TMath::ATan2(TMath::Sin( phi1 - phi2), TMath::Cos(phi1-phi2)),2) );
}


// function to make an event list object for MET filtering
std::map<unsigned, std::set<unsigned> > readEventList(char const* _fileName) {
  std::map<unsigned, std::set<unsigned> > list;
  ifstream listFile(_fileName);
  if (!listFile.is_open())
    throw std::runtime_error(_fileName);

  unsigned iL(0);
  std::string line;
  while (true) {
    std::getline(listFile, line);
    if (!listFile.good())
      break;

    if (line.find(":") == std::string::npos || line.find(":") == line.rfind(":"))
      continue;

    unsigned run(std::atoi(line.substr(0, line.find(":")).c_str()));
    unsigned event(std::atoi(line.substr(line.rfind(":") + 1).c_str()));

    list[run].insert(event);

    ++iL;
  }

  std::cout << "Loaded " << iL << " events" << std::endl;

  return list;
}
