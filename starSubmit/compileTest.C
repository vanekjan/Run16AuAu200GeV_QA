#include "TString.h"
#include "TROOT.h"

void compileTest(const char* name = "runPicoDpmAnaMaker.C"){ //orig. runPicoHFMyAnaMaker.C
  Long_t ret= gROOT->ProcessLine(".L StRoot/macros/loadSharedHFLibraries.C");
  cout << ret << endl;
  loadSharedHFLibraries();
  ret = gROOT->ProcessLine(Form(".L StRoot/macros/%s++", name));
}
