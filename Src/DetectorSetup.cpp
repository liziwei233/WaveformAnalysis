#include "DetectorSetup.h"

DetectorSetup::~DetectorSetup()
{
    for(int i = 0 ; i < Detectors.size(); ++i)
        delete Detectors.at(i);
}

void DetectorSetup::CreateMCP()
{
    Detector* MCP = new Detector();
    MCP->SetMCP();
    Detectors.push_back(MCP);
}
void DetectorSetup::CreateTR()
{
    Detector* TR = new Detector();
    TR->SetTR();
    Detectors.push_back(TR);
}

void DetectorSetup::Fill_Tree()
{
    OutTree->Fill();
}

void DetectorSetup::Finalize_Tree(const char* outrootfile)
{
    OutFile = new TFile(outrootfile,"recreate");
    OutTree->Write();
    OutFile->Write();
    OutFile->Close();
    delete OutFile;
    delete OutTree;
}

void DetectorSetup::SetDetectorWaveform(int i, std::vector<float> x, std::vector<float> y)
{
    Detectors.at(i)->SetWaveX(x);
    Detectors.at(i)->SetWaveY(y);
}

