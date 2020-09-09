#ifndef __detectorsetup_h__
#define __detectorsetup_h__

#include <string>
#include <vector>

#include <TFile.h>
#include <TVectorD.h>
#include <TObject.h>
#include <TMath.h>
#include <TTree.h>

#include "Detector.h"
class Detector;

class DetectorSetup
{
    public:
        ~DetectorSetup();
        void CreateMCP();
        void CreateTR();
        void SetDetectorWaveform(int i, std::vector<float> x, std::vector<float> y);
        //void DumpOld();
        void Fill_Tree();
        void Finalize_Tree(const char* outrootfile);
    protected:
        std::vector<Detector*> Detectors;
        int baseline_region_end;
        int max_region_end;
        int NofDetectors;

        TFile* OutFile;//("outrootfile.root","recreate");
        TTree* OutTree;//("pico","Analysis Output");
};
#endif
