#include "ExpSetup.h"
#include <iostream>
#include "TLine.h"
#include "TVirtualPad.h"

void ExpSetup::Analysis()
{
    //ScaleAndShiftTimes();
    for (int i = 0; i < NofDetectors; ++i)
    {
        Detectors.at(i)->InvertY();
        if(Detectors.at(i)->type>0) Detectors.at(i)->ConvertToTR();
    }
    //Dynamically find the end for the baseline calculation region
    //for(int i = 0; i < NofDetectors; ++i)
    //{
    //}
    //===
    for (int i = 0; i < NofDetectors; ++i)
    {
        baseline_region_end = gNsample;
        Detectors.at(i)->FindGlobalMaximum(0, Detectors.at(i)->waveform_y.size() - 1);
        if (Detectors.at(i)->global_maximum.position < baseline_region_end)
            baseline_region_end = Detectors.at(i)->global_maximum.position;
        //std::cout<<"NofDet="<<i<<",NofPoint="<<Detectors.at(i)->waveform_y.size()<<",baseline_region_end="<<baseline_region_end-200<<std::endl;

        baseline_region_end -= 15./ghorizontal_interval; // assume risetime<5ns
        if (baseline_region_end < 0)
        {

            baseline_region_end = 100.;
            if (baseline_region_end > Detectors.at(i)->global_maximum.position)
            {
                baseline_region_end = 50.;
                if (baseline_region_end > Detectors.at(i)->global_maximum.position)
                    baseline_region_end = 25.;
            }
        }
        record_blregion_end[i] = baseline_region_end;
        max_region_end = 2000 + baseline_region_end;

        Detectors.at(i)->SubstractBaseline(baseline_region_end);

        Detectors.at(i)->FindGlobalMaximum(baseline_region_end, max_region_end);
        Detectors.at(i)->FindStartPoint(baseline_region_end);
        Detectors.at(i)->FindEndPoint(max_region_end);
        Detectors.at(i)->FindInvertMaximum(baseline_region_end, Detectors.at(i)->global_maximum.position);
        Detectors.at(i)->CalculateCharges();
        Detectors.at(i)->FindNaiveTiming();
        Detectors.at(i)->FindRiseTime(); //0.2-0.8 Leading Edge

        Detectors.at(i)->TimeInformation();
    }
}

void ExpSetup::Dump(int id)
//void TestBeamSetup::Dump(int id)
{
    TFile dumpfile("CheckYourWaveform.root", "update");
    //TTree *newtree = OutTree->CloneTree(0);
    //newtree->Fill();

    TObject integer;

    integer.SetUniqueID(Detectors.size());
    integer.Write("n");

    TGraph gr;
    TGraph gMP;  // main peak
    TGraph gCTP; //crosstalk peak
    for (int i = 0; i < Detectors.size(); ++i)
    {
        char str1[20];
        sprintf(str1, "CH%dgraph_%d", Channel_IDs.at(i), id);
        gr = TGraph(Detectors.at(i)->waveform_x.size(), &Detectors.at(i)->waveform_x[0], &Detectors.at(i)->waveform_y[0]);
        gr.Write(str1);

        gMP.SetPoint(0, Detectors.at(i)->global_maximum.x, Detectors.at(i)->global_maximum.y);
        gCTP.SetPoint(0, Detectors.at(i)->invert_maximum.x, Detectors.at(i)->invert_maximum.y);
        //sprintf(str1,"graph_%d",i,id);
        gMP.SetMarkerSize(2);
        gMP.SetMarkerStyle(32);
        gMP.SetMarkerColor(3);
        gCTP.SetMarkerSize(2);
        gCTP.SetMarkerStyle(46);
        gCTP.SetMarkerColor(2);

       

        gr.Draw("AL");
        gMP.Draw("Psame");
        gCTP.Draw("Psame");
        TLine *linebl = new TLine(Detectors.at(i)->waveform_x[0], 0, Detectors.at(i)->waveform_x.at(record_blregion_end[i]), 0);
        linebl->SetLineColor(2);
        linebl->SetLineWidth(2);
        linebl->Draw();
        TLine *lineQ = new TLine(Detectors.at(i)->start_point.x, Detectors.at(i)->global_maximum.y, Detectors.at(i)->end_point.x, Detectors.at(i)->global_maximum.y);
        lineQ->SetLineColor(2);
        lineQ->SetLineWidth(2);
        lineQ->Draw();
        char str2[20];
        sprintf(str2, "CH%dC_%d", Channel_IDs.at(i), id);
        gPad->Write(str2);
    }

    //Detectors.at(0)->pre_filter_backup->Sigmoid.fit_func.Write("rawsigmoid");
    //newtree->Write();
    dumpfile.Close();
}

void ExpSetup::init(std::vector<int> channel_ids)
{
    Channel_IDs = channel_ids;
    NofDetectors = Detectors.size();
    //max_region_end = 20001;
    init_tree();
}

void ExpSetup::init_tree()
{
    //OutFile = new TFile("outrootfile.root","recreate");
    OutTree = new TTree("Pico", "Analysis Output");
    int mcp_no = 0;
    int tr_no = 0;
    std::cout << "===List of different types of characteristics===" << std::endl;
    OutTree->Branch("Baseline_Window", &baseline_region_end, "Baseline_Window/I");
    for (int i = 0; i < Detectors.size(); ++i)
    {
        Detector *det;
        int type = Detectors.at(i)->type;
        std::string typestr;
        if(type>0){
            tr_no = Channel_IDs.at(i);
            char str[20];
            sprintf(str,"TR%d_",tr_no);
            typestr = str;
        }
        else{

        mcp_no = Channel_IDs.at(i);
        char str[20];
        sprintf(str, "MCP%d_", mcp_no);
        typestr = str;
        }

        det = Detectors.at(i);
        //std::cout << "mcp" << i << std::endl;

        std::string varname;
        std::string leafname;

        varname = typestr + "baseline_level";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->baseline_level, leafname.c_str());

        varname = typestr + "baseline_rms";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->baseline_rms, leafname.c_str());

        varname = typestr + "global_maximum_y";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->global_maximum.y, leafname.c_str());

        varname = typestr + "global_maximum_x";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->global_maximum.x, leafname.c_str());

        varname = typestr + "start_x";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->start_point.x, leafname.c_str());

        varname = typestr + "end_x";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->end_point.x, leafname.c_str());

        varname = typestr + "invert_maximum_x";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->invert_maximum.x, leafname.c_str());

        varname = typestr + "invert_maximum_y";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->invert_maximum.y, leafname.c_str());

        varname = typestr + "all_charge";
        leafname = varname + "[4]/D";
        OutTree->Branch(varname.c_str(), det->charge_all, leafname.c_str());

        varname = typestr + "rise_time";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->rise_time, leafname.c_str());

        varname = typestr + "CFDtime";
        leafname = varname + "[8]/D";
        OutTree->Branch(varname.c_str(), det->CFDtime, leafname.c_str());

        varname = typestr + "CFDfrac";
        leafname = varname + "[8]/D";
        OutTree->Branch(varname.c_str(), det->CFDfrac, leafname.c_str());

        varname = typestr + "CFDfailed";
        leafname = varname + "[8]/O";
        OutTree->Branch(varname.c_str(), det->CFDfailed, leafname.c_str());

        varname = typestr + "LEDtime";
        leafname = varname + "[14]/D";
        OutTree->Branch(varname.c_str(), det->LEDtime, leafname.c_str());

        varname = typestr + "LEDthrd";
        leafname = varname + "[14]/D";
        OutTree->Branch(varname.c_str(), det->LEDthrd, leafname.c_str());

        varname = typestr + "LEDfailed";
        leafname = varname + "[14]/O";
        OutTree->Branch(varname.c_str(), det->LEDfailed, leafname.c_str());

        std::cout << typestr << std::endl;
    }
}
