#include "ExpSetup.h"
#include <iostream>
#include "TLine.h"
#include "TVirtualPad.h"

ExpSetup::~ExpSetup()
{
    if (!avers.empty())
    {

        for (int i = 0; i < avers.size(); ++i)
            delete avers.at(i);
    }
}

void ExpSetup::Analysis()
{
    //ScaleAndShiftTimes();
    for (int i = 0; i < NofDetectors; ++i)
    {
        Detectors.at(i)->InvertY();
        //if (Detectors.at(i)->type == 0)
        //    Detectors.at(i)->ConvertToTR();
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

        baseline_region_end -= 15. / ghorizontal_interval; // assume risetime<5ns
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
        //Detectors.at(i)->FindGlobalMaximum(baseline_region_end, max_region_end);
        //Detectors.at(i)->FindStartPoint(baseline_region_end);
        //Detectors.at(i)->FindEndPoint(max_region_end);
        //Detectors.at(i)->FindRiseTime(); //0.2-0.8 Leading Edge

        Detectors.at(i)->FindFirstPeak(baseline_region_end, max_region_end);
        Detectors.at(i)->ConvertFirstPeak2GlobalMaximum();
        Detectors.at(i)->FindStartPoint(baseline_region_end);
        Detectors.at(i)->FindEndPoint(max_region_end);
        Detectors.at(i)->CalculateCharges();
        //Detectors.at(i)->FindNaiveTiming();
        //Detectors.at(i)->FindeightypercentTiming();
        Detectors.at(i)->FindRiseTime(); //0.2-0.8 Leading Edge
        Detectors.at(i)->FindWidth();
        Detectors.at(i)->FindInvertMaximum(baseline_region_end, Detectors.at(i)->global_maximum.position);
        Detectors.at(i)->FindSecondInvertPeak(Detectors.at(i)->global_maximum.position);
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
    TGraph gRP;  // ringing peak
    for (int i = 0; i < Detectors.size(); ++i)
    {
        char str1[20];
        sprintf(str1, "CH%dgraph_%d", Channel_IDs.at(i), id);
        gr = TGraph(Detectors.at(i)->waveform_x.size(), &Detectors.at(i)->waveform_x[0], &Detectors.at(i)->waveform_y[0]);
        gr.Write(str1);

        gMP.SetPoint(0, Detectors.at(i)->global_maximum.x, Detectors.at(i)->global_maximum.y);
        gCTP.SetPoint(0, Detectors.at(i)->invert_maximum.x, Detectors.at(i)->invert_maximum.y);
        gRP.SetPoint(0, Detectors.at(i)->SecondInvertPeak.x, Detectors.at(i)->SecondInvertPeak.y);
        //sprintf(str1,"graph_%d",i,id);
        gMP.SetMarkerSize(2);
        gMP.SetMarkerStyle(32);
        gMP.SetMarkerColor(3);
        gCTP.SetMarkerSize(2);
        gCTP.SetMarkerStyle(46);
        gCTP.SetMarkerColor(2);
        gRP.SetMarkerSize(2);
        gRP.SetMarkerStyle(41);
        gRP.SetMarkerColor(2);

        gr.Draw("AL");
        gMP.Draw("Psame");
        gCTP.Draw("Psame");
        gRP.Draw("Psame");
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
    //dumpfile->Delete();
}

void ExpSetup::init(std::vector<int> channel_ids)
{
    //TFile create("CheckYourWaveform.root","recreate");
    //create.Close();
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
        if (type == 0)
        {
            tr_no = Channel_IDs.at(i);
            char str[20];
            sprintf(str, "TR%d_", tr_no);
            typestr = str;
        }
        else
        {

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

        varname = typestr + "secondinvertpeak_x";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->SecondInvertPeak.x, leafname.c_str());

        varname = typestr + "secondinvertpeak_y";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->SecondInvertPeak.y, leafname.c_str());

        varname = typestr + "all_charge";
        leafname = varname + "[4]/D";
        OutTree->Branch(varname.c_str(), det->charge_all, leafname.c_str());

        varname = typestr + "rise_time";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->rise_time, leafname.c_str());

        varname = typestr + "width";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(), &det->width, leafname.c_str());

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
void ExpSetup::CreateAverageTools()
{
    for (int i = 0; i < NofDetectors; ++i)
    {
        avers.push_back(new AverageTool());
    }
}

void ExpSetup::SetWaveformToAverage()
{
    //ScaleAndShiftTimes();
    for (int i = 0; i < NofDetectors; ++i)
    {
        Detectors.at(i)->InvertY();
    }
    int j = 0;
    double ref_time[50] = {0};
    double n = 0;
    double normalization = 0;
    float yfactor = 0;
    double cut[10] = {9e5, 9e5, 9e5, 9e5, 9e5};
    double qcut[10] = {-9e5, -9e5, -9e5, -9e5};
    int ntrue = 0;
    for (int i = 0; i < NofDetectors; ++i)
    {
        baseline_region_end = gNsample;
        normalization = 1;
        Detectors.at(i)->FindGlobalMaximum(0, Detectors.at(i)->waveform_y.size() - 1);
        if (Detectors.at(i)->global_maximum.position < baseline_region_end)
            baseline_region_end = Detectors.at(i)->global_maximum.position;

        baseline_region_end -= 15. / ghorizontal_interval;
        if (baseline_region_end <= 0)
        {

            baseline_region_end = 100.;
            if (baseline_region_end > Detectors.at(i)->global_maximum.position)
            {
                baseline_region_end = 50.;
                if (baseline_region_end > Detectors.at(i)->global_maximum.position)
                    baseline_region_end = 25.;
            }
        }
        max_region_end = 2000 + baseline_region_end;
        record_blregion_end[i] = baseline_region_end;
        //===
        //

        Detectors.at(i)->SubstractBaseline(baseline_region_end);
        if (Detectors.at(i)->type == 0)
            Detectors.at(i)->FindGlobalMaximum(baseline_region_end, max_region_end);
        else
        {
            Detectors.at(i)->FindFirstPeak(baseline_region_end, max_region_end);
            Detectors.at(i)->ConvertFirstPeak2GlobalMaximum();
        }

        Detectors.at(i)->FindStartPoint(baseline_region_end);
        Detectors.at(i)->FindEndPoint(max_region_end);
        //Detectors.at(i)->FindElectronPeakEndPoint();
        Detectors.at(i)->CalculateCharges();
        Detectors.at(i)->FindNaiveTiming();
        //Detectors.at(i)->FindRiseTime();
        //Detectors.at(i)->FindFirstPeak();
        //Detectors.at(i)->FindMaxDerivative();

        //Detectors.at(i)->TimeTwentyPercent();
        //Detectors.at(i)->TimeInflection();
        //ref_time += Detectors.at(i)->Inflection.timing;
        
        ref_time[i] = Detectors.at(i)->naive_point.x;
        if (ref_time[i] <= 0 || ref_time[i] > 4e3)
        ref_time[i] = 15.;
        if (Detectors.at(i)->global_maximum.y < cut[i] && Detectors.at(i)->charge_all[0] > qcut[i])
        {
            ntrue++;
        }
    }
    std::cout << "the delta ref_time= " << ref_time[1]- ref_time[0]<< std::endl;
    if (ntrue > NofDetectors - 1)
    {
        for (int i = 0; i < NofDetectors; ++i)
        {

            //if (std::abs(Detectors.at(i)->global_maximum.y) > 0.5e-3)
            //    normalization = 1. / Detectors.at(i)->global_maximum.y;
            //avers.at(i)->initial();

            //std::cout<<" global maximum y &cut ="<<Detectors.at(i)->global_maximum.y<<"\t"<<cut[Channel_IDs.at(i)]<<std::endl;
            if (Detectors.at(i)->type == 0)
                yfactor = gvertical_gain_tr[Channel_IDs.at(i) - 100];
            else
                yfactor = gvertical_gain_ch[Channel_IDs.at(i)];
            avers.at(i)->SetWaveform(Detectors.at(i)->waveform_x, Detectors.at(i)->waveform_y, ref_time[i], yfactor, baseline_region_end);
            avers.at(i)->StandardAverage();
        }
    }

    //ref_time/= n;
    //if (ref_time > 1.e3)
    //    ref_time = 15.;

    //double normalization = 1. / Detectors.at(j)->global_maximum.y;
    //for(int i = 0; i < Detectors.at(j)->waveform_y.size(); ++i)
    //{
    //    normalization += Detectors.at(j)->waveform_y.at(i);
    //}
    //normalization /= Detectors.at(j)->waveform_y.size();

    //aver.SetWaveform(Detectors.at(j)->waveform_x, Detectors.at(j)->waveform_y, ref_time, normalization, baseline_region_end);
}
void ExpSetup::Finalize_AverageTools(const char *outfile_name)
{
    //char str1[200];
    TFile *f = new TFile(outfile_name, "recreate");
    f->Close();
    delete f;
    char ch_name[1024];
    for (int i = 0; i < avers.size(); ++i)
    {
        //sprintf(str1, "CH%d%s.root", Channel_IDs.at(i),outfile_name);
        avers.at(i)->Finalize();
        if (Detectors.at(i)->type == 0)
            sprintf(ch_name, "TR%d", Channel_IDs.at(i));
        else
            sprintf(ch_name, "CH%d", Channel_IDs.at(i));
        avers.at(i)->Write(outfile_name, ch_name);
    }
}