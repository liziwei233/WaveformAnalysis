#ifndef __detector_h__
#define __detector_h__

#include <vector>

#include <TMath.h>
#include <TVirtualFFT.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "global.h"
#include "DetectorSetup.h"

struct WaveformPoint
{
    double y;
    double x;
    int position;
};

struct TimingInfo
{
    double y;
    double x;
    double slope;
    double intersect;
    double timing;
    double parameters[4];//pol3 parameters
    bool failed = 0;
};



struct SigmoidInfo
{
    TF1 fit_func;
    double parameters[4];
    bool failed;
    double chisquare;
    int degrees_freedom;
};

class Detector
{
    friend class DetectorSetup;
    //friend class MuonSetup;
    //friend class LaserSetup;
    //friend class CalibrationSetup;
    //friend class SimulationSetup;
    friend class ExpSetup;

    public:
    Detector();
    ~Detector();
    void SetWaveY(std::vector<float> wave);
    void SetWaveX(std::vector<float> wave);
    void SetMCP();
    void SetTR();
    void InvertY();
    void ConvertToTR();
    void SubstractBaseline(int base_region_end);
    void FindGlobalMaximum(int start, int end);
    void FindInvertMaximum(int start, int end);
    void FindSecondInvertPeak(int start);
    void FindFirstPeak(int start, int end);
    void ConvertFirstPeak2GlobalMaximum();
    void FindStartPoint(int start);
    void FindEndPoint(int start);
    void CalculateCharges();
    void FindNaiveTiming();
    void FindeightypercentTiming();
    double linear_interpolation(double x1, double x2, double y1, double y2, double y);
    bool FitPol3(double* x, double* y, double* fit_parameters);
    void LineFitLeastSquares(double *data_x, double *data_y, int data_n, std::vector<double> &vResult);
    void TimeInformation();
    TimingInfo Time(double fac,int partype);
    TimingInfo Time_linear(double fac,int partype,int Npoint);
    WaveformPoint FindTimingPoint(double fac,int partype);
    void TimeSigmoid();
    void FindRiseTime();
    void FindWidth();

    private:
    int type;

    std::vector<float> waveform_y;
    std::vector<float> waveform_x;


    double baseline_level;
    double baseline_rms;

    WaveformPoint global_maximum;
    WaveformPoint firstpeak;
    WaveformPoint invert_maximum;
    WaveformPoint SecondInvertPeak;
    WaveformPoint start_point;
    WaveformPoint end_point;
    WaveformPoint naive_point;
    WaveformPoint eightypercent_naive;

    double charge_all[4];

    double naive_time;//ns
    double eightypercent_naive_time;//ns
    double rise_time;// defined from 20% to 80% height of the pulse
    double width;

    SigmoidInfo Sigmoid;
    
    double CFDtime[8];
    double CFDfrac[8];
    bool CFDfailed[8];

    double LEDtime[14];
    double LEDthrd[14];
    bool LEDfailed[14];

    TimingInfo CFD;
    TimingInfo LED;

};

double fermi_dirac(double *x, double *par);

#endif
