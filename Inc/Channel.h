#ifndef __channel_h
#define __channel_h__

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>

#include <string.h>

#include "global.h"
#include "Detector.h"
#include "DetectorSetup.h"

#pragma pack(push, 1)

#pragma pack(pop)

class Channel
{
    public:
        Channel(const char* dir, const char* file_midname,int ID);
        ~Channel();
        //Channel(const Channel&) = delete;
        std::vector<float> GetWaveformX(); 
        std::vector<float> GetWaveformY(); 
        int OpenFile();
        bool GetNextEvent();
        
    private:
        std::vector<float> wave;

        std::vector<float> waveform_x;
        std::vector<float> waveform_y;

        std::string directory_path;
        float amp;
        int event_number;
        int ReadBinary(std::ifstream &infile,std::vector<float> * wave);
        std::string filename;
        std::string full_path;

};



#endif
