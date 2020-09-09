#ifndef __filereader_h__
#define __filereader_h__

#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <string.h>

#include "Detector.h"
#include "DetectorSetup.h"
#include "Channel.h"

class FileReader
{
    public:
        FileReader(std::vector<int> channel_ids, const char* dir_path, const char* file_midname);
        bool GetNextEvent();
        void SetDetectorSetup(DetectorSetup &det);
        void OpenTriggerChannel(const char* dir_path, const char* file_midname,int ch_id);
    private:
        std::vector<Channel> Channels;
        bool trigger_channel_is_open = 0;
        Channel *TriggerChannel;
        DetectorSetup* Detectors;
};
#endif




