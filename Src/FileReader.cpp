#include "FileReader.h"

FileReader::FileReader(std::vector<int> channel_ids, const char *dir_path, const char *file_midname)
{
    for (int i = 0; i < channel_ids.size(); ++i)
    {
        if (channel_ids.at(i) >= 100 && channel_ids.at(i) < 200)
        {
            Channel C = Channel(dir_path, "TR_0", channel_ids.at(i) - 100,1);
            Channels.push_back(C);
        }
        else if (channel_ids.at(i) >= 200)
        {
            Channel C = Channel(dir_path, "TR_1", channel_ids.at(i) - 200,1);
            Channels.push_back(C);
        }
        else
        {

            Channel C = Channel(dir_path, file_midname, channel_ids.at(i),0);
            Channels.push_back(C);
        }
    }
}

void FileReader::OpenTriggerChannel(const char *dir_path, const char *file_midname, int ch_id)
{
    trigger_channel_is_open = 1;
    TriggerChannel = new Channel(dir_path, file_midname, ch_id,1);
}

bool FileReader::GetNextEvent()
{
    for (int i = 0; i < Channels.size(); ++i)
    {
        if (Channels.at(i).GetNextEvent())
        {
            std::cout << "End of events" << std::endl;
            return 0;
        }
        else if (trigger_channel_is_open)
        {
            TriggerChannel->GetNextEvent();
        }

        Detectors->SetDetectorWaveform(i, Channels.at(i).GetWaveformX(), Channels.at(i).GetWaveformY());
    }
    return 1;
}
void FileReader::SetDetectorSetup(DetectorSetup &det)
{
    Detectors = &det;
}