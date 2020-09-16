#include "Channel.h"

Channel::Channel(const char *dir, const char *file_midname, int ID)
{
    directory_path = std::string(dir);
    if (directory_path.back() != '/')
        directory_path.push_back('/');
    if (gBinaryFlag)
        filename = std::string(file_midname) + std::string("_") + std::to_string(ID) + std::string(".dat");
    else
        filename = std::string(file_midname) + std::string("_") + std::to_string(ID) + std::string(".txt");

    wave.reserve(1e9);
    //waveform_x.reserve(4002);
    //waveform_y.reserve(4002);

    event_number = 0;

    full_path = directory_path + filename;
    OpenFile();
    std::cout <<"\t \t "<< full_path << std::endl;
}

Channel::~Channel()
{
}

int Channel::OpenFile()
{
    wave.clear();
    std::ifstream infile;
    infile.open(full_path.c_str(), std::ifstream::binary);
    if (!infile)
    {
        std::cerr << "Unknown format " << std::endl
                  << std::endl;
        return 1;
    }
    else if(ReadBinary(infile, &wave)){
        std::cout << "Read file is over: " << std::endl;
            return 0;
    }
    infile.close();
    return 0;
}
int Channel::ReadBinary(std::ifstream &infile, std::vector<float> *wave)
{

    amp = 0;

    while (1)
    {
        if (infile.eof())
        {
            return 1;
        }
        else
        {

            infile.read((char *)&amp, sizeof(float));
            if(amp!=0&&amp<1e-15) amp=0;
            wave->push_back(amp);
            //std::cout << ", waveform_y= " << wave->back() << std::endl;
        }
    }

    return 0;
}

std::vector<float> Channel::GetWaveformX() { return waveform_x; }
std::vector<float> Channel::GetWaveformY() { return waveform_y; }

bool Channel::GetNextEvent()
{
    int waveform_size = gNsample;
    if(event_number >= wave.size()/waveform_size) 
    {
            return 1;
    }

    waveform_x.clear();
    waveform_y.clear();

    for (int i = 0; i < waveform_size; ++i)
    {
        waveform_x.push_back(i * ghorizontal_interval + horizontal_offset);
        waveform_y.push_back(wave.at(i+event_number*waveform_size) * gvertical_gain_ch + vertical_offset);
    }

    event_number++;
    return 0;
}