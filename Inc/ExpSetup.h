#ifndef __expsetup_h__
#define __expsetup_h__

#include "DetectorSetup.h"
#include "global.h"
class ExpSetup : public DetectorSetup
{
    public:
        ~ExpSetup();
        void Analysis();
        void CreateAverageTools();
        void SetWaveformToAverage();
        void Finalize_AverageTools(const char* outfile_name);
        void Dump(int id);
        void init(std::vector<int> channel_ids);
        void init_tree();
    protected:
        int record_blregion_end[100];
        std::vector<int> Channel_IDs;
        std::vector<AverageTool*> avers;
};
#endif
