#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>

#include <TBenchmark.h>

#include "FileReader.h"
#include "Detector.h"
#include "ExpSetup.h"


int main(int argc, char *argv[])
{
    //1st input argument is is directory path with files
    //2nd input argument is number of files for each channel in the directory
    //3rd input argument is name the output root file will have
    //4th input argument (optional), the number of the first channel
    //5th input argument (optional), if it is 1 then dumping will be enabled

    TBenchmark bench;
    bench.Start("full");

    int start_id = 0;
    bool enable_dumping = false;
    int dumping_id = 0;
    if (argc > 4)
    {
        start_id = atoi(argv[4]);
    }

    if (argc > 5)
    {
        //if(atoi(argv[4]) == 1)
        enable_dumping = true;
        dumping_id = atoi(argv[5]);
    }
    //Specify the channels written in the files
    int Nch = atoi(argv[2]);
    ExpSetup mysetup;
    std::vector<int> channel_IDs;
    for(int i = start_id; i<Nch+start_id;i++)
    {
    channel_IDs.push_back(i); 
    mysetup.CreateMCP();
    }
    channel_IDs.push_back(100); 
    mysetup.CreateTR();
    //channel_IDs.push_back(101); 
    //mysetup.CreateTR();
    //return 0;
    FileReader myfile(channel_IDs, argv[1], "wave");
    
    mysetup.init(channel_IDs);

    myfile.SetDetectorSetup(mysetup);

    const char *output_rootfile_name = argv[3];

    ////////////////////////////////////
    std::cout << std::endl;
    int i = 0;

    while (myfile.GetNextEvent())
    {

        printf("\rEvent is: %.5d ===============================================\n", i); //magic

        mysetup.Analysis();

        if ( enable_dumping && i  == dumping_id) 
            mysetup.Dump(i);
     
        mysetup.Fill_Tree();
        i++;
        //if(i==10) break;
    }
    mysetup.Finalize_Tree(output_rootfile_name);

    bench.Show("full");
    //system("pause");
    return 0;
}
/*
void analyze()
{
    //1st input argument is is directory path with files
    //2nd input argument is number of files for each channel in the directory
    //3rd input argument is name the output root file will have
    //4th input argument (optional), if it is 1 then dumping will be enabled
    int argc=3;
    char argv[5][100]={"","../binarytest/","1","test.root",""};
    TBenchmark bench;
    bench.Start("full");

    bool enable_dumping = false;
    int dumping_id = 0;

    if (argc > 4)
    {
        //if(atoi(argv[4]) == 1)
        enable_dumping = true;
        dumping_id = atoi(argv[4]);
    }
    //Specify the channels written in the files
    int Nch = atoi(argv[2]);
    ExpSetup mysetup;
    std::vector<int> channel_IDs;
    for(int i = 1; i<Nch;i++)
    {
        std::cout<<"wave:"<< i<<std::endl;
    channel_IDs.push_back(i); 
    mysetup.CreateMCP();
    }
    
    FileReader myfile(channel_IDs, argv[1], "wave");

    
    mysetup.init();

    myfile.SetDetectorSetup(mysetup);

    const char *output_rootfile_name = argv[3];

    ////////////////////////////////////
    std::cout << std::endl;
    int i = 0;

    while (myfile.GetNextEvent())
    {

        printf("\rEvent is: %.5d ===============================================", i); //magic

        mysetup.Analysis();

        if ( enable_dumping && i  == dumping_id) 
            mysetup.Dump(i);
     
        mysetup.Fill_Tree();
        i++;
        if(i==10) break;
    }
    mysetup.Finalize_Tree(output_rootfile_name);

    bench.Show("full");
}
*/