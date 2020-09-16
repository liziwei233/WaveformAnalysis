    char path[1024] = "/mnt/f/XOPtest/4Anode/crosstalk/A3";
void checkdata(const char* name=""){
    ifstream infile;
    char buff[1024];
    float wave[1024];
    sprintf(buff,"%s/%s",path,name);
    infile.open(buff, ios::binary);
    if (!infile)
    {
        std::cerr << "Unknown format " << std::endl
                  << std::endl;
        return ;
    }
    
    float amp = 0;
    int waveNum=0;
    while (1)
    {
        if (infile.eof())
        {
            std::cout << waveNum << std::endl;
            return;
        }
        else
        {
            for(int i=0;i<1024;i++)
            {

                
            infile.read((char *)&amp, sizeof(float));
            wave[i]=amp;
            if(amp!=0&& amp<1e-15) {
                std::cout << waveNum << std::endl;
            return;
            }
            if(waveNum==9873) cout<<i<<"\t"<<amp<<endl;
            }
            waveNum++;
        }
    }
    std::cout << waveNum << std::endl;
    
    infile.close();
}