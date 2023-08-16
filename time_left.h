#include <iostream>
#include <ctime>
using namespace std;

void DisplayTimeToCalculate(Int_t EvNum, Int_t total_entries, time_t start_time)
{
    if (EvNum%10000 == 0) 
    {
        std::cout<< u8"\033[2J\033[1;1H"; 
        std::cout << (Float_t)EvNum/(Float_t)total_entries*100 << "%" << std::endl;
        time_t time_left = (time(NULL)-start_time)*(float)(total_entries-EvNum)/(float)(EvNum);
        std::cout << " time left: ";
        if (time_left/3600 > 0) cout << time_left/3600 <<"h ";
        if ((time_left%3600)/60 > 0 || time_left/3600 == 0) cout << (time_left%3600)/60 << "m ";
        cout << (time_left%3600)%60<< "s " <<std::endl;
    }        
}

