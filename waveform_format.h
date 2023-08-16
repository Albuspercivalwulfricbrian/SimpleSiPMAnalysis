
#ifndef WF_FORMAT
#define WF_FORMAT
//#include <iostream>
#include <string>
#include <sstream>
#include "constants.h"



class DataFileReader
{
    private:
    string fileName;
    std::ifstream file;
    public:
    Short_t wf[TOTAL_CHANNELS][SNAPSHOT_LENGTH];
    DataFileReader(string Name)
    {
        fileName = Name;
        file.open(Name,std::ifstream::binary);
    }
    ~DataFileReader(){};

    void MoveToEventWf(Int_t eventnum)
    {
        file.seekg(HEADER_SIZE+(POINT_SIZE_IN_BYTES*(TOTAL_CHANNELS)*SNAPSHOT_LENGTH+SNAPSHOT_PREFIX)*eventnum+SNAPSHOT_PREFIX);
    }

    void DemonstrateWf(Int_t eventnum)
    {
        char b;
        file.seekg(HEADER_SIZE+(TOTAL_CHANNELS*SNAPSHOT_LENGTH+SNAPSHOT_PREFIX)*eventnum);
        for (int i = 0; i < 30; i++)
        {

            file.read(reinterpret_cast<char *>(&b), sizeof(b));
            cout << b;
        } 

        cout <<endl; 
        short a =-10009;

        for (int cc = 0; cc < 3;cc++)
        {
        MoveToEventWf(eventnum+cc);

            for (int i = 0; i < TOTAL_CHANNELS*SNAPSHOT_LENGTH; i++)
            {
                file.read(reinterpret_cast<char *>(&a), sizeof(a));
                cout << a << " ";
            }     
            cout <<endl; 
            for (int i = 0; i < 30; i++)
            {

                file.read(reinterpret_cast<char *>(&b), sizeof(b));
                cout << b;
            }    
            cout <<endl; 
        }

    }

    void GetEntry(Int_t eventnum)
    {
        MoveToEventWf(eventnum);

        Short_t a = -10009;
        for (int i = 0; i < TOTAL_CHANNELS*SNAPSHOT_LENGTH; i++)
        {
            file.read(reinterpret_cast<char *>(&a), sizeof(a));
            if (i%2==0) wf[0][i/2] = a;
            else wf[1][i/2] = a;
        }  
    }
    short GetTotalEvents()
    {
        short events = 0;
        return events;
    }
};

#endif WF_FORMAT