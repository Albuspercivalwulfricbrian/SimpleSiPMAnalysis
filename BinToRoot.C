#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "waveform_format.h"
#include "ChannelEntry.h"
const Int_t H_SIZE = 7168;
const Int_t channels_number = 2;

void BinToRoot(TString source_file = "/home/admin/HOLODILNIK/HAMAMATSU/LED/T26/HV42_led2.5V")
{
    // string pathName = "/media/doc/DATA/SiPM_low_energy_detector/2_days_cooling with good reflection/";
    //TString source_file = gSystem->GetFromPipe("yad --file-selection");
    //TString source_file = "/home/admin/Bochka_cross_tolk_kitayskiy/T145_32V_1024";
    DataFileReader *dataFile = new DataFileReader((string)source_file);
    ChannelEntry wf_data[TOTAL_CHANNELS];

    TFile* RootDataFile = new TFile ((TString)source_file + ".root", "RECREATE");
    TTree* RootDataTree = new TTree ("adc64_data","adc64_data");
    for (int i = 0; i < channels_number; i++) wf_data[i].Initialize();
    for (int i = 0; i < channels_number; i++) RootDataTree->Branch((Form("channel_%i", i)),&wf_data[i], "wf_size/I:wf[wf_size]/S");

    const int calculate_events = 1000000; 
    for (int i = 0; i < calculate_events; i++)
    {
        if (i >  calculate_events) break;
        //dataFile->MoveToEventWf(i);
        dataFile->GetEntry(i);
        for (int j = 0; j < channels_number; j++)
        {
            wf_data[j].Initialize();
            if (sizeof(dataFile->wf[j])!=0)
            {
                wf_data[j].FillWf(dataFile->wf[j]);
            }
            // cout << wf_data[j].wf_size << endl;
            // cout << endl;
        }
        RootDataTree->Fill();
        if (i % 1000 == 0)cout << (float)i/(float)calculate_events*100 <<"%" << endl;
        if (wf_data[0].wf[0]==-10009 && wf_data[1].wf[0]==-10009 ) break;
    }
    RootDataTree->Write();
    RootDataFile->Close();
}
        



