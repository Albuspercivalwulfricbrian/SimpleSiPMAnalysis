#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TObjString.h>
#include "ChannelEntry.h"
#include "Coeffs.h"
#include <stdlib.h>
#include "time_left.h"
const int GATE_BEG = 200;
///////////////////////Define fit constants range
void GAGG_shit(TString source_file = "/home/admin/HOLODILNIK/HAMAMATSU/LED/T26/HV42_led2.5V.root")
{
    TString prefix = "_analyzed";
    // TString source_file = gSystem->GetFromPipe("yad --file-selection");

    // std::filesystem::path fsp = source_file.Data();
    string fullname = source_file.Data();
    size_t lastindex = (fullname).find_last_of("."); 
    TString rawname = (TString)fullname.substr(0, lastindex); 

    const Int_t TOTAL_CHANNELS = 2;
    ChannelEntry channel_info[TOTAL_CHANNELS];
    short_energy_ChannelEntry short_channel_info[TOTAL_CHANNELS];
    diff_short_energy_ChannelEntry diff_short_channel_info[TOTAL_CHANNELS];
    Coeffs_struct cal_coeff[TOTAL_CHANNELS];
    TFile *combined_root = new TFile (rawname+prefix+".root", "RECREATE");
    TTree *combined_tree = new TTree ("adc64_data","adc64_data");
    TTree *diff_tree = new TTree ("diff_tree","diff_tree");
    for (Int_t channel_number = 0; channel_number < TOTAL_CHANNELS; channel_number++)
        combined_tree->Branch(TString::Format("channel_%i", channel_number),&short_channel_info[channel_number]);
    TChain *PMT_tree = new TChain;
    PMT_tree->AddFile( (source_file + "/adc64_data").Data() );
    for (Int_t channel_number = 0; channel_number < TOTAL_CHANNELS; channel_number++)
        // (channel_info[channel_number]).SetBranch(PMT_tree,channel_number);
        PMT_tree->SetBranchAddress(Form("channel_%i",channel_number),&channel_info[channel_number]);

    Int_t total_entries = PMT_tree->GetEntries();
    time_t start_time = time(NULL);
    for (Long_t entry_num = 0;  entry_num< total_entries; entry_num++)
    {
        PMT_tree->GetEntry(entry_num);
        for (Int_t channel_number = 0; channel_number < TOTAL_CHANNELS; channel_number++)
        {
            channel_info[channel_number].SplineWf();
            short_channel_info[channel_number].Initialize();
            channel_info[channel_number].Set_Zero_Level_Area(GATE_BEG);
            //Int_t zero_level = channel_info[channel_number].Get_Zero_Level();
            int zero_level = 32000; channel_info[channel_number].Set_Zero_Level(zero_level);
            //short_channel_info[channel_number].zl_rms = channel_info[channel_number].Get_Zero_Level_RMS();
            
            //channel_info[channel_number].SplineWf();
            Int_t pp = channel_info[channel_number].Get_time();
            //channel_info[channel_number].AssumeSmartScope(); 
            channel_info[channel_number].SetBoarders(610,660); 
            short_channel_info[channel_number].amp = channel_info[channel_number].Get_Amplitude();
            short_channel_info[channel_number].charge = channel_info[channel_number].Get_Charge();
            short_channel_info[channel_number].time = channel_info[channel_number].Get_time_gauss();
            short_channel_info[channel_number].zl = zero_level;
            short_channel_info[channel_number].II = channel_info[channel_number].GetIntegralInfo();
        }
        combined_tree->Fill();
        DisplayTimeToCalculate(entry_num,total_entries, start_time);
        if ((entry_num + 1) % 100000 == 0) cout << (float)entry_num/total_entries*100<<"%"<<endl;
        for (int ch = 0; ch < TOTAL_CHANNELS; ch++) channel_info[ch].Initialize();
    }

    delete PMT_tree;
    combined_tree->Write();
    combined_root->Close();
}
        



