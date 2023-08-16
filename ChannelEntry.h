#ifndef CHANNEL_ENTRY_H
#define CHANNEL_ENTRY_H
#include<TTree.h>
#include"constants.h"
using namespace std;
// const int MAX_N_SAMPLES = 2048;

struct IntegralInfo
{
    Short_t signal_length = 0;
    Short_t npeaks = 0;
    Int_t end_amplitude = 0;
    void Initialize();
};

struct short_energy_ChannelEntry
{
    Float_t charge;
    Float_t time;
    UShort_t amp;
    Float_t zl;
    Float_t zl_rms;
    IntegralInfo II;
    static TString GetChName(Int_t channel_num);
    TBranch* CreateBranch(TTree *tree, Int_t channel_num);
    Int_t SetBranch(TTree *tree, Int_t channel_num);
    void Initialize();
};

struct diff_short_energy_ChannelEntry
{
    Short_t min_diff;
    Short_t min_diff_time;
    Short_t max_diff;
    Short_t max_diff_time;
    static TString GetChName(Int_t channel_num);
    TBranch* CreateBranch(TTree *tree, Int_t channel_num);
    void Initialize();
};

class ChannelEntry {

    public:
    Int_t wf_size;
    Short_t wf[MAX_N_SAMPLES];

    private:
    Int_t fZlLeft = 0;
    Int_t fZlRight = 20;
    Float_t zl;
    IntegralInfo II;
    Int_t amp = 0;
    Short_t peak_position;
    Int_t fGATE_BEG = 1000000;
    Int_t fGATE_END = -1000000;
    public:
    static TString GetChName(Int_t channel_num);
    Int_t SetBranch(TTree *tree, Int_t channel_num);
    void Initialize();
    void SplineWf();
    void DiffWf();
    void AssumeSmartScope();
    void SetBoarders(Int_t,Int_t);
    void FindDiffWfPars(Short_t &min_diff, Short_t &min_time, Short_t &max_diff, Short_t &max_time);
    void Set_Zero_Level(int);
    void Set_Zero_Level_Area(Int_t i);
    Int_t Get_Zero_Level();
    Float_t Get_Zero_Level_RMS();
    Float_t Get_Charge();
    Short_t Get_time();
    Float_t Get_time_gauss();
    UShort_t Get_Amplitude();
    IntegralInfo GetIntegralInfo();
    void FillWf(Short_t *Ewf);

};

//##################
//################
struct mini_tree_nrg
{
    Float_t EdepIntermediate0;
    Float_t EdepIntermediate1;
    Float_t EdepScat0;
    Float_t EdepScat1;
    Float_t EdepDet0;
    Float_t EdepDet1;  
    Short_t DetNum0;
    Short_t DetNum1;
    Short_t EventType; 
    Int_t SetBranch(TTree *tree);
    Int_t Initialize();
};

struct mini_tree_time
{
    Float_t TimeIntermediate0;
    Float_t TimeIntermediate1;
    Float_t TimeScat0;
    Float_t TimeScat1;
    Float_t TimeDet0;
    Float_t TimeDet1;
    Int_t SetBranch(TTree *tree);
    Int_t Initialize();
};

#endif CHANNEL_ENTRY_H
