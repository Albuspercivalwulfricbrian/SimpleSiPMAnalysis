#ifndef CHANNEL_ENTRY_H
#define CHANNEL_ENTRY_H
#include <TTree.h>
#include "constants.h"
#include "PronyFitter.h"
using namespace std;
// const int MAX_N_SAMPLES = 2048;

struct IntegralInfo
{
    Short_t signal_length = 0;
    Short_t npeaks = 0;
    Int_t end_amplitude = 0;
    void Initialize();
};

struct FitResults
{
    Float_t fit_charge;
    Float_t fit_amp;
    Float_t fit_time;
    
};

struct HarmonicsStruct
{
    Float_t first_harmonic, second_harmonic, fit_chi2, fit_R2, fit_charge, charge;

    void Initialize()
    {
        first_harmonic = 0.;
        second_harmonic = 0.;
        fit_chi2 = 999;
        fit_R2 = 999;
        charge = 0.;
        fit_charge = 0.;
    }
};

struct short_energy_ChannelEntry
{
    Float_t charge;
    Float_t time;
    UShort_t amp;
    Float_t zl;
    Float_t zl_rms;
    IntegralInfo II;
    FitResults fit;
    Short_t nCoincidencePeaks;
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
    Short_t dwf[MAX_N_SAMPLES] = {0};
    Int_t fZlLeft = 0;
    Int_t fZlRight = 200;
    Float_t zl = 0;
    IntegralInfo II;
    FitResults fit;
    Int_t amp = 0;
    Short_t peak_position = 0;
    Int_t fGATE_BEG = 1000000;
    Int_t fGATE_END = -1000000;
    public:
    static TString GetChName(Int_t channel_num);
    Int_t SetBranch(TTree *tree, Int_t channel_num);
    void Initialize();
    void SplineWf();
    void CalculateDiffWf();
    void AssumeSmartScope();
    void SetBoarders(Int_t,Int_t);
    void FindDiffWfPars(Short_t &min_diff, Short_t &min_time, Short_t &max_diff, Short_t &max_time);
    void Set_Zero_Level(int);
    void Set_Zero_Level_Area(Int_t i);
    Short_t CountCoincidencePeaks(Int_t, Int_t);
    Int_t PointAmpl(Int_t);
    Float_t CalculateZlwithNoisePeaks(int);
    Int_t Get_Zero_Level();
    Float_t Get_Zero_Level_RMS();
    Float_t Get_Charge();
    Short_t Get_time();
    Float_t Get_time_gauss();
    UShort_t Get_Amplitude();
    IntegralInfo GetIntegralInfo();
    FitResults GetFitResults();
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




























// #include <TTree.h>
// #include "ChannelEntry.h"
#include <iostream>

using namespace std;

    void IntegralInfo::Initialize()
    {
        signal_length = 0;
        npeaks = 0;
        end_amplitude = 0;    
    }

    void FitResults::Initialize()
    {
        fit_amp = 0;
        fit_charge = 0;
        fit_time = 0;
    }

    void short_energy_ChannelEntry::Initialize()
    {
        charge = 0.;
        time = 0.;
        amp = 0; 
        zl_rms = 0.;
        zl = 0.;
        nCoincidencePeaks = 0;
        II.Initialize();
    }


    TString short_energy_ChannelEntry::GetChName(Int_t channel_num)
    {
	    return TString::Format("channel_%i", channel_num);
    }

    Int_t short_energy_ChannelEntry::SetBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
    


    TString diff_short_energy_ChannelEntry::GetChName(Int_t channel_num)
    {
	    return TString::Format("diff_channel_%i", channel_num);
    }
    TBranch* diff_short_energy_ChannelEntry::CreateBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->Branch(GetChName(channel_num).Data(), this, "min_diff/S:min_diff_time/S:max_diff/S:max_diff_time/S");
    }

    void diff_short_energy_ChannelEntry::Initialize()
    {
        min_diff = 0;
        min_diff_time = 0;
        max_diff = 0;
        max_diff_time = 0;
    }    


    TString ChannelEntry::GetChName(Int_t channel_num)
    {
	    return TString::Format("channel_%i", channel_num);
    }

    Int_t ChannelEntry::SetBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
    
    void ChannelEntry::Initialize()
    {
        for (int i = 0; i < sizeof(wf)/sizeof(wf[0]); i++) {wf[i] = 0; dwf[i] = 0;}
        wf_size = 0;
    }

    void ChannelEntry::SplineWf()
    {
        Float_t wf1[MAX_N_SAMPLES] = {0};
        const Int_t SplineWidth = 2;
        for (Int_t i = 0; i < wf_size; i++)
        {
            Int_t il=i-SplineWidth; Int_t ir=i+SplineWidth;
            if (il<0) il=0;
            if (ir>wf_size-1) ir=wf_size-1;
            Float_t counter = 0;
            for (Int_t in = il; in <=ir; in++) {wf1[i]+=wf[in];counter++;}
            wf1[i]/=counter;
        }
        for (Int_t i = 0; i < wf_size; i++) wf[i] = wf1[i];
    }

    Int_t ChannelEntry::PointAmpl(Int_t i)
    {
        Int_t v = zl - (Int_t)wf[i];
        return v;

    }
    Short_t ChannelEntry::CountCoincidencePeaks(Int_t threshold1, Int_t threshold2)
    {
        Short_t nCoincidencePeaks=0;
        Int_t epsilon1 = 10;
        Int_t epsilon2 = 100;
        // Short_t wf1[MAX_N_SAMPLES] = {0};
        for (Int_t i = fGATE_BEG; i < fGATE_END; i++)
        {
            Int_t v = PointAmpl(i);
            Int_t vl = PointAmpl(i-1);
            Int_t vll = PointAmpl(i-2);
            Int_t vr = PointAmpl(i+1);
            Int_t vrr = PointAmpl(i+2);

            if ( (v-vr > epsilon1) && (v-vl > epsilon1) && (v-vrr > epsilon2) && (v-vll > epsilon2) && ( v > threshold1 && v < threshold2)) nCoincidencePeaks++;

        }
        // for (Int_t i = 0; i < wf_size; i++) wf[i] = dwf[i];
        return nCoincidencePeaks;
    } 

    void ChannelEntry::CalculateDiffWf()
    {
        const Float_t Diff_window = 4;
        // Short_t wf1[MAX_N_SAMPLES] = {0};
        for (Int_t i = 0; i < wf_size; i++)
        {
            Int_t il=i-Diff_window; Int_t ir=i+Diff_window;
            if (il<0) il=0;
            if (ir>wf_size-1) ir=wf_size-1;
            dwf[i]=(Short_t)((Float_t)(wf[ir]-wf[il])/(Float_t)(ir-il));
        }
        // for (Int_t i = 0; i < wf_size; i++) wf[i] = dwf[i];
    } 

    Float_t ChannelEntry::CalculateZlwithNoisePeaks(int a)
    {
        CalculateDiffWf();
        Float_t sum = 0;
        Float_t counter = 0;
        for (int s=fZlLeft+1; s < fZlRight; ++s) {

            if (abs(dwf[s]) < a && abs(dwf[s-1]) < a && abs(dwf[s+1]) < a) {sum += wf[s]; counter++;}
        }
        zl = sum/counter;
    return zl;
    }

    void ChannelEntry::AssumeSmartScope()
    {
        fGATE_BEG = peak_position;
        fGATE_END = peak_position;
        
        while (1)
        {
            fGATE_BEG--;
            if (fGATE_BEG < 0) {fGATE_BEG++; break;} 
            if (wf[fGATE_BEG] > zl) break;
        }
        while (1)
        {
            fGATE_END++;
            if (fGATE_END >= wf_size) {fGATE_END--; break;} 
            // if (wf[fGATE_END] > zl) break;
            if (wf[fGATE_END] > zl) break;
        }
    }

    void ChannelEntry::SetBoarders(Int_t BEG, Int_t END)
    {
        fGATE_BEG = BEG;
        fGATE_END = END;
    }
    void ChannelEntry::FindDiffWfPars(Short_t &min_diff, Short_t &min_time, Short_t &max_diff, Short_t &max_time)
    {
        for (Short_t s=fGATE_BEG; s < fGATE_END; ++s) {
            Short_t v = wf[s];
            if (v < min_diff) 
            {
                min_diff = v;
                min_time = 16*s;
            }
            if (v > max_diff) 
            {
                max_diff = v;
                max_time = 16*s;
            }      
        }
    }

    void ChannelEntry::Set_Zero_Level_Area(Int_t i)
    {
        fZlLeft = 0;
        fZlRight = i;
    }
    
    void ChannelEntry::Set_Zero_Level(int EZL)
    {
      zl = EZL;
    }

    Int_t ChannelEntry::Get_Zero_Level()
    {
        const Int_t interv_num = 1;
        int zero_lvl = 0;
        int best_spread = -1;
        for (int i=0; i < interv_num; ++i) {
            int vmin = numeric_limits<int>::max();
            int vmax = numeric_limits<int>::min();
            int sum = 0;
            for (int s=fZlRight/interv_num * i; s < fZlRight/interv_num * (i+1); ++s) {
                int v = wf[s]; 
                sum += v;
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
            }
            int spread = vmax - vmin;
            if (best_spread < 0) best_spread = spread;
            if (spread <= best_spread) {
                best_spread = spread;
                zero_lvl = sum / (fZlRight/interv_num);
            }
            zl = zero_lvl;
        }
        return zero_lvl;
    }
    Float_t ChannelEntry::Get_Zero_Level_RMS()
    {
        const Int_t interv_num = 1;
        Float_t best_spread = -1;
        Float_t rms_zl = -1;
        for (Int_t i=0; i < interv_num; ++i) {
            Int_t vmin = numeric_limits<int>::max();
            Int_t vmax = numeric_limits<int>::min();
            Float_t sum = 0; Float_t sumsquare = 0; Float_t sum_counter = 0;
            for (Int_t s=fZlRight/interv_num * i; s < fZlRight/interv_num * (i+1); ++s) {
                Int_t v = wf[s]; 
                sum += (Float_t)v;
                sum_counter++;
            }
            sum /=sum_counter;
            sumsquare = 0.;
            
            for (Int_t s=fZlRight/interv_num * i; s < fZlRight/interv_num * (i+1); ++s) {
                sumsquare += (Float_t)(wf[s] - sum)*(wf[s] - sum)/sum_counter;
            }            
            rms_zl = sqrt(sumsquare);
            if (best_spread < 0) best_spread = rms_zl;

        }
        return best_spread;
    }

    Float_t ChannelEntry::Get_Charge()
    {
        Float_t gateInteg = 0;
        {
            //for (int s=180; s <= 680; ++s) {

             for (int s=fGATE_BEG; s <= fGATE_END; ++s) {
                //if (zl < (Float_t)wf[s] && s > peak_position+7) {II.signal_length = s - fGATE_BEG; II.end_amplitude = (float)zl - (float)wf[s];break;}
                gateInteg +=  ((Float_t)zl - (Float_t)wf[s]) ;
            }                
        }

        return gateInteg;
    }
    
    IntegralInfo ChannelEntry::GetIntegralInfo()
    {
        return II;
    }

    FitResults ChannelEntry::GetFitResults()
    {
        return fit;
    }

    Short_t ChannelEntry::Get_time()
    {
        amp = 0;
        peak_position = 0;
        if (wf_size == 0) 
        {
            amp = 0;
            peak_position = 0;
            return peak_position;
        }
        for (int s=fGATE_BEG; s < fGATE_END; ++s) {
            Int_t v = zl - wf[s];
            if (v > amp) {
                amp = v;
                peak_position = s;
            }
        }
        // if ( amp > 0) cout << amp << endl;
        return peak_position;
    }
    Float_t ChannelEntry::Get_time_gauss()
    {
        if (wf_size == 0) return 0;
        Float_t peak_search = 0.;
        Float_t ampl_sum = 0;
        for (Int_t s= fGATE_BEG; s < fGATE_END; ++s) {
            Int_t v = zl - wf[s];
            if (v > amp*0.1)
            {
                ampl_sum += (Float_t) v;
                peak_search+= (Float_t) v*s;
            }
        }
        peak_search /= ampl_sum;
    return 16.0*peak_search;
    }
////////////
    UShort_t ChannelEntry::Get_Amplitude()
    {
        return (UShort_t)amp;
    }

    void ChannelEntry::FillWf(Short_t *Ewf)
    {
        for (int i = 0; i < SNAPSHOT_LENGTH; i++)
        {
            wf[i] = Ewf[i];
            wf_size++;
            // cout << wf[i] << endl;
        }
    }
//##################

    Int_t mini_tree_nrg::Initialize()
    {
        EdepIntermediate0 = 0;
        EdepIntermediate1 = 0;

        EdepScat0 = 0;
        EdepScat1 = 0;
        EdepDet0 = 0;
        EdepDet1 = 0;  
        DetNum0 = 0;
        DetNum1 = 0; 
        EventType = -10;
        return 0;
    }

    Int_t mini_tree_time::Initialize()
    {
        TimeScat0 = 0;
        TimeScat1 = 0;
        TimeDet0 = 0;
        TimeDet1 = 0;
        TimeIntermediate0 = 0;
        TimeIntermediate1 = 0;

        return 0;
    }



















#endif CHANNEL_ENTRY_H
