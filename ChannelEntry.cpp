
#include <TTree.h>
#include "ChannelEntry.h"
#include <iostream>

using namespace std;

    void IntegralInfo::Initialize()
    {
        signal_length = 0;
        npeaks = 0;
        end_amplitude = 0;    
    }
    void short_energy_ChannelEntry::Initialize()
    {
        charge = 0.;
        time = 0.;
        amp = 0; 
        zl_rms = 0.;
        zl = 0.;
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

             for (int s=fGATE_BEG; s <= fGATE_END; ++s) {
                if (zl < (Float_t)wf[s] && s > peak_position+7) {II.signal_length = s - fGATE_BEG; II.end_amplitude = (float)zl - (float)wf[s];break;}
                gateInteg +=  ((Float_t)zl - (Float_t)wf[s]) ;
            }                
        }

        return gateInteg;
    }
    
    IntegralInfo ChannelEntry::GetIntegralInfo()
    {
        return II;
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
