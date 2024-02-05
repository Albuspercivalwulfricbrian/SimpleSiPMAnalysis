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

void Calculate_fit_harmonics(event_fit_struct &result_fit_event, Int_t &event_counter, Int_t event_num, Int_t ch_num, Int_t gate_beg, Int_t gate_end, Int_t gate_maximum_beg, Int_t gate_maximum_end, TString source_path, TString run_name, Float_t ampl_scale, TObjArray *check_fit_arr)
{

    Int_t wf_size = ChannelEntry[channel].wf_size;

    Int_t ch_iter = ch_num;

    result_fit_event.reset();
    if(samples_data == NULL) return;

    
    Float_t *scale_samples = new Float_t[wf_size];
    for(Int_t isaml = 0; isaml < wf_size; isaml++)
	scale_samples[isaml] = ampl_scale * (float)samples_data[ch_iter][isaml];

    //Zero level calculation
    const int n_gates = 3;
    int gate_npoints = (int)floor((gate_beg-2.)/n_gates);


    Float_t gates_mean[n_gates], gates_rms[n_gates];
    for(Int_t igate = 0; igate < n_gates; igate++)
	MeanRMScalc(scale_samples, gates_mean+igate, gates_rms+igate, igate*gate_npoints, (igate+1)*gate_npoints);

    int best_gate = 0;
    for(Int_t igate = 0; igate < n_gates; igate++)
	if(gates_rms[igate] < gates_rms[best_gate]) best_gate = igate;

    Float_t zero_level_ = gates_mean[best_gate];
    Short_t MAX_in_gate_ = -32760;
    Int_t time_max_in_gate_ = 0;

    for(UShort_t sample_curr = 0; sample_curr < wf_size; sample_curr++)
    {
	Float_t val_sample = scale_samples[sample_curr];	
	if((sample_curr >= gate_maximum_beg) && (sample_curr <= gate_maximum_end))
	{
	    if(val_sample > MAX_in_gate_)
	    {
		MAX_in_gate_ = (Short_t)val_sample;
		time_max_in_gate_ = sample_curr;
	    }
	}
    }

    Double_t charge_ = 0.;

    for(UShort_t sample_curr = 0; sample_curr < wf_size; sample_curr++)
    {
	Float_t val_sample = scale_samples[sample_curr];
	
	if((sample_curr >= gate_beg) && (sample_curr <= gate_end)) //in of gate
	{
	    charge_ += (val_sample-zero_level_);
	}
	}
    
    if(MAX_in_gate_ - zero_level_ > 0)
    {
	Bool_t IsDebug = false;
	if(IsDebug) printf("event for harmonic %i event counter %i\n", event_num, event_counter);

	// Calculating timing
	Float_t Amplitude = MAX_in_gate_ - zero_level_;
	Float_t trsh_03 = zero_level_ + Amplitude*0.3;
	Float_t trsh_09 = zero_level_ + Amplitude*0.9;

	Int_t point = time_max_in_gate_;
	Float_t front_time_beg_03 = GoToLevel(scale_samples, trsh_03, &point, -1, wf_size);

	point = time_max_in_gate_;
	Float_t front_time_end = GoToLevel(scale_samples, trsh_09, &point, -1, wf_size);

	//Fitting


	PronyFitter Pfitter(2, 2, gate_beg, gate_end, wf_size);
	//Pfitter.SetDebugMode(1);
	Pfitter.SetWaveform(scale_samples, zero_level_);
	int SignalBeg = Pfitter.CalcSignalBegin(front_time_beg_03, front_time_end);
	Int_t best_signal_begin = Pfitter.ChooseBestSignalBeginHarmonics(SignalBeg-1, SignalBeg+2);
	double *harmonics;
	if(best_signal_begin>gate_beg)
	{
		Pfitter.SetSignalBegin(best_signal_begin);
		Pfitter.CalculateFitHarmonics();
		Pfitter.CalculateFitAmplitudes();
		harmonics = Pfitter.GetHarmonics();
		Float_t fit_charge = Pfitter.GetIntegral(gate_beg, gate_end);
		Float_t fit_chi2 = Pfitter.GetChiSquare(gate_beg, gate_end, time_max_in_gate_);
		Float_t fit_R2 = Pfitter.GetRSquare(gate_beg, gate_end);

		if(fit_R2<0.02)
		{
			event_counter++;
			result_fit_event.first_harmonic = (Float_t)harmonics[1];
			result_fit_event.second_harmonic = (Float_t)harmonics[2];
			result_fit_event.fit_chi2 = fit_chi2;
			result_fit_event.fit_R2 = fit_R2;
			result_fit_event.charge = (Float_t)charge_;
			result_fit_event.fit_charge = fit_charge;
			TString signal_name = Form("ch_num %i fit_charge %.0f integral %.0f chi2 %.1f fit_R2 %.3f", 
								ch_iter, fit_charge, charge_, fit_chi2, fit_R2);
			if(event_counter<5) Pfitter.DrawFit(check_fit_arr, signal_name);
		}


	}
    }

    delete[] scale_samples;

}