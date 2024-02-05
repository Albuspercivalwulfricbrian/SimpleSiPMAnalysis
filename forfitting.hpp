void BinDataReader::Calculate_waveform(event_data_struct &result_event, Int_t ch_num, Int_t gate_beg, Int_t gate_end, Int_t gate_maximum_beg, Int_t gate_maximum_end, Float_t ampl_scale, Bool_t IsFIT, Double_t first_fit_harmonic, Double_t second_fit_harmonic, TString FIT_QA_mode, TObjArray *check_fit_arr, Float_t *fitQA_arg_arr)
{

  if(readdebug) printf("Calculating wave form -------------------------------------------\n");
  if(readdebug) printf("h_num: %i; gate_beg: %i; gate_end %i; ampl_scale %i; IsFIT %i\n", ch_num, gate_beg, gate_end, ampl_scale, IsFIT);

    Int_t ch_iter = ch_num;

    result_event.reset();
    if(samples_data == NULL) return;

    
    Float_t *scale_samples = new Float_t[sample_total];
    for(Int_t isaml = 0; isaml < sample_total; isaml++)
	scale_samples[isaml] = ampl_scale * (float)samples_data[ch_iter][isaml];
    
    Int_t samples_in_gate = 0;
    Int_t samples_out_gate = 0;

    //Zero level calculation
    const int n_gates = 3;
    int gate_npoints = (int)floor((gate_beg-2.)/n_gates);


    Float_t gates_mean[n_gates], gates_rms[n_gates];
    for(Int_t igate = 0; igate < n_gates; igate++)
	MeanRMScalc(scale_samples, gates_mean+igate, gates_rms+igate, igate*gate_npoints, (igate+1)*gate_npoints);

    int best_gate = 0;
    for(Int_t igate = 0; igate < n_gates; igate++)
	if(gates_rms[igate] < gates_rms[best_gate]) best_gate = igate;

    result_event.zero_level = gates_mean[best_gate];
    result_event.zero_level_RMS = gates_rms[best_gate];

    
    //PASS 1
    result_event.mean_in_gate = 0.;
    result_event.mean_out_gate = 0.;
    result_event.MAX_in_gate = -32760;
    result_event.MIN_in_gate = 32767;
    result_event.MAX_out_gate = -32760;
    result_event.MIN_out_gate = 32767;


    
    for(UShort_t sample_curr = 0; sample_curr < sample_total; sample_curr++)
    {
	Float_t val_sample = scale_samples[sample_curr];
	
	//if(sample_curr < gate_beg) //out of gate
	if((sample_curr < gate_beg) || (sample_curr > gate_end)) //out of gate
	{
	    result_event.mean_out_gate += val_sample;
	    
	    if((Short_t)val_sample < result_event.MIN_out_gate)result_event.MIN_out_gate =(Short_t) val_sample;
	    if((Short_t)val_sample > result_event.MAX_out_gate)result_event.MAX_out_gate = (Short_t) val_sample;
	    
	    samples_out_gate++;
	}
	
	if((sample_curr >= gate_beg) && (sample_curr <= gate_end)) //in of gate
	{
	    result_event.mean_in_gate += val_sample;
	    
	    if(val_sample < result_event.MIN_in_gate)
	    {
		result_event.MIN_in_gate = (Short_t)val_sample;
		result_event.time_min_in_gate = sample_curr;
	    }

	    if((sample_curr >= gate_maximum_beg) && (sample_curr <= gate_maximum_end))
	    {
		if(val_sample > result_event.MAX_in_gate)
		{
		    result_event.MAX_in_gate = (Short_t)val_sample;
		    result_event.time_max_in_gate = sample_curr;
		}
	    }
	    samples_in_gate++;
	}
    }
    result_event.mean_in_gate /= (float)samples_in_gate;
    result_event.mean_out_gate /= (float)samples_out_gate;
    //result_event.integral_in_gate /= (float)samples_in_gate;
    
    
    //PASS 2
    Double_t integral_in_gate_ = 0.;
	Double_t integral_in_gate_noises_2 = 0.;
	Double_t integral_in_gate_noises = 0.;
	Double_t integral_in_gate_noises_1 = 0.;
	Double_t gate_left = 120.;
	Double_t gate_right = 200.;
    result_event.RMS_in_gate = 0.;
    result_event.RMS_out_gate = 0.;
    for(UShort_t sample_curr = 0; sample_curr < sample_total; sample_curr++)
    {
	Float_t val_sample = scale_samples[sample_curr];
	//if(sample_curr < gate_beg) //out of gate
	if((sample_curr < gate_beg) || (sample_curr > gate_end)) //out of gate
	{
	    result_event.RMS_out_gate +=
	            (val_sample-result_event.mean_out_gate) * (val_sample-result_event.mean_out_gate);
	}
	
	if((sample_curr >= gate_beg) && (sample_curr <= gate_end)) //in of gate
	{
	    integral_in_gate_ += (val_sample-result_event.zero_level);
	    //printf("%f \n", integral_in_gate_);
	    result_event.RMS_in_gate +=
	            (val_sample-result_event.mean_in_gate) * (val_sample-result_event.mean_in_gate);
	}
	if ((sample_curr > gate_left) && (sample_curr < (gate_right+gate_left)/2)) integral_in_gate_noises_1 += (val_sample);
	if ((sample_curr > (gate_right+gate_left)/2) && (sample_curr < gate_right)) integral_in_gate_noises_2 += (val_sample);
//	integral_in_gate_noises = 2*(gate_end-gate_beg)/(gate_right-gate_left)*(integral_in_gate_noises_2-integral_in_gate_noises_1);
	integral_in_gate_noises = (integral_in_gate_noises_2-integral_in_gate_noises_1);
	
	if (integral_in_gate_noises==0) integral_in_gate_noises=10000;

    }
    
    
    result_event.integral_in_gate = (Int_t)integral_in_gate_;
    result_event.integral_in_gate_noises = (Int_t)integral_in_gate_noises;

    result_event.RMS_in_gate = TMath::Sqrt( result_event.RMS_in_gate / (float)samples_in_gate );
    result_event.RMS_out_gate = TMath::Sqrt( result_event.RMS_out_gate / (float)samples_out_gate );
    //printf("res %i \n", result_event.integral_in_gate);
    
    // Calculating timing
    Float_t Amplitude = result_event.MAX_in_gate - result_event.zero_level;
    Float_t trsh_01 = result_event.zero_level + Amplitude*0.1;
    Float_t trsh_03 = result_event.zero_level + Amplitude*0.3;
    Float_t trsh_05 = result_event.zero_level + Amplitude*0.5;
    Float_t trsh_09 = result_event.zero_level + Amplitude*0.9;

    Int_t point = result_event.time_max_in_gate;
    Float_t front_time_beg = GoToLevel(scale_samples, trsh_01, &point, -1, sample_total);

    point = result_event.time_max_in_gate;
    Float_t front_time_beg_03 = GoToLevel(scale_samples, trsh_03, &point, -1, sample_total);

    point = result_event.time_max_in_gate;
    Float_t time_front = GoToLevel(scale_samples, trsh_05, &point, -1, sample_total);

    point = result_event.time_max_in_gate;
    Float_t front_time_end = GoToLevel(scale_samples, trsh_09, &point, -1, sample_total);

    Float_t time = (front_time_end + front_time_beg)*0.5;
    //printf("begin: %f; end: %f; time_front: %f; time: %f\n", front_time_beg, front_time_end, time_front, time);
    result_event.time_cross = time_front;
    result_event.time_half = time;
    
    result_event.timestamp = TimeStamp_date*1e9 + (TimeStamp_ns >> 2);
	
	
    //Fitting
    if(IsFIT)
    {
	if(readdebug) printf("\n\nFit signal procedure --- \n");
        if(result_event.MAX_in_gate - result_event.zero_level>0) //cosmic ZS selection
        {
	    PronyFitter Pfitter(2, 2, gate_beg, gate_end, sample_total);
	    //Pfitter.SetDebugMode(1);
	    Pfitter.SetWaveform(scale_samples, result_event.zero_level);
	    //Pfitter.MakePileUpRejection(result_event.time_max_in_gate+1);  
	    Pfitter.SetExternalHarmonics(first_fit_harmonic, second_fit_harmonic);
	    if(!isfinite(front_time_end)) front_time_end = result_event.time_max_in_gate - 1;
	    int SignalBeg = round(front_time_end -
						(log(second_fit_harmonic) - log(first_fit_harmonic)) / (second_fit_harmonic - first_fit_harmonic));
	    //Pfitter.SetSignalBegin(SignalBeg);
	    Int_t best_signal_begin = Pfitter.ChooseBestSignalBegin(SignalBeg-1, SignalBeg+1);
	    Pfitter.SetSignalBegin(best_signal_begin);
	    Pfitter.CalculateFitAmplitudes();
	    result_event.FIT_integral_in_gate = Pfitter.GetIntegral(gate_beg, gate_end);
	    result_event.FIT_MAX_in_gate = Pfitter.GetMaxAmplitude();
	    //result_event.FIT_MAX_in_gate = Pfitter.GetFitValue((Int_t)result_event.time_max_in_gate);
	    result_event.FIT_Chi2 = Pfitter.GetChiSquare(gate_beg, gate_end, result_event.time_max_in_gate);
	    result_event.FIT_R2_gate = Pfitter.GetRSquare(gate_beg, gate_end);
	    result_event.FIT_R2_signal = Pfitter.GetRSquareSignal();
	    result_event.FIT_zero_level = Pfitter.GetZeroLevel();
	    result_event.FIT_first_harmonic = (Float_t)first_fit_harmonic;
	    result_event.FIT_second_harmonic = (Float_t)second_fit_harmonic;
	    Float_t FIT_amplitude = result_event.FIT_MAX_in_gate - result_event.FIT_zero_level;
	    result_event.FIT_time_half = Pfitter.GetX(result_event.FIT_zero_level + FIT_amplitude/2., SignalBeg-1, round(front_time_end), 0.25);
	    result_event.FIT_time_max_in_gate = (Pfitter.GetSignalBeginFromPhase() > 0)? Pfitter.GetSignalBeginFromPhase() + 
						(log(second_fit_harmonic) - log(first_fit_harmonic)) / (second_fit_harmonic - first_fit_harmonic) : 0;

	    //FIT QA
	    if(!FIT_QA_mode.Contains("false"))
	    {
		if(check_fit_arr->GetLast()+1 < (Int_t)fitQA_arg_arr[0])
		{
		    Bool_t FIT_QA = true;
		    if(FIT_QA_mode.Contains("default")) FIT_QA *= true;
		    if(FIT_QA_mode.Contains("neg_fit_integral")) FIT_QA *= result_event.FIT_integral_in_gate < 0;
		    if(FIT_QA_mode.Contains("diff_fit_int_and_wf_int")) FIT_QA *= 
			abs(result_event.FIT_integral_in_gate - result_event.integral_in_gate) > (fitQA_arg_arr[2])/100.*abs(result_event.FIT_integral_in_gate);
		    if(FIT_QA_mode.Contains("diff_fit_ampl_and_wf_ampl")) FIT_QA *= 
			abs(result_event.FIT_MAX_in_gate-result_event.zero_level - result_event.MAX_in_gate+result_event.zero_level) > 
										(fitQA_arg_arr[3])/100.*abs(result_event.MAX_in_gate-result_event.zero_level);
		    if(FIT_QA_mode.Contains("saturation")) FIT_QA *= result_event.integral_in_gate > fitQA_arg_arr[4];
		    if(FIT_QA_mode.Contains("integral_region")) FIT_QA *= (result_event.integral_in_gate > fitQA_arg_arr[5])&&
										(result_event.integral_in_gate < fitQA_arg_arr[6]);
		    if(FIT_QA_mode.Contains("amlitude_region")) FIT_QA *= (result_event.MAX_in_gate-result_event.zero_level > fitQA_arg_arr[7])&&
										(result_event.MAX_in_gate-result_event.zero_level < fitQA_arg_arr[8]);
		    if(FIT_QA_mode.Contains("chi_square")) FIT_QA *= (result_event.FIT_Chi2 > fitQA_arg_arr[9])&&
													(result_event.FIT_Chi2 < fitQA_arg_arr[10]);
		    if(FIT_QA_mode.Contains("R_square")) FIT_QA *= (result_event.FIT_R2_gate > fitQA_arg_arr[11])&&
													(result_event.FIT_R2_gate < fitQA_arg_arr[12]);
		    if(FIT_QA_mode.Contains("R_square_signal")) FIT_QA *= (result_event.FIT_R2_signal > fitQA_arg_arr[13])&&
													(result_event.FIT_R2_signal < fitQA_arg_arr[14]);
		    if(FIT_QA_mode.Contains("time_half")) FIT_QA *= (result_event.FIT_time_half > fitQA_arg_arr[15])&&
													(result_event.FIT_time_half < fitQA_arg_arr[16]);
		    if(FIT_QA_mode.Contains("time_max_in_gate")) FIT_QA *= (result_event.FIT_time_max_in_gate > fitQA_arg_arr[17])&&
													(result_event.FIT_time_max_in_gate < fitQA_arg_arr[18]);
		    Bool_t Selected_channel = true;
		    if(fitQA_arg_arr[1] != -1) if(ch_iter != (Int_t)fitQA_arg_arr[1]) Selected_channel = false;

		    if(FIT_QA&&Selected_channel)
		    {
			TString signal_name = Form("FIT QA mode '%s' ch_num %i fit_integral %i integral %i chi2 %.1f R2 %.2f", 
				FIT_QA_mode.Data(), ch_iter, (Int_t)result_event.FIT_integral_in_gate, result_event.integral_in_gate, 
												result_event.FIT_Chi2, result_event.FIT_R2_gate);
			Pfitter.DrawFit(check_fit_arr, signal_name);
		    }
		}
	    }
	}

        if(readdebug) printf("------------------------ \n");
    }
}




void BinDataReader::Calculate_fit_harmonics(event_fit_data_struct &result_fit_event, Int_t &event_counter, Int_t event_num, Int_t ch_num, Int_t gate_beg, Int_t gate_end, Int_t gate_maximum_beg, Int_t gate_maximum_end, TString source_path, TString run_name, Float_t ampl_scale, TObjArray *check_fit_arr)
{


  if(readdebug) printf("Calculating fit harmonics -------------------------------------------\n");
  if(readdebug) printf("h_num: %i; gate_beg: %i; gate_end %i; ampl_scale %i\n", ch_num, gate_beg, gate_end, ampl_scale);

    Int_t ch_iter = ch_num;

    result_fit_event.reset();
    if(samples_data == NULL) return;

    
    Float_t *scale_samples = new Float_t[sample_total];
    for(Int_t isaml = 0; isaml < sample_total; isaml++)
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

    for(UShort_t sample_curr = 0; sample_curr < sample_total; sample_curr++)
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

    Double_t integral_in_gate_ = 0.;

    for(UShort_t sample_curr = 0; sample_curr < sample_total; sample_curr++)
    {
	Float_t val_sample = scale_samples[sample_curr];
	
	if((sample_curr >= gate_beg) && (sample_curr <= gate_end)) //in of gate
	{
	    integral_in_gate_ += (val_sample-zero_level_);
	    //printf("%f \n", integral_in_gate_);
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
	Float_t front_time_beg_03 = GoToLevel(scale_samples, trsh_03, &point, -1, sample_total);

	point = time_max_in_gate_;
	Float_t front_time_end = GoToLevel(scale_samples, trsh_09, &point, -1, sample_total);

	//Fitting

	if(readdebug) printf("\n\nFit signal procedure --- \n");

	PronyFitter Pfitter(2, 2, gate_beg, gate_end, sample_total);
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
		Float_t fit_integral = Pfitter.GetIntegral(gate_beg, gate_end);
		Float_t fit_chi2 = Pfitter.GetChiSquare(gate_beg, gate_end, time_max_in_gate_);
		Float_t fit_R2 = Pfitter.GetRSquare(gate_beg, gate_end);
		if(IsDebug) printf("fit integral %.0f integral %.0f chi2 %.1f R2 %.3f\n", fit_integral, integral_in_gate_, fit_chi2, fit_R2);

		if(fit_R2<0.02)
		{
			event_counter++;
			result_fit_event.first_harmonic = (Float_t)harmonics[1];
			result_fit_event.second_harmonic = (Float_t)harmonics[2];
			result_fit_event.fit_chi2 = fit_chi2;
			result_fit_event.fit_R2 = fit_R2;
			result_fit_event.integral_in_gate = (Float_t)integral_in_gate_;
			result_fit_event.fit_integral_in_gate = fit_integral;
			TString signal_name = Form("ch_num %i fit_integral %.0f integral %.0f chi2 %.1f fit_R2 %.3f", 
								ch_iter, fit_integral, integral_in_gate_, fit_chi2, fit_R2);
			if(event_counter<5) Pfitter.DrawFit(check_fit_arr, signal_name);
		}


		if(readdebug) printf("------------------------ \n");
	}
    }

    delete[] scale_samples;

    if(readdebug) printf("\n-----------------------------------------------------------------\n");
}


























#if IsUseNewFitHarmonics
    TString result_fit_file_name;
    result_fit_file_name =  save_at_same_folder ?
                source_path + run_name + "/" + run_name + result_file_prefix + "_fit_harmonics.root":
                result_path + run_name + result_file_prefix + "_fit_harmonics.root";
    TFile *result_fit_file = new TFile(result_fit_file_name.Data(), "RECREATE"); 
    TTree *data_fit_tree = new TTree("adc64_fit_data", "adc64_fit_data");
    event_fit_data_struct *event_fit_struct = new event_fit_data_struct[total_channels];
    for(Int_t ch = 0; ch < total_channels; ch++)
	(event_fit_struct+ch)->CreateBranch(data_fit_tree, ch);
    TObjArray *check_fit_harmonics_arr = new TObjArray;
    Int_t *event_counters = new Int_t[total_channels];
    for(Int_t i = 0; i < total_channels; i++)
	event_counters[i] = 0;
#else
    //create result root file, create tree
    TFile *result_file = new TFile(result_file_name.Data(), "RECREATE");

    if(IsFIT)
    {
	TFile *fit_harmonics_file = new TFile(fit_harmonics_path.Data(), "READONLY");
	TTree *data_fit_tree = dynamic_cast<TTree*>(fit_harmonics_file->FindObjectAny("adc64_fit_data"));
	event_fit_data_struct *event_fit_struct = new event_fit_data_struct[total_channels];
	for(Int_t ch = 0; ch < total_channels; ch++)
	{
	    data_fit_tree->SetBranchAddress((event_fit_struct+ch)->GetChName(ch).Data(), (event_fit_struct+ch));
	}
	Int_t *harmonic_counter = new Int_t[total_channels];
	for(Int_t i = 0; i < total_channels; i++)
	    harmonic_counter[i] = 0;

	for(Int_t entry = 0; entry <= data_fit_tree->GetEntries(); entry++)
	{
	    data_fit_tree->GetEntry(entry);
	    for(Int_t ch_iter = 0; ch_iter < total_channels; ch_iter++)
	    {
		if( (event_fit_struct[ch_iter].first_harmonic > 0.1)&&(event_fit_struct[ch_iter].first_harmonic < 1.) && 
			(event_fit_struct[ch_iter].second_harmonic > 0.1)&&(event_fit_struct[ch_iter].second_harmonic < 1.) )
		{	
		    first_fit_harmonic[ch_iter] += event_fit_struct[ch_iter].first_harmonic;
		    second_fit_harmonic[ch_iter] += event_fit_struct[ch_iter].second_harmonic;
		    harmonic_counter[ch_iter] ++;
		}
	    }
	}

	for(Int_t ch_iter = 0; ch_iter < total_channels; ch_iter++)
	{
	    if(harmonic_counter[ch_iter] != 0) first_fit_harmonic[ch_iter] /= harmonic_counter[ch_iter];
	    else first_fit_harmonic[ch_iter] = 0;

	    if(harmonic_counter[ch_iter] != 0) second_fit_harmonic[ch_iter] /= harmonic_counter[ch_iter];	    
	    else second_fit_harmonic[ch_iter] = 0;
	}

	fit_harmonics_file->Close();

#if IsUseCommonHarmonics

	Double_t first_fit_harmonic_common = 0.;
	Double_t second_fit_harmonic_common = 0.;
	Int_t non_zero_counter = 0;
	for(Int_t ch_iter = 0; ch_iter < total_channels; ch_iter++)
	{
	    if((first_fit_harmonic[ch_iter] > 0.)&&(second_fit_harmonic[ch_iter] > 0.))
	    {
		first_fit_harmonic_common += first_fit_harmonic[ch_iter];
		second_fit_harmonic_common += second_fit_harmonic[ch_iter];
		non_zero_counter++;
	    }
	}

	for(Int_t ch_iter = 0; ch_iter < total_channels; ch_iter++)
	{
	    first_fit_harmonic[ch_iter] = first_fit_harmonic_common/non_zero_counter;
	    second_fit_harmonic[ch_iter] = second_fit_harmonic_common/non_zero_counter;
	}
		
#endif

    }
#endif









#if IsUseNewFitHarmonicss
	    if(readdebug)   printf("\n\nCall fit harmonics calculator: Event %i channel %i board %i ch_board %i\n", nevent, channel, ch_board, ch_num);

	    bin_reader_arr[ch_board].Calculate_fit_harmonics(event_fit_struct[channel], event_counters[channel], nevent, ch_num, gate_beg, gate_end, gate_maximum_beg, gate_maximum_end, source_path, run_name, ampl_scale, check_fit_harmonics_arr);

#else
	    if(readdebug)   printf("\n\nCall waveform calculator: Event %i channel %i board %i ch_board %i\n", nevent, channel, ch_board, ch_num);

	    bin_reader_arr[ch_board].Calculate_waveform(event_struct[channel], ch_num, gate_beg, gate_end, gate_maximum_beg, gate_maximum_end, ampl_scale, IsFIT, first_fit_harmonic[channel], second_fit_harmonic[channel], FIT_QA_mode, check_fit_arr, fitQA_arg_arr);

#endif

















#if IsUseNewFitHarmonics
    response = Form("%i events written in %s file\nlast readed event from file(0): %i; events in dig: %i",
                    events_to_process, result_fit_file_name.Data(), last_readed, events_in_dig);

    result_fit_file->cd();
    data_fit_tree->Write();
    TIter nx_iter((TCollection*)(check_fit_harmonics_arr));
    TObject* obj_ptr;
    while ( obj_ptr=(TObjArray*)nx_iter() )
    {
	//if(obj_ptr == 0) continue; 
	obj_ptr->Write();
    }
    result_fit_file->Close();
    delete check_fit_harmonics_arr;
#else
    //copy info file
    if(!save_at_same_folder)gROOT->ProcessLine( Form(".cp %s %s.txt",readme_reader_ptr->Get_info_file_path().Data(), (result_path+run_name+result_file_prefix).Data()) );


    response = Form("%i events writed in %s file\nlast readed event from file(0): %i; events in dig: %i",
                    events_to_process, result_file_name.Data(), last_readed, events_in_dig);

#if !IsUseChunkStructure
    result_file->cd();
    data_tree->Write();
    result_file->Close();
#endif

#endif