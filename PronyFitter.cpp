#include "PolinomRoots.h"
#include "PronyFitter.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "../my_macros/like_ivashkin_wants_it.h"

PronyFitter::PronyFitter(Int_t model_order, Int_t exponent_number, Int_t gate_beg, Int_t gate_end, Int_t sample_total)
{
	 Initialize(model_order, exponent_number, gate_beg, gate_end, sample_total);
}

void PronyFitter::SetDebugMode(Bool_t IsDebug)
{
	fIsDebug = IsDebug;
}

void PronyFitter::Initialize(Int_t model_order, Int_t exponent_number, Int_t gate_beg, Int_t gate_end, Int_t sample_total)
{
	fModelOrder = model_order;
	fExpNumber = exponent_number;
	fGateBeg = gate_beg;
	fGateEnd = gate_end;
	fPileUpBeg = gate_end;
	fSampleTotal = sample_total;
	AllocData();
}

void PronyFitter::AllocData()
{
	fValSampleArr = new Float_t[fSampleTotal];
	fFITValSampleArr = new Float_t[fSampleTotal];
	for(Int_t i = 0; i < fSampleTotal; i++)
	{
		fValSampleArr[i] = 0.;
		fFITValSampleArr[i] = 0.;
	}

	fz = new Double_t[fExpNumber+1];
	fh = new Double_t[fExpNumber+1];
	for(Int_t i = 0; i < fExpNumber+1; i++)
	{
		fz[i] = 0.;
		fh[i] = 0.;
	}
}

void PronyFitter::SetWaveform(Float_t *val_sample_arr, Float_t zero_level)
{
	memcpy(fValSampleArr, val_sample_arr, (fSampleTotal)*sizeof(float));
	fZeroLevel = zero_level;
}

Int_t PronyFitter::CalcSignalBegin(Float_t front_time_beg_03, Float_t front_time_end)
{
	return ceil((3*front_time_beg_03 - front_time_end)/2.) ;
}

void PronyFitter::SetSignalBegin(Int_t SignalBeg)
{
	fSignalBegin = SignalBeg;
	if(fIsDebug) printf("\nsignal begin %i  zero level %.0f\n", fSignalBegin, fZeroLevel);
}

void PronyFitter::CalculateFitHarmonics()
{
	Double_t **Rkm_arr = new Double_t*[fModelOrder];
	for(Int_t i =0; i < fModelOrder; i++)
	{
		Rkm_arr[i] = new Double_t[fModelOrder];
		for(Int_t j = 0; j < fModelOrder; j++) Rkm_arr[i][j] = 0.;
	}

	Double_t *R0k_arr = new Double_t[fModelOrder];
	for(Int_t i =0; i < fModelOrder; i++)
		R0k_arr[i] = 0.;

	Double_t *a_arr = new Double_t[fModelOrder+1];
	for(Int_t i =0; i <= fModelOrder; i++)
		a_arr[i] = 0.;

	for(Int_t i = 1; i <= fModelOrder; i++)
	    for(Int_t j = 1; j <= fModelOrder; j++)
		for(UShort_t sample_curr = fSignalBegin + fModelOrder; sample_curr < fGateEnd; sample_curr++)
			Rkm_arr[i-1][j-1] += (Double_t) (fValSampleArr[sample_curr-i]-fZeroLevel)*(fValSampleArr[sample_curr-j]-fZeroLevel);

	for(Int_t i = 1; i <= fModelOrder; i++)
	    for(UShort_t sample_curr = fSignalBegin + fModelOrder; sample_curr < fGateEnd; sample_curr++)
		R0k_arr[i-1] -= (Double_t) (fValSampleArr[sample_curr]-fZeroLevel)*(fValSampleArr[sample_curr-i]-fZeroLevel);

	if(fIsDebug)
	{
		printf("system forward\n");
		for(Int_t i =0; i < fModelOrder; i++)
		{
			for(Int_t j = 0; j < fModelOrder; j++) printf("%e ", Rkm_arr[i][j]);
			printf("%e\n", R0k_arr[i]);
		}
	}

	Double_t *a = new Double_t[fModelOrder];
	for(Int_t i = 0; i < fModelOrder; i++)
	    a[i] = 0.;

	SolveSLEGauss(a, Rkm_arr, R0k_arr, fModelOrder);
	if(fIsDebug)
	{
		printf("SLE roots ");
		for(Int_t i =0; i < fModelOrder; i++) printf(" %e ", a[i]);
		printf("\n");
	}	

	for(Int_t i = 0; i < fModelOrder; i++) a_arr[i+1] = a[i];
	a_arr[0] = 1.;

	Double_t *z = new Double_t[fModelOrder];
	for(Int_t i = 0; i < fModelOrder; i++)
	    z[i] = 0.;

	Int_t total_roots;
	polynomRealRoots(z, fModelOrder, a_arr, total_roots);
	if(fIsDebug)
	{
		printf("forward polinom roots ");
		for(Int_t i = 0; i < fModelOrder; i++)
		printf("%.5f ", z[i]);
		printf("\n");

		Float_t *omega = new Float_t[fModelOrder];
		for(Int_t i = 0; i < fModelOrder; i++) omega[i] = log(z[i]);
		printf("forward freqs ");
		for(Int_t i = 0; i < fModelOrder; i++)
		printf("%.5f ", omega[i]);
		printf("\n");
		delete[] omega;
	}

	Double_t *z_arr = new Double_t[fExpNumber+1];
	for(Int_t i =0; i <= fExpNumber; i++)
		z_arr[i] = 0.;

	for(Int_t i = 0; i < fExpNumber; i++) z_arr[i+1] = z[i];
	z_arr[0] = 1.;
	SetHarmonics(z_arr);
	fTotalPolRoots = total_roots;

	for(Int_t i =0; i < fModelOrder; i++) delete[] Rkm_arr[i];
	delete[] Rkm_arr;
	delete[] R0k_arr;
	delete[] a;
	delete[] a_arr;
	delete[] z;
	delete[] z_arr;
}

void PronyFitter::SetHarmonics(Double_t *z)
{
	memcpy(fz, z, (fExpNumber+1)*sizeof(double));	
}

void PronyFitter::SetExternalHarmonics(Double_t z1, Double_t z2)
{
	Double_t *z_arr = new Double_t[fExpNumber+1];
	for(Int_t i =0; i <= fExpNumber; i++)
		z_arr[i] = 0.;
	z_arr[0] = 1.;
	z_arr[1] = z1;
	z_arr[2] = z2;
	SetHarmonics(z_arr);
	delete[] z_arr;	
}

Double_t *PronyFitter::GetHarmonics()
{
	return fz;	
}

Int_t PronyFitter::GetNumberPolRoots()
{
	return fTotalPolRoots;	
}

void PronyFitter::MakePileUpRejection(Int_t time_max)
{
	Int_t smoothing_gate = 5;
	Float_t *ValSampleArrSmoothed = new Float_t[fSampleTotal];
	for(Int_t sample_iter = 0; sample_iter < fSampleTotal; sample_iter++)
	    ValSampleArrSmoothed[sample_iter] = 0.;

	for(Int_t sample_iter = time_max; sample_iter < fGateEnd; sample_iter++)
	{
	    for(Int_t smooth_iter = 0; smooth_iter < smoothing_gate; smooth_iter++) ValSampleArrSmoothed[sample_iter] += fValSampleArr[sample_iter+smooth_iter];
	    ValSampleArrSmoothed[sample_iter] /= smoothing_gate;
	}

	for(Int_t sample_iter = time_max+1; sample_iter < fGateEnd; sample_iter++)
	    if(ValSampleArrSmoothed[sample_iter] >= ValSampleArrSmoothed[sample_iter-1]) {fPileUpBeg = sample_iter; break;}

	delete[] ValSampleArrSmoothed;	
}

void PronyFitter::CalculateFitAmplitudes()
{

	Double_t **Zik_arr = new Double_t*[fExpNumber+1];
	for(Int_t i = 0; i < fExpNumber+1; i++)
	{
		Zik_arr[i] = new Double_t[fExpNumber+1];
		for(Int_t j = 0; j < fExpNumber+1; j++) Zik_arr[i][j] = 0.;
	}

	Double_t *Zyk_arr = new Double_t[fExpNumber+1];
	for(Int_t i = 0; i < fExpNumber+1; i++)
		Zyk_arr[i] = 0.;

	Double_t *z_power = new Double_t[fExpNumber+1];
	for(Int_t i = 0; i < fExpNumber+1; i++)	z_power[i] = (Double_t) pow(fz[i], fSignalBegin);

	for(Int_t i = 0; i <= fExpNumber; i++)
	    for(Int_t j = 0; j <= fExpNumber; j++)
		if(fz[i]*fz[j] != 1) Zik_arr[i][j] = (Double_t) (pow(fz[i]*fz[j], fPileUpBeg) - pow(fz[i]*fz[j], fSignalBegin))/(fz[i]*fz[j] - 1.);
		else Zik_arr[i][j] = (Double_t) fPileUpBeg-fSignalBegin;

	for(Int_t i = 0; i <= fExpNumber; i++)
	{
	    for(UShort_t sample_curr = fSignalBegin; sample_curr < fPileUpBeg; sample_curr++)
	    {
		Zyk_arr[i] += (Double_t) ( fValSampleArr[sample_curr] * z_power[i] );
		z_power[i] *= fz[i];
	    }
	}

	if(fIsDebug)
	{
		printf("\nampl calculation\n");
		for(Int_t i = 0; i <= fExpNumber; i++)
		{
			for(Int_t j = 0; j <= fExpNumber; j++) printf("%e ", Zik_arr[i][j]);
			printf("%e\n", Zyk_arr[i]);
		}
	}

	Double_t *h = new Double_t[fExpNumber+1];
	for(Int_t i = 0; i < fExpNumber+1; i++)
	    h[i] = 0.;

	SolveSLEGauss(h, Zik_arr, Zyk_arr, fExpNumber+1);
	//SolveSLEOverDetermined_3x2_(h, Zik_arr, Zyk_arr);
	//SolveSLECholesky(h, Zik_arr, Zyk_arr, fExpNumber+1);
	memcpy(fh, h, (fExpNumber+1)*sizeof(double));

	if(fIsDebug)
	{
		printf("amplitudes\n%.0f ", fh[0]);
		for(Int_t i = 1; i < fExpNumber+1; i++)
		printf("%e ", fh[i]);
		printf("\n\n");
	}

	Float_t *FITValSampleArr = new Float_t[fSampleTotal];
	for(Int_t i = 0; i < fSampleTotal; i++)
	    FITValSampleArr[i] = 0.;

	for(Int_t i = 0; i < fExpNumber+1; i++)	z_power[i] = (Double_t) pow(fz[i], fSignalBegin);

	for(UShort_t sample_curr = 0; sample_curr < fSampleTotal; sample_curr++)
	{
	    if((sample_curr>=fSignalBegin)&&(sample_curr<=fGateEnd))
	    {
		for(Int_t i = 0; i < fExpNumber+1; i++)
		{ 
		    FITValSampleArr[sample_curr] += (Float_t) fh[i] * z_power[i]; 
		    z_power[i] *= fz[i];
		}
	    }
	    else FITValSampleArr[sample_curr] = (Float_t) fh[0] * 1.;
	}
	memcpy(fFITValSampleArr, FITValSampleArr, (fSampleTotal)*sizeof(float));

	for(Int_t i =0; i < fExpNumber+1; i++) delete[] Zik_arr[i];
	delete[] Zik_arr;
	delete[] Zyk_arr;
	delete[] h;
	delete[] FITValSampleArr;
	delete[] z_power;
}

Double_t *PronyFitter::GetAmplitudes()
{
	return fh;	
}

Float_t PronyFitter::GetIntegral(Int_t gate_beg, Int_t gate_end)
{
	Float_t integral = 0.;
	for(UShort_t sample_curr = gate_beg; sample_curr < gate_end; sample_curr++)
		integral += fFITValSampleArr[sample_curr] - fh[0];

	if(isfinite(integral)) return integral;
	return 0;
}

Float_t PronyFitter::GetFitValue(Int_t sample_number)
{
	if(isfinite(fFITValSampleArr[sample_number])) return fFITValSampleArr[sample_number];
	return 0;
}

Float_t PronyFitter::GetZeroLevel()
{
	return (Float_t)fh[0];
}

Float_t PronyFitter::GetSignalBeginFromPhase()
{
	if(-fh[2]/fh[1] > 0) return log(-fh[2]/fh[1])/log(fz[1]/fz[2]);
	return 0;
}

Float_t PronyFitter::GetMaxAmplitude()
{
	Float_t amplitude = 0.;
	if(GetSignalBeginFromPhase() > 0)
	{
		Float_t signal_max_time = GetSignalBeginFromPhase() + (log(fz[2])-log(fz[1]))/(fz[2]-fz[1]);
		for(Int_t i = 0; i < fExpNumber+1; i++)
			amplitude += (Float_t) fh[i] * pow(fz[i], signal_max_time);

		if(isfinite(amplitude)) return amplitude;
	}
	return 0;
}

Float_t PronyFitter::GetX(Float_t level, Int_t first_sample, Int_t last_sample)
{
	Int_t step = 0;
	if(first_sample<last_sample) step = 1;
	else step = -1;
	Float_t result_sample = 0.;
	Int_t sample_to_check = first_sample;
	Float_t amplitude = 0.;
	Float_t amplitude_prev = fFITValSampleArr[sample_to_check-step];
	while( (first_sample-sample_to_check)*(last_sample-sample_to_check) <= 0 )
	{	
		amplitude = fFITValSampleArr[sample_to_check];
		if( (level - amplitude)*(level - amplitude_prev) <= 0 )
		{
			result_sample = LevelBy2Points(sample_to_check, amplitude, sample_to_check-step, amplitude_prev, level);
			return result_sample;
		}
		amplitude_prev = amplitude;
		sample_to_check += step;
	}

	return 0;
}

Float_t PronyFitter::GetX(Float_t level, Int_t first_sample, Int_t last_sample, Float_t step)
{
	Float_t result_sample = 0.;
	Float_t sample_to_check = first_sample;
	Float_t amplitude = 0.;
	Float_t amplitude_prev = 0.;
	for(Int_t i = 0; i < fExpNumber+1; i++)
		amplitude_prev += (Float_t) (fh[i] * pow(fz[i], sample_to_check-step));
	while( (first_sample-sample_to_check)*(last_sample-sample_to_check) <= 0 )
	{	
		amplitude = 0.;
		for(Int_t i = 0; i < fExpNumber+1; i++)
			amplitude += (Float_t) (fh[i] * pow(fz[i], sample_to_check));		
		if( (level - amplitude)*(level - amplitude_prev) <= 0 )
		{
			if(amplitude != amplitude_prev)
				result_sample = LevelBy2Points(sample_to_check, amplitude, sample_to_check-step, amplitude_prev, level);
			return result_sample;
		}
		amplitude_prev = amplitude;
		sample_to_check += step;
	}

	return 0;	
}

Float_t PronyFitter::LevelBy2Points(Float_t X1, Float_t Y1, Float_t X2, Float_t Y2, Float_t Y0)
{
	return (X1*Y0 - X1*Y2 - X2*Y0 + X2*Y1) / (Y1-Y2);
}

Float_t PronyFitter::GetRSquare(Int_t gate_beg, Int_t gate_end)
{
	Float_t R2 = 0.;
	Float_t RSS = 0.;
	Float_t TSS = 0.;	
	Int_t m = gate_end-gate_beg;
	Int_t params_number = 1+2*fModelOrder;
	if(m<=params_number) return 999;
	Float_t average = 0.;
	for(UShort_t sample_curr = gate_beg; sample_curr < gate_end; sample_curr++) average += fValSampleArr[sample_curr];
	average /= m;

	for(UShort_t sample_curr = gate_beg; sample_curr < gate_end; sample_curr++)
	{
	    RSS += (fFITValSampleArr[sample_curr]-fValSampleArr[sample_curr])*(fFITValSampleArr[sample_curr]-fValSampleArr[sample_curr]);
	    TSS += (fValSampleArr[sample_curr]-average)*(fValSampleArr[sample_curr]-average);
	}
	if(TSS == 0) return 999;
	R2 = RSS/TSS; // correct definition is R2=1.-RSS/TSS, but R2=RSS/TSS is more convenient

	Float_t R2_adj = R2*(m-1)/(m-params_number);
	return R2_adj;
}

Float_t PronyFitter::GetRSquareSignal()
{
	return GetRSquare(fSignalBegin, fPileUpBeg);
}

Float_t PronyFitter::GetChiSquare(Int_t gate_beg, Int_t gate_end, Int_t time_max)
{
	Float_t chi2 = 0.;
	Int_t freedom_counter = 0;
	Int_t regions_number = 10;
	Float_t amplitude_max = abs(fValSampleArr[time_max] - fZeroLevel);
	if(amplitude_max == 0) return 999;

	Int_t *probability_exp = new Int_t[regions_number];
	Int_t *probability_theor = new Int_t[regions_number];
	Float_t *amplitude_regions = new Float_t[regions_number+1];
	amplitude_regions[0] = 0.;
	for(Int_t i = 0; i < regions_number; i++)
	{
		probability_exp[i] = 0;
		probability_theor[i] = 0;
		amplitude_regions[i+1] = (i+1)*amplitude_max/regions_number;
	}

	for(UShort_t sample_curr = gate_beg; sample_curr < gate_end; sample_curr++)
	{
	    for(Int_t i = 0; i < regions_number; i++)
	    {
		if((abs(fValSampleArr[sample_curr]-fZeroLevel) > amplitude_regions[i])&&(abs(fValSampleArr[sample_curr]-fZeroLevel) <= amplitude_regions[i+1]))
		    probability_exp[i]++;
		if((abs(fFITValSampleArr[sample_curr]-fh[0]) > amplitude_regions[i])&&(abs(fFITValSampleArr[sample_curr]-fh[0]) <= amplitude_regions[i+1]))
		    probability_theor[i]++;
	    }
	}

	for(Int_t i = 0; i < regions_number; i++)
	{
	    if(probability_exp[i] > 0)
	    {
		chi2 += pow(probability_exp[i] - probability_theor[i], 2.) / (probability_exp[i]);
		freedom_counter++;
	    }
	}

	if(freedom_counter>0) chi2 /= freedom_counter;
	delete[] probability_exp;
	delete[] probability_theor;
	delete[] amplitude_regions;

	return chi2;
}

Float_t PronyFitter::GetDeltaInSample(Int_t sample)
{
	return fFITValSampleArr[sample] - fValSampleArr[sample];
}

void PronyFitter::DrawFit(TObjArray *check_fit_arr, TString hist_title)
{
	Float_t *sample_arr = new Float_t[fSampleTotal];
	for(Int_t i = 0; i < fSampleTotal; i++)
	    sample_arr[i] = (Float_t) i;

	TGraph* tgr_ptr = new TGraph( fSampleTotal, sample_arr, fValSampleArr);
	TGraph* tgr_ptr_fit = new TGraph( fSampleTotal, sample_arr, fFITValSampleArr);
	TCanvas *canv_ptr = new TCanvas(hist_title.Data());
	
	tgr_ptr->SetTitle(hist_title.Data());
	
	tgr_ptr->Draw();

	tgr_ptr_fit->SetLineColor(kRed);
	tgr_ptr_fit->SetLineWidth(2);

	tgr_ptr_fit->Draw("same");

	check_fit_arr->Add(canv_ptr);

	delete[] sample_arr;
}

void PronyFitter::DrawFit_for_presentation(TCanvas *canv_ptr, TString source_path, TString run_name, TString fit_quality,TString hist_title)
{
	Float_t *sample_arr = new Float_t[fSampleTotal];
	Float_t *yarr = new Float_t[fSampleTotal];
	Float_t *yarrfit = new Float_t[fSampleTotal];

	for(Int_t i = 0; i < fSampleTotal; i++)
	{

	    sample_arr[i] = (Float_t) i;
		yarr[i] = fValSampleArr[i]-fFITValSampleArr[0];
		yarrfit[i] = fFITValSampleArr[i]-fFITValSampleArr[0];
	}


	TGraph* tgr_ptr = new TGraph( fSampleTotal, sample_arr, yarr);
	TGraph* tgr_ptr_fit = new TGraph( fSampleTotal, sample_arr, yarrfit);
	//TCanvas *canv_ptr = new TCanvas(run_name);
	canv_ptr->cd();
	tgr_ptr->SetTitle("");
	Int_t color = 1;
	if(fit_quality =="good_fit") color = 1;
	if(fit_quality=="bad_fit") color = 4;
	graph_like_ivashkin_wants_it(tgr_ptr, color, "Time [a.u.]", "Amplitude [ADC ch]");
	tgr_ptr->Draw("");

	if(fit_quality =="good_fit") 	tgr_ptr_fit->SetLineColor(kRed);
	if(fit_quality=="bad_fit") 	tgr_ptr_fit->SetLineColor(kGreen);


	tgr_ptr_fit->SetLineWidth(2);

	tgr_ptr_fit->Draw("same");

	canv_ptr->SaveAs ((source_path + run_name+ "/"+run_name+"_"+fit_quality + ".pdf)").Data ());

	delete[] sample_arr;
}

Int_t PronyFitter::ChooseBestSignalBeginHarmonics(Int_t first_sample, Int_t last_sample) 
{
	Float_t best_R2 = 0.;
	Int_t best_signal_begin = 0;
	Bool_t IsPositiveRoot;
	Bool_t IsGoodFit = false;
	Int_t good_fit_counter = 0;

	for(Int_t signal_begin = first_sample; signal_begin <= last_sample; signal_begin++)
	{
	    SetSignalBegin(signal_begin);
	    CalculateFitHarmonics();
  	    IsPositiveRoot = true;
  	    for(Int_t j = 0; j < fExpNumber; j++) IsPositiveRoot *= (fz[j+1] > 1e-2)&&(fz[j+1] < 1e1);
	    IsGoodFit = (fTotalPolRoots>0)&&(IsPositiveRoot);

	    if(IsGoodFit)
	    {
		good_fit_counter++;
		CalculateFitAmplitudes();
		Float_t R2 = GetRSquare(fGateBeg, fPileUpBeg);
		if(good_fit_counter == 1) {best_R2 = R2; best_signal_begin = signal_begin;}
		if(R2<best_R2)
		{
		    best_R2 = R2;
		    best_signal_begin = signal_begin;
		}		
	    }
	}

	return best_signal_begin;
}

Int_t PronyFitter::ChooseBestSignalBegin(Int_t first_sample, Int_t last_sample) 
{
	Float_t best_R2 = 0.;
	Int_t best_signal_begin = first_sample;

	for(Int_t signal_begin = first_sample; signal_begin <= last_sample; signal_begin++)
	{
	    SetSignalBegin(signal_begin);
	    CalculateFitAmplitudes();
	    Float_t R2 = GetRSquare(fGateBeg, fPileUpBeg);
	    if(signal_begin == first_sample) best_R2 = R2;
	    if(R2<best_R2)
	    {
		best_R2 = R2;
		best_signal_begin = signal_begin;
	    }
	}

	return best_signal_begin;
}

void PronyFitter::SolveSLEGauss(Double_t *x, Double_t **r, Double_t *b, Int_t n) 
{
	Bool_t solvable = true;
	Int_t maxRow;
	Double_t maxEl, tmp, c;
	Double_t **a = new Double_t *[n];
	for(Int_t i = 0; i < n; i++)
	{
		a[i] = new Double_t [n+1];
		for(Int_t j = 0; j < n+1; j++) a[i][j] = 0.;
	}

	for(Int_t i = 0; i < n; i++)
	{
		for(Int_t j = 0; j < n; j++) a[i][j] = r[i][j];
		a[i][n] = b[i];
	}
 
	for(Int_t i = 0; i < n; i++)
	{
		maxEl = abs(a[i][i]);
		maxRow = i;
		for(Int_t k = i+1; k < n; k++)
			if(abs(a[k][i]) > maxEl)
			{
				maxEl = abs(a[k][i]);
				maxRow = k;
			}

		if(maxEl == 0) { solvable = false; if(fIsDebug) printf("SLE has no solution\n"); } 

		for(Int_t k = i; k < n+1; k++)
		{
			tmp = a[maxRow][k];
			a[maxRow][k] = a[i][k];
			a[i][k] = tmp;
		}

		for(Int_t k = i+1; k < n; k++)
		{
			c = -a[k][i]/a[i][i];
			for(Int_t j = i; j < n+1; j++)
			{
				if(i==j) a[k][j] = 0.;
				else a[k][j] += c * a[i][j];
			}
		}
	}

	for(Int_t i = n-1; i >= 0; i--)
	{
		x[i] = a[i][n]/a[i][i];
		for(Int_t k = i-1; k >= 0; k--)
			a[k][n] -= a[k][i] * x[i];
	}

	if(!solvable) { for(Int_t i = n-1; i >= 0; i--) x[i] = 0.;}

	for(Int_t i=0;i<n;i++) delete[] a[i];
	delete[] a;
}

void PronyFitter::SolveSLECholesky(Double_t *x, Double_t **a, Double_t *b, Int_t n) 
{
	Double_t temp;
	Double_t **u = new Double_t *[n];
	for(Int_t i = 0; i < n; i++)
	{
		u[i] = new Double_t [n];
		for(Int_t j = 0; j < n; j++) u[i][j] = 0.;
	}

	Double_t *y = new Double_t[n];
	for(Int_t i = 0; i < n; i++) y[i] = 0.;
 
	for(Int_t i = 0; i < n; i++)
	{
		temp = 0.;
		for(Int_t k = 0; k < i; k++)
			temp = temp + u[k][i] * u[k][i];
		u[i][i] = sqrt(a[i][i] - temp);
		for(Int_t j = i; j < n; j++)
		{
			temp = 0.;
			for(Int_t k = 0; k < i; k++)
				temp = temp + u[k][i] * u[k][j];
			u[i][j] = (a[i][j] - temp) / u[i][i];
		}
	}

	for(Int_t i = 0; i < n; i++)
	{
		temp = 0.;
		for(Int_t k = 0; k < i; k++)
			temp = temp + u[k][i] * y[k];
		y[i] = (b[i] - temp) / u[i][i];
	}

	for(Int_t i = n - 1; i >= 0; i--)
	{
		temp = 0.;
		for(Int_t k = i + 1; k < n; k++)
			temp = temp + u[i][k] * x[k];
		x[i] = (y[i] - temp) / u[i][i];
	}

	for(Int_t i = 0; i < n; i++) delete[] u[i];
	delete[] u;
	delete[] y;
}

void PronyFitter::SolveSLEOverDetermined_3x2_(Double_t *x, Double_t **ai, Double_t *b) 
{
	Int_t n = 3;
	Double_t **a = new Double_t *[n];
	for(Int_t i = 0; i < n; i++)
	{
		a[i] = new Double_t [n-1];
		for(Int_t j = 0; j < n-1; j++) a[i][j] = 0.;
	}

	for(Int_t i = 0; i < n; i++)
	{
		a[i][0] = ai[i][0];
		a[i][1] = ai[i][1] - ai[i][2] * pow(fz[1]/fz[2], fSignalBegin);
	}

	Double_t **at = new Double_t *[n-1];
	for(Int_t i = 0; i < n-1; i++)
	{
		at[i] = new Double_t [n];
		for(Int_t j = 0; j < n; j++) at[i][j] = a[j][i];
	}

	Double_t **ata = new Double_t *[n-1];
	for(Int_t i = 0; i < n-1; i++)
	{
		ata[i] = new Double_t [n-1];
		for(Int_t j = 0; j < n-1; j++) ata[i][j] = 0.;
	}

	Double_t tmpf;
	for(Int_t i = 0; i < n-1; i++)
	{
		for(Int_t j = 0; j < n-1; j++)
		{
			tmpf = 0.;
			for(Int_t k = 0; k < n; k++)
				tmpf += at[i][k]*a[k][j];
			ata[i][j] = tmpf;
		}
	}	

	Double_t **ata_inv = new Double_t *[n-1];
	for(Int_t i = 0; i < n-1; i++)
	{
		ata_inv[i] = new Double_t [n-1];
		for(Int_t j = 0; j < n-1; j++) ata_inv[i][j] = 0.;
	}

	Double_t det = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
	if(det != 0.)
	{
		ata_inv[0][0] = ata[1][1] / det;
		ata_inv[0][1] = -ata[0][1] / det;
		ata_inv[1][0] = -ata[1][0] / det;
		ata_inv[1][1] = ata[0][0] / det;			
	}	

	Double_t *atb = new Double_t[n-1];
	for(Int_t i = 0; i < n-1; i++) atb[i] = 0.;

	for(Int_t i = 0; i < n-1; i++)
	{
		tmpf = 0.;
		for(Int_t k = 0; k < n; k++)
			tmpf += at[i][k]*b[k];
		atb[i] = tmpf;
	}

	for(Int_t i = 0; i < n-1; i++)
	{
		tmpf = 0.;
		for(Int_t k = 0; k < n-1; k++)
			tmpf += ata_inv[i][k]*atb[k];
		x[i] = tmpf;
	}

	x[2] = -x[1] * pow(fz[1]/fz[2], fSignalBegin);

	if(fIsDebug)
	{
		printf("\n       AtA_inv                                At                                b\n");
		for(Int_t i =0; i < n-1; i++)
		{
			for(Int_t j = 0; j < n-1; j++)
			{
				printf("%e ", ata_inv[i][j]);
			}
			printf("     ");
			for(Int_t j = 0; j < n; j++)
			{
				printf("%e ", at[i][j]);
			}
			printf("     ");
			if(i == 0) for(Int_t j = 0; j < n; j++) {printf("%e ", b[j]);}
			printf("\n");
		}
	}

	for(Int_t i = 0; i < n; i++) delete[] a[i];
	delete[] a;
	for(Int_t i = 0; i < n-1; i++)
	{
		delete[] at[i];
		delete[] ata[i];
		delete[] ata_inv[i];
	}
	delete[] at;
	delete[] ata;
	delete[] ata_inv;
	delete[] atb;
}

void PronyFitter::DeleteData()
{
	delete[] fValSampleArr;
	delete[] fFITValSampleArr;
	delete[] fz;
	delete[] fh;
}

void PronyFitter::Clear()
{    
	fModelOrder = 0;
	fExpNumber = 0;
	fGateBeg = 0;
	fGateEnd = 0;
	fSampleTotal = 0;
	fZeroLevel = 0.;
	fSignalBegin = 0;
	fTotalPolRoots = 0;
	DeleteData();
}


