/** @file   PronyFitter.h
    @class  PronyFitter
    @author Nikolay Karpushkin (nkarpushkin@mail.ru)
    @brief  Class to fit waveform using Prony least squares method
*/

#ifndef PronyFitter_H
#define PronyFitter_H 1

class PronyFitter
{
        
    public:
        
        /**   Default constructor   **/
	PronyFitter() {};
	PronyFitter(Int_t, Int_t, Int_t, Int_t, Int_t);
	
        /**   Default destructor   **/
        ~PronyFitter() {Clear();};

	Int_t CalcSignalBegin(Float_t, Float_t);
	Int_t ChooseBestSignalBeginHarmonics(Int_t, Int_t);
	Int_t ChooseBestSignalBegin(Int_t, Int_t);
	void MakePileUpRejection(Int_t);       
        void CalculateFitHarmonics();
        void CalculateFitAmplitudes();
        void SolveSLEGauss(Double_t *, Double_t **, Double_t *, Int_t);
        void SolveSLECholesky(Double_t *, Double_t **, Double_t *, Int_t);
	void SolveSLEOverDetermined_3x2_(Double_t *, Double_t **, Double_t *);
	void DrawFit(TObjArray *, TString);
	void DrawFit_for_presentation(TCanvas *, TString, TString, TString, TString);

	Float_t LevelBy2Points(Float_t, Float_t, Float_t, Float_t, Float_t);
//         
//         Setters
//       
	void SetDebugMode(Bool_t);
        void SetWaveform(Float_t *, Float_t);
        void SetSignalBegin(Int_t);
        void SetHarmonics(Double_t *);
        void SetExternalHarmonics(Double_t, Double_t);
//         
//         Getters
//         
        Double_t *GetHarmonics();
	Int_t GetNumberPolRoots();
	Double_t *GetAmplitudes();
	Float_t GetIntegral(Int_t, Int_t);
	Float_t GetFitValue(Int_t);
	Float_t GetZeroLevel();
	Float_t GetX(Float_t, Int_t, Int_t);
	Float_t GetX(Float_t, Int_t, Int_t, Float_t);
	Float_t GetRSquare(Int_t, Int_t);
	Float_t GetRSquareSignal();
	Float_t GetChiSquare(Int_t, Int_t, Int_t);
	Float_t GetDeltaInSample(Int_t);
	Float_t GetSignalBeginFromPhase();
	Float_t GetMaxAmplitude();

    private:
        void Initialize(Int_t, Int_t, Int_t, Int_t, Int_t);
	void AllocData();
	void DeleteData();
	void Clear();

	Bool_t fIsDebug = false;
	Int_t fModelOrder;
	Int_t fExpNumber;
	Int_t fGateBeg;
	Int_t fGateEnd;
	Int_t fPileUpBeg;
	Int_t fSampleTotal;

	Float_t *fValSampleArr;
	Float_t fZeroLevel;
	Int_t fSignalBegin;
	Int_t fTotalPolRoots;

	Double_t *fz;
	Double_t *fh;
	Float_t *fFITValSampleArr;
};

#endif

