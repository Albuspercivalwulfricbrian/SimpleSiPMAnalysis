
#include<TTree.h>

struct Coeffs_struct
{
    Float_t peak_for_calibration;
    Float_t time_in_hours;

    static TString GetChName(Int_t channel_num)
    {
	    return TString::Format("channel_%i", channel_num);
    }

    TBranch* CreateBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->Branch(GetChName(channel_num).Data(), this, "peak_for_calibration/F:time_in_hours/F");
    }

    Int_t SetBranch(TTree *tree, Int_t channel_num)
    {
        return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
    
};
