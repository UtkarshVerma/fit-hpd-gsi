#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TList.h>

// extern TFile* _file0;
// extern TList* l;

void fit(unsigned int index) {
    gErrorIgnoreLevel = kError;  // Disable warnings temporarily

    auto l = (TList*)_file0->Get("H");
    auto h = (TH1D*)l->At(index);
    h->Draw();

    auto r = h->Fit("gaus", "SQ", "", 40, 80);
    if (r == -1) {
        fprintf(stderr, "error: first fit failed\n");
        return;
    }

    auto amp0 = r->Parameter(0);
    auto mean0 = r->Parameter(1);
    auto sigma0 = r->Parameter(2);
    auto f2 = new TF1("f2", "gaus + gaus(3) + gaus(6)");
    f2->SetParameter(0, amp0);
    f2->SetParameter(1, mean0 - 0.05);
    f2->SetParameter(2, 0.75 * sigma0);
    f2->SetParameter(3, amp0 / 5);
    f2->SetParameter(4, mean0 + 2 * sigma0);
    f2->SetParameter(5, 1.35 * sigma0);
    f2->SetParameter(6, amp0 / 14);
    f2->SetParameter(7, mean0 + 6 * sigma0);
    f2->SetParameter(8, 2.7 * sigma0);

    auto r2 = h->Fit("f2", "SQ", "", mean0 - 3 * sigma0, mean0 + 15 * sigma0);
    if (r2 == -1) {
        fprintf(stderr, "error: second fit failed\n");
    }

    auto mean1 = r2->Parameter(1);
    auto sigma1 = r2->Parameter(2);
    auto mean3 = r2->Parameter(7);
    auto sigma3 = r2->Parameter(8);

    h->SetAxisRange(mean1 - 3 * sigma1, mean3 + 3 * sigma3);
    gErrorIgnoreLevel = kWarning;
}
