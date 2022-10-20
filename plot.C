#include <DllImport.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TList.h>
#include <TStyle.h>
#include <TTree.h>

#include <cstdio>

#ifndef __CLING__
extern TFile *_file0;
extern TList *l;
#endif

void plot() {
    TCanvas *c1 = new TCanvas("c1", "multipads");
    const unsigned int rows = 4, cols = 3;
    c1->Divide(cols, rows, 0.01, 0.01);

    auto D = (TTree *)_file0->Get("D");

    c1->cd(1);
    D->Draw("sigmas[0]>>h1(10,0,500)");

    c1->cd(2);
    D->Draw("sigmas[1]>>h2(10,0,500)");

    c1->cd(3);
    D->Draw("sigmas[2]>>h3(10,0,500)");

    c1->cd(4);
    D->Draw("means[0]>>h4(10,0,100)");

    c1->cd(5);
    D->Draw("means[1]>>h5(10,0,100)");

    c1->cd(6);
    D->Draw("means[2]>>h6(10,0,100)");

    c1->cd(7);
    D->Draw("amplitudes[0]>>h7(10,0,200)");

    c1->cd(8);
    D->Draw("amplitudes[1]>>h8(10,0,200)");

    c1->cd(9);
    D->Draw("amplitudes[2]>>h9(10,0,200)");

    c1->cd(cols * (rows - 1) + 1);
    D->Draw("rms>>h10(100,0,2000)");
    return;

    // const char *branches[] = {"sigmas", "means", "amplitudes", "rms"};
    //
    // for (unsigned int i = 0; i < rows - 1; i++) {
    //     for (unsigned int j = 0; j < cols; j++) {
    //         c1->cd(i * cols + j + 1);
    //         auto branch = TString::Format("%s[%d]", branches[i], j);
    //         D->Draw(branch);
    //
    //         if (i == 0) {
    //             gPad->Range(0, 0, 100, 10);
    //             gPad->Draw();
    //         }
    //     }
    // }
}
