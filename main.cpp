#include <TArrayD.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH2D.h>
#include <TLeaf.h>
#include <TList.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

using std::string;

#define LEN(x) (sizeof(x) / sizeof(x[0]))
#define DIRICH_CHAN_COUNT 32
#define TRIGGER_CHAN 994

const unsigned int dirichStartChans[] = {1024, 1056};

void writeToFile(string filename, TFitResultPtr fitResults[][DIRICH_CHAN_COUNT],
                 double cuts[][DIRICH_CHAN_COUNT][2]) {
    FILE *file = fopen(filename.data(), "w");
    if (file == NULL) {
        fprintf(stderr, "\nerror: could not open file '%s'\n", filename.data());
    }

    fprintf(file,
            "Channel\tToT Cut[ns]\tMean[ns]\tSigma[ps]\tSigma Error[ps]\n");

    for (unsigned int i = 0; i < LEN(dirichStartChans); i++) {
        for (unsigned int j = 0; j < DIRICH_CHAN_COUNT; j++) {
            int channel = dirichStartChans[i] + j;
            auto fitResult = fitResults[i][j];
            if (fitResult == -1) continue;

            fprintf(file, "%d\t%f\t%f\t%f\t%f\n", channel, cuts[i][j][0],
                    fitResult->Parameter(1), fitResult->Parameter(2) * 1e3,
                    fitResult->ParError(2) * 1e3);
        }
    }

    fclose(file);
}

void scanDiRICH(unsigned int dirichIndex, double *fLeadEdgeBase,
                double *fTotBase, TH2D *hists[][DIRICH_CHAN_COUNT]) {
    unsigned int startChan = dirichStartChans[dirichIndex];
    double fLeadEdgeChan, fLeadEdgeTrigger, fTotChan;

    for (unsigned int i = 0; i < LEN(hists[0]); i++) {
        int chan = startChan + i;
        if (chan == TRIGGER_CHAN) continue;

        fLeadEdgeChan = *(fLeadEdgeBase + chan);
        fLeadEdgeTrigger = *(fLeadEdgeBase + TRIGGER_CHAN);
        fTotChan = *(fTotBase + chan);

        hists[dirichIndex][i]->Fill(fTotChan, fLeadEdgeChan - fLeadEdgeTrigger);
    }
}

int main(int argc, const char **argv) {
    if (argc != 2) {
        printf("Usage:\n\tfits pulser_data.root\n");
        return 1;
    }

    string ext = ".A.root";
    string inputFilename = argv[1];
    string inputStem =
        inputFilename.substr(0, inputFilename.length() - ext.length());

    // Disable buffering for STDOUT
    setbuf(stdout, NULL);

    // Disable referencing for histograms globally in favour of variables
    // Note that all histogram classes extend TH1
    TH1::AddDirectory(kFALSE);

    printf("Loading data");
    auto file = new TFile(inputFilename.data());
    if (file == NULL) {
        printf("error: could not open '%s'\n", inputFilename.data());
        return 1;
    }
    printf(" ...done\n");

    printf("Populating the 2D histogram");
    TH2D *hists[LEN(dirichStartChans)][DIRICH_CHAN_COUNT];
    for (unsigned int i = 0; i < LEN(hists); i++) {
        for (unsigned int j = 0; j < LEN(hists[0]); j++) {
            hists[i][j] = new TH2D("", "", 2000, 20, 50, 2000, 20, 60);
        }
    }

    auto tree = (TTree *)file->Get("A");
    auto leadEdgeBranch = tree->GetBranch("fLeadingEdge");
    auto leadEdgeLeaf = leadEdgeBranch->GetLeaf("fLeadingEdge");
    auto totBranch = tree->GetBranch("fTot");
    auto *totLeaf = totBranch->GetLeaf("fTot");

    double *leadEdgeBasePointer, *totBasePointer;
    long long nEntries = tree->GetEntriesFast();
    for (long long i = 0; i < nEntries; i++) {
        leadEdgeBranch->GetEntry(i);
        totBranch->GetEntry(i);
        leadEdgeBasePointer = (double *)leadEdgeLeaf->GetValuePointer();
        totBasePointer = (double *)totLeaf->GetValuePointer();

        for (unsigned int j = 0; j < LEN(dirichStartChans); j++) {
            scanDiRICH(j, leadEdgeBasePointer, totBasePointer, hists);
        }
    }
    delete tree;
    printf(" ...done\n");

    printf("Finding the cuts");
    double cuts[LEN(hists)][LEN(hists[0])][2];
    int maxBins[LEN(hists)][LEN(hists[0])][2];
    for (unsigned int i = 0; i < LEN(cuts); i++) {
        for (unsigned int j = 0; j < LEN(cuts[0]); j++) {
            auto hist = hists[i][j];

            auto xAxis = hist->GetXaxis();
            auto yAxis = hist->GetYaxis();
            int maxBin = hist->GetMaximumBin();
            int binX, binY, binZ;
            hist->GetBinXYZ(maxBin, binX, binY, binZ);

            cuts[i][j][0] = xAxis->GetBinCenter(binX);
            cuts[i][j][1] = yAxis->GetBinCenter(binY);
            maxBins[i][j][0] = binX;
            maxBins[i][j][1] = binY;
        }
    }
    printf(" ...done\n");

    printf("Making the cuts");
    TH1D *projHists[LEN(hists)][LEN(hists[0])];
    for (unsigned int i = 0; i < LEN(projHists); i++) {
        for (unsigned int j = 0; j < LEN(projHists[0]); j++) {
            int maxBinX = maxBins[i][j][0];
            const int nBinsCut = 2;
            projHists[i][j] = hists[i][j]->ProjectionY(
                "", maxBinX - nBinsCut / 2, maxBinX + nBinsCut / 2);
        }
    }
    printf(" ...done\n");

    printf("Fitting the histograms");
    // Only show fatal errors
    gErrorIgnoreLevel = kFatal;

    auto fitFunc = new TF1("fitFunc", "gaus");
    bool needsNewline = true;
    TFitResultPtr fitResults[LEN(projHists)][LEN(projHists[0])];
    for (unsigned int i = 0; i < LEN(fitResults); i++) {
        for (unsigned int j = 0; j < LEN(fitResults[0]); j++) {
            double mean = cuts[i][j][1];
            auto projHist = projHists[i][j];

            // Make the cut
            const double cutWidth = 2;
            projHist->SetAxisRange(mean - cutWidth / 2, mean + cutWidth / 2);
            fitFunc->SetParameter(1, mean);
            fitFunc->SetParameter(2, 0.1);
            fitResults[i][j] = projHist->Fit("fitFunc", "SNQ", "");
            if (fitResults[i][j] == -1) {
                if (needsNewline) {
                    fprintf(stderr, "\n");
                    needsNewline = false;
                }
                fprintf(stderr, "warning: fit failed for channel %d\n",
                        dirichStartChans[i] + j);
            }
        }
    }

    // Show warning messages
    gErrorIgnoreLevel = kWarning;
    printf(" ...done\n");

    printf("Writing computed data to ROOT file");
    string outRootFilename = inputStem;
    outRootFilename.append(".D.root");

    // const unsigned int arrSize = LEN(fitResults) * LEN(fitResults[0]);
    // double sigma[arrSize];
    auto outRootFile = new TFile(TString(outRootFilename), "RECREATE");
    // auto dataTree = new TTree("D", "Computed data");
    // dataTree->Branch("sigma", sigma, "sigma[64]/D");
    // dataTree->Branch("rms", rms, "rms[64]/D");

    // // Fill sigmas
    // for (unsigned int i = 0; i < LEN(fitResults); i++) {
    //     for (unsigned int j = 0; j < LEN(fitResults[0]); j++) {
    //         auto fitResult = fitResults[i][j];
    //         unsigned int index = LEN(fitResults[0]) * i + j;
    //         if (fitResult == -1) {
    //             continue;
    //         }
    //
    //         sigma[index] = fitResult->Parameter(2) * 1e3;
    //     }
    // }
    //
    //
    auto projHistList = new TList();
    for (unsigned int i = 0; i < LEN(hists); i++) {
        for (unsigned int j = 0; j < LEN(hists[0]); j++) {
            projHistList->Add(projHists[i][j]);
        }
    }

    // dataTree->Fill();
    // dataTree->Write("D");
    projHistList->Write("histList", TObject::kSingleKey);
    outRootFile->Close();
    printf(" ...done\n");

    printf("Writing fit data to TSV file");
    string outDataFilename = inputStem;
    outDataFilename.append(".tsv");
    writeToFile(outDataFilename, fitResults, cuts);
    printf(" ...done\n");

    file->Close();
    return 0;
}
