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

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/range/detail/implementation_help.hpp>
#include <complex>
#include <cstdio>

using std::string;
using std::vector;
namespace po = boost::program_options;

#define LEN(x) (sizeof(x) / sizeof(x[0]))
#define DIRICH_CHAN_COUNT 32

vector<unsigned int> dirichStartChans;
unsigned int triggerChan;

void writeToTSVFile(TFitResultPtr fitResults[][DIRICH_CHAN_COUNT],
                    double cuts[][DIRICH_CHAN_COUNT][2], string filename) {
    FILE *file = fopen(filename.data(), "w");
    if (file == NULL) {
        fprintf(stderr, "\nerror: could not open file '%s'\n", filename.data());
    }

    fprintf(file,
            "Channel\tToT Cut[ns]\t"
            "Mean1[ns]\tSigma1[ps]\tSigmaErr1[ps]\t"
            "Mean2[ns]\tSigma2[ps]\tSigmaErr2[ps]\t"
            "Mean3[ns]\tSigma3[ps]\tSigmaErr3[ps]\n");

    for (unsigned int i = 0; i < dirichStartChans.size(); i++) {
        for (unsigned int j = 0; j < DIRICH_CHAN_COUNT; j++) {
            unsigned int chan = dirichStartChans[i] + j;
            auto r = fitResults[i][j];
            if (r == -1) continue;

            fprintf(file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", chan,
                    cuts[i][j][0], r->Parameter(1), r->Parameter(2) * 1e3,
                    r->ParError(2) * 1e3, r->Parameter(4),
                    r->Parameter(5) * 1e3, r->ParError(5) * 1e3,
                    r->Parameter(7), r->Parameter(8) * 1e3,
                    r->ParError(8) * 1e3);
        }
    }

    fclose(file);
}

void writeToRootFile(TFitResultPtr fitResults[][DIRICH_CHAN_COUNT],
                     string outFile) {
    const unsigned int arrSize = dirichStartChans.size() * LEN(fitResults[0]);
    double sigmaErrs[3][arrSize], sigmas[3][arrSize], means[3][arrSize];

    auto outRootFile = new TFile(TString(outFile), "RECREATE");
    auto dataTree = new TTree("D", "Computed data");
    dataTree->Branch("sigmas", sigmas,
                     TString::Format("sigmas[3][%d]/D", arrSize));
    dataTree->Branch("means", means,
                     TString::Format("means[3][%d]/D", arrSize));

    // Fill sigmas
    for (unsigned int i = 0; i < dirichStartChans.size(); i++) {
        for (unsigned int j = 0; j < LEN(fitResults[0]); j++) {
            auto fitResult = fitResults[i][j];
            unsigned int index = LEN(fitResults[0]) * i + j;
            if (fitResult == -1) continue;

            for (unsigned int j = 0; j < 3; j++) {
                means[j][index] = fitResult->Parameter(3 * j + 1);
                sigmas[j][index] = fitResult->Parameter(3 * j + 2) * 1e3;
                sigmaErrs[j][index] = fitResult->ParError(3 * j + 2) * 1e3;
            }
        }
    }

    // auto projHistList = new TList();
    // for (unsigned int i = 0; i < LEN(hists); i++) {
    //     for (unsigned int j = 0; j < LEN(hists[0]); j++) {
    //         projHistList->Add(projHists[i][j]);
    //     }
    // }

    dataTree->Fill();
    dataTree->Write("D");
    // projHistList->Write("histList", TObject::kSingleKey);
    outRootFile->Close();
}

void fill2DHistograms(TFile *file, TH2D *hists[][DIRICH_CHAN_COUNT]) {
    for (unsigned int i = 0; i < dirichStartChans.size(); i++) {
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
    for (long long k = 0; k < nEntries; k++) {
        leadEdgeBranch->GetEntry(k);
        totBranch->GetEntry(k);
        leadEdgeBasePointer = (double *)leadEdgeLeaf->GetValuePointer();
        totBasePointer = (double *)totLeaf->GetValuePointer();

        for (unsigned int i = 0; i < dirichStartChans.size(); i++) {
            unsigned int startChan = dirichStartChans[i];
            double fLeadEdgeChan, fLeadEdgeTrigger, fTotChan;

            for (unsigned int j = 0; j < DIRICH_CHAN_COUNT; j++) {
                unsigned int chan = startChan + j;
                if (chan == triggerChan) continue;

                fLeadEdgeChan = *(leadEdgeBasePointer + chan);
                fLeadEdgeTrigger = *(leadEdgeBasePointer + triggerChan);
                fTotChan = *(totBasePointer + chan);

                hists[i][j]->Fill(fTotChan, fLeadEdgeChan - fLeadEdgeTrigger);
            }
        }
    }
    delete tree;
}

void findCuts(TH2D *hists[][DIRICH_CHAN_COUNT],
              double cuts[][DIRICH_CHAN_COUNT][2],
              int maxBins[][DIRICH_CHAN_COUNT][2]) {
    for (unsigned int i = 0; i < dirichStartChans.size(); i++) {
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
}

void cutHistograms(TH2D *hists[][DIRICH_CHAN_COUNT],
                   TH1D *projHists[][DIRICH_CHAN_COUNT],
                   int maxBins[][DIRICH_CHAN_COUNT][2]) {
    for (unsigned int i = 0; i < dirichStartChans.size(); i++) {
        for (unsigned int j = 0; j < LEN(projHists[0]); j++) {
            int maxBinX = maxBins[i][j][0];
            const int nBinsCut = 2;
            projHists[i][j] = hists[i][j]->ProjectionY(
                "", maxBinX - nBinsCut / 2, maxBinX + nBinsCut / 2);
        }
    }
}

void fitHistograms(TH1D *projHists[][DIRICH_CHAN_COUNT],
                   double cuts[][DIRICH_CHAN_COUNT][2],
                   TFitResultPtr fitResults[][DIRICH_CHAN_COUNT]) {
    gErrorIgnoreLevel = kError;  // Disable warnings temporarily

    bool needsNewline = true;
    for (unsigned int i = 0; i < dirichStartChans.size(); i++) {
        for (unsigned int j = 0; j < LEN(fitResults[0]); j++) {
            double mean = cuts[i][j][1];
            auto projHist = projHists[i][j];

            const double cutWidth = 4;
            projHist->SetAxisRange(mean - cutWidth / 2, mean + cutWidth / 2);
            auto fitResult = projHist->Fit("gaus", "SNQ", "");
            if (fitResult == -1) {
                if (needsNewline) {
                    fprintf(stderr, "\n");
                    needsNewline = false;
                }
                fprintf(stderr, "warning: initial fit failed for channel %d\n",
                        dirichStartChans[i] + j);
                continue;
            }

            // Avoiding fit mean guess because it is skewed rightwards
            double amp1 = fitResult->Parameter(0);
            double sigma1 = fitResult->Parameter(2);
            auto fitFunc =
                new TF1("fitFunc", "gaus + gaus(3) + gaus(6)", 0, 100);
            fitFunc->SetParameter(0, amp1);
            fitFunc->SetParameter(1, mean);
            fitFunc->SetParameter(2, sigma1);
            fitFunc->SetParameter(3, amp1 / 4);
            fitFunc->SetParameter(4, mean + sigma1);
            fitFunc->SetParameter(5, sigma1);
            fitFunc->SetParameter(6, amp1 / 8);
            fitFunc->SetParameter(7, mean + 3 * sigma1);
            fitFunc->SetParameter(8, sigma1);
            fitResult = projHist->Fit("fitFunc", "SNQ", "");
            if (fitResult == -1) {
                if (needsNewline) {
                    fprintf(stderr, "\n");
                    needsNewline = false;
                }
                fprintf(stderr, "warning: final fit failed for channel %d\n",
                        dirichStartChans[i] + j);
            }
            fitResults[i][j] = fitResult;
        }
    }

    gErrorIgnoreLevel = kWarning;  // Show warning messages
}

void printUsage() {
    printf(
        "This program takes in the PMT data from the DiRICH boards as a ROOT "
        "file. It then makes cut on the (fLeadingEdge - fLeadingEdge[trigger]) "
        "vs fTot banana. The resulting 1D histogram is then fit using three "
        "Gaussians.\n"
        "The fit results are saved in two forms - a ROOT file and a TSV "
        "file.\n\n"
        "Usage:\n"
        "\t" BINARY
        " -i <input ROOT file> -t <trigger channel> -d <dirich start "
        "channels>\n\n"
        "For example, the command below reads data from data.A.root and runs "
        "the fit routine for channels 1000-1032 and 1040-1072 with trigger "
        "being at channel 990. The final results are saved as data.tsv and "
        "data.D.root.\n"
        "\t" BINARY " -i data.A.root -t 990 -d 1000 1040\n\n");
}

int main(int argc, char **argv) {
    // Parse CLI arguments
    string inputFilename;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Display help message")(
        "input-file,i", po::value(&inputFilename)->required(),
        "Input ROOT file")("trigger-channel,t",
                           po::value(&triggerChan)->required(),
                           "Channel number of the trigger")(
        "dirich-start-channels,d",
        po::value<std::vector<std::string>>()->multitoken()->required(),
        "Starting channels of each DiRICH board being used");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    if (argc == 1 || vm.count("help")) {
        printUsage();
        desc.print(std::cout);
        return 0;
    }

    try {
        po::notify(vm);
    } catch (po::error &e) {
        fprintf(stderr, "error: %s\n", e.what());
        desc.print(std::cout);
        return 1;
    }

    if (vm.count("dirich-start-channels")) {
        for (auto &startChans :
             vm["dirich-start-channels"].as<vector<string>>()) {
            dirichStartChans.push_back(std::stoi(startChans));
        }
    }

    // Disable buffering for STDOUT
    setbuf(stdout, NULL);

    // Disable referencing for histograms globally in favour of variables
    // Note that all histogram classes extend TH1
    TH1::AddDirectory(kFALSE);

    printf("Loading data");
    auto file = new TFile(inputFilename.data());
    if (file == NULL) {
        printf("\nerror: could not open '%s'\n", inputFilename.data());
        return 1;
    }
    printf(" ...done\n");

    printf("Populating the 2D histogram");
    TH2D *hists[dirichStartChans.size()][DIRICH_CHAN_COUNT];
    fill2DHistograms(file, hists);
    printf(" ...done\n");

    printf("Finding the cuts");
    double cuts[LEN(hists)][LEN(hists[0])][2];
    int maxBins[LEN(hists)][LEN(hists[0])][2];
    findCuts(hists, cuts, maxBins);
    printf(" ...done\n");

    printf("Making the cuts");
    TH1D *projHists[LEN(hists)][LEN(hists[0])];
    cutHistograms(hists, projHists, maxBins);
    printf(" ...done\n");

    printf("Fitting the histograms");
    TFitResultPtr fitResults[LEN(projHists)][LEN(projHists[0])];
    fitHistograms(projHists, cuts, fitResults);
    printf(" ...done\n");

    string ext = ".A.root";
    string inputStem =
        inputFilename.substr(0, inputFilename.length() - ext.length());

    printf("Writing computed data to ROOT file");
    string outRootFilename = inputStem;
    outRootFilename.append(".D.root");
    writeToRootFile(fitResults, outRootFilename);
    printf(" ...done\n");

    printf("Writing fit data to TSV file");
    string outTSVFilename = inputStem;
    outTSVFilename.append(".tsv");
    writeToTSVFile(fitResults, cuts, outTSVFilename);
    printf(" ...done\n");

    file->Close();
    return 0;
}
