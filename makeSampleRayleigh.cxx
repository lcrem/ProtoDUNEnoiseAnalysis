#include "AveragePowerSpectrum.h"

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "lardataobj/RawData/RawDigit.h"
#include "AveragePowerSpectrum.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TColor.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace art;
using namespace std;
using namespace std::chrono;


string folder="/dune/data/users/lcremone/philUncompressed/";///dune/data/users/rodriges/protodune-noise/";
string filename="np04_signal_run004643_0001_dl1_evt0-5.txt";

int main(int argc, char *argv[]){

  Int_t run;
  string inFileName;
  Int_t maxChannelsToProcess=10;

  if((argc!=3)&&(argc!=4)){
    std::cerr << "Usage 1: " << argv[0] << " [irun] [maxChannels] (infile)" << std::endl;
    return 1;
  } else {
    run = atoi(argv[1]);
    maxChannelsToProcess = atoi(argv[2]);
    if (argc==4) inFileName += argv[3];
    else  inFileName +=  folder + "/" + filename;
  }

  cout << " Processing " << inFileName << endl;

  int offlineChannel;

  double x[6000];
  //  double y[6000];
  
  double dt = NOMINAL_SAMPLING_DELTAT;
  double deltaf = 2./NUM_SAMPLES;

  for (int i=0; i<NUM_SAMPLES; i++) x[i] = i*dt;


  AveragePowerSpectrum *avgPowSpec = new AveragePowerSpectrum("avgPowSpec", "ProtoDUNE average power spectrum");


  TGraph *sampleGraph;
  bool doneGraph=false;


  // while (getline(infile,line) && countlines<maxChannelsToProcess){

  //   stringstream ss (line.c_str());
  //   ss >> eventNumber >> offlineChannel;
    
  //   int countpoints=0;
  //   int adc=0;
  //   double mean=0;
  //   while ( ss >> adc ){
  //     y[countpoints] = adc;
  //     countpoints++;
  //     mean += adc;
  //   }

  //   // Remove offset
  //   for (int i=0; i<NUM_SAMPLES; i++) y[i]-=mean;
  //   TGraph *gtemp = new TGraph(NUM_SAMPLES, x, y);

  //   avgPowSpec->add(gtemp);
  //   cout << "Event and channel " << eventNumber << " " << offlineChannel << " " << countpoints << " " << NUM_SAMPLES << endl;

  //   if (doneGraph==false){
  //     sampleGraph = new TGraph(NUM_SAMPLES, x, y);
  //     doneGraph=true;
  //   }

  //   delete gtemp;

  //   countlines++;
  // }

  // infile.close();


  std::string filenameList="../filenameList.txt";
 
  int nevents=1;

  InputTag const daq_tag{ "tpcrawdecoder", "daq", "DecoderandReco" }; //RunRawDecoder" };
  //  InputTag daq_tag{ "daq" };

  // Create a vector of length 1, containing the given filename.                                                                                   
  double samples[6000];

  ifstream FileList(filenameList.c_str());
 std:string filename;
  vector<string> filenames;
  while(FileList >> filename) {
    filenames.push_back(filename);
  }
  
  std::cout << "Have got a vector of filenames with length " << filenames.size() << "\n";

  int iev=0;
  double mean=0;
  int countsamples=0;

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    
    //    if(iev>=nevents) break;
    std::cout << "Event " << iev << std::endl;
    
    if (iev<nevents){
    // Look at the digits                                                                                                                          

    auto& digits =
        *ev.getValidHandle<std::vector<raw::RawDigit>>(daq_tag);

    for(auto&& digit: digits){
      if(digit.Compression()!=0){
        std::cout << "Compression type " << digit.Compression() << std::endl;
      }
      // std::cout << iev << " " << digit.Channel() <<  " ";
      offlineChannel=digit.Channel();
      mean=0;
      countsamples=0;
      for(auto&& sample: digit.ADCs()){
	// std::cout << sample << " ";
	samples[countsamples]=(sample);
	mean+=sample;
	countsamples++;
      }
      mean/=(countsamples*1.);      
      
      for (int i=0; i<NUM_SAMPLES; i++) samples[i]-=mean;
      TGraph *gtemp = new TGraph(NUM_SAMPLES, x, samples);
      
      avgPowSpec->add(gtemp);
      cout << "Event and channel " << iev << " " << offlineChannel << " " << mean << " "<< countsamples << " " << NUM_SAMPLES << endl;
      
      if (doneGraph==false){
	sampleGraph = new TGraph(NUM_SAMPLES, x, samples);
	doneGraph=true;
      }
      
      delete gtemp;

    } // end loop over digits (=?channels)                                                                                                        
    

    }

    ++iev;
  } // end loop over events                                                                                                                        



  avgPowSpec->rebinAllRayleighHistograms(2); // emperically determined... for now.                                                                      
  avgPowSpec->fitAllRayleighHistograms();

  double powX[NUM_FREQS];
  Double_t rAmplitudes[NUM_FREQS];
  Double_t summedPowSpec[NUM_FREQS];
  Double_t rChiSquares[NUM_FREQS];
  Double_t rChiSquaresOverNdf[NUM_FREQS];
  Int_t rNdf[NUM_FREQS];
  Double_t rChiSquaresFullRange[NUM_FREQS];
  Int_t rNdfFullRange[NUM_FREQS];
  Double_t xHigh[NUM_FREQS];
  Double_t numOutliers[NUM_FREQS];

  memset(powX, 0, sizeof(powX));
  memset(rAmplitudes, 0, sizeof(rAmplitudes));
  memset(rChiSquares, 0, sizeof(rChiSquares));
  memset(rNdf, 0, sizeof(rNdf));
  memset(rChiSquaresFullRange, 0, sizeof(rChiSquaresFullRange));
  memset(rNdfFullRange, 0, sizeof(rNdfFullRange));


  for (int ifreq=0; ifreq<NUM_FREQS; ifreq++){

    powX[ifreq] = ifreq*deltaf;

    summedPowSpec[ifreq] = avgPowSpec->summedPowSpec[ifreq];
    rAmplitudes[ifreq] = avgPowSpec->rayleighAmplitudes[ifreq];
    rChiSquares[ifreq] = avgPowSpec->rayleighFitChiSquares[ifreq];
    rNdf[ifreq] = avgPowSpec->rayleighNdf[ifreq];
    rChiSquaresOverNdf[ifreq] = rChiSquares[ifreq]/rNdf[ifreq];
    rChiSquaresFullRange[ifreq] = avgPowSpec->rayleighFitChiSquaresFullRange[ifreq];
    rNdfFullRange[ifreq] = avgPowSpec->rayleighNdfFullRange[ifreq];
     
    xHigh[ifreq] = avgPowSpec->xHigh[ifreq];
    numOutliers[ifreq] = avgPowSpec->numOutliers[ifreq];
    cout << " What number is this " << powX[ifreq]  << " " << summedPowSpec[ifreq]  << " " << rAmplitudes[ifreq] << endl;
  } 


  gStyle->SetOptFit(1);
  
  TFile *fout = new TFile("temp.root", "recreate");
  sampleGraph->Write("sampleGraph");


  TGraph *gAvgPower = new TGraph(NUM_FREQS, powX, summedPowSpec);
  gAvgPower->SetTitle("Average power spectrum;Frequency [MHz];Amplitude");
  gAvgPower->Write("gAvgPowSpectrum");

  TGraph *gRayleighAmpl = new TGraph(NUM_FREQS, powX, rAmplitudes);
  gRayleighAmpl->SetTitle("Rayleigh Amplitude;Frequency [MHz];Amplitude [ADC/MHz]");
  gRayleighAmpl->Write("gRayleighAmpl");

  
  TGraph *gChiSquare = new TGraph(NUM_FREQS, powX, rChiSquaresOverNdf);
  gChiSquare->SetTitle("Rayleigh fit reduced chi square;Frequency [MHz];Chi square / NDF");
  gChiSquare->Write("gChiSquare");


  TDirectory *adir = gDirectory->mkdir("RayleighDistributions");
  adir->cd();

  for (int i=0; i<NUM_FREQS; i++){
    avgPowSpec->getRayleighHistogram(i)->Write();
    avgPowSpec->getRayleighHistogramFit(i)->Write();
  }

  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << std::endl;
  fout->Close();
  delete avgPowSpec;
  delete gAvgPower;
  delete gRayleighAmpl;
  delete gChiSquare;
  delete sampleGraph;
  delete fout;

  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << std::endl;

  return 0;

}

