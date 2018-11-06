#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "lardataobj/RawData/RawDigit.h"
#include "AveragePowerSpectrum.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TAxis.h"

using namespace art;
using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]){


  std::string filenameList="../filenameList.txt";
 
  int nevents=1;

  InputTag daq_tag{ "tpcrawdecoder", "daq", "RunRawDecoder" };

  // Create a vector of length 1, containing the given filename.                                                                                   
  int eventNum;
  int channel;
  std::vector<int> samples;

  ifstream FileList(filenameList.c_str());
 std:string filename;
  vector<string> filenames;
  while(FileList >> filename) {
    filenames.push_back(filename);
  }
  
  std::cout << "Have got a vector of filenames with length " << filenames.size() << "\n";

  int iev=0;
  
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    
    if(iev>=nevents) break;
    std::cout << "Event " << iev << std::endl;
    
    // Look at the digits                                                                                                                          

    auto& digits =
        *ev.getValidHandle<std::vector<raw::RawDigit>>(daq_tag);

    for(auto&& digit: digits){
      if(digit.Compression()!=0){
        std::cout << "Compression type " << digit.Compression() << std::endl;
      }
      // std::cout << iev << " " << digit.Channel() <<  " ";
      channel=digit.Channel();
      samples.clear();
      for(auto&& sample: digit.ADCs()){
	// std::cout << sample << " ";
	samples.push_back(sample);
      }
      

    } // end loop over digits (=?channels)                                                                                                        
    
    ++iev;
  } // end loop over events                                                                                                                        
  
  
}

 


