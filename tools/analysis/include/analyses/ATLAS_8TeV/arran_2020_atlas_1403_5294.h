#ifndef ARRAN_2020_ATLAS_1403_5294_H_
#define ARRAN_2020_ATLAS_1403_5294_H_
// AUTHOR: Arran Freegard
//  EMAIL: acf1g14@soton.ac.uk
#include "AnalysisBase.h"

class Arran_2020_atlas_1403_5294 : public AnalysisBase {
  public:
    Arran_2020_atlas_1403_5294() : AnalysisBase()  {}               
    ~Arran_2020_atlas_1403_5294() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
};

#endif
