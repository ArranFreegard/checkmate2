#ifndef ARRAN_2020_ATLAS_HIGG_2013_03_H_
#define ARRAN_2020_ATLAS_HIGG_2013_03_H_
// AUTHOR: Arran Freegard
//  EMAIL: acf1g14@soton.ac.uk
#include "AnalysisBase.h"

class Arran_2020_atlas_higg_2013_03 : public AnalysisBase {
  public:
    Arran_2020_atlas_higg_2013_03() : AnalysisBase()  {}               
    ~Arran_2020_atlas_higg_2013_03() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
};

#endif
