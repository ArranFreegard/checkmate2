#ifndef ATLAS_2020_1403_5294_H_
#define ATLAS_2020_1403_5294_H_
// AUTHOR: Arran Freegard
//  EMAIL: acf1g14@soton.ac.uk
#include "AnalysisBase.h"

class Atlas_2020_1403_5294 : public AnalysisBase {
  public:
    Atlas_2020_1403_5294() : AnalysisBase()  {}               
    ~Atlas_2020_1403_5294() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
};

#endif
