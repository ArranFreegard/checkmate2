#ifndef ATLAS_1403_5294_v2_H_
#define ATLAS_1403_5294_v2_H_
// AUTHOR: Arran Freegard
//  EMAIL: acf1g14@soton.ac.uk
#include "AnalysisBase.h"

class Atlas_1403_5294_v2 : public AnalysisBase {
  public:
    Atlas_1403_5294_v2() : AnalysisBase()  {}               
    ~Atlas_1403_5294_v2() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
};

#endif
