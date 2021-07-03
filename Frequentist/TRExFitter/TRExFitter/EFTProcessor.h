#ifndef EFTPROCESSOR_H
#define EFTPROCESSOR_H

/// Framework includes
/// ROOT includes
#include "Rtypes.h"

/// c++ includes
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

/// Forward class declaration
class NormFactor;
class Region;

class EFTProcessor {
public:


    explicit EFTProcessor();
    explicit EFTProcessor(std::string param, std::string title);
    ~EFTProcessor();


    class EFTExtrap {
    public:
      explicit EFTExtrap();
      int nBins;
      std::vector<double> p0s;
      std::vector<double> p1s;
      std::vector<double> p2s;
    };

    void Print() const;
    void DrawEFTInputs(std::vector < Region* > Regions) const;
    void FitEFTInputs(std::vector < Region* > Regions, const std::string& fileName);
    void ApplyMuFactExpressions(std::vector < Region* > Regions,
                                std::vector < std::shared_ptr<NormFactor> > &NormFactors);
    void ReadEFTFitResults(std::vector < Region* > Regions, const std::string& fileName);

    std::string fName;
    std::string fTitle;
    std::map<std::string,std::vector<std::string> > fRefMap;
    std::map<std::string,std::unique_ptr<EFTExtrap> > fExtrapMap;
};

#endif


