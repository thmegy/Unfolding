#ifndef LIKELIHOOHSVANMANAGER_H
#define LIKELIHOOHSVANMANAGER_H

#include <memory>
#include <string>
#include <vector>

class RooDataSet;
class RooWorkspace;

class LikelihoodScanManager {
public:

    struct Result2D {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<std::vector<double> > z;
    };

    explicit LikelihoodScanManager();
    ~LikelihoodScanManager();

    LikelihoodScanManager(const LikelihoodScanManager& m) = delete;
    LikelihoodScanManager(LikelihoodScanManager&& m) = delete;
    LikelihoodScanManager& operator=(const LikelihoodScanManager& m) = delete;
    LikelihoodScanManager& operator=(LikelihoodScanManager&& m) = delete;

    inline void SetScanParamsX(const double min, const double max, const int steps) {
        fScanMinX = min;
        fScanMaxX = max;
        fStepsX = steps;
    }
    
    inline void SetScanParamsY(const double min, const double max, const int steps) {
        fScanMinY = min;
        fScanMaxY = max;
        fStepsY = steps;
    }

    inline void SetNCPU(const int n){fCPU = n;}
    
    inline void SetOffSet(const bool flag){fUseOffset = flag;}
    
    inline void SetParalel2Dstep(const int step) {
        fParal2D = true;
        fParal2Dstep = step;
    }
    
    inline void SetBlindedParameters(const std::vector<std::string>& par){fBlindedParameters = par;}

    inline void SetUseNll(const bool flag) {fUseNllInLHscan = flag;}

    inline float GetMinValX() const {return fMinValX;}
    inline float GetMinValY() const {return fMinValY;}
    inline float GetMaxValX() const {return fMaxValX;}
    inline float GetMaxValY() const {return fMaxValY;}

    typedef std::pair<std::vector<double>, std::vector<double> > scanResult1D;

    scanResult1D Run1DScan(const RooWorkspace* ws,
                           const std::string& varName,
                           RooDataSet* data) const;

    Result2D Run2DScan(const RooWorkspace* ws,
                       const std::pair<std::string,std::string>& varNames,
                       RooDataSet* data);

private:
    double fScanMinX;
    double fScanMinY;
    int fStepsX;
    double fScanMaxX;
    double fScanMaxY;
    int fStepsY;
    bool fUseOffset;
    int fCPU;
    bool fParal2D;
    int fParal2Dstep;
    bool fUseNllInLHscan;
    std::vector<std::string> fBlindedParameters;
    float fMinValX;
    float fMinValY;
    float fMaxValX;
    float fMaxValY;
};

#endif
