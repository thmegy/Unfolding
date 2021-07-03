#ifndef FITUTILS_H
#define FITUTILS_H

#include <map>
#include <memory>
#include <vector>

class FittingTool;
class NormFactor;
class RooArgSet;
class RooWorkspace;
class RooSimultaneous;
class TDirectory;
class TFile;

namespace RooStats {
    class ModelConfig;
}

namespace FitUtils {

/**
 * Implementing external constraints on the workspace
 * @param workspace
 * @param Fitting tool
 * @param the PDF
 * @param vector of norma factors
 */
void ApplyExternalConstraints(RooWorkspace* ws,
                              FittingTool* fitTool,
                              RooSimultaneous* simPdf,
                              const std::vector<std::shared_ptr<NormFactor> >& normFactors,
                              int type = 0);

/**
 * A helper function to set BinnedLikelihoodOptimisation
 * @param workspace
 */
void SetBinnedLikelihoodOptimisation(RooWorkspace* ws);

/**
 * A helper function to injects values of NPs
 * @param workspace
 * @param map of NP values
 */
void InjectGlobalObservables(RooWorkspace* ws, const std::map< std::string, double >& npValues);

/**
 * A helper function to set saturated model shapefactors to const
 * @param workspace
 */
void DisableSaturatedModel(RooWorkspace* ws);

/**
 * A helper function to set POI in a file
 * @param path to the file
 * @param WS name
 * @param POI
 */
void SetPOIinFile(const std::string& path, const std::string& wsName, const std::string& poi);

/**
 * A helper function to set POI in a file
 * @param ModelConfig
 * @return parameters
 */
std::vector<std::string> GetAllParameters(const RooStats::ModelConfig* mc);

/**
 * A helper function to get number of free parameters
 * @param ModelConfig
 * @return number of free parameter
 */
std::size_t NumberOfFreeParameters(const RooStats::ModelConfig* mc);

/**
 * A helper function to set all paramet to const except the specified ones
 * @param ModelConfig
 * @param names of the exceptions
 */
void FixAllParameters(RooStats::ModelConfig* mc, const std::vector<std::string>& except);

/**
 * A helper function to set all paramet to floating
 * @param ModelConfig
 */
void FloatAllParameters(RooStats::ModelConfig* mc, const bool saturatedModel);

/**
 * A helper function to propagate Expression using Roofit
 * @param Workspace
 * @param File with roofit results
 * @param name of the parameter to be calculated
 * @return Returns a vector of doubles with a following convention:
 *         - element 0: nominal
 *         - element 1: up
 *         - element 2: down
 *         - element 3: error
 * Return 0 size vector if parameter cannot be found
 */
std::vector<double> CalculateExpressionRoofit(RooWorkspace* ws,
                                              TFile* fitResultFile,
                                              const std::string& name);

/**
 * A helper function to disable all saturated model parameters in a WS 
 * The WS will be replaced with one without the saturated model parameters.
 * Also set the POI
 * @param path to the ROOT file
 * @param name of the input WS
 * @param name of the output WS
 * @param name of ther POI
 */
void DisableSaturatedModelFileAndSetPOI(const std::string& path,
                                        const std::string& wsName,
                                        const std::string& wsOutName,
                                        const std::string& poi);

/**
 * A helper function to get RooArgSet without some parameter 
 * @param MC
 * @param names of the morph parameters
 */
RooArgSet GetPOIsWithout(const RooStats::ModelConfig* mc, const std::vector<std::string>& names);
}

#endif
