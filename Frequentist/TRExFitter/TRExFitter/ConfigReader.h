#ifndef CONFIGREADER_H
#define CONFIGREADER_H

/// c++ includes
#include <memory>
#include <string>
#include <vector>

///Forward class declaration
class ConfigParser;
class ConfigSet;
class TRExFit;
class Region;
class Sample;
class Systematic;

/**
 * \class ConfigReader
 * \brief Class that reads config parameters
 * \author Tomas Dado
 */

class ConfigReader {
    public:

        /**
          * The default constructor.
          * @param fitter A pointer to TRExFit class
          */
        explicit ConfigReader(TRExFit* fitter);

        /**
          * The default destructor
          */
        ~ConfigReader() = default;

        /**
          * Deleted constructors and operators
          */
        ConfigReader() = delete;
        ConfigReader(const ConfigReader& c) = delete;
        ConfigReader(ConfigReader&& c) = delete;
        ConfigReader& operator=(const ConfigReader& c) = delete;
        ConfigReader& operator=(ConfigReader&& c) = delete; 

        /**
          * Reads the config and passes parameters to TRExFit
          * @param string Config path
          * @param string Running options
          * @param string Additional options
          * @return int status code
          */
        int ReadFullConfig(const std::string& fileName, const std::string& opt, const std::string& option);

    private:

        /**
          * Helper function to read Command line settings
          * @param string config options
          * @return int status code
          */
        int ReadCommandLineOptions(const std::string& option);

        /**
          * Helper function to read JOB settings
          * @return int status code
          */
        int ReadJobOptions();

        /**
          * Helper function to read plotting settings in JOB
          * @param ConfigSet A pointer needed to parse the config
          * @return int status code
          */
        int SetJobPlot(ConfigSet *confSet);

        /**
          * Helper function to read plotting settings in Options
          * @return int status code
          */
        int ReadGeneralOptions();

        /**
          * Helper function to read Fit settings
          * @return int status code
          */
        int ReadFitOptions();

        /**
          * Helper function to read Limit settings
          * @return int status code
          */
        int ReadLimitOptions();

        /**
          * Helper function to read Significance settings
          * @return int status code
          */
        int ReadSignificanceOptions();

        /**
          * Helper function to read Region settings
          * @param Running options
          * @return int status code
          */
        int ReadRegionOptions(const std::string& opt);

        /**
          * Helper function to read Region settings based on input type
          * @param Region A pointer needed to parse the input
          * @param ConfigSet A pointer neede to parse the input
          * @return int status code
          */
        int SetRegionHIST(Region* reg, ConfigSet *confSet);

        /**
          * Helper function to read Region settings based on input type
          * @param Region A pointer needed to parse the input
          * @param ConfigSet A pointer neede to parse the input
          * @return int status code
          */
        int SetRegionNTUP(Region* reg, ConfigSet *confSet);

        /**
          * Helper function to check if config has settings for NTUP
          * @param ConfigSet A pointer needed to parse the input
          * @return bool
          */
        bool ConfigHasNTUP(ConfigSet* confSet);

        /**
          * Helper function to check if config has settings for HIST
          * @param ConfigSet A pointer needed to parse the input
          * @return bool
          */
        bool ConfigHasHIST(ConfigSet* confSet);

        /**
          * Helper function to read Sample settings
          * @return int status code
          */
        int ReadSampleOptions();

        /**
          * Helper function to read NormFactor settings
          * @return int status code
          */
        int ReadNormFactorOptions();

        /**
          * Helper function to read ShapeFactor settings
          * @return int status code
          */
        int ReadShapeFactorOptions();

        /**
          * Helper function to read Systematic settings
          * @return int status code
          */
        int ReadSystOptions();

        /**
          * A helper function to read Unfolding settings
          * @return int status code
          */ 
        int ReadUnfoldingOptions();
        
        /**
          * A helper function to read TruthSample
          * @return int status code
          */ 
        int ReadTruthSamples();

        /**
          * A helper function to read UnfoldingSample
          * @return int status code
          */ 
        int ReadUnfoldingSampleOptions();

        /**
          * A helper function to read UnfoldingSystematic
          * @return int status code
          */ 
        int ReadUnfoldingSystematicOptions();

        /**
          * Helper function to read Part of Syst config
          * @param COnfigSet A pointer needed for reading
          * @param Systematic A pointer to syst that is being set
          * @param vector of strings Needed for this setting
          * @param vector of strings Needed for this setting
          * @return int status code
          */
        int SetSystNoDecorelate(ConfigSet *confSet, std::shared_ptr<Systematic> sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude);

        /**
          * Helper function to read Part of Syst config
          * @param COnfigSet A pointer needed for reading
          * @param Systematic A pointer to syst that is being set
          * @param vector of strings Needed for this setting
          * @param vector of strings Needed for this setting
          * @param vector of strings Needed for this setting
          * @param int Flag needed for for this setting
          * @return int status code
          */
        int SetSystRegionDecorelate(ConfigSet *confSet,
                                    std::shared_ptr<Systematic> sys,
                                    const std::vector<std::string>& samples,
                                    const std::vector<std::string>& exclude,
                                    const std::vector<std::string>& regions,
                                    int type);

        /**
          * Helper function to read Part of Syst config
          * @param COnfigSet A pointer needed for reading
          * @param Systematic A pointer to syst that is being set
          * @return int status code
          */
        int SetSystSampleDecorelate(ConfigSet *confSet, std::shared_ptr<Systematic> sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude);

        /**
          * Helper function to read Part of Syst config
          * @param COnfigSet A pointer needed for reading
          * @param Systematic A pointer to syst that is being set
          * @return int status code
          */
        int SetSystShapeDecorelate(ConfigSet *confSet, std::shared_ptr<Systematic> sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude);

        /**
          * Helper function that propagates samples and systematics when Unfolding is used
          * @return int status code
          */
        int UnfoldingCorrections();

        /**
          * Helper function that is run after config is read
          * @param Running options
          * @return int status code
          */
        int PostConfig(const std::string& opt);

        /**
          * Helper function to check if elements of one vector are present in another
          * @param vector of parameters to check
          * @param vector of paramaeters to check to
          * @return True if all exist, False if at least one does not exist
          */
        bool CheckPresence(const std::vector<std::string> &v1, const std::vector<std::string> &v2);

        /**
          * Helper function to check if elements of one vector are present in another
          * @param vector of parameters to check
          * @param vector of paramaeters to check to
          * @param vector of paramaeters to check to
          * @return True if all exist, False if at least one does not exist
          */
        bool CheckPresence(const std::vector<std::string> &v1, const std::vector<std::string> &v2, const std::vector<std::string> &v3);

        /**
          * A helper function to check of the same uses ghost samples with the same name
          * @param Sample
          * return true if ok
          */ 
        bool SampleIsOk(const Sample* sample) const;

        /**
          * Helper function to check if the regions from command line exist
          * @return A verctor of region names
          */
        std::vector<std::string> GetAvailableRegions();

        /**
          * Helper function to check if the samples from command line exist
          * @return A verctor of sample names
          */
        std::vector<std::string> GetAvailableSamples();

        /**
          * Helper function to check if the systematics from command line exist
          * @return A verctor of systematic names
          */
        std::vector<std::string> GetAvailableSysts();

        /**
          * Helper function to check if the name of the systematic is problematic
          * @param name of the syst
          * @return flag if the systematic is problematic
          */
        bool SystHasProblematicName(const std::string& name);

        /**
          * Helper function to convert UnfoldingSample to Samples
          * @return status code
          */
        int ProcessUnfoldingSamples();
        
        /**
          * Helper function to convert UnfoldingSystematics to Systematics
          * @return status code
          */
        int ProcessUnfoldingSystematics();
        
        /**
          * Helper function to add NormFactors for truth bins when running unfolding
          * @return status code
          */
        int AddUnfoldingNormFactors();

        /**
          * A helper function to fix reference samples at the end of the config
          * @param Given systematic
          */
        void FixReferenceSamples(std::shared_ptr<Systematic> sys) const;

        /**
          * A helepr function to check if POIs are properly set
          */ 
        int CheckPOIs() const;

        /**
          * Pointer to TRExFit class, set during initialization
          */
        TRExFit *fFitter;

        /**
          * Pointer to ConfigParser used to parse the text
          */
        std::unique_ptr<ConfigParser> fParser;

        /**
          * flag to control if wrong samples/regions are ok
          */
        bool fAllowWrongRegionSample;

        /**
          * flag to control if other than ghost sampels have been set already
          */
        bool fNonGhostIsSet;

        /**
          * vector of strings, one for each sample, needed for cross-checks
          */
        std::vector< std::string > fSamples;

        /**
          * vector of strings, one for each sample, needed for cross-checks
          */
        std::vector< std::string > fGhostSamples;

        /**
          * vector of strings, one for each sample, needed for cross-checks
          */
        std::vector< std::string > fEFTSamples;

        /**
          * vector of strings, one for each sample defined in config file
          */
        std::vector< std::string > fAvailableSamples;

        /**
          * vector of strings, one for each region, needed for cross-checks
          */
        std::vector< std::string > fRegions;

        /**
          * vector of strings, one for each region defined in config file
          */
        std::vector< std::string > fAvailableRegions;

        /**
          * vector of strings, one for each region
          */
        std::vector< std::string > fOnlyRegions;

        /**
          * vector of strings, one for each sample
          */
        std::vector< std::string > fOnlySamples;

        /**
          * vector of strings, one for each systematics
          */
        std::vector< std::string > fOnlySystematics;

        /**
          * vector of strings, one for each exclude region
          */
        std::vector< std::string > fToExclude;

        /**
          * vector of strings, one for each exclude region sample
          */
        std::vector< std::string > fExcludeRegionSample;
        /**
          * vector of names, one for each region
          */
        std::vector<std::string> fRegNames;

        /**
          *  vector of strings for signal only
          */
        std::vector<std::string> fOnlySignals;

        /**
          *  string for LH scan values from command line
          */
        std::string fOnlyLHscan;

        /**
          * bool to check if there is at least one valid region
          */
        bool fHasAtLeastOneValidRegion;

        /**
          * bool to check if there is at least one valid sample for the fit
          */
        bool fHasAtLeastOneValidSample;

        /**
          * A container for Tau parameters
          */
        std::vector<std::pair<int, double> > fTaus;

        /**
          * A container for expressions (used to set relationships between unflding parameters)
          */
        std::vector<std::pair<std::string, std::string> > fExpressions;
};

#endif
