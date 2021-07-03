from array import array
import os
import hjson
from utils import writejson, getBinnedMean, makeResponseMatrix, responseTimesTruth
import numpy as np
import sys




def calculate_bin_edges(binning, granu):

    bounds = []
    bounds.append(0)
    
    for j in range(len(binning)-1):
        if j==0:
            bounds.append(abs(binning[j]-binning[j+1]) / granu)
        else:
            bounds.append(bounds[-1] + abs(binning[j]-binning[j+1]) / granu)

    ## Ugly workaround to remove pathologic float values                                                                                
    for ir,r in enumerate(bounds):
        if int(r)!=r:
            if abs(int(r)-r) < 0.000001:
                bounds[ir] = int(r)
            elif abs(int(r)+1-r) < 0.000001:
                bounds[ir] = int(r)+1
    
    return [int(b) for b in bounds]




def rebin_hist(hist, binning, name):
    import ROOT

    granu = float( (hist.GetXaxis().GetXmax()-hist.GetXaxis().GetXmin()) / hist.GetNbinsX() )
    bin_edges = calculate_bin_edges(binning, granu) # Get bin edges in histogram to rebin
        
    newh = ROOT.TH1D('newh_'+name, 'newh', len(binning)-1, array('d', binning))

    for i in range(len(binning)-1):
        error = ROOT.Double()
        content = hist.IntegralAndError(1+bin_edges[i], bin_edges[i+1], error)
        newh.SetBinContent(i+1, content)
        newh.SetBinError(i+1, error)
            
    return newh




def rebin_2Dhist(hist, binning_truth, binning_reco, name):
    import ROOT

    granu_truth = float( (hist.GetYaxis().GetXmax()-hist.GetYaxis().GetXmin()) / hist.GetNbinsY() )
    bin_edges_truth = calculate_bin_edges(binning_truth, granu_truth) # Get bin edges in histogram to rebin

    granu_reco = float( (hist.GetXaxis().GetXmax()-hist.GetXaxis().GetXmin()) / hist.GetNbinsX() )
    bin_edges_reco = calculate_bin_edges(binning_reco, granu_reco) # Get bin edges in histogram to rebin

    newh = ROOT.TH2D('newh_'+name, 'newh', len(binning_reco)-1, array('d', binning_reco), len(binning_truth)-1, array('d', binning_truth))

    for i in range(len(binning_reco)-1):
        for j in range(len(binning_truth)-1):
            error = ROOT.Double()
            content = hist.IntegralAndError(1+bin_edges_reco[i], bin_edges_reco[i+1], 1+bin_edges_truth[j], bin_edges_truth[j+1], error)
            newh.SetBinContent(i+1, j+1, content)
            newh.SetBinError(i+1, j+1, error)

    return newh




def getTruthParameter(config, outdir, outfile, hist_list, var, suff, binning_list):
    import ROOT

    truth_parameters = {}
    truth_parameters_unbinned = {}

    factor = config['factor'][var]
    rfile = ROOT.TFile(outfile, 'READ')

    for h in hist_list:
        hist = rfile.Get(h)
        truth_bins = np.array( [hist.GetBinContent(i+1) for i in range(hist.GetNbinsX())] )
        passed_bins = 0

        for i, binning in enumerate(binning_list):
            nbin = len(binning)-1

            param_unbinned = (factor * hist.GetMean(), abs(factor * hist.GetMeanError()))
            param_binned = factor * getBinnedMean(truth_bins[passed_bins:passed_bins+nbin], binning)

            if i == 0:
                truth_parameters_unbinned[h] = [param_unbinned]
                truth_parameters[h] = [param_binned]
            else:
                truth_parameters_unbinned[h].append(param_unbinned)
                truth_parameters[h].append( param_binned )

            passed_bins += nbin 
    rfile.Close()

    writejson('{}/truth_{}{}_binned.json'.format(outdir, var, suff), truth_parameters)
    print('\nThe binned truth values of {0} for the nominal and reweighted signals has been written to {1}/truth_{0}{2}_binned.json'.format(var, outdir, suff))
    writejson('{}/truth_{}{}_unbinned.json'.format(outdir, var, suff), truth_parameters_unbinned)
    print('The unbinned truth values of {0} for the nominal and reweighted signals has been written to {1}/truth_{0}{2}_unbinned.json'.format(var, outdir, suff))


        

def main(args):
    import ROOT

    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))

    if rebin_config['differential_measurement']:
        # flatten binning list and translate to match input range (0 to nbin_input*ndiff)
        binning_truth = [ (val+2*ir+1)*rebin_config['nbin_input']/2 for ir,r in enumerate(rebin_config['binning']) for val in r ]
        binning_truth = np.unique(np.array(binning_truth))
        suff = f'_{rebin_config["diff_variable"]}_diff'
        # check if provided binning is compatible with provided number of differential bins
        if len(rebin_config['diff_bins']) != len(rebin_config['binning']):
            print('\nThe provided binning is incompatible with the differential binning. Please inspect your rebin_config file!\n')
            sys.exit()
    else:
        binning_truth = rebin_config['binning'][0]
        suff = ''
        if len(rebin_config['binning']) != 1:
            print('\nThe provided binning is incompatible with the differential binning. Please inspect your rebin_config file!\n')
            sys.exit()

    syst_list_tot = []
        
    for camp in args.campaign:
        for chan in args.channel:
            for btag in config['btag_regions'][chan]:
                outdir = '{}/{}/{}/{}'.format(args.outpath, camp, chan, btag)
                os.makedirs(outdir, exist_ok=True)

                # check if channel is used as control region
                if rebin_config['is_control']==True and chan in rebin_config['control_channels'].keys():
                    is_control = True
                    bounds = rebin_config['control_variable_bounds'][rebin_config['control_channels'][chan]]
                    binning_reco = np.linspace(bounds[0], bounds[1], rebin_config['nbin_control']+1).astype(np.int)
                else:
                    is_control = False
                    binning_reco = binning_truth

                for var in args.variable:

                    print('\nRebinning file {}/{}/{}/{}{}.root'.format(camp, chan, btag, var, suff))
                    rfile = ROOT.TFile('{}/{}/{}/{}/{}{}.root'.format(args.inpath, camp, chan, btag, var, suff), 'READ')

                    # rebin bootstraps if any
                    bootstrap_list = []
                    if rebin_config['bootstrap']:
                        for ib in range(rebin_config['n_bootstrap']):
                            bootstrap_list.append('{}_{}'.format(rebin_config['bootstrap_name'], ib))
                
                    newh_nom = []
                    for h in rebin_config['histo_1D']+config['background']+bootstrap_list:
                        try:
                            hist = rfile.Get(h)
                            newh_nom.append( (rebin_hist(hist, binning_reco, h), h) )
                        except:
                            print(f'Issue with {h} histogram!!')

                    # get truth distribution and rebin
                    truth = rfile.Get('truth')
                    truth_rebin = rebin_hist(truth, binning_truth, 'truth')
                    newh_nom.append( (truth_rebin, 'truth') )
                    # get migration matrix and rebin
                    migmat = rfile.Get(config['migration_matrix_name'])
                    migmat_rebin = rebin_2Dhist(migmat, binning_truth, binning_reco, 'nom')
                    newh_nom.append( (migmat_rebin, config['migration_matrix_name']) )
                    migmat_norm = migmat_rebin.Integral()
                    # compute response matrix
                    resmat = makeResponseMatrix(truth_rebin, migmat_rebin)
                    newh_nom.append( (resmat, config['response_matrix_name']) )


                    newh_syst = {}
                    for proc in ['all']+config['background']+['ttbar']:
                    
                        syst_list = []
                        # get two-sided systematics
                        try:
                            syst_list += [ s[0] for s in config['twosided_syst'][proc] ]
                            syst_list += [ s[1] for s in config['twosided_syst'][proc] ]
                        except:
                            print('No two-sided systematics for {}'.format(proc))
                        
                        # get one-sided systematics
                        try:
                            syst_list += config['onesided_syst'][proc]
                        except:
                            print('No one-sided systematics for {}'.format(proc))

                        # get reference samples
                        try:
                            syst_list += list( set( config['reference'][proc].values() ) )
                        except:
                            print('No reference sample for {}'.format(proc))

                        if len(syst_list) > 0:
                            syst_list_tot.append( (proc, syst_list) )
                            
                        for syst in syst_list:
                            newh = []
                            print(f'\n{syst}')

                            if proc == 'all' or proc == 'ttbar':
                                # varied migration matrix
                                migmat = rfile.Get('{}/{}'.format(syst, config['migration_matrix_name']))
                                migmat_rebin = rebin_2Dhist(migmat, binning_truth, binning_reco, syst) 
                                # varied truth
                                try:
                                    truth_syst = rfile.Get(f'{syst}/truth')
                                    truth_syst_rebin = rebin_hist(truth_syst, binning_truth, f'truth_{syst}')
                                except:
                                    print('Using nominal truth!')
                                    truth_syst_rebin = truth_rebin
                            
                                # compute varied response matrix
                                # then normalise it in order for nominal and systematically varied reco distributions to have the same normalisation (if shape-only syst requested)
                                resmat = makeResponseMatrix(truth_syst_rebin, migmat_rebin)
                                if args.syst == 'shape':
                                    resmat.Scale( migmat_norm / responseTimesTruth(truth_rebin, resmat) )
                                newh.append( (resmat, config['response_matrix_name']) )

                                if proc == 'all':
                                    for h in config['background']:
                                        print(h)
                                        hist = rfile.Get('{}/{}'.format(syst, h))
                                        newh.append( (rebin_hist(hist, binning_reco, f'{h}_{syst}'), h) )

                            else:
                                hist = rfile.Get('{}/{}'.format(syst, proc))
                                newh.append( (rebin_hist(hist, binning_reco, f'{proc}_{syst}'), h) )

                            newh_syst[syst] = newh

                    outfile_name = '{}/{}{}.root'.format(outdir, var, suff)
                    outfile = ROOT.TFile(outfile_name, 'RECREATE')
                    outfile.cd()
                    for h in newh_nom:
                        h[0].Write(h[1])
                    for syst, hlist in newh_syst.items():
                        tdir = ROOT.TDirectoryFile(syst, syst)
                        tdir.cd()
                        for h in hlist:
                            h[0].Write(h[1])
                        tdir.SaveSelf()
                        outfile.cd()

                    outfile.Close()
                    rfile.Close()

                    print('The rebinned input file has been written to {}/{}{}.root'.format(outdir, var, suff))


                    # Extract binned and unbinned truth parameter (for nominal and reweighted signal)
                    hist_list = [h for h in rebin_config['histo_1D'] if 'truth' in h] + ['truth']
                    getTruthParameter(config, outdir, outfile_name, hist_list, var, suff, rebin_config['binning'])


                
   
if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter

    # get path for root files
    parser = ArgumentParser(description='Rebin input histograms', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inpath', type=str, help='Path to input', required=True)
    parser.add_argument('-o', '--outpath', type=str, help='Path to input', default='unfoldingInput/base/')
    parser.add_argument('--config', type=str, help='Config file', default='config/config.json')
    parser.add_argument('-v', '--variable', nargs='*', default=['CosThetaKplus', 'CosThetaNplus', 'CosThetaRplus', 'CosThetaKminus', 'CosThetaNminus', 'CosThetaRminus', 'CorrKK', 'CorrNN', 'CorrRR', 'CorrNKplus', 'CorrNKminus', 'CorrNRplus', 'CorrNRminus', 'CorrRKplus', 'CorrRKminus'], help='Observables to prepare the inputs for')
    parser.add_argument('-c', '--channel', nargs='*', default=['em', 'ee', 'mm'], help='Channels to prepare the inputs for')
    parser.add_argument('--campaign', nargs='*', default=['mc16a', 'mc16d', 'mc16e', 'all'], help='Campaigns to prepare the inputs for')
    parser.add_argument('-s', '--syst', choices=['full', 'shape'], help='Shape-only systematics or shape+syst (full)', default='full')
    args = parser.parse_args()

    main(args)
