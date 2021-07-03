import os
import hjson
from utils import getSpinParameter, plotResmat, plotUnfoldedDistribution, plotSpinParameter, rootToNumpy, getBinnedMean, getNPsFromFitResults, writejson, plotBrasilian
import numpy as np
import sys
import re


    
def writeConfig(args, i):
    '''
    Write TRExFitter configuration file.
    '''

    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    binning = rebin_config['binning']   # Get binning
    n_bin = sum([len(b)-1 for b in binning])
    labels = hjson.load(open(config['labels'], 'r'))

    factor = config['factor'][args.variable]  # factor to multiply the mean by

    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])
        diff_bins = rebin_config['diff_bins']
    else:
        diff_bins = ['inclusive']
    
    config_path = args.config.replace( re.findall('/\w+.json', args.config)[0], '' )
    cfgname = '{}/trex/{}_{}_{}_{}.config'.format(config_path, args.campaign, '_'.join(args.channel), file_name, args.data[i])
    os.makedirs(f'{config_path}/trex/', exist_ok=True)
    path = '{}/{}/{}/'.format(os.getcwd(), args.inpath, args.campaign)
    
    with open(cfgname, 'w') as f:
        f.write('Job: {}\n'.format(args.data[i]))
        f.write('  OutputDir: {}/{}_{}_{}\n'.format(args.outpath, args.campaign, '_'.join(args.channel), file_name))
        f.write('  Label: {}_{}_{}_{}\n'.format(args.campaign, '_'.join(args.channel), file_name, args.data[i]))
        if 'all' in args.campaign:
            lumi_campaign = args.campaign[:3]
        else:
            lumi_campaign = args.campaign[:5]
        f.write('  LumiLabel: {} fb^{{-1}}\n'.format(config['lumi'][lumi_campaign][0]))
        f.write('  HistoPath: {}\n'.format(path))
        f.write('  ResponseMatrixPath: {}\n'.format(path))
        writeConfigBlock(config, 'job', f)
        
        f.write('Fit: myFit\n')
        writeConfigBlock(config, 'fit', f)
        
        f.write('Unfolding: Unfolding\n')
        truth_path = '{}/{}/{}/'.format(path, args.channel[0], config['btag_regions'][args.channel[0]][0])
        f.write('  TruthDistributionPath: {}\n'.format(truth_path))
        f.write('  TruthDistributionFile: {}\n'.format(file_name))
        f.write('  NumberOfTruthBins: {}\n'.format( n_bin ))
        f.write('  TitleX: {}\n'.format(labels['observables'][args.variable].replace('$', '')))
        if args.measurement == 'param':
            import ROOT
            rf = ROOT.TFile('{}/{}.root'.format(truth_path, file_name), 'READ')
            truth = rootToNumpy( rf.Get('truth') )
            writeExpression(binning, truth, f, factor)
        writeConfigBlock(config, 'unfolding', f)

        # add norm factors to measure the unfolded mean if likelihood is re-parametrised
        if args.measurement == 'param':
            iglob = 0
            for ib, b in enumerate(binning):
                if len(b) > 2:
                    f.write('NormFactor: \"param_{}\"\n'.format(ib))
                    f.write( '  Nominal: {}\n'.format(factor * getBinnedMean(truth[iglob:iglob+len(b)-1], b)) )
                    f.write('  Min: -2\n')
                    f.write('  Max: 2\n')
                    f.write('  Samples: none\n')
                    f.write('\n')
                
                iglob += len(b)-1

        
        f.write('TruthSample: PP8\n')
        writeConfigBlock(config, 'truthsample', f)

        region_list = []
        for channel in args.channel:
            for btag in config['btag_regions'][channel]:
                f.write('Region: {}_{}\n'.format(channel, btag))
                region_list.append('{}_{}'.format(channel, btag))
                f.write('  ResponseMatrixPathSuff: {}/{}\n'.format( channel, btag ))
                f.write('  HistoPathSuff: {}/{}\n'.format( channel, btag ))
                if rebin_config['is_control']==True and channel in rebin_config['control_channels'].keys():
                    f.write('  NumberOfRecoBins: {}\n'.format( rebin_config['nbin_control'] ))
                else:
                    f.write('  NumberOfRecoBins: {}\n'.format( n_bin ))
                f.write('  Label: {} {}\n'.format(labels['regions'][channel], labels['regions'][btag]))
                f.write('  ShortLabel: {}\n'.format(channel))
                writeConfigBlock(config, 'region', f)
        
        f.write('Sample: Data\n')
        f.write('  HistoFile: {}\n'.format( file_name ))
        f.write('  HistoName: {}\n'.format(args.data[i]))
        f.write('  Regions: {}\n'.format(', '.join(region_list)))
        writeConfigBlock(config, 'datasample', f)
        
        f.write('UnfoldingSample: ttbar\n')
        f.write('  ResponseMatrixFile: {}\n'.format( file_name ))
        f.write('  ResponseMatrixName: {}\n'.format(config['response_matrix_name']))
        f.write('  Regions: {}\n'.format(', '.join(region_list)))
        writeConfigBlock(config, 'unfoldingsample', f)
        
        for bkg in config['background']:
            f.write('Sample: {}\n'.format(bkg))
            f.write('  Title: {}\n'.format(bkg.capitalize()))
            f.write('  HistoFile: {}\n'.format( file_name ))
            f.write('  HistoName: {}\n'.format(bkg))
            f.write('  FillColor: {}\n'.format(config['bkgcolour'][bkg]))
            f.write('  Regions: {}\n'.format(', '.join(region_list)))
            writeConfigBlock(config, 'bkgsample', f)

        # write down ghost samples, to be used as reference for some systematics
        for proc in config['reference'].keys():
            for ref in list( set( config['reference'][proc].values() ) ):
                if proc == 'ttbar':
                    f.write('UnfoldingSample: {}\n'.format(ref))
                    f.write('  Type: GHOST\n')
                    f.write('  ResponseMatrixFile: {}\n'.format( file_name ))
                    f.write('  ResponseMatrixName: {}/{}\n'.format(ref, config['response_matrix_name']))
                    f.write('  Regions: {}\n'.format(', '.join(args.channel)))
                    f.write('\n')
                else:
                    f.write('Sample: {}\n'.format(ref))
                    f.write('  Type: GHOST\n')
                    f.write('  HistoFile: {}\n'.format( file_name ))
                    f.write('  HistoName: {}/{}\n'.format(ref, bkg))
                    f.write('  Regions: {}\n'.format(', '.join(args.channel)))
                    f.write('\n')

        for proc in config['background']+['ttbar']:
            writeSystBlock(config, proc, config['twosided_syst'], f, False)
            writeSystBlock(config, proc, config['onesided_syst'], f, True)

    return config, rebin_config, labels, cfgname, file_name, binning, diff_bins, lumi_campaign, factor


            
def writeConfigBlock(config, block, outfile):
    opt = config[block]
    for key, value in opt.items():
        outfile.write('  {}: {}\n'.format(key, value))
    outfile.write('\n')



def writeSystBlock(config, proc, syst_dict, outfile, onesided=False):
    # get all systs for process
    syst_list = []
    if 'all' in syst_dict.keys():
        syst_list += syst_dict['all']
    if proc in syst_dict.keys():
        syst_list += syst_dict[proc]


    # check if some systs use reference sample different from nominal
    try:
        syst_list_ref = config['reference'][proc].keys()
    except:
        syst_list_ref = []


    # write config for systematics
    # one block for each systematic with a reference sample
    for s in syst_list_ref:
        if s in syst_list:
            syst_list.remove(s)  # remove syst from list if it uses a reference sample

            syst_name = s
            if onesided:
                if proc is 'ttbar':
                    syst_path = '{}/{}'.format(s, config['response_matrix_name'])
                else:
                    syst_path = '{}/{}'.format(s, proc)
            else:
                syst_name += syst[0].replace('_up', '').replace('_UP', '') + ' ; '
                if proc is 'ttbar':
                    syst_path_up += '{}/{}'.format(s[0], config['response_matrix_name'])
                    syst_path_down += '{}/{}'.format(s[1], config['response_matrix_name'])
                else:
                    syst_path_up += '{}/{}'.format(s[0], proc)
                    syst_path_down += '{}/{}'.format(s[1], proc)

            if proc is 'ttbar':
                outfile.write('UnfoldingSystematic: {}\n'.format(syst_name))
                if onesided:
                    outfile.write('  Symmetrisation: ONESIDED\n')
                    outfile.write('  ResponseMatrixNameUp: {}\n'.format(syst_path))
                else:
                    outfile.write('  ResponseMatrixNameDown: {}\n'.format(syst_path_down))
                    outfile.write('  ResponseMatrixNameUp: {}\n'.format(syst_path_up))                
                    outfile.write('  Symmetrisation: TWOSIDED\n')
            else:
                outfile.write('Systematic: {}\n'.format(syst_name))
                if onesided:
                    outfile.write('  Symmetrisation: ONESIDED\n')
                    outfile.write('  HistoNameUp: {}\n'.format(syst_path))
                else:
                    outfile.write('  HistoNameDown: {}\n'.format(syst_path_down))
                    outfile.write('  HistoNameUp: {}\n'.format(syst_path_up))
                    outfile.write('  Symmetrisation: TWOSIDED\n')

            outfile.write('  NuisanceParameter: {}\n'.format(syst_name))
            outfile.write('  Title: {}\n'.format(syst_name))
            outfile.write('  Type: HISTO\n')
            outfile.write('  Samples: {}\n'.format(proc))
            outfile.write('  ReferenceSample: {}\n'.format(config['reference'][proc][s]))
            outfile.write('\n')

    # one common block for all systematics without reference sample
    if len(syst_list) > 0:
        syst_name = ''
        if onesided:
            syst_path = ''
            for syst in syst_list:
                syst_name += syst + ' ; '
                if proc is 'ttbar':
                    syst_path += '{}/{} ; '.format(syst, config['response_matrix_name'])
                else:
                    syst_path += '{}/{} ; '.format(syst, proc)
        else:
            syst_path_up = ''
            syst_path_down = ''
            for syst in syst_list:
                syst_name += syst[0].replace('_up', '').replace('_UP', '').replace('__1up', '') + ' ; '
                if proc is 'ttbar':
                    syst_path_up += '{}/{} ; '.format(syst[0], config['response_matrix_name'])
                    syst_path_down += '{}/{} ; '.format(syst[1], config['response_matrix_name'])
                else:
                    syst_path_up += '{}/{} ; '.format(syst[0], proc)
                    syst_path_down += '{}/{} ; '.format(syst[1], proc)

        if proc is 'ttbar':
            outfile.write('UnfoldingSystematic: {}\n'.format(syst_name))
            if onesided:
                outfile.write('  Symmetrisation: ONESIDED\n')
                outfile.write('  ResponseMatrixNameUp: {}\n'.format(syst_path))
            else:
                outfile.write('  ResponseMatrixNameDown: {}\n'.format(syst_path_down))
                outfile.write('  ResponseMatrixNameUp: {}\n'.format(syst_path_up))                
                outfile.write('  Symmetrisation: TWOSIDED\n')
        else:
            outfile.write('Systematic: {}\n'.format(syst_name))
            if onesided:
                outfile.write('  Symmetrisation: ONESIDED\n')
                outfile.write('  HistoNameUp: {}\n'.format(syst_path))
            else:
                outfile.write('  HistoNameDown: {}\n'.format(syst_path_down))
                outfile.write('  HistoNameUp: {}\n'.format(syst_path_up))
                outfile.write('  Symmetrisation: TWOSIDED\n')

        outfile.write('  NuisanceParameter: {}\n'.format(syst_name))
        outfile.write('  Title: {}\n'.format(syst_name))
        outfile.write('  Type: HISTO\n')
        outfile.write('  Samples: {}\n'.format(proc))

        outfile.write('\n')





def writeExpression(binning, truth, outfile, factor):
    # write expression to re-parametrise the likelihood in order to measure the mean of the unfolded distribution
    expr_full = '  Expressions: '
    iglob = 0
    for ib, b in enumerate(binning):
        if len(b) > 2:
            bin_center = [ (b[i]+b[i+1])/2 for i in range(len(b)-1) ]
            expr = ''
            dep = 'param_{}[-2,2]'.format(ib)
            for ibc, bc in enumerate(bin_center[1:]):
                if ibc > 0:
                    expr += '+'
                if (iglob+ibc+2) < 10:
                    expr += '({}-param_{}/{})*Bin_00{}_mu*{}'.format(bc, ib, factor, iglob+ibc+2, truth[iglob+ibc+1])
                    dep += ',Bin_00{}_mu[0,2]'.format(iglob+ibc+2)
                else:
                    expr += '({}-param_{}/{})*Bin_0{}_mu*{}'.format(bc, ib, factor, iglob+ibc+2, truth[iglob+ibc+1])
                    dep += ',Bin_0{}_mu[0,2]'.format(iglob+ibc+2)
            expr = '({})/((param_{}/{}-{})*{})'.format(expr, ib, factor, bin_center[0], truth[iglob])
        
            if ib > 0:
                expr_full += ','
            if (iglob+1) < 10:
                expr_full += '\"Bin_00{}_mu\"=\"{}\":\"{}\"'.format(iglob+1, expr, dep)
            else:
                expr_full += '\"Bin_0{}_mu\"=\"{}\":\"{}\"'.format(iglob+1, expr, dep)

        iglob += len(b)-1

    outfile.write('{}\n'.format(expr_full))




def main(args):
    for i in range(len(args.data)):
        # write TRExFitter config file
        config, rebin_config, labels, trex_config, file_name, binning, diff_bins, lumi_campaign, factor = writeConfig(args, i)
        path = '{}/{}_{}_{}/{}'.format(args.outpath, args.campaign, '_'.join(args.channel), file_name, args.data[i])


        # run TRExFitter
        os.system( 'trex-fitter uhwf {}'.format(trex_config) )


        # run unfolding for null hypothesis if significance calculation requested
        if args.variable in config['significance'].keys() and args.measurement == 'param':
            with open(trex_config, 'r') as f: # open original config file
                sig_cfg = f.read()
                sig_cfg = sig_cfg.replace('param_0[-2,2]', f'param_0[{config["significance"][args.variable]}]') # fix param value for null hypothesis
                path_sig = path.replace(f'/{args.data[i]}', '') + '_sig'  # output directory for null hypothesis unfolding
                sig_cfg = sig_cfg.replace(path.replace(f'/{args.data[i]}', ''), path_sig)
            trex_config_sig = trex_config.replace('.config', '_sig.config')
            with open(trex_config_sig, 'w') as f: # save modified config
                f.write(sig_cfg)
                
            # run null hypothesis unfolding
            print(f'\n\n\nRunning unfolding for {factor} * <{args.variable}> = {config["significance"][args.variable]}\n\n')
            os.system( 'trex-fitter uhwf {}'.format(trex_config_sig) )

            # get measured minimum NLLs
            with open(f'{path}/Fits/{args.data[i]}.txt', 'r') as f:
                lf = f.readlines()
                NLL_unconstrained = float(lf[len(lf)-2].replace('\n', ''))
            with open(f'{path_sig}/{args.data[i]}/Fits/{args.data[i]}.txt', 'r') as f:
                lf = f.readlines()
                NLL_nullhyp = float(lf[len(lf)-2].replace('\n', ''))

            sig = np.sqrt( 2*(NLL_nullhyp-NLL_unconstrained) )  # using q0 formula from asymptotic approximation

            

        # Save results to json file
        result_dict = {}
        result_dict['observable'] = args.variable
        result_dict['binning'] = rebin_config['binning']
        if rebin_config['differential_measurement']:
            result_dict['differential_observable'] = rebin_config['diff_variable']
            result_dict['differential_observable_unit'] = rebin_config['diff_variable_unit']
            result_dict['differential_binning'] = rebin_config['diff_bins']
        else:
            result_dict['differential_binning'] = 'inclusive'

        if args.measurement == 'param':
            params = getNPsFromFitResults('{}/Fits/{}.root'.format(path, args.data[i]), get_param=True)

            mean_list = []
            uncert_up_list = []
            uncert_down_list = []

            for name, p in params.items():
                print( '{} * <{}> = {:.4f} + {:.4f} - {:.4f}'.format(factor, args.variable, p[0], p[1], abs(p[2])) )
                mean_list.append( p[0] )
                uncert_up_list.append( p[1] )
                uncert_down_list.append( abs(p[2]) )

            result_dict['mean'] = mean_list
            result_dict['uncert_up'] = uncert_up_list
            result_dict['uncert_down'] = uncert_down_list

            if args.variable in config['significance'].keys():
                result_dict['significance'] = sig
                print(f'significance = {sig:.4f} (H0: {factor} * <{args.variable}> = {config["significance"][args.variable]:.4f})')

            writejson('{}/{}_unfolded_parameter.json'.format(path, args.variable), result_dict)
            print('\nThe unfolded parameters are saved to: {}/{}_unfolded_parameter.json\n'.format(path, args.variable))

        elif args.measurement == 'distrib':
            import yaml
            ## Get unfolded bins content and uncertainties
            with open('{}/UnfoldingData.yaml'.format(path)) as f:
                bins = yaml.full_load(f)
            unfolded_data = np.array( [ bins[i]['mean'] for i in range(len(bins)) ] )
            unfolded_err_up = np.array( [  bins[i]['uncertaintyUp'] for i in range(len(bins)) ] )
            unfolded_err_down = np.array( [  abs(bins[i]['uncertaintyDown']) for i in range(len(bins)) ] )
            
            if len(unfolded_data) == len(binning)*len(binning[0]): ## only if all diff bin have same binning
                unfolded_data_2d = unfolded_data.reshape( (len(rebin_config['binning']), len(rebin_config['binning'][0])-1) )
                unfolded_err_up_2d = unfolded_err_up.reshape( (len(rebin_config['binning']), len(rebin_config['binning'][0])-1) )
                unfolded_err_down_2d = unfolded_err_down.reshape( (len(rebin_config['binning']), len(rebin_config['binning'][0])-1) )
                
                result_dict['xsec'] = np.sum( unfolded_data_2d / config['lumi'][lumi_campaign][1], axis=1)
                
                result_dict['mean'] = ( unfolded_data_2d.T / np.sum(unfolded_data_2d, axis=1) ).T
                result_dict['uncert_up'] = ( unfolded_err_up_2d.T / np.sum(unfolded_data_2d, axis=1) ).T
                result_dict['uncert_down'] = ( unfolded_err_down_2d.T / np.sum(unfolded_data_2d, axis=1) ).T

                ## Get correlation matrix
                with open('{}/CorrelationMatrix.yaml'.format(path)) as f:
                    corrmat = np.array( yaml.full_load(f)[1]['correlation_rows'] )
                    ## Keep only correlations for truth bins
                    nbins = len(bins)
                    corrmat = corrmat[0:nbins, 0:nbins]
                    result_dict['correlation_matrix'] = corrmat
                
                writejson('{}/{}_unfolded_distribution.json'.format(path, args.variable), result_dict)
                print('\nThe unfolded parameters are saved to: {}/{}_unfolded_distribution.json\n'.format(path, args.variable))


        # plot unfolded distribution and response matrix if requested
        if args.plot:
            import ROOT
            rfile = ROOT.TFile('{}/UnfoldingHistograms/FoldedHistograms.root'.format(path), 'read')

            # plot response matrix
            for chan in args.channel:
                for btag in config['btag_regions'][chan]:
                    mat = '{}_{}_ttbar_response'.format(chan, btag)
                    h_resmat = rfile.Get(mat)
                    resmat = rootToNumpy(h_resmat, is_2d=True)
    
                    # create response matrix normalised to 1 per truth bin
                    resmat_norm = ( resmat.T / np.sum(resmat, axis=1) ).T

                    # plot response matrices (not normalised and normalised)
                    plotResmat('response_matrix', resmat, args.config, path, chan, btag, args.variable, lumi_campaign, diff_bins, binning, docmap=True)
                    plotResmat('response_matrix_normalised', resmat_norm, args.config, path, chan, btag, args.variable, lumi_campaign, diff_bins, binning)

            # plot spin parameter vs differential variable
            if args.measurement == 'param':
                if len(mean_list) > 1:
                    truth_param = hjson.load( open(f'{args.inpath}/{args.campaign}/{args.channel[0]}/{config["btag_regions"][args.channel[0]][0]}/truth_{args.variable}_{rebin_config["diff_variable"]}_diff_binned.json', 'r') )['truth']
                    plotSpinParameter(path, mean_list, uncert_up_list, uncert_down_list, truth_param, args.variable, lumi_campaign, args.config, diff_bins)

            # plot unfolded distribution
            elif args.measurement == 'distrib':
                truth_hist = rfile.Get('truth_distribution')
                truth = [ truth_hist.GetBinContent(i+1) for i in range(truth_hist.GetNbinsX()) ]
                plotUnfoldedDistribution(path, unfolded_data, unfolded_err_up, unfolded_err_down, truth, args.variable, lumi_campaign, args.config, diff_bins, binning)

            # plot NP pulls and constraints
            NP_dict = getNPsFromFitResults('{}/Fits/{}.root'.format(path, args.data[i]))
            plotBrasilian(path, NP_dict, labels)



    


if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter

    # get path for root files
    parser = ArgumentParser(description='Prepare config file TRExFitter.', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inpath', type=str, help='Path to input', required=True)
    parser.add_argument('-o', '--outpath', type=str, default='results/base/', help='Path to input')
    parser.add_argument('--config', type=str, help='Config file', default='config/config.json')
    parser.add_argument('-v', '--variable', type=str, help='Spin observable to run')
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels to include')
    parser.add_argument('--campaign', type=str, help='MC campaign to run')
    parser.add_argument('-d', '--data', nargs='*', default=['mcdata'], help='Dataset to run on')
    parser.add_argument('-m', '--measurement', choices=['param', 'distrib'], default='param', help='Measure unfolded distribution via original likelihood or spin parameter via re-parametrisation of likelihood.')
    parser.add_argument('-p', '--plot', action='store_true', help='plot unfolded distribution and response matrix')
    args = parser.parse_args()

    main(args)
