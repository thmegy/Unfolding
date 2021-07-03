import os, sys
import json, hjson
from utils import writejson, submit_job
import numpy as np
from datetime import datetime



def getBootstraps(config):
    rebin_config = hjson.load(open(config['rebin']))
    n_bootstrap = rebin_config['n_bootstrap']
    bootstrap_name = rebin_config['bootstrap_name']
    return n_bootstrap, bootstrap_name    



def unfoldBootstraps(args):
    ## Make output directory
    if args.output == 'observable_correlations/':
        timestamp = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
        outdir = '{}/{}/'.format(args.output, timestamp)
    else:
        outdir = args.output

    print('\nCreating the output directory: {}'.format(outdir))

    if len(args.variable) < 2:
        print('You need at leat two observables to measure correlations! Exiting..')
        sys.exit()

    # get name and number of bootstraps
    config = hjson.load(open(args.config))
    n_bootstrap, bootstrap_name = getBootstraps(config)
    config['job']['DebugLevel'] = 0
    config['job']['StatOnly'] = 'TRUE'
    config['onesided_syst'] = {}
    config['twosided_syst'] = {}
    try:
        del config['fit']['UseMinos']
    except:
        None
    try:
        del config['fit']['doLHscan']
    except:
        None
    os.makedirs('{}/config/'.format(outdir), exist_ok=True)
    writejson('{}/config/config.json'.format(outdir), config)

    ## Loop over tested numbers of bins
    for var in args.variable:
        print('\nLaunching jobs for the {} observable...'.format(var))
        # loop over bootstraps
        for ib in range(n_bootstrap):
            ## Launch job on batch system: unfolding
            outdir_output_tmp = '{}/results/'.format(outdir)
            os.makedirs(outdir_output_tmp, exist_ok=True)
            command = 'python scripts/runTRExFitter.py -i {} -o {} --config {}/config/config.json -v {} -c {} --campaign {} -d {}_{}\n'.format(args.inpath, outdir_output_tmp, outdir, var, ' '.join(args.channel), args.campaign, bootstrap_name, ib)
            result_path = '{0}/{1}_{2}_{3}/{4}_{5}/'.format(outdir_output_tmp, args.campaign, '_'.join(args.channel), var, bootstrap_name, ib)
            command += 'rm -rf {0}/Histograms/ {0}/Systematics/ {0}/RooStats/ {0}/Fits/ {0}/UnfoldingHistograms/ {0}/Systematics/ {0}/*pdf'.format(result_path)
            submit_job('unfold_{}_{}'.format(bootstrap_name, ib), var, outdir, command, mem='4gb')




def getCorrelations(args):
    import pandas as pd
    import matplotlib as mpl
    mpl.use('Agg') # load Agg backed which does not require X11 connection to make plots
    import matplotlib.pyplot as plt

    # get name and number of bootstraps
    config = hjson.load(open(args.config))
    n_bootstrap, bootstrap_name = getBootstraps(config)
    labels = hjson.load(open(config['labels'], 'r'))

    # check if unfolding has been succesful for all bootstraps of all variables
    excluded_bootstraps = set()
    for var in args.variable:
        for ib in range(n_bootstrap):
            result_path = '{}/results/{}_{}_{}/{}_{}/spin_parameter.json'.format(args.inpath, args.campaign, '_'.join(args.channel), var, bootstrap_name, ib)
            if os.path.isfile(result_path):
                continue
            else:
                print('The output for bootstrap {} of {} does not exist...'.format(ib, var))
                excluded_bootstraps.add(ib)

    print('\nThe following bootstraps are excluded: {}\n'.format(excluded_bootstraps))

    # construct panda dataframe
    result_array = []
    for ib in range(n_bootstrap):
        if ib in excluded_bootstraps:
            continue
        else:
            result_array_bootstrap = []
            for var in args.variable:
                result_path = '{}/results/{}_{}_{}/{}_{}/spin_parameter.json'.format(args.inpath, args.campaign, ' '.join(args.channel), var, bootstrap_name, ib)
                param = json.load(open(result_path, 'r'))
                result_array_bootstrap.append( param['mean'] )
            result_array.append( result_array_bootstrap )

    names = [ labels['parameters'][var] for var in args.variable ]
    result = pd.DataFrame(data=result_array, columns=names)
    plt.rcParams['axes.labelsize'] = 20

    # plot scatter matrix
    axes = pd.plotting.scatter_matrix(result, figsize=(20, 20), diagonal='hist', rasterized=True)
    for ax in axes.flatten():
        ax.xaxis.label.set_rotation(90)
        ax.yaxis.label.set_rotation(0)
        ax.yaxis.label.set_ha('right')
    outpath = 'results/observable_correlations/'
    os.makedirs(outpath, exist_ok=True)
    outpath += 'scatter_matrix_{}_{}.pdf'.format(args.campaign, '_'.join(args.channel), '_'.join(args.channel))
    plt.savefig(outpath, bbox_inches='tight')
    print('Saving {}'.format(outpath))

    # plot correlation matrix
    corr = result.corr()
    fig, ax = plt.subplots( figsize=(10, 8) )
    opts = {'cmap': 'RdBu_r', 'vmin': -1, 'vmax': +1}
    heatmap = ax.pcolor(corr, **opts)
    ax.set_yticks(np.arange(0.5, len(names), 1))
    ax.set_xticks(np.arange(0.5, len(names), 1))
    ax.set_xticklabels(labels=names, fontsize=20, rotation=90, ha='center')
    ax.set_yticklabels(labels=names, fontsize=20)
    cbar = plt.colorbar(heatmap)
    outpath = 'results/observable_correlations/correlation_matrix_{}_{}.pdf'.format(args.campaign, '_'.join(args.channel), '_'.join(args.channel))
    fig.savefig(outpath, bbox_inches='tight')
    print('Saving {}'.format(outpath))
    
    with open(outpath.replace('.pdf', '.tex'), 'w') as f:
        f.write(corr.to_latex())
        print('Saving {}'.format(outpath.replace('.pdf', '.tex')))
    




if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter
	
    parser = ArgumentParser('Measure correlations between spin observables by runnig unfolding with bootstraps', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--inpath', help='Path to input histograms', required=True)
    parser.add_argument('-o','--output', default='observable_correlations/', help='Path to output')
    parser.add_argument('-v', '--variable', nargs='*', default=['CosThetaKplus', 'CosThetaNplus', 'CosThetaRplus', 'CosThetaKminus', 'CosThetaNminus', 'CosThetaRminus', 'CorrKK', 'CorrNN', 'CorrRR', 'CorrNKplus', 'CorrNKminus', 'CorrNRplus', 'CorrNRminus', 'CorrRKplus', 'CorrRKminus'], help='Observables to prepare the inputs for')
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels to run', required=True)
    parser.add_argument('--campaign', type=str, help='MC campaign to run', required=True)
    parser.add_argument('-g', '--config', type=str, default='config/config.json', help='Path to config file')

    parser.add_argument('--plot-correlation', action='store_true', help='Plot correlations between spin parameters')    
    
    args = parser.parse_args()

    if not args.plot_correlation:
        unfoldBootstraps(args)
    else:
        getCorrelations(args)
