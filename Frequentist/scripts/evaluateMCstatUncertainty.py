import numpy as np
import os, hjson, json
from utils import writejson, submit_job, makeToys, plotEnsembleTestHist
from datetime import datetime




def make_config(config, outdir, migmat_name):
    """
    Create config with toy migmat 
    """

    config['job']['DebugLevel'] = 0
    config['job']['StatOnly'] = 'TRUE'
    config['onesided_syst'] = {}
    config['twosided_syst'] = {}

    config['migration_matrix_name'] = migmat_name

    try:
        del config['fit']['UseMinos']
    except:
        None
    try:
        del config['fit']['doLHscan']
    except:
        None
    
    os.makedirs('{}/'.format(outdir), exist_ok=True)
    writejson('{}/config.json'.format(outdir), config)




def plotEnsembleTest(args):
    outpath = 'results/MCstat/{}_{}/'.format(args.campaign, '_'.join(args.channel))
    os.makedirs(outpath, exist_ok=True)

    config = hjson.load(open(args.config, 'r'))
    labels = hjson.load(open(config['labels'], 'r'))

    ## use correct path if differential measurement
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    Ndiff = len(rebin_config['binning'])
    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])
        diff_bins = rebin_config['diff_bins']
    else:
        diff_bins = ['']
    
    ensemble_pull = []
    ensemble_constraint = []
    
    # plot result of ensemble test as histogram for each tested binning
    param_mean = []
    for itoy in range(args.ntoys):
        try:
            param = json.load( open('{}/results/toymigmat_{}/{}_{}_{}/mcdata/spin_parameter.json'.format(args.inpath, itoy, args.campaign, '_'.join(args.channel), file_name) ) )
            param_mean.append( param['mean'] )

        except Exception as e:
            print(e)

    # ensemble test for parameter
    param_mean = np.array(param_mean)
    label_xaxis = labels['parameters'][args.variable]
    for idiffbin in range(Ndiff):
        outpath_tmp = '{}/ensemble_test_{}_param.pdf'.format(outpath, file_name.replace('diff', '_'.join(diff_bins[idiffbin])))
        plotEnsembleTestHist(param_mean[:,idiffbin], label_xaxis, outpath_tmp)




def runToys(args):
    timestamp = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
    outdir = '{}/{}/'.format(args.output, timestamp)
    config = hjson.load(open(args.config, 'r'))

    print('\nRunning ensemble test to determine systematic uncertainty due to MC stat. uncertainty')
    print('Will produce corresponding inputs and outputs under directory {}'.format(outdir))

    # generate toy migration matrices to perform pseudoexperiments
    for channel in args.channel:
        inpath_new = '{}/unfoldingInput/{}/{}/'.format(outdir, args.campaign, channel)
        os.makedirs(inpath_new, exist_ok=True)
        os.system('cp {}/{}/{}/{}.root {}'.format(args.inpath, args.campaign, channel, args.variable, inpath_new))
        makeToys('{}/{}.root'.format(inpath_new, args.variable), 'migmat', args.ntoys, is_2d=True, pdf=args.method)


    for itoy in range(args.ntoys):
        ## Create config file
        outdir_cfg = '{}/config/toymigmat_{}'.format(outdir, itoy)
        make_config(config, outdir_cfg, 'toymigmat_{}'.format(itoy))

        ## Launch job on batch system: unfolding
        outdir_output_tmp = '{}/results/toymigmat_{}/'.format(outdir, itoy)
        os.makedirs(outdir_output_tmp, exist_ok=True)
        dep_unfolding = []
        command = 'python scripts/runTRExFitter.py -i {}/unfoldingInput/ -o {} --config {}/config.json -v {} -c {} --campaign {} \n'.format(outdir, outdir_output_tmp, outdir_cfg, args.variable, ' '.join(args.channel), args.campaign)
        command += 'rm -rf {0}/{1}_{2}_{3}/Histograms/ {0}/{1}_{2}_{3}/Systematics/ {0}/{1}_{2}_{3}/RooStats/'.format(outdir_output_tmp, args.campaign, '_'.join(args.channel), args.variable)
        dep_unfolding.append( submit_job('unfold', 'toymigmat_{}'.format(itoy), outdir, command) )

    ## Launch job to rank impacts and draw summary plot
    command = 'python scripts/evaluateMCstatUncertainty.py -i {} -v {} -c {} --campaign {} --ntoys {} --plot-ensemble\n'.format( outdir, args.variable, ' '.join(args.channel), args.campaign, args.ntoys )
    submit_job('plotensemble', 'MCstat', outdir, command, dependency=':'.join(dep_unfolding))

    print( '\nThe summary plot will be saved to results/MCstat/{}_{}/'.format(args.campaign, '_'.join(args.channel)) )
        




if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description='Determine systematic uncertainty due to limited MC statistics, from ensemble test.', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--inpath', help='Path to input histograms', required=True)
    parser.add_argument('-o','--output', default='MCstat/', help='Path to output')
    parser.add_argument('-v', '--variable', type=str, help='Spin observable to run', required=True)
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels to run', required=True)
    parser.add_argument('--campaign', type=str, help='MC campaign to run', required=True)
    parser.add_argument ('-d','--data', type=str, help='dataset to run on', default='mcdata')
    parser.add_argument('-g', '--config', type=str, default='config/config.json', help='Path to config file')
    parser.add_argument('--ntoys', type=int, help='number of pseudoexperiments', default=200)
    parser.add_argument('-m', '--method', choices=['Gaus', 'MCstat'], help='Use gaussian (Gaus) or Poisson (MCstat) to draw toy migration matrices.', default='MCstat')

    parser.add_argument('--plot-ensemble', action='store_true', help='Plot ensemble test and get systematic from MC stat. uncertainty.')

    args = parser.parse_args()
    
    if not args.plot_ensemble:
        runToys(args)
    else:
        plotEnsembleTest(args)
