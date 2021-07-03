from utils import submit_job, writejson, plotEnsembleTestPull, plotEnsembleTestHist, makeToys
import os
import hjson, json
import numpy as np




def make_config(config, cfgdir):
    """
    Create config with new generated binning 
    """
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
    
    writejson('{}'.format( cfgdir.replace('.json', '_stat_only.json') ), config)

        


def plotEnsembleTest(args):
    config = hjson.load(open(args.config, 'r'))
    labels = hjson.load(open(config['labels'], 'r'))

    ## use correct path if differential measurement
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])
        diff_bins = rebin_config['diff_bins']
    else:
        diff_bins = ['']

    Ndiff = len( [b for b in rebin_config['binning'] if len(b)>2] ) # number of diff bins with more than one bin

    outpath = '{}/{}_{}_{}/pull_test/'.format(args.outpath, args.campaign, '_'.join(args.channel), file_name)
    os.makedirs(outpath, exist_ok=True)
    
    # plot result of ensemble test as histogram
    param_mean = []
    for itoy in range(args.ntoys):
        try:
            param = json.load( open('{}/{}_{}_{}/toymcdata_{}/{}_unfolded_parameter.json'.format(args.outpath, args.campaign, '_'.join(args.channel), file_name, itoy, args.variable) ) )
            param_mean.append( param['mean'] )                
        except Exception as e:
            print(e)

    # get nominal mean and std to get ensemble pull and constraint
    nom_param = json.load( open('{}/stat_only/{}_{}_{}/mcdata/{}_unfolded_parameter.json'.format(args.outpath, args.campaign, '_'.join(args.channel), file_name, args.variable) ) )

    # ensemble test for parameter
    param_mean = np.array(param_mean)
    label_xaxis_hist = labels['parameters'][args.variable]
    label_xaxis_pull = '$\\left({0}{1} - {0}{2}\\right) / \\Delta {0}$'.format(labels['parameters'][args.variable].replace('$', ''), '_{ensemble}', '_{asimov}')
    for idiffbin in range(Ndiff):
        pull_up = (np.mean(param_mean[:,idiffbin]) - nom_param['mean'][idiffbin]) / nom_param['uncert_up'][idiffbin]   # just pull up because using only stat. uncertainties (i.e. err_up = err_down)
        constraint = np.std(param_mean[:,idiffbin]) / np.abs(nom_param['uncert_up'][idiffbin])
        outpath_tmp = '{}/ensemble_test_{}_param.pdf'.format(outpath, file_name.replace('diff', '_'.join(diff_bins[idiffbin])))
        plotEnsembleTestHist(param_mean[:,idiffbin], label_xaxis_hist, pull_up, constraint, label_xaxis_pull, outpath_tmp)





def main(args):
    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))

    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])

    make_config(config, args.config)

    os.system('rm -rf pull_test')

    # generate toy datasets to perform pseudoexperiments
    for channel in args.channel:
        for btag in config['btag_regions'][channel]:
            makeToys('{}/{}/{}/{}/{}.root'.format(args.inpath, args.campaign, channel, btag, file_name), 'mcdata', args.ntoys, is_2d=False)

    # run nominal stat-onlu unfolding
    command = 'python scripts/runTRExFitter.py -i {} -o {}/stat_only/ --config {} -v {} -c {} --campaign {} -d mcdata\n'.format(args.inpath, args.outpath, args.config.replace('.json', '_stat_only.json'),
                                                                                                                     args.variable, ' '.join(args.channel), args.campaign)
    submit_job('unfold_stat_only_nom', args.variable, 'pull_test', command)

    # run unfolding for all toy datasets
    for itoy in range(args.ntoys):
        command = 'python scripts/runTRExFitter.py -i {} -o {} --config {} -v {} -c {} --campaign {} -d toymcdata_{}\n'.format(args.inpath, args.outpath, args.config.replace('.json', '_stat_only.json'),
                                                                                                                               args.variable, ' '.join(args.channel), args.campaign, itoy)
        command += 'rm -rf {0}/{1}_{2}_{3}/toymcdata{4}/Histograms/ {0}/{1}_{2}_{3}/toymcdata{4}/Systematics/ {0}/{1}_{2}_{3}/toymcdata{4}/RooStats/'.format(args.outpath, args.campaign, 
                                                                                                                                                             '_'.join(args.channel), file_name, itoy)
        submit_job('unfold_toymcdata_{}'.format(itoy), args.variable, 'pull_test', command)




if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter

    parser = ArgumentParser('Find best binning configurations based on linearity test', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--inpath', type=str, help='path to rebinned input')
    parser.add_argument('-o','--outpath', type=str, help='path to output')
    parser.add_argument('-v','--variable', type=str, help='name of variable')
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels', required=True)
    parser.add_argument('--campaign', type=str, help='MC campaign to run', required=True)
    parser.add_argument('-g', '--config', type=str, default='config/config.json', help='Path to config file')
    parser.add_argument('--ntoys', type=int, help='number of pseudoexperiments', default=1000)
    parser.add_argument('--plot-ensemble', action='store_true', help='plot results of ensemble tests')
    args = parser.parse_args()

    if args.plot_ensemble == False:
        main(args)
    else:
        plotEnsembleTest(args)
