import matplotlib as mpl
mpl.use('Agg') # load Agg backed which does not require X11 connection to make plots
import matplotlib.pyplot as plt
import os
import math
import json, hjson, yaml
from utils import submit_job, writejson, plotEnsembleTestPull, plotEnsembleTestHist, plotScatter, normaliseFeature, makeToys
import numpy as np
import re


def find_binning_with_the_smallest_parameter(distance_dict):
    previous_binning_parameter = 1000
    best_binning = ''
    
    for key, value in distance_dict.items():
        figure_of_merit = value
        current_binning_parameter = float(figure_of_merit)
        if current_binning_parameter < previous_binning_parameter:
            previous_binning_parameter = current_binning_parameter
            best_binning = key

    return best_binning, previous_binning_parameter




def find_other_binnings_with_smallest_parameter(previous_best_binning, distance_dict):
    previous_binning_parameter = 1000
    next_best_binning = ''

    for key, value in distance_dict.items():
        figure_of_merit = value
        current_binning_parameter = float(figure_of_merit)
        if current_binning_parameter > previous_best_binning and current_binning_parameter < previous_binning_parameter:
            previous_binning_parameter = current_binning_parameter
            next_best_binning = key

    return next_best_binning, previous_binning_parameter

        


def plotEnsembleTest(args):
    outpath = 'results/binning_optimisation/{}_{}/'.format(args.campaign, '_'.join(args.channel))
    os.makedirs(outpath, exist_ok=True)

    options = hjson.load(open(args.config, 'r'))
    labels = hjson.load(open(options['labels'], 'r'))
    best_binnings = hjson.load(open('{}/config/best_binnings.json'.format(args.inpath), 'r'))

    ## use correct path if differential measurement
    rebin_config = hjson.load(open(options['rebin'], 'r'))
    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])
        diff_bins = rebin_config['diff_bins']
    else:
        diff_bins = ['']
    Ndiff = len(diff_bins)
    
    ensemble_pull = []
    ensemble_constraint = []
    
    # plot result of ensemble test as histogram for each tested binning
    for binning in best_binnings.keys():

        param_mean = []
        for itoy in range(args.ntoys):
            try:
                param = json.load( open('{}/results/{}/{}_{}_{}/toymcdata_{}/spin_parameter.json'.format(args.inpath, binning, args.campaign, '_'.join(args.channel), file_name, itoy) ) )
                param_mean.append( param['mean'] )

            except Exception as e:
                print(e)

        # ensemble test for parameter
        param_mean = np.array(param_mean)
        label_xaxis = labels['parameters'][args.variable]
        for idiffbin in range(Ndiff):
            outpath_tmp = '{}/ensemble_test_{}_{}_param.pdf'.format(outpath, binning, file_name.replace('diff', '_'.join(diff_bins[idiffbin])))
            plotEnsembleTestHist(param_mean[:,idiffbin], label_xaxis, outpath_tmp)

        # get nominal mean and std to get ensemble pull and constraint
        nom_param = json.load( open('{}/results/{}/{}_{}_{}/mcdata/spin_parameter.json'.format(args.inpath, binning, args.campaign, '_'.join(args.channel), file_name) ) )
        pull_up = (np.mean(param_mean, axis=0) - np.array(nom_param['mean'])) / np.array(nom_param['uncert_up'])   # just pull up because using only stat. uncertainties (i.e. err_up = err_down)
        constraint = np.std(param_mean, axis=0) / np.abs(np.array(nom_param['uncert_up']))
        ensemble_pull.append( np.mean(pull) )
        ensemble_constraint.append( np.mean(constraint) )

        
    # summary plot: spin parameter pull and constraint for all tested binnings
    label_yaxis = '$\\left({0}{1} - {0}{2}\\right) / \\Delta {0}$'.format(labels['parameters'][args.variable].replace('$', ''), '_{ensemble}', '_{asimov}')
    outpath_tmp = '{}/ensemble_test_{}_pulls.pdf'.format(outpath, file_name)
    plotEnsembleTestPull(ensemble_pull, ensemble_constraint, label_yaxis, best_binnings, outpath_tmp)

    # normalise pulls and constraints and compute ensmble distance
    ensemble_pull = normaliseFeature( np.abs( np.array(ensemble_pull) ) )
    ensemble_constraint = normaliseFeature( np.abs( 1-np.array(ensemble_constraint) ) )
    distance_param = np.sqrt(ensemble_pull**2 + ensemble_constraint**2)

    # summary plot: linearity vs ensemble test --> spin parameter
    distance_param = np.array(distance_param)
    linearity = np.array([best_binnings[b][1] for b in best_binnings.keys()])
    label_scatter = [best_binnings[b][0] for b in best_binnings.keys()]
    x_label = '$D_{linearity}$'
    y_label = '$D{}^{}$'.format('_{ensemble}', '{'+labels['parameters'][args.variable].replace('$', '')+'}')
    outpath_tmp = '{}/linearity_vs_ensemble_test_{}_param.pdf'.format(outpath, file_name)
    plotScatter(linearity, distance_param, label_scatter, x_label, y_label, outpath_tmp, legendonframe=False)

    # save results
    results = {b:{'binning':best_binnings[b][0],
                  'linearity':linearity[i],
                  'ensemble test_param':distance_param[i],
                  'combined_param':np.sqrt(linearity[i]**2+distance_param[i]**2)}
               for i,b in enumerate(best_binnings.keys())}
    outpath_tmp = '{}/final_results_{}.json'.format(outpath, file_name)
    writejson(outpath_tmp, results)
    print('\nWriting results of binning optimisation to {}\n'.format(outpath_tmp))

    
                       


def pickBinnings(args):

    outpath = 'results/binning_optimisation/{}_{}/'.format(args.campaign, '_'.join(args.channel))
    os.makedirs(outpath, exist_ok=True)

    binnings = json.load( open('{}/config/binnings.json'.format(args.inpath), 'r') )
    ## use correct path if differential measurement
    options = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(options['rebin'], 'r'))
    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])


    ########################################################################
    ## draw scatter plots with results of linearity test for all binnings ##
    ########################################################################

    nbin = list( {re.match('\d+bin', i).group().replace('bin', '') for i in binnings.keys()} )  # retrieve how many bins are used in the binnings
    nbin.sort()

    key_array = []

    # lists needed for scatter plots
    slope_array = []
    offset_array = []
    distance_array = []
    uncert_array = []

    # lists needed for figure of merit
    slope_flat = []
    offset_flat = []
    uncert_flat = []

        
    fig, ax = plt.subplots()

    for ibin in nbin:
        if len(ibin)==1:
            keys = [key for key in binnings.keys() if ibin == key[0]]
        elif len(ibin)==2:
            keys = [key for key in binnings.keys() if ibin == key[0:2]]
               
        slope = []
        slope_err = []
        offset = []
        offset_err = []
        
        uncert = []

        distance = []
        bin_size = []
        key_list = []

        for key in keys:
            # get output of linearity test
            infile = '{}/results/{}/{}_{}_{}/linearity_test/linearity_test.json'.format(args.inpath, key, args.campaign, '_'.join(args.channel), file_name)
            if not os.path.isfile( infile ):
                print('The output of the linearity test does not exist for binning {} !'.format(key))
                continue
            else:
                key_list.append(key)
                linearity = json.load( open(infile, 'r') )
                slope.append(linearity['slope'])
                slope_err.append(linearity['slope_err'])
                offset.append(linearity['offset'])
                offset_err.append(linearity['offset_err'])
                
                # compute linearity distance
                distance_diff = []
                for idiffbin in range(len(linearity['slope'])):
                    distance_diff.append( math.sqrt( (1-linearity['slope'][idiffbin])**2 + linearity['offset'][idiffbin]**2 ) )
                distance.append( np.mean( np.array(distance_diff) ) )
                                
                central_bin_index = int( len(binnings[key])/2 ) + 1
                bin_size.append( binnings[key][central_bin_index] )

            # get uncertainty on unfolded parameter
            with open('{}/results/{}/{}_{}_{}/mcdata/{}_unfolded_parameter.json'.format(args.inpath, key, args.campaign, '_'.join(args.channel), file_name, args.variable)) as f:
                param = json.load(f)
                uncert_up = np.array( param['uncert_up'] )
                uncert_down = np.array( param['uncert_down'] )
                uncert.append( (np.mean(uncert_up)+np.mean(uncert_down)) / 2 )
                print(uncert_down)

        key_array += key_list

                
        ## plot figure of merit vs size of bin closest to 0
        fig2, ax2 = plt.subplots()
        ax2.scatter(bin_size, distance)
        ax2.text(0.05, 0.65, '{}-bin binnings'.format(ibin), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, weight='bold')
        ax2.set_xlabel('Central-bin size')
        ax2.set_ylabel('Distance to perfect linearity')
        fig2.savefig('{}/{}_linearity_vs_binsize_{}bins.pdf'.format(outpath, file_name, ibin), bbox_inches='tight')
        print('Saving {}/{}_linearity_vs_binsize_{}bins.pdf\n'.format(outpath, file_name, ibin))

        slope_array.append( np.mean(np.array(slope), axis=1) )
        slope_flat += slope
        offset_array.append( np.mean(np.array(offset), axis=1) )
        offset_flat += offset
        distance_array.append( np.array(distance) )

        uncert_array.append(np.array(uncert))
        uncert_flat += uncert

    # plot linearity offset vs linearity slope
    label_scatter = [ibin+' bins' for ibin in nbin]
    x_label = 'linearity slope'
    y_label = 'linearity offset'
    outpath_tmp = '{}/{}_linearity_scatter.pdf'.format(outpath, file_name)
#    plotScatter(slope_array, offset_array, label_scatter, x_label, y_label, outpath_tmp, size=5, xlim=(1.02,1.2), ylim=(-0.075, -0.01))
    plotScatter(slope_array, offset_array, label_scatter, x_label, y_label, outpath_tmp, size=5)


    # plot linearity distance vs uncertainty on unfolded parameter
    x_label = 'Unfolded uncertainty'
    y_label = 'linearity distance'
    outpath_tmp = '{}/{}_uncert_vs_linearity_scatter.pdf'.format(outpath, file_name)
    plotScatter(uncert_array, distance_array, label_scatter, x_label, y_label, outpath_tmp, size=5)


    ##########################################################################
    ## find best binnings: select best binnings for each tested bin numbers ##
    ##########################################################################

    # normalise features used to determine the best binnings and build dictionary
    slope_flat = np.mean(np.array(slope_flat), axis=1)
    offset_flat = np.mean(np.array(offset_flat), axis=1)

    slope = normaliseFeature(np.abs(1-np.array(slope_flat)))
    offset = normaliseFeature(np.abs(offset_flat))
    uncert = normaliseFeature(uncert_flat)

    distance = np.sqrt( (slope**2 + offset**2 + uncert**2 ) )
    distance_dict = { key_array[i] : distance[i] for i in range(len(key_array)) }
                    

    print('Searching for the binnings achieving the best linearities...')
    best_binnings = []

    for ibin in nbin:
        distance_dict_bin = { key:value for key, value in distance_dict.items() if ibin == key[0:2] or ibin == key[0] }
        best_binning, best_parameter = find_binning_with_the_smallest_parameter(distance_dict_bin)
        print('The best binning for the {}-bin configurations is {}'.format(ibin, best_binning))
        best_binnings.append(best_binning)

        next_parameter = best_parameter
        number_of_binnings_passed = 0
        if len(distance_dict_bin) < args.nbinnings:
            nbinnings = len(distance_dict_bin)
        else:
            nbinnings = args.nbinnings

        while number_of_binnings_passed < nbinnings-1:

            next_binning, new_parameter = find_other_binnings_with_smallest_parameter(next_parameter, distance_dict_bin)
            next_parameter = new_parameter

            pass_to_file = True

            for item_in_list in range(len(best_binnings)):
                if best_binnings[item_in_list] == next_binning:
                    pass_to_file = False

            if pass_to_file == True:
                number_of_binnings_passed += 1
                print('next Binning = ', next_binning)

            best_binnings.append(next_binning)	
        
    print('\nThe selected binnings are: {}'.format(best_binnings))
    print('Saving them to {}/config/best_binnings.json'.format(args.inpath))
    writejson( '{}/config/best_binnings.json'.format(args.inpath), { b:(binnings[b], distance_dict[b]) for b in best_binnings } )


    ##############################################
    ## run pseudo-experiments for best binnings ##
    ##############################################

    if args.pseudoexp:
        print('\nRunning pseudo-experiments for the selected best binnings')

        for binning in best_binnings:
            print('\nLaunching pseudo-experiments for binning {}'.format(binning))

            # generate toy datasets to perform pseudoexperiments
            for channel in args.channel:
                makeToys('{}/unfoldingInput/{}/{}/{}/{}.root'.format(args.inpath, binning, args.campaign, channel, file_name), 'mcdata', args.ntoys, is_2d=False)

            for itoy in range(args.ntoys):
                ## Launch job on batch system: unfolding
                outdir_input_tmp = '{}/unfoldingInput/{}/'.format(args.inpath, binning)
                outdir_output_tmp = '{}/results/{}/'.format(args.inpath, binning)
                outdir_tmp = '{}/config/{}/'.format(args.inpath, binning)
                os.makedirs(outdir_output_tmp, exist_ok=True)
                command = 'python scripts/runTRExFitter.py -i {} -o {} --config {}/config.json -v {} -c {} --campaign {} -d toymcdata_{}\n'.format(outdir_input_tmp, outdir_output_tmp, outdir_tmp, args.variable, ' '.join(args.channel), args.campaign, itoy)
                command += 'rm -rf {0}/{1}_{2}_{3}/toymcdata{4}/Histograms/ {0}/{1}_{2}_{3}/toymcdata{4}/Systematics/ {0}/{1}_{2}_{3}/toymcdata{4}/RooStats/'.format(outdir_output_tmp, args.campaign, '_'.join(args.channel), file_name, itoy)
                submit_job('unfold_toymcdata_{}'.format(itoy), binning, args.inpath, command)




if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter

    parser = ArgumentParser('Find best binning configurations based on linearity test', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--inpath', type=str, help='path to input, i.e. output of binning-optimisation linearity tests')
    parser.add_argument('-v','--variable', type=str, help='name of variable')
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels', required=True)
    parser.add_argument('--campaign', type=str, help='MC campaign to run', required=True)
    parser.add_argument('-g', '--config', type=str, default='config/config.json', help='Path to config file')
    parser.add_argument('-n','--nbinnings', type=int, help='number of binnings picked', default=2)
    parser.add_argument('-p','--pseudoexp', action='store_true', help='run pseudo-experiments for picked best binnings')
    parser.add_argument('--ntoys', type=int, help='number of pseudoexperiments', default=1000)
    parser.add_argument('--plot-ensemble', action='store_true', help='plot results of ensemble tests')
    args = parser.parse_args()

    if args.plot_ensemble == False:
        pickBinnings(args)
    else:
        plotEnsembleTest(args)
