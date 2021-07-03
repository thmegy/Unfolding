import matplotlib as mpl
mpl.use('Agg') # load Agg backed which does not require X11 connection to make plots
import numpy as np
import json, hjson
import math
from utils import writejson, getDiffBinLabel
import matplotlib.pyplot as plt
from iminuit import Minuit
import sys, os



def fitandplot(args, outputdir, plot_name, list1, error_up, error_down, list1name, list2, list2name, diff_name, config, rebin_config, labels, subplot=False):

    if subplot == True:
        fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        ax[0].errorbar(list2, list1, yerr=[error_down, error_up], barsabove=True, fmt='o', ecolor='g', color='g')    
    else:
        fig, ax = plt.subplots()
        ax.errorbar(list2, list1, yerr=[error_down, error_up], barsabove=True, fmt='o', ecolor='g', color='g')    

        
    def f(x, a, b): return b + a*x
    
    def chi2(a, b):
        c2 = 0.
        for x, y, yerr in zip(list2, list1, (error_down+error_up)/2):
            c2 += (f(x, a, b) - y)**2 / yerr**2
        return c2

    m = Minuit(chi2)
    m.migrad()
    m.hesse()
    a = m.values['a']
    b = m.values['b']
    a_err = m.errors['a']
    b_err = m.errors['b']

    toret = [a, a_err, b, b_err]
    
    list2_sorted = sorted(list2)
    x = np.linspace(list2_sorted[0]*1.1, list2_sorted[-1]*1.1, 10)
    y = f(x,a,b)

    if subplot ==True:

        ax[0].plot(x, y, 'r')
        ax[1].errorbar(x, y-x, barsabove=True, fmt='+',markersize=4,capsize=3,ecolor='g',label='ensemble')
        #create plots with slope and offset displayed, numbers rounded on 5 decimal places
        ax[0].text(0.55, 0.15, r'slope = %.5f +/- %.5f'%(a,a_err),horizontalalignment='left', transform=ax[0].transAxes)
        ax[0].text(0.55, 0.05, r'offset = %.5f +/- %.5f'%(b,b_err),horizontalalignment='left', transform=ax[0].transAxes)   
        ax[0].grid(True)
        ax[0].set_ylabel(list1name + ' ' + labels['parameters'][args.variable])
        ax[1].grid(True)
        ax[1].set_xlabel(list2name + ' ' + labels['parameters'][args.variable])
        ax[1].set_ylabel('Fit - Model')

    else:
        ax.plot(x, y, 'r')
        #create plots with slope and offset displayed, numbers rounded on 5 decimal places
        ax.text(0.55, 0.13, r'slope = %.5f +/- %.5f'%(a,a_err),horizontalalignment='left', transform=ax.transAxes)
        ax.text(0.55, 0.05, r'offset = %.5f +/- %.5f'%(b,b_err),horizontalalignment='left', transform=ax.transAxes)
        ax.grid(True)
        ax.set_xlabel(list2name + ' ' + labels['parameters'][args.variable])       
        ax.set_ylabel(list1name + ' ' + labels['parameters'][args.variable]) 
        ax.text(0.04, 0.96, r'ATLAS Internal', transform=ax.transAxes, horizontalalignment='left', verticalalignment='top')
        t = ax.text(0.04, 0.885, r'$\sqrt{s}=13$ TeV, $\mathcal{L}$ = ' + str(config['lumi'][args.campaign][0]) + r' fb$^{-1}$', transform=ax.transAxes, horizontalalignment='left', verticalalignment='top')
        t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
        ax.text(0.04, 0.81, diff_name , transform=ax.transAxes, horizontalalignment='left', verticalalignment='top')


    fig.savefig('{}/{}.pdf'.format(outputdir, plot_name))
    print('\nSaving {}/{}.pdf'.format(outputdir, plot_name))
    
    return toret




def main(args):
    # check if several differential bins
    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    labels = hjson.load(open(config['labels'], 'r'))

    file_name = args.variable
    plot_name = ''
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])
        diff_bins = rebin_config['diff_bins']

    Ndiff = len( [b for b in rebin_config['binning'] if len(b)>2] ) # number of diff bins with more than one bin

    outpath = '{}/{}_{}_{}/{}/'.format(args.inpath, args.campaign, '_'.join(args.channel), file_name, args.outpath)
    os.makedirs(outpath, exist_ok=True)


    # load truth spin parameter values (dictionary for variables) from json (computed from histograms)
    try:
        truth = json.load( open( '{}/{}/{}/{}/truth_{}_binned.json'.format(args.truthdir, args.campaign, args.channel[0], config['btag_regions'][args.channel[0]][0], file_name) ) )
    except KeyError:
        print('ERROR: The truth values are not available for the {} observable!'.format(args.variable))
        sys.exit(3)

    # load unbinned truth spin-parameter values (i.e. true values)
    try:
        truth_unbinned = json.load( open( '{}/{}/{}/{}/truth_{}_unbinned.json'.format(args.truthdir, args.campaign, args.channel[0], config['btag_regions'][args.channel[0]][0], file_name) ) )
    except KeyError:
        print('ERROR: The unbinned truth values are not available for the {} observable!'.format(args.variable))
        sys.exit(3)

    print('Considering the following reweightings: {}'.format(args.reweightlist))


    unfold_mean = []
    unfold_uncert_up = []
    unfold_uncert_down = []
    truth_binned_mean = []
    truth_unbinned_mean = []
    truth_unbinned_std = []

    for i in args.reweightlist+['mcdata']:
        unfold = json.load( open('{}/{}_{}_{}/{}/{}_unfolded_parameter.json'.format(args.inpath, args.campaign, '_'.join(args.channel), file_name, i, args.variable)) )

        unfold_mean.append( unfold['mean'] )
        unfold_uncert_up.append( unfold['uncert_up'] )
        unfold_uncert_down.append( unfold['uncert_down'] )
        truth_binned_mean.append( truth[i.replace('mcdata', 'truth')] )
 #       truth_unbinned_mean.append( truth_unbinned[i.replace('mcdata', 'truth')][0] )
 #       truth_unbinned_std.append( truth_unbinned[i.replace('mcdata', 'truth')][1] )
        
    # plot linearity test: unfolded value vs truth value (computed from histogram)
    slope = []
    slope_err = []
    offset = []
    offset_err = []

    for idiff in range(Ndiff): # loop over differential bins
        if rebin_config['differential_measurement']:
            plot_name = '_{}_{}'.format(rebin_config['diff_variable'], '_'.join(diff_bins[idiff]))
            diff_name = getDiffBinLabel(labels, rebin_config['diff_variable'], rebin_config['diff_variable_unit'], diff_bins[idiff])
        else:
            diff_name = 'inclusive'

        result = fitandplot(args, outpath, 'linearity_test{}'.format(plot_name), 
                            np.array(unfold_mean)[:,idiff], np.array(unfold_uncert_up)[:,idiff], np.array(unfold_uncert_down)[:,idiff], 'Unfolded', np.array(truth_binned_mean)[:,idiff], 'True', 
                            diff_name, config, rebin_config, labels, args.subplots)

        print('\n\n unfolded = ({}+-{}) * truth + ({}+-{})\n\n'.format(result[0], result[1], result[2], result[3]))
        slope.append( result[0] )
        slope_err.append( result[1] )
        offset.append( result[2] )
        offset_err.append( result[3] )

    ParamData = {}
    ParamData['slope'] = slope
    ParamData['slope_err'] = slope_err
    ParamData['offset'] = offset
    ParamData['offset_err'] = offset_err
    
    writejson('{}/linearity_test.json'.format(outpath), ParamData)
    print('Saving results to {}/linearity_test.json\n'.format(outpath))


    # plot mean calibration test: binned truth value (i.e. computed from histogram) vs true value
#    for idiff in range(Ndiff): # loop over differential bins
#        if rebin_config['differential_measurement']:
#            plot_name = '_{}_{}'.format(rebin_config['diff_variable'], '_'.join(diff_bins[idiff]))
#        result = fitandplot(outpath, 'mean_calibration{}'.format(plot_name), np.array(truth_unbinned_mean), np.array(truth_unbinned_std), np.array(truth_unbinned_std), 'True', np.array(truth_binned_mean)[:,idiff], 'Binned', args.variable, args.subplots)
#
#        print('\n\n binned = ({}+-{}) * true + ({}+-{})\n\n'.format(result[0], result[1], result[2], result[3]))




if __name__=="__main__":

    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter\
        
    parser = ArgumentParser(description='Perform linearity test and mean calibration', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--reweightlist', nargs='*', default=['mcdata_rw0', 'mcdata_rw1', 'mcdata_rw2', 'mcdata_rw3', 'mcdata_rw4'], help='identifiers of signal reweightings')
    parser.add_argument('--config', type=str, help='Config file', default='config/config.json')
    parser.add_argument('-v', '--variable', type=str, help='name of spin observable')
    parser.add_argument('-c', '--channel', nargs='*', help='name of dileptonic channel')
    parser.add_argument('--campaign', type=str, help='name of MC campaign')
    parser.add_argument('-i', '--inpath',  help='path to fit output', default='results/')
    parser.add_argument('-o', '--outpath',  help='path to linearity test output', default='linearity_test/')
    parser.add_argument('--truthdir',  help='path to truth parameters', default='unfoldingInput/')
    parser.add_argument('-d','--subplots', help='Set if subplot with difference between fit and ideal y=x function should be included.', action='store_true')
    args = parser.parse_args()

    main(args)
