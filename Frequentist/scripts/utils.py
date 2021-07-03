import json
import hjson
import numpy as np
import matplotlib as mpl
mpl.use('Agg') # load Agg backed which does not require X11 connection to make plots
import matplotlib.pyplot as plt



def submit_job(job_type, job_name, outdir, command, dependency=None, mem = '2gb'):
    import subprocess
    from os import makedirs
    
    #Create bash script to run with torque
    queue = 'prod2C7'
    makedirs('{}/batchSubmission/logs/'.format(outdir), exist_ok=True)
    scriptname = '{}/batchSubmission/{}_job_{}.sh'.format(outdir, job_type, job_name)
    
    script = open(scriptname, 'w')
    script.write('#!/bin/bash\n')
    script.write('#PBS -j oe\n')
    script.write('#PBS -o {}/batchSubmission/logs/\n'.format(outdir))
    script.write('#PBS -N {}_job_{}\n'.format(job_type, job_name))
    script.write('cd\n')
    script.write('cd home2/PL_unfolding\n')
    script.write('source init.sh\n')
    script.write(command)
    script.close()

    job_options = '-q {}@clratlserv04'.format(queue)
    if dependency is not None:
        job_options += ' -W depend=afterok:{}'.format(dependency)
    job_options += ' -l mem={0} -l vmem={0}'.format(mem)

    output = subprocess.getoutput('qsub {} {}'.format(job_options, scriptname))
    print(output)
    return output




def writejson(fname, datatowrite):
    """
    Write data into JSON file. Standard data structures such as dictionaries are
    supported. In addition, writing of numpy arrays is also supported.
    Parameters:
        fname -- Name of the output JSON file. Can ommit .json suffix
        datatowrite -- the object to write into JSON
    """
    if len(fname) < 5 or fname[-5:] != '.json':
        fname += '.json'
    with open(fname, 'w') as outfile:
        outfile.write(json.dumps(datatowrite, cls=MyJSONEncoder, indent=4))

        


class NumpyEncoder(json.JSONEncoder):
    """
    For converting numpy arrays into lists that can be written to JSON file
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

    

    
class MyJSONEncoder(json.JSONEncoder):
    """
    For pretty-printing JSON in a more human-readable form
    """
    def __init__(self, *args, **kwargs):
        super(MyJSONEncoder, self).__init__(*args, **kwargs)
        self.current_indent = 0
        self.current_indent_str = ""

    def encode(self, o):
        from numpy import ndarray
        # Special Processing for lists
        if isinstance(o, (list, tuple, ndarray)):
            primitives_only = True
            for item in o:
                if isinstance(item, (list, tuple, dict, ndarray)):
                    primitives_only = False
                    break
            output = []
            if primitives_only:
                for item in o:
                    output.append(json.dumps(item, cls=NumpyEncoder))
                return "[ " + ", ".join(output) + " ]"
            else:
                self.current_indent += self.indent
                self.current_indent_str = "".join([" " for x in range(self.current_indent)])
                for item in o:
                    output.append(self.current_indent_str + self.encode(item))
                self.current_indent -= self.indent
                self.current_indent_str = "".join([" " for x in range(self.current_indent)])
                return "[\n" + ",\n".join(output) + "\n" + self.current_indent_str + "]"
        elif isinstance(o, dict):
            output = []
            self.current_indent += self.indent
            self.current_indent_str = "".join([" " for x in range(self.current_indent)])
            for key, value in o.items():
                output.append(self.current_indent_str + json.dumps(key, cls=NumpyEncoder) + ": " + self.encode(value))
            self.current_indent -= self.indent
            self.current_indent_str = "".join([" " for x in range(self.current_indent)])
            return "{\n" + ",\n".join(output) + "\n" + self.current_indent_str + "}"
        else:
            return json.dumps(o, cls=NumpyEncoder)

        


class Axisoffset():
    """
    Workaround class for lack of axis exponent location setting in matplotlib.
    Parameters:
        ax -- instance of matplotlib Axes
        axis -- string, either 'x' or 'y'
        pos -- tuple of (x,y) position of the exponent location
            (coordinates relative [0, 1])
        label, label_args (optional) -- if instead want to append the axis
            exponent string to the axis label as in "label [x10^N]"
    """
    def __init__(self, ax, axis, pos, label=None, label_args=None):
        self.axis = axis
        self.axis_dict = {"y":ax.yaxis, "x":ax.xaxis}
        self.axes = ax
        self.pos = pos
        self.label = label
        if label_args is None:
            self.label_args = {}
        else:
            self.label_args = label_args
        ax.callbacks.connect(axis+'lim_changed', self.update)
        ax.figure.canvas.draw()
        self.update()

    def update(self):
        fmt = self.axis_dict[self.axis].get_major_formatter()
        self.axis_dict[self.axis].offsetText.set_visible(False)
        if self.label is None:
            self.axes.text(self.pos[0], self.pos[1], fmt.get_offset(), transform=self.axes.transAxes,
                    horizontalalignment='right', verticalalignment='bottom')
        else:
            set_label = {'x': self.axes.set_xlabel, 'y': self.axes.set_ylabel}
            if fmt.get_offset() != '':
                offset_str = ' [' + fmt.get_offset() + ']'
            else:
                offset_str = ''
            set_label[self.axis](self.label + offset_str, **self.label_args)
        



def getSpinParameter(config, fitoutput, variable, binning, diff_bins):
    '''
    Measure spin parameter and its uncertainty using MC method to propagate errors on the unfolded bins.
    '''
    import yaml

    Ndiff = len(binning) # how many differential bins
    
    factor = config['factor'][variable]

    ## Get correlation matrix for truth bins
    with open('{}/CorrelationMatrix.yaml'.format(fitoutput)) as f:
        corrmat = np.array( yaml.full_load(f)[1]['correlation_rows'] )

    ## Get truth bins content and uncertainties
    with open('{}/UnfoldingData.yaml'.format(fitoutput)) as f:
        bins = yaml.full_load(f)

    ## Keep only correlations for truth bins
    nbins = len(bins)
    if nbins != corrmat.shape[0]:
        corrmat = corrmat[0:nbins, 0:nbins]

    ## Compute mean and mean error
    result_dict = {}
    mean_list = []
    uncert_up_list = []
    uncert_down_list = []
    iglob = 0
    
    for ibin in range(Ndiff):
        if diff_bins[ibin] == 'inclusive':
            diffbin_name = diff_bins[ibin]
        else:
            diffbin_name = ' '.join(diff_bins[ibin])

        x = np.array( [ (binning[ibin][i] + binning[ibin][i+1]) / 2 for i in range(len(binning[ibin])-1) ] )
        corrmat_tmp = corrmat[iglob:iglob+len(binning[ibin])-1, iglob:iglob+len(binning[ibin])-1] # get correlation matrix only for unfolded bins of differential bin
        
        import scipy.stats
        unfolded_hist = np.array( [ bins[iglob+i]['mean'] for i in range(len(binning[ibin])-1) ] )
        unfolded_err_up = np.array( [  bins[iglob+i]['uncertaintyUp'] for i in range(len(binning[ibin])-1) ] )
        unfolded_err_down = np.array( [  bins[iglob+i]['uncertaintyDown'] for i in range(len(binning[ibin])-1) ] )
        x = np.array(x)
        covmat_up = (corrmat_tmp.T * unfolded_err_up).T * unfolded_err_up
        covmat_down = (corrmat_tmp.T * unfolded_err_down).T * unfolded_err_down
        
        pdf_up = scipy.stats.multivariate_normal(mean=unfolded_hist, cov=covmat_up)
        mc_up = pdf_up.rvs(1000000)
        mean_mc_up = np.sum(mc_up * x, axis=1) / np.sum(mc_up, axis=1)
        pdf_down = scipy.stats.multivariate_normal(mean=unfolded_hist, cov=covmat_down)
        mc_down = pdf_down.rvs(1000000)
        mean_mc_down = np.sum(mc_down * x, axis=1) / np.sum(mc_down, axis=1)
        print('\nDifferential bin {}'.format(diffbin_name))
        print( '{:.4f} + {:.4f} - {:.4f}'.format(factor*np.sum(unfolded_hist*x)/np.sum(unfolded_hist), abs(factor*np.std(mean_mc_up)), abs(factor*np.std(mean_mc_down)) ))

        mean_list.append(factor*np.mean( mean_mc_up ))
        uncert_up_list.append(abs(factor)*np.std( mean_mc_up ))
        uncert_down_list.append(abs(factor)*np.std( mean_mc_down ))

        iglob += len(binning[ibin])-1
            
    result_dict['mean'] = mean_list
    result_dict['uncert_up'] = uncert_up_list
    result_dict['uncert_down'] = uncert_down_list

    writejson('{}/spin_parameter.json'.format(fitoutput), result_dict)
    print('\nThe unfolded parameters are saved to: {}/spin_parameter.json\n'.format(fitoutput))




def getNPsFromFitResults(filepath, get_param=False):
    import ROOT

    rfile = ROOT.TFile(filepath, 'read')
    try:
        fitres = rfile.Get('nll_simPdf_newasimovData_with_constr')
        fitres.floatParsFinal()
    except:
        fitres = rfile.Get('nll_simPdf_newasimovData')

    NP_postfit = {}
    for param in fitres.floatParsFinal():
        if get_param:
            if 'param' in param.GetName():
                NP_postfit[param.GetName()] = [ param.getValV(),
                                                param.getErrorHi(),
                                                param.getErrorLo() ]
        else:
            if 'alpha' in param.GetName() or 'gamma' in param.GetName():
                NP_postfit[param.GetName()] = [ param.getValV(),
                                                param.getErrorHi(),
                                                param.getErrorLo() ]
            
    rfile.Close()
    return NP_postfit




def makeResponseMatrix(truth, migmat):
    # make response matrix from migration matrix and truth distribution
    resmat = migmat.Clone()
    for tbin in range(truth.GetNbinsX()):
        for rbin in range(migmat.GetNbinsX()):
            resmat.SetBinContent( rbin+1, tbin+1, float(migmat.GetBinContent(rbin+1, tbin+1) / truth.GetBinContent(tbin+1)) )
            if migmat.GetBinContent(rbin+1, tbin+1) > 0 and truth.GetBinContent(tbin+1) > 0:
                error = (resmat.GetBinContent(rbin+1, tbin+1) * 
                         np.sqrt( (migmat.GetBinError(rbin+1, tbin+1)/migmat.GetBinContent(rbin+1, tbin+1))**2 + (truth.GetBinError(tbin+1)/truth.GetBinContent(tbin+1))**2 ))
            else: 
                error = 0.
            resmat.SetBinError( rbin+1, tbin+1, error )
    return resmat




def responseTimesTruth(truth, resmat):
    weight = 0
    for tbin in range(truth.GetNbinsX()):
        for rbin in range(resmat.GetNbinsX()):
            weight += resmat.GetBinContent(rbin+1, tbin+1) * truth.GetBinContent(tbin+1)
    return weight




def getBinnedMean(hist, binning):
    bins_center = np.array( [ (binning[i] + binning[i+1]) / 2 for i in range( len(binning)-1 ) ] )  # Get center of bins
    mean = np.sum(hist*bins_center) / np.sum(hist)
    return mean




def rootToNumpy(hist, is_2d=False):
    # Turn root histogram into numpy array
    new_hist = []
    if not is_2d:
        for i in range(hist.GetNbinsX()):
            new_hist.append(hist.GetBinContent(i+1))
    else:
        for j in range(hist.GetNbinsY()):
            row = []
            for i in range(hist.GetNbinsX()):
                row.append(hist.GetBinContent(i+1, j+1))
            new_hist.append(row)

    return np.array(new_hist)




def normaliseFeature(array):
    '''
    Normalise feature so that all values are between 0 and 1.
    '''
    array = np.array(array)
    array_norm = ( array - np.min(array) ) / ( np.max(array) - np.min(array) )

    return array_norm




def makeToys(file_name, hist_name, ntoys, is_2d=False, pdf='Poisson'):
    import ROOT
    r = ROOT.TRandom3(0)

    rfile = ROOT.TFile(file_name, 'update')
    hist = rfile.Get(hist_name)
    
    for itoy in range(ntoys):

        hist_new = hist.Clone('mcdata_toy{}'.format(itoy))

        if type(hist) == ROOT.TH2D:
            for i in range( 1, hist_new.GetNbinsX()+1 ):
                for j in range( 1, hist_new.GetNbinsY()+1 ):
                    if pdf=='Poisson':
                        hist_new.SetBinContent( i, j, r.Poisson(hist_new.GetBinContent(i,j)) )
                    elif pdf=="Gaus":
                        hist_new.SetBinContent( i, j, r.Gaus(hist_new.GetBinContent(i,j), hist_new.GetBinError(i,j)) )
                    elif pdf=="MCstat":
                        Neff = (hist_new.GetBinContent(i,j) / hist_new.GetBinError(i,j))**2 # effective number of MC events
                        hist_new.SetBinContent( i, j, r.Poisson(Neff) * hist_new.GetBinContent(i,j) / Neff )
        else:
            for i in range( 1, hist_new.GetNbinsX()+1 ):
                if pdf=='Poisson':
                    hist_new.SetBinContent( i, r.Poisson(hist_new.GetBinContent(i)) )
                elif pdf=="Gaus":
                    hist_new.SetBinContent( i, r.Gaus(hist_new.GetBinContent(i), hist_new.GetBinError(i)) )
                elif pdf=="MCstat":
                    Neff = (hist_new.GetBinContent(i) / hist_new.GetBinError(i))**2 # effective number of MC events
                    hist_new.SetBinContent( i, r.Poisson(Neff) * hist_new.GetBinContent(i) / Neff )

        rfile.cd()
        hist_new.Write('toy{}_{}'.format(hist_name, itoy))

    rfile.Close()




def getGammaNPLabel(np_name):
    import re
    bin_number = re.findall('\d', np_name)
    channel = re.findall('stat_..', np_name)[0].replace('stat_', '').replace('m', '\mu')
    if len(bin_number) > 1:
        label = '$\gamma_{{s,{}}}^{{{}}}$'.format(channel, ''.join(bin_number))
    else:
        label = '$\gamma_{{b,{}}}^{{{}}}$'.format(channel, bin_number[0])
    return label
        

    



#####################################################
####################   Plotting  ####################
#####################################################




def getDiffBinLabel(labels, diff_variable, unit, diff_bin):
    if len(diff_bin) == 1:
        diffbin_label = '{} > {} {}'.format(labels['differential_variables'][diff_variable], diff_bin[0], unit)
    elif len(diff_bin) == 2:
        diffbin_label = '{} < {} < {} {}'.format(diff_bin[0], labels['differential_variables'][diff_variable], diff_bin[1], unit)
    else:
        print('Unexpected number of bin edges!')

    return diffbin_label




def plotUnfoldedDistribution(path, unfolded_data, unfolded_err_up, unfolded_err_down, truth, var, campaign, config, diff_bins, binning):
    config = hjson.load( open(config, 'r') )
    labels = hjson.load( open(config['labels'], 'r') )
    rebin = hjson.load( open(config['rebin'], 'r') )
    
    Ndiff = len(diff_bins)

    truth = np.array(truth)
    unfolded_data = np.array(unfolded_data)
    unfolded_err_down = np.array(unfolded_err_down)
    unfolded_err_up = np.array(unfolded_err_up)
    
    fig, axes = plt.subplots(2, Ndiff, sharex='col', sharey='row', gridspec_kw={'height_ratios': [3,1]})


    iglob = 0
    for idiff in range(Ndiff):
        b = binning[idiff]
        Nbin = len(b)-1
        bins = [ (b[ib]+b[ib+1]) / 2 for ib in range(Nbin) ]
        binwidth = [ abs((b[ib]-b[ib+1])) / 2 for ib in range(Nbin) ]

        data = unfolded_data[iglob:iglob+Nbin]
        truth_bin = truth[iglob:iglob+Nbin]
        err_up = unfolded_err_up[iglob:iglob+Nbin]
        err_down = unfolded_err_down[iglob:iglob+Nbin]

        if Ndiff == 1:
            ax1 = axes[0]
            ax2 = axes[1]
        else:
            ax1 = axes[0, idiff]
            ax2 = axes[1, idiff]
            # write diff bin label
            ax1.text(0.04, 0.76, getDiffBinLabel(labels, rebin['diff_variable'], rebin['diff_variable_unit'], diff_bins[idiff]) , transform=ax1.transAxes, horizontalalignment='left', verticalalignment='top')

        ax1.errorbar(bins, truth_bin, xerr=binwidth, yerr=0., fmt='o', label='Truth', color='red')
        ax1.errorbar(bins, data, yerr=[err_down, err_up], fmt='o', label='Unfolded data', color='black')

        ax1.xaxis.set_tick_params(which='major', length=0, width=0)
        ax1.xaxis.set_tick_params(which='minor', length=0, width=0)

        iglob += Nbin
        
        ## draw ratio plot
        ax2.errorbar(bins, data / truth_bin, xerr=binwidth, yerr=[err_down/truth_bin, err_up/truth_bin], fmt='o', label='Truth', color='black')
        ax2.axhline(y=1, linestyle='dashed', color='red', linewidth=1)

        # draw atlas label only on leftmost subplot
        if idiff == 0:
            ax1.set_ylim(0, 1.5*max(unfolded_data))
            ax1.set_ylabel('Events')
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,3))

            ax1.text(0.04, 0.96, r'ATLAS Internal', transform=ax1.transAxes, horizontalalignment='left', verticalalignment='top')
            t = ax1.text(0.04, 0.885, r'$\sqrt{s}=13$ TeV, $\mathcal{L}$ = ' + str(config['lumi'][campaign][0]) + r' fb$^{-1}$', transform=ax1.transAxes, horizontalalignment='left', verticalalignment='top')
            t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))

            ax2.set_ylabel('Data / Truth')
        else:
            ax1.yaxis.set_tick_params(which='major', length=0, width=0)
            ax1.yaxis.set_tick_params(which='minor', length=0, width=0)
            ax2.yaxis.set_tick_params(which='major', length=0, width=0)
            ax2.yaxis.set_tick_params(which='minor', length=0, width=0)

        # draw legend only on rightmost subplot
        if idiff == Ndiff-1:
            ax1.legend()
            ax2.set_xlabel(labels['observables'][var], loc='right')
        
    fig.savefig('{}/unfolded_data.pdf'.format(path), bbox_inches='tight')
    print('Saving {}/unfolded_data.pdf\n'.format(path))




def plotSpinParameter(path, unfolded_param, unfolded_err_up, unfolded_err_down, truth_param, var, campaign, config, diff_bins):
    config = hjson.load( open(config, 'r') )
    labels = hjson.load( open(config['labels'], 'r') )
    rebin = hjson.load( open(config['rebin'], 'r') )
    
    # add maximum if differential variable has no upper bound
    xmax = None
    if len(diff_bins[-1])==1:
        xmax = 2*int(diff_bins[-1][0])

    # get center of differential bins
    bins = [ int(val) for db in diff_bins for val in db ]
    if xmax is not None:
        bins.append(xmax)
    bins = np.unique( np.array( bins ) )
    x = [ (bins[i+1] + bins[i]) / 2 for i in range(len(bins)-1) ]
    binwidth = [ (bins[i+1] - bins[i]) / 2 for i in range(len(bins)-1) ]    

    # convert data to numpy
    unfolded_param = np.array(unfolded_param)
    unfolded_err_up = np.array(unfolded_err_up)
    unfolded_err_down = np.array(unfolded_err_down)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3,1]})
    ax1.errorbar(x, truth_param, xerr=binwidth, yerr=0., fmt='o', label='Truth', color='red')
    ax1.errorbar(x, unfolded_param, yerr=[unfolded_err_down, unfolded_err_up], fmt='o', label='Unfolded data', color='black')

    if xmax is not None:
        ax1.set_xlim(int(rebin['diff_bins'][0][0]), xmax)
#    ax1.set_ylim( np.min(unfolded_param-unfolded_err_down), 1.4*np.max(unfolded_param+unfolded_err_up) )
    ax1.set_ylabel(labels['parameters'][var])

    ax1.text(0.02, 0.96, r'ATLAS Internal', transform=ax1.transAxes, horizontalalignment='left', verticalalignment='top')
    t = ax1.text(0.02, 0.885, r'$\sqrt{s}=13$ TeV, $\mathcal{L}$ = ' + str(config['lumi'][campaign][0]) + r' fb$^{-1}$', transform=ax1.transAxes, horizontalalignment='left', verticalalignment='top')
    t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))

    ax1.legend(loc=4)

    ## draw ratio plot
    ax2.errorbar(x, unfolded_param / truth_param, xerr=binwidth, yerr=[unfolded_err_down/truth_param, unfolded_err_up/truth_param], fmt='o', label='Truth', color='black')
    ax2.axhline(y=1, linestyle='dashed', color='red', linewidth=1)
    ax2.set_ylabel('Data / Truth')
    xlabel = labels['differential_variables'][rebin['diff_variable']]
    if rebin['diff_variable_unit'] != '':
        xlabel += f' [{rebin["diff_variable_unit"]}]'
    ax2.set_xlabel(xlabel, loc='right')
    
    fig.savefig('{}/unfolded_param.pdf'.format(path), bbox_inches='tight')
    print('Saving {}/unfolded_param.pdf\n'.format(path))




def plotResmat(name, resmat, config, outpath, channel, btag, var, campaign, diff_bins, binning, docmap=False):
    config = hjson.load( open(config, 'r') )
    labels = hjson.load( open(config['labels'], 'r') )
    rebin = hjson.load( open(config['rebin'], 'r') )
    fontProperties = {'family':'sans-serif'}

    # check if channel is used as control region
    if rebin['is_control'] and channel in rebin['control_channels'].keys():
        is_control = True
    else:
        is_control = False
    
    if docmap:
        fig, ax = plt.subplots( figsize=(resmat.shape[1]*1.25, resmat.shape[0]) )
        opts = {'cmap': 'Greens'}
        heatmap = ax.pcolor(resmat, **opts)
        cbar = plt.colorbar(heatmap, ax=ax)
        #cbar.ax.set_yticklabels( [i.get_text().strip('$') for i in cbar.ax.get_yticklabels()], **fontProperties )
    else:
        fig, ax = plt.subplots( figsize=(resmat.shape[1], resmat.shape[0]) )
        opts = {'cmap': 'Greens', 'vmin': 0, 'vmax': +1}
        heatmap = ax.pcolor(resmat, **opts)
        for irow in range(resmat.shape[0]):
            for icol in range(resmat.shape[1]):
                ax.text(icol+0.5, irow+0.5, '{:.3f}'.format(resmat[irow][icol]), ha="center", va="center", color="black")

    ylabels = [ '[{},{}]'.format(b[i], b[i+1]) for b in binning for i in range(len(b)-1) if len(b)>2 ]
    ax.set_yticks(np.arange(len(ylabels))+0.5, minor=False)
    ax.set_yticklabels(labels=ylabels, fontdict=fontProperties, fontsize=8, minor=False)
    if is_control:
        bounds = rebin['control_variable_bounds'][rebin['control_channels'][channel]]
        xlabels = np.linspace(bounds[0], bounds[1], 11).astype(np.int)
        ax.set_xticks(np.linspace(0, rebin['nbin_control'], len(xlabels)), minor=False)
        ax.set_xticklabels(labels=xlabels, fontdict=fontProperties, fontsize=8, minor=False)
    else:
        xlabels = ylabels
        ax.set_xticks(np.arange(len(xlabels))+0.5, minor=False)
        ax.set_xticklabels(labels=xlabels, fontdict=fontProperties, fontsize=8, minor=False, rotation=45, ha='right')
        
    if len(diff_bins) > 1:
        labels_minor = [ getDiffBinLabel(labels, rebin['diff_variable'], rebin['diff_variable_unit'], db) for db in diff_bins ]
        xticks_minor = []
        yticks_minor = []
        lines = []
        nglob = 0
        for b in binning:
            xticks_minor.append( (len(b)-1)/2 + nglob )
            yticks_minor.append( (len(b)-1)/2 + nglob + 0.4 )
            nglob += len(b)-1
            lines.append(nglob)
        if not is_control:
            ax.set_xticks(xticks_minor, minor=True)
            ax.set_xticklabels(labels=labels_minor, fontdict=fontProperties, fontsize=8, minor=True)
            ax.xaxis.set_tick_params(which='minor', pad=45, length=0, width=0)
        ax.set_yticks(yticks_minor, minor=True)
        ax.set_yticklabels(labels=labels_minor, fontdict=fontProperties, fontsize=8, rotation=90, minor=True)
        ax.yaxis.set_tick_params(which='minor', pad=45, length=0, width=0)

        ## draw vertical and horizontal lines to separate differential bins
        for i in range(len(lines)-1):
            if not is_control:
                ax.axvline(x=lines[i], linestyle='dashed', color='black')
            ax.axhline(y=lines[i], linestyle='dashed', color='black')

    ax.set_ylabel(f'True {labels["observables"][var]}')
    if is_control:
        ax.set_xlabel(f'Reconstructed {labels["differential_variables"][rebin["control_channels"][channel]]} [{rebin["control_variable_unit"]}]')
    else:
        ax.set_xlabel(f'Reconstructed {labels["observables"][var]}')
    ax.set_title('ATLAS simulation\n{} = 13 TeV, {}={}{}\n{}, {}'.format('$\\sqrt{s}$', '$\mathcal{L}$', 
                                                                        config['lumi'][campaign][0], 'fb$^{-1}$', labels['regions'][channel], labels['regions'][btag]), loc='left')

    
    fig.savefig('{}/{}_{}_{}.pdf'.format(outpath, name, channel, btag), bbox_inches='tight')        
    print('Saving {}/{}_{}_{}.pdf'.format(outpath, name, channel, btag))




def plotBrasilian(outpath, NP_dict, labels):
    NP_labels = [ labels['NP_names'][syst.replace('alpha_', '')] for syst in NP_dict if 'alpha' in syst ]
    NP_values = np.array( [ val for key, val in NP_dict.items() if 'alpha' in key ] )
    
    if len(NP_labels) > 0:
        y_pos = np.array([i for i in range(len(NP_values))], dtype=float)

        fig, ax = plt.subplots(figsize=(8, len(y_pos)*0.35))
        
        ax.fill_betweenx([y_pos[0] - 1, y_pos[-1] + 1], x1=-1, x2=1, color='#33cc33', zorder=1)
        # 2-sigma yellow band
        ax.fill_betweenx([y_pos[0] - 1, y_pos[-1] + 1], x1=-2, x2=-1, color='#ffed40', zorder=1)
        ax.fill_betweenx([y_pos[0] - 1, y_pos[-1] + 1], x1=1, x2=2, color='#ffed40', zorder=1)
        
        ax.errorbar(NP_values[:,0], y_pos, xerr=[-NP_values[:,2],NP_values[:,1]],
                    color='k', marker='o', capsize=0, elinewidth=3, markersize=12, linestyle='')
        
        ax.axvline(x=0, color='k', linestyle='dashed')
        ax.set_xlabel(r'$\hat{\theta}$')
        ax.set_xlim(-3,3)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(NP_labels)
        ax.set_ylim(y_pos[0] - 1, y_pos[-1] + 1)
        ax.margins(y=0.001)
        
        outname = '{}/NP_pulls.pdf'.format(outpath)
        fig.savefig(outname, bbox_inches='tight')
        print('\nSaving {}'.format(outname))
        



def plotEnsembleTestPull(pull_array, constraint_array, label_yaxis, label_xticks, outpath):
    fig, ax = plt.subplots()
    x_pos = np.arange(len(pull_array))
    ax.errorbar(x_pos, pull_array, yerr=constraint_array, color='k', marker='o', capsize=0, markersize=4, linestyle='', label='Ensemble test pull')
    ax.set_ylim(-1.7, 1.7)
    ax.set_ylabel(label_yaxis)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(label_xticks, rotation=45, ha='right')
    ax.set_yticks([-1, 0, 1], minor=False)
    ax.set_yticks([-1.5, -0.5, 0.5, 1.5], minor=True)
    ax.yaxis.grid(True, which='major', ls='--')
    ax.legend()
    fig.savefig(outpath, bbox_inches='tight')
    print('Saving {}'.format(outpath))




def plotEnsembleTestHist(mean_array, label_xaxis_hist, pull, constraint, label_xaxis_pull, outpath):
    fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4,1]})

    # ensemble histogram
    axes[0].hist(mean_array, label='unfolded parameter')
    axes[0].axvline(x=np.mean(mean_array), linestyle='solid', color='C1', label='{}$_{{ensemble}}$={:.3f}\n$\\Delta${}={:.3f}'.format(label_xaxis_hist, np.mean(mean_array), label_xaxis_hist, np.std(mean_array)))
    axes[0].set_xlabel(label_xaxis_hist)
    axes[0].set_ylabel('Pseudo-experiments')
    axes[0].legend()

    # pull plot
    axes[1].errorbar(pull, 1, xerr=constraint, color='k', marker='o', capsize=0, markersize=4, linestyle='', label='Ensemble test pull')
    axes[1].set_xlabel(label_xaxis_pull)
    axes[1].set_xticks([-1, 0, 1], minor=False)
    axes[1].set_xlim([-1.25, 1.25])
    axes[1].xaxis.grid(True, which='major', ls='--')
    axes[1].yaxis.set_visible(False)

    fig.tight_layout()
    fig.savefig(outpath, bbox_inches='tight')
    print('Saving {}\n'.format(outpath))




def plotScatter(x_array, y_array, label_array, x_label, y_label, outpath, legendonframe=True, xlim=None, ylim=None, size=20):
    fig, ax = plt.subplots()
    for x, y, label in zip(x_array, y_array, label_array):
        sc = ax.scatter(x, y, label=label, s=size)
        
        # draw arrow if point is outside of defined range
        if (ylim is not None and xlim is None):
            if not isinstance(y, list) and not isinstance(y, np.ndarray):
                if y > ylim[1]:
                    ax.annotate( '', xy=(x, 0.92*ylim[1]), xytext=(0, 24), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                elif y < ylim[0]:
                    ax.annotate( '', xy=(x, 0.08*ylim[1]), xytext=(0, -24), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
            else:
                for yp, xp in zip(y, x):
                    if yp > ylim[1]:
                        ax.annotate( '', xy=(xp, 0.92*ylim[1]), xytext=(0, 24), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                    elif yp < ylim[0]:
                        ax.annotate( '', xy=(xp, 0.08*ylim[1]), xytext=(0, -24), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                        
        elif (xlim is not None and ylim is None):
            if not isinstance(x, list) and not isinstance(x, np.ndarray):
                if x > xlim[1]:
                    ax.annotate( '', xy=(0.92*xlim[1], y), xytext=(24, 0), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                elif x < xlim[0]:
                    ax.annotate( '', xy=(0.08*xlim[1], y), xytext=(-24, 0), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
            else:
                for yp, xp in zip(y, x):
                    if xp > xlim[1]:
                        ax.annotate( '', xy=(0.92*xlim[1], yp), xytext=(24, 0), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                    elif xp < xlim[0]:
                        ax.annotate( '', xy=(0.08*xlim[1], yp), xytext=(-24, 0), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                        
        elif (ylim is not None and xlim is not None):
            if not isinstance(x, list) and not isinstance(x, np.ndarray):
                if not (y > ylim[1] or y < ylim[0]) and (x > xlim[1] or x < xlim[0]):
                    if x > xlim[1]:
                        ax.annotate( '', xy=(0.92*xlim[1], y), xytext=(24, 0), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                    elif x < xlim[0]:
                        ax.annotate( '', xy=(0.08*xlim[1], y), xytext=(-24, 0), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                elif (y > ylim[1] or y < ylim[0]) and not (x > xlim[1] or x < xlim[0]):
                    if y > ylim[1]:
                        ax.annotate( '', xy=(x, 0.92*ylim[1]), xytext=(0, 24), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                    elif y < ylim[0]:
                        ax.annotate( '', xy=(x, 0.08*ylim[1]), xytext=(0, -24), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
            else:
                for yp, xp in zip(y, x):
                    if not (yp > ylim[1] or yp < ylim[0]) and (xp > xlim[1] or xp < xlim[0]):
                        if xp > xlim[1]:
                            ax.annotate( '', xy=(0.92*xlim[1], yp), xytext=(24, 0), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                        elif xp < xlim[0]:
                            ax.annotate( '', xy=(0.08*xlim[1], yp), xytext=(-24, 0), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                    elif (yp > ylim[1] or yp < ylim[0]) and not (xp > xlim[1] or xp < xlim[0]):
                        if yp > ylim[1]:
                            ax.annotate( '', xy=(xp, 0.92*ylim[1]), xytext=(0, 24), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
                        elif yp < ylim[0]:
                            ax.annotate( '', xy=(xp, 0.08*ylim[1]), xytext=(0, -24), textcoords='offset points', arrowprops=dict(arrowstyle="<|-", color=sc.get_facecolors()[0].tolist()) )
            
                        
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    if legendonframe:
        ax.legend()
    else:
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1), frameon=False)
    fig.savefig(outpath, bbox_inches='tight')
    print('Saving {}\n'.format(outpath))

