import os
import hjson
from utils import writejson, submit_job, rootToNumpy
import numpy as np
from datetime import datetime
import sys



def make_config(config, rebin_config, outdir, binning, data):
    """
    Create config with new generated binning 
    """

    rebin_config['binning'] = binning
    truth_reweight = [ d.replace('mcdata', 'truth') for d in data if 'mcdata_rw' in d ]
    rebin_config['histo_1D'] = data + truth_reweight

    config['rebin'] = '{}/rebin_config.json'.format(outdir)
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
    
    os.makedirs('{}/'.format(outdir), exist_ok=True)
    writejson('{}/config.json'.format(outdir), config)
    writejson('{}/rebin_config.json'.format(outdir), rebin_config)




def main(args):
    ## Make output directory
    timestamp = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
    outdir = '{}/{}/'.format(args.output, timestamp)
    print('\nCreating directory: {}\n'.format(outdir))
            
    binnings_dict = {}   # dictionary of binnings and their names

    ## Get number of differential bins
    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    if rebin_config['differential_measurement']:
        Ndiff = len(rebin_config['diff_bins'])
    else:
        Ndiff = 1
    nbin_input = rebin_config['nbin_input']

    ## Get resolution as a function of unfolded observable
    resolution = hjson.load(open('config/resolution/resolution_mc16a_em_{}.json'.format(args.variable), 'r'))

    ## Loop over tested numbers of bins
    for ibin in args.nbins:

        # Generate all possible binnings
        ibin = int(ibin)
        if args.binning_phi:
            # for Phi, two bins below 0, the rest above
            n_bin_edges = ibin - 2 - 1       # number of bin edges to generate
            n_poss = int(nbin_input / 2 - 1)  # number of possible positions for bin edges
            ncomb_up = np.math.factorial( n_poss ) / (np.math.factorial(n_bin_edges) * np.math.factorial(n_poss-n_bin_edges)) # number of possible binnings
            ncomb_down = np.math.factorial(19) / np.math.factorial(18)
            ncomb = ncomb_up * ncomb_down
        else:
            n_bin_edges = int(ibin / 2 - 1)   # number of bin edges to generate
            n_poss = int(nbin_input / 2 - 1)  # number of possible positions for bin edges
            ncomb = np.math.factorial( n_poss ) / (np.math.factorial(n_bin_edges) * np.math.factorial(n_poss-n_bin_edges)) # number of possible binnings
        n_discarded = 0                   # count number of binnings discarded, to keep continuous numbering scheme
        binnings = []
        binnings_filtered = []
        
        while len(binnings) < ncomb:
            if len(binnings)%1000 == 0:
                print(len(binnings))
            ## generate random binning until all possibilities are covered
            if args.binning_phi:
                bin_edges_up = np.random.randint(1,20, n_bin_edges)/20
                bin_edges_down = np.random.randint(-19,0)/20
                bin_edges = np.append(bin_edges_up, bin_edges_down)
                bin_edges = np.append(bin_edges,[1, -1, 0])
            else:
                bin_edges = np.array([np.random.randint(1,20)/20 for i in range(n_bin_edges)])
                bin_edges = np.append(bin_edges,[1])
                bin_edges = np.append(bin_edges, bin_edges*-1)
                bin_edges = np.append(bin_edges, [0])
            
            if np.unique(bin_edges).shape == bin_edges.shape:
                binning = list(np.sort(bin_edges))
                if binning in binnings:
                    continue
                else:
                    binnings.append(binning)

                    # Check that bin width is not smaller than resolution (-15% to release the constraint)
                    bin_center = [ (binning[i]+binning[i+1])/2 for i in range(len(binning)-1) ]
                    bin_width = [ binning[i+1]-binning[i] for i in range(len(binning)-1) ]
                    res = []
                    for bc in bin_center:
                        for ires in range(len(resolution[0])-1):
                            if bc > resolution[0][ires] and bc <= resolution[0][ires+1]:
                                res.append( resolution[1][ires] )
                    if np.all(np.array(res)-0.15*np.array(res)<np.array(bin_width)):
                        binnings_filtered.append(binning)
            else:
                continue

        for igen, fbinning in enumerate(binnings_filtered):
            binning_name = '{}bin{}'.format(ibin, igen)
            binnings_dict[binning_name] = fbinning
            print('\n{}'.format(fbinning))

            # Use same binning for all differential bins
            fbinning = [ fbinning for i in range(Ndiff) ]
            #fbinning = [fbinning, [-1,1]]

            
            ## generate config
            outdir_tmp = '{}/config/{}/'.format(outdir, binning_name)
            make_config(config, rebin_config, outdir_tmp, fbinning, args.data)

            
            ## Launch job on batch system: rebin         
            outdir_input_tmp = '{}/unfoldingInput/{}/'.format(outdir, binning_name)
            os.makedirs(outdir_input_tmp, exist_ok=True)
            command = 'python scripts/rebin.py -i {} -o {} --config {}/config.json -v {} -c {} --campaign {}'.format(args.inpath, outdir_input_tmp, outdir_tmp, args.variable, ' '.join(args.channel), args.campaign)
            dep_input = submit_job('rebin', binning_name, outdir, command)


            ## Launch job on batch system: unfolding
            outdir_output_tmp = '{}/results/{}/'.format(outdir, binning_name)
            os.makedirs(outdir_output_tmp, exist_ok=True)
            dep_unfolding = []
            for d in args.data:
                command = 'python scripts/runTRExFitter.py -i {} -o {} --config {}/config.json -v {} -c {} --campaign {} -d {}\n'.format(outdir_input_tmp, outdir_output_tmp, outdir_tmp, args.variable, ' '.join(args.channel), args.campaign, d)
                command += 'rm -rf {0}/{1}_{2}_{3}/{4}/Histograms/ {0}/{1}_{2}_{3}/{4}/Systematics/ {0}/{1}_{2}_{3}/{4}/RooStats/'.format(outdir_output_tmp, args.campaign, '_'.join(args.channel), args.variable, d)
                dep_unfolding.append( submit_job('unfold_'+d, binning_name, outdir, command, dependency=dep_input) )

                
            ## Launch job on batch system: linearity test
            command = 'python scripts/linearity_test.py -i {} --truthdir {} -r {} -v {} -c {} --campaign {} --config {}/config.json\n'.format(outdir_output_tmp, outdir_input_tmp, ' '.join(args.data), args.variable, ' '.join(args.channel), args.campaign, outdir_tmp)
            submit_job('linearity', binning_name, outdir, command, dependency=':'.join(dep_unfolding))


        writejson('{}/config/binnings.json'.format(outdir), binnings_dict)




if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter
	
    parser = ArgumentParser('Generate binnings for binning optimisation', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--inpath', help='Path to input histograms', required=True)
    parser.add_argument('-o','--output', default='binning_optimisation/', help='Path to output')
    parser.add_argument('-v', '--variable', type=str, help='Spin observable to run', required=True)
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels to run', required=True)
    parser.add_argument('--campaign', type=str, help='MC campaign to run', required=True)
    parser.add_argument('-d','--data', nargs='*', help='list of reweighted asimov data', default=['mcdata', 'mcdata_rw0', 'mcdata_rw1', 'mcdata_rw2', 'mcdata_rw3', 'mcdata_rw4'])
    parser.add_argument('-g', '--config', type=str, default='config/config.json', help='Path to config file')
    parser.add_argument('-n','--nbins', nargs='*', help='Tested number of bins', default=[4, 6, 8, 10])
    parser.add_argument('--binning-phi', action='store_true', help='Generates binning for phi observable (non-symmetric)')
    args = parser.parse_args()

    main(args)
