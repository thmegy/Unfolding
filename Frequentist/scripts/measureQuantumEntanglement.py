import os
import hjson
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg') # load Agg backed which does not require X11 connection to make plots




def plotSummary(args):
    var = 'DeltaPhi'

    mean = []
    uncert_up = []
    uncert_down = []
    sig = []

    cuts = [float(c) for c in args.list_cut]

    for cut in args.list_cut:
        with open('{}/{}/{}_{}_{}_Mttbar_diff/{}/DeltaPhi_unfolded_parameter.json'.format(args.output, cut, args.campaign, '_'.join(args.channel), var, args.data), 'r') as f:
            res = hjson.load(f)
            mean.append( res['mean'][0] )
            uncert_up.append( res['uncert_up'][0] )
            uncert_down.append( res['uncert_down'][0] )
            sig.append( res['significance'] )

    fig, axes = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [1,3]})
    
    # plot significance
    axes[0].plot(cuts, sig, 'ro', linestyle=':')
    axes[0].grid(axis='y', linestyle='dashed')
    axes[0].set_ylabel('$\\sqrt{q_0}$')
    axes[0].set_ylim([0.8*min(sig), 1.1*max(sig)])

    # plot D vs mttbar cut
    axes[1].errorbar(cuts, mean, yerr=[uncert_down, uncert_up], color='black', fmt='o', label='Unfolded parameter')
    axes[1].axhline(y=-1/3, linestyle='dashed', color='red', linewidth=1)
    axes[1].set_xlabel('$m_{t\\bar{t}}$ threshold [GeV]')
    axes[1].set_ylabel('$D$')

    # find minimum to set y-range
    low_array = np.array(mean) - np.array(uncert_down)
    low_range = 1.1 * np.min(low_array)
    axes[1].set_ylim([low_range, -0.25])

    axes[1].set_yticks([-1/3], minor=True)
    axes[1].set_yticklabels(labels=['-1/3'], color='red', fontsize=12, minor=True)

    axes[1].legend()

    fig.savefig('{}/quantum_entanglement/summary_{}.pdf'.format(args.output, args.data), bbox_inches='tight')
    print('Saving {}/quantum_entanglement/summary_{}.pdf\n'.format(args.output, args.data))

    # also write results in latex table
    with open('{}/quantum_entanglement/summary_{}.tex'.format(args.output, args.data), 'w') as f:
        f.write( '\\begin{tabular}{cc}\n' )
        f.write( '\\toprule\n' )
        f.write( '$m_{t\\bar{t}}$ cut & D\\\\\n' )
        for l, m, up, do in zip(args.list_cut, mean, uncert_up, uncert_down):
            f.write( '\\hline\n' )
            f.write( f'{l} & {m:.3f}$^{{+{up:.3f}}}_{{-{do:.3f}}}$\\\\\n' )
        f.write( '\\bottomrule\n' )
        f.write( '\\end{tabular}' )
    print('Saving {}/quantum_entanglement/summary_{}.tex'.format(args.output, args.data))





def main(args):
    var = 'DeltaPhi'
    os.makedirs(f'{args.output}/quantum_entanglement/', exist_ok=True)

    for cut in args.list_cut:
        command = 'python scripts/rebin.py -i {}/{} -o {}/{} --config {} -v {} -c {} --campaign {} -s {}'.format(args.inpath, cut, args.output_rebin, cut, args.config, var, ' '.join(args.channel), args.campaign, args.syst)
        os.system(command)
        command = 'python scripts/runTRExFitter.py -i {}/{} -o {}/{} --config {} -v {} -c {} --campaign {} -d {}\n'.format(args.output_rebin, cut, args.output, cut, args.config, var, ' '.join(args.channel), args.campaign, args.data)
        os.system(command)

    plotSummary(args)




if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter
	
    parser = ArgumentParser('Measure quantum entanglement and its significance', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--inpath', help='Path to input histograms', required=True)
    parser.add_argument('-o','--output', default='results/base/', help='Path to output')
    parser.add_argument('-or','--output-rebin', default='unfoldingInput/base/', help='Path to rebin output')
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels to run', required=True)
    parser.add_argument('--campaign', type=str, help='MC campaign to run', default='all')
    parser.add_argument('-d', '--data', type=str, default='mcdata', help='Dataset to run on')
    parser.add_argument('-g', '--config', type=str, default='config/DeltaPhi/config.json', help='Path to config file')
    parser.add_argument('-l', '--list-cut', nargs='*', required=True, help='List of Mttbar cuts used')
    parser.add_argument('-p', '--plot', action='store_true', help='plot summary')
    parser.add_argument('-s', '--syst', choices=['full', 'shape'], help='Shape-only systematics or shape+syst (full)', default='full')
    args = parser.parse_args()

    if not args.plot:
        main(args)
    else:
        plotSummary(args)
