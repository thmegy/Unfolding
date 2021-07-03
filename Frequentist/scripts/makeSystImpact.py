import numpy as np
import os, hjson, json
from utils import writejson, submit_job, getNPsFromFitResults, Axisoffset, getDiffBinLabel, getGammaNPLabel
from datetime import datetime
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker




def makeConfig(outdir, config_baseline, NP, NP_dict, is_up=True, is_postfit=True):
    
    config = hjson.load(open(config_baseline))

    if is_postfit:
        if is_up:
            fixed_NP = NP_dict[NP][0]+NP_dict[NP][1]
        else:
            fixed_NP = NP_dict[NP][0]+NP_dict[NP][2]
    else:
        if is_up:
            fixed_NP = 1
        else:
            fixed_NP = -1
        
    config['fit']['FixNPs'] = '{}:{:3f}'.format(NP, fixed_NP)
    config['job']['DebugLevel'] = 0
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



    
def rankImpact(args):
    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    labels_config = hjson.load(open(config['labels'], 'r'))

    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])
        diff_bins = rebin_config['diff_bins']
    else:
        diff_bins = ['']

    # Perform ranking based on postfit impacts
    NP_dict = getNPsFromFitResults('{}/{}_{}_{}/{}/Fits/{}.root'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), file_name, args.data, args.data))
    print('{}/{}_{}_{}/{}/Fits/{}.root'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), file_name, args.data, args.data))

    nom_list = hjson.load( open('{}/{}_{}_{}/{}/{}_unfolded_parameter.json'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)) )['mean']
        
    # Loop over differential bins
    for ibin, db in enumerate(diff_bins):

        try:
            nom = nom_list[ibin]
        except Exception as e:
            print(f'\nResults for differential bin {ibin} not found!')
            print(e)
            continue

        ## For each systematics, retrieve the +1 and -1 sigma impact and store them
        impact_post_plus = []
        impact_post_minus = []
        impact_pre_plus = []
        impact_pre_minus = []
        impact = []
        NP_list = []
        for NP in NP_dict.keys():
            try:
                sigplus_post = hjson.load( open('{}/postfit/{}_up/results/{}_{}_{}/{}/{}_unfolded_parameter.json'.format(args.inpath, NP, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)) )['mean'][ibin]
                sigminus_post = hjson.load( open('{}/postfit/{}_down/results/{}_{}_{}/{}/{}_unfolded_parameter.json'.format(args.inpath, NP, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)) )['mean'][ibin]

                sigplus_pre = hjson.load( open('{}/prefit/{}_up/results/{}_{}_{}/{}/{}_unfolded_parameter.json'.format(args.inpath, NP, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)) )['mean'][ibin]
                sigminus_pre = hjson.load( open('{}/prefit/{}_down/results/{}_{}_{}/{}/{}_unfolded_parameter.json'.format(args.inpath, NP, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)) )['mean'][ibin]
            except Exception as e:
                print(e)
                continue
            impact_post_plus.append( (sigplus_post-nom)/nom )
            impact_post_minus.append( (sigminus_post-nom)/nom )
            impact.append( abs(sigplus_post-sigminus_post) )
            impact_pre_plus.append( (sigplus_pre-nom)/nom )
            impact_pre_minus.append( (sigminus_pre-nom)/nom )
            NP_list.append(NP)

        ## The ranking is happening here
        impact = np.array(impact)
        sorted_args = impact.argsort()

        impact_plus_sel = []
        impact_minus_sel = []
        systs_sel = []
        n = 1

        if len(impact) < args.nrank:
            nrank = len(impact)
        else:
            nrank = args.nrank

        while n <= nrank:
            arg = sorted_args[-n]
            systs_sel.append(NP_list[arg])
            if 'gamma' not in NP_list[arg]:
                impact_plus_sel.append( (impact_post_plus[arg], impact_pre_plus[arg]) )
                impact_minus_sel.append( (impact_post_minus[arg], impact_pre_minus[arg]) )
            else:
                # no pre-fit for gamma NPs
                impact_plus_sel.append( (impact_post_plus[arg], impact_post_plus[arg]) )
                impact_minus_sel.append( (impact_post_minus[arg], impact_post_minus[arg]) )
            n += 1

        # Get pulls and constraints corresponding to the selected systematics
        pull = { NP:NP_dict[NP][0] for NP in NP_dict }
        constraint_plus = { NP:NP_dict[NP][1] for NP in NP_dict }
        constraint_minus = { NP:NP_dict[NP][2] for NP in NP_dict }
        pull_sel = [ pull[i] for i in systs_sel ]
        constraint_plus_sel = [ constraint_plus[i] for i in systs_sel ]
        constraint_minus_sel = [ abs(constraint_minus[i]) for i in systs_sel ]

        # Save full ranking to .tex file
        outpath = '{}/syst_impact/'.format(args.unconstrained_output)
        os.makedirs(outpath, exist_ok=True)
        file_name = args.variable
        if db is not '':
            file_name = '{}_{}_{}'.format(args.variable, rebin_config['diff_variable'], '_'.join(db))

        with open('{}/syst_impact_ranking_{}_{}_{}_{}.tex'.format(outpath, args.campaign, '_'.join(args.channel), file_name, args.data), 'w') as f:

            f.write( '\\begin{tabular}{cccc}\n' )
            f.write( '\\toprule\n' )
            f.write( 'Rank & Uncertainty source & Impact & Constraint\\\\\n' )
            f.write( '\\hline\n' )

            for n in range(len(sorted_args)):
                arg = sorted_args[-(n+1)]
                NP = NP_list[arg]
                if 'alpha' in NP:
                    label = labels_config['NP_names'][NP.replace('alpha_', '')]
                else:
                    label = getGammaNPLabel(NP)
                f.write( '{} &{} & {:.3f} / {:.3f} & +{:.3f} / {:.3f}\\\\\n'.format(n+1, label, impact_post_plus[arg], impact_post_minus[arg], constraint_plus[NP], constraint_minus[NP] ) )
            f.write( '\\bottomrule\n' )
            f.write( '\\end{tabular}' )
            print( 'Saving {}/syst_impact_ranking_{}_{}_{}_{}.tex'.format(outpath, args.campaign, '_'.join(args.channel), file_name, args.data) )


        # Plot only post-fit impacts
        plotImpact(args, config, rebin_config, labels_config, file_name, db, systs_sel, np.array(impact_plus_sel), np.array(impact_minus_sel), True)        
        # Plot pre- and post-fit impacts
        plotImpact(args, config, rebin_config, labels_config, file_name, db, systs_sel, np.array(impact_plus_sel), np.array(impact_minus_sel), False,
                   pull=pull_sel, constraint=[constraint_minus_sel, constraint_plus_sel])





def plotImpact(args, config, rebin_config, labels_config, file_name, diff_bin, systs, impact_plus, impact_minus, onlypost, pull=None, constraint=None):
    # Retrieve NP labels
    y_pos = np.arange(len(systs))
    y_labels = [ labels_config['NP_names'][syst.replace('alpha_', '')] if 'alpha' in syst else getGammaNPLabel(syst) for syst in systs ]

    fig, ax = plt.subplots(figsize=(4,len(systs)*0.75))

    if onlypost:
        ax.barh(y_pos, impact_plus[:,0], color='#ea9293', height=0.8, label=r'$+1\sigma$ impact')
        ax.barh(y_pos, impact_minus[:,0], color='#8ebad9', height=0.8, label=r'$-1\sigma$ impact')
        ax.grid(which='major', axis='x', ls='--', c='k')
    else:
        ax.barh(y_pos, impact_plus[:,0], color='#ea9293', height=0.8, label=r'Post-fit $+1\sigma$ impact')
        ax.barh(y_pos, impact_minus[:,0], color='#8ebad9', height=0.8, label=r'Post-fit $-1\sigma$ impact')
        ax.barh(y_pos, impact_plus[:,1], color='#00000000', edgecolor='#d62728', height=0.8, label=r'Pre-fit $+1\sigma$ impact')
        ax.barh(y_pos, impact_minus[:,1], color='#00000000', edgecolor='#1f77b4', height=0.8, label=r'Pre-fit $-1\sigma$ impact')

    if pull is not None and constraint is not None:
        # Add NP pulls and constraints
        ax2 = ax.twiny()
        ax2.errorbar(pull, y_pos, xerr=constraint, color='k', marker='o', capsize=0, markersize=4, linestyle='', label='NP pull')
        ax2.set_xlim(-2, 2)
        ax2.set_xlabel(r'$\hat{\theta}$')
        ax2.xaxis.set_label_position('bottom')
        ax2.xaxis.tick_bottom()
        ax2.grid(which='major', axis='x', ls='--', c='k')

    xmax = np.max(np.abs(impact_plus))
    if np.max(np.abs(impact_minus)) > xmax:
        xmax = np.max(np.abs(impact_minus))
    ax.set_xlim(-xmax*1.2, xmax*1.2)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-2,0))
    ax.invert_yaxis()
    ax.set_xlabel('$\\Delta$ {0} / {0}'.format(labels_config['parameters'][args.variable]))
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    formatter = ticker.ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((-4, 4))
    ax.xaxis.set_major_formatter(formatter)
    axoffset = Axisoffset(ax, 'x', (1.0, 1.0))
    ax.margins(y=0.001)

    if pull is not None and constraint is not None:
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines_aligned = lines + lines2
        labels_aligned = labels + labels2
        ax.legend(lines_aligned, labels_aligned, loc='lower right', bbox_to_anchor=(-0.05, 1.1), frameon=False)
    else:
        ax.legend(loc='lower right', bbox_to_anchor=(-0.05, 1.1), frameon=False)
    xmax = ax.get_xlim()
    ax.text(xmax[1], -2.5, r'ATLAS Internal', horizontalalignment='right', verticalalignment='bottom', fontsize=12)
    if 'all' in args.campaign:
        lumi_campaign = args.campaign[:3]
    else:
        lumi_campaign = args.campaign[:5]
    ax.text(xmax[1], -2, f'$\sqrt{{s}} = 13$ TeV, {config["lumi"][lumi_campaign][0]} fb$^{{-1}}$', horizontalalignment='right', verticalalignment='bottom', fontsize=12)

    os.makedirs('{}/syst_impact/'.format(args.unconstrained_output), exist_ok=True)
    plot_name = '{}/syst_impact/syst_impact_ranking_{}_{}_{}_{}'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), file_name, args.data)
    if onlypost:
        plot_name += '_postfit.pdf'
    else:
        plot_name += '_full.pdf'
    fig.savefig(plot_name, bbox_inches='tight')
    print('Saving {}'.format(plot_name))




def runImpact(args):

    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])
    
    # get post-fit value of NPs
    NP_dict = getNPsFromFitResults('{}/{}_{}_{}/{}/Fits/{}.root'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), file_name, args.data, args.data))
    
    timestamp = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
    outdir = '{}/{}/'.format(args.output, timestamp)
    os.makedirs('{}/prefit/'.format(outdir), exist_ok=True)
    os.makedirs('{}/postfit/'.format(outdir), exist_ok=True)

    print('Determining systematics impact...')
    print('Will produce corresponding inputs and outputs under directory {}'.format(outdir))

    dep_unfolding = []

    for NP, values in NP_dict.items():
        for fitstage in ['pre', 'post']:
            for syst in ['up', 'down']:
                outdir_tmp = '{}/{}fit/{}_{}/'.format(outdir, fitstage, NP, syst)
        
                # Create config files with NP frozen to +1sigma or -1sigma
                is_postfit=True
                if fitstage == 'pre':
                    is_postfit=False
                is_up=True
                if syst == 'down':
                    is_up=False
                makeConfig(outdir_tmp, args.config, NP, NP_dict, is_postfit=is_postfit, is_up=is_up)     

                # Run unfolding on batch
                outdir_output_tmp = '{}/results/'.format(outdir_tmp)
                os.makedirs(outdir_output_tmp, exist_ok=True)

                command = 'python scripts/runTRExFitter.py -i {} -o {} --config {}/config/config.json -v {} -c {} --campaign {} -d {}\n'.format(args.inpath, outdir_output_tmp, outdir_tmp, args.variable, ' '.join(args.channel), args.campaign, args.data)
                command += 'rm -rf {0}/{1}_{2}_{3}/{4}/Histograms/ {0}/{1}_{2}_{3}/{4}/Systematics/ {0}/{1}_{2}_{3}/{4}/RooStats/'.format(outdir_output_tmp, args.campaign, '_'.join(args.channel), file_name, args.data)
                dep_unfolding.append( submit_job('unfold', '{}_{}_{}'.format(NP, syst, fitstage), outdir, command, mem='3gb') )


    ## Launch job to rank impacts and draw summary plot
    command = 'python scripts/makeSystImpact.py -i {} -v {} -c {} --campaign {} -d {} -u {} --rank-plot -n {} --config {}\n'.format(
        outdir, args.variable, ' '.join(args.channel), args.campaign, args.data, args.unconstrained_output, args.nrank, args.config )
    submit_job('rank', 'systs', outdir, command, dependency=':'.join(dep_unfolding))

    print('\nThe summary plots will be saved to {}/syst_impact/\n'.format(args.unconstrained_output))



if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description='Determine impact of systematic uncertainties on measured spin parameter.', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--inpath', help='Path to input histograms', required=True)
    parser.add_argument('-o','--output', default='syst_impact/', help='Path to output')
    parser.add_argument('-v', '--variable', type=str, help='Spin observable to run', required=True)
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels to run', required=True)
    parser.add_argument('--campaign', type=str, help='MC campaign to run', required=True)
    parser.add_argument ('-d','--data', type=str, help='dataset to run on', default='mcdata')
    parser.add_argument('-g', '--config', type=str, default='config/config.json', help='Path to config file')
    parser.add_argument('-u', '--unconstrained-output', type=str, default='results/base/', help='Path to output of unconstrained fit')

    parser.add_argument('-r', '--rank-plot', action='store_true', help='Rank impacts and plot summary.')
    parser.add_argument('-n', '--nrank', type=int, default=8, help='Plot n systematics with the largest impact. Only use with -r argument.')

    args = parser.parse_args()
    
    if not args.rank_plot:
        runImpact(args)
    else:
        rankImpact(args)
