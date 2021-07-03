import numpy as np
import os, hjson, json
from utils import writejson, submit_job, getNPsFromFitResults, getGammaNPLabel
from datetime import datetime




def makeConfig(outdir, config_baseline, NP_group, NP_dict):
    config = hjson.load(open(config_baseline))

    fixed_NPs = [ '{}:{:3f}'.format(NP, NP_dict[NP][0]) for NP in NP_group ]        
        
    config['fit']['FixNPs'] = ','.join(fixed_NPs)
    config['job']['DebugLevel'] = 0

    os.makedirs('{}/config/'.format(outdir), exist_ok=True)
    writejson('{}/config/config.json'.format(outdir), config)



    
def makeSystGroups(NP_dict, breakdown_config):
    syst_groups_pattern = hjson.load( open(breakdown_config, 'r') )['syst_groups']
    syst_groups = {}
    for group, syst_patterns in syst_groups_pattern.items():
        systs = []
        for pat in syst_patterns:
            systs += [ NP for NP in NP_dict.keys() if pat in NP ]
        syst_groups[group] = systs

    return syst_groups




def writeSystGroups(args, syst_groups, labels):
    outpath = '{}/breakdown/'.format(args.unconstrained_output)
    os.makedirs(outpath, exist_ok=True)
    with open('{}/NP_groups.tex'.format(outpath), 'w') as f:
        f.write( '\\begin{tabular}{cc}\n' )
        f.write( '\\toprule\n' )
        f.write( 'Uncertainty source & Nuisance Parameters\\\\\n' )
        for group, np_list in syst_groups.items():
            f.write( '\\hline\n' )
            for i, NP in enumerate(np_list):
                if 'alpha' in NP:
                    np_name = labels['NP_names'][NP.replace('alpha_', '')]
                elif 'gamma' in NP:
                    np_name = getGammaNPLabel(NP)
                if i==0:
                    f.write( '\\multirow{{{}}}{{*}}{{{}}} & {}\\\\\n'.format(len(np_list), labels['group_breakdown_names'][group], np_name) )
                else:
                    f.write( '& {}\\\\\n'.format(np_name) )
        f.write( '\\bottomrule\n' )
        f.write( '\\end{tabular}' )
    print('Saving {}/NP_groups.tex'.format(outpath))
    



def makeBreakdown(args):
    NP_groups = hjson.load( open('{}/NP_groups.json'.format(args.inpath), 'r') )
    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    labels = hjson.load(open(config['labels'], 'r'))
            
    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])
        diff_bins = rebin_config['diff_bins']
    else:
        diff_bins = ['']

    print(file_name)

    if args.measurement == 'param':
        nom_name = '{}/{}_{}_{}/{}/{}_unfolded_parameter.json'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)
        nom_res = hjson.load( open(nom_name) )

        syst_list_up = []
        syst_list_down = []
        stat_list = []

        for ibin, db in enumerate(diff_bins):
            uncert_tot_up = nom_res['uncert_up'][ibin]
            uncert_tot_down = nom_res['uncert_down'][ibin]
            
            breakdown_up = []
            breakdown_down = []
            group_label = []
            uncert_up_tmp = uncert_tot_up
            uncert_down_tmp = uncert_tot_down
            for group in NP_groups:
                # Get uncertainty group label
                group_label.append( labels['group_breakdown_names'][group] )
                # Get breakdown
                res = hjson.load( open('{}/{}/results/{}_{}_{}/{}/{}_unfolded_parameter.json'.format(args.inpath, group, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)) )
                uncert_up = res['uncert_up'][ibin]
                uncert_down = res['uncert_down'][ibin]
                breakdown_up.append( np.sqrt(uncert_up_tmp**2 - uncert_up**2) )
                breakdown_down.append( np.sqrt(uncert_down_tmp**2 - uncert_down**2)  )
                uncert_up_tmp, uncert_down_tmp = uncert_up, uncert_down
        

            syst_list_up.append(breakdown_up)
            syst_list_down.append(breakdown_down)
            stat_list.append(uncert_up_tmp)

            # Save breakdown to .tex file
            outpath = '{}/breakdown/'.format(args.unconstrained_output)
            os.makedirs(outpath, exist_ok=True)
            with open('{}/breakdown_{}_{}_{}_{}.tex'.format(outpath, args.campaign, '_'.join(args.channel), file_name.replace('diff', '_'.join(db)), args.data), 'w') as f:
                f.write( '\\begin{tabular}{cc}\n' )
                f.write( '\\toprule\n' )
                f.write( 'Uncertainty source & $\\Delta$ {}\\\\\n'.format(labels['parameters'][args.variable]) )
                f.write( '\\hline\n' )
                for group, bdu, bdd in zip(group_label, breakdown_up, breakdown_down):
                    f.write( '{} & +{:.3f} / -{:.3f}\\\\\n'.format(group, bdu, bdd) )
                f.write( '\\hline\n' )
                f.write( 'Total systematic uncert. & +{:.3f} / -{:.3f}\\\\\n'.format( np.sqrt(uncert_tot_up**2-uncert_up_tmp**2), np.sqrt(uncert_tot_down**2-uncert_down_tmp**2) ) )
                f.write( '{} & $\\pm${:.3f}\\\\\n'.format(labels['group_breakdown_names']['stat'], uncert_up_tmp) )
                f.write( 'Total & +{:.3f} / -{:.3f}\\\\\n'.format(uncert_tot_up, uncert_tot_down) )
                f.write( '\\bottomrule\n' )
                f.write( '\\end{tabular}' )
                print('Saving {}/breakdown_{}_{}_{}_{}.tex'.format(outpath, args.campaign, '_'.join(args.channel), file_name.replace('diff', '_'.join(db)), args.data))

        # Transpose to get diff bins on columns and systematic groups on rows
        syst_list_up = np.array(syst_list_up).T
        syst_list_down = np.array(syst_list_down).T

        # Add syst/stat breakdown to .json result file
        nom_res['syst_groups'] = NP_groups
        nom_res['uncert_syst_up'] = syst_list_up
        nom_res['uncert_syst_down'] = syst_list_down
        nom_res['uncert_stat'] = stat_list
        writejson(nom_name, nom_res)
        print('\nThe systematics/statistics breakdown has been added to {}\n'.format(nom_name))


    elif args.measurement == 'distrib':
        nom_name = '{}/{}_{}_{}/{}/{}_unfolded_distribution.json'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)
        nom_res = hjson.load( open(nom_name) )
        uncert_tot_up = np.array(nom_res['uncert_up'])
        uncert_tot_down = np.array(nom_res['uncert_down'])
        
        breakdown_up = []
        breakdown_down = []
        uncert_up_tmp = uncert_tot_up
        uncert_down_tmp = uncert_tot_down
        # Get breakdown
        for group in NP_groups:
            res = hjson.load( open('{}/{}/results/{}_{}_{}/{}/{}_unfolded_distribution.json'.format(args.inpath, group, args.campaign, '_'.join(args.channel), file_name, args.data, args.variable)) )
            uncert_up = np.array(res['uncert_up'])
            uncert_down = np.array(res['uncert_down'])
            breakdown_up.append( np.sqrt(uncert_up_tmp**2 - uncert_up**2) )
            breakdown_down.append( np.sqrt(uncert_down_tmp**2 - uncert_down**2) )
            uncert_up_tmp, uncert_down_tmp = uncert_up, uncert_down
        
        # Add syst/stat breakdown to .json result file
        nom_res['syst_groups'] = NP_groups
        nom_res['uncert_syst_up'] = breakdown_up
        nom_res['uncert_syst_down'] = breakdown_down
        nom_res['uncert_stat'] = uncert_up
        writejson(nom_name, nom_res)
        print('\nThe systematics/statistics breakdown has been added to {}\n'.format(nom_name))

            




def runBreakdown(args):
    timestamp = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
    outdir = '{}/{}/'.format(args.output, timestamp)

    print('Determining uncertainties breakdown...')
    print('Will produce corresponding inputs and outputs under directory {}'.format(outdir))

    config = hjson.load(open(args.config, 'r'))
    rebin_config = hjson.load(open(config['rebin'], 'r'))
    file_name = args.variable
    if rebin_config['differential_measurement']:
        file_name = '{}_{}_diff'.format(args.variable, rebin_config['diff_variable'])

    # get post-fit value of NPs
    NP_dict = getNPsFromFitResults('{}/{}_{}_{}/{}/Fits/{}.root'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), file_name, args.data, args.data))
    # make groups of systematics
    syst_groups = makeSystGroups(NP_dict, args.breakdown_config)
    syst_groups = { g:s for g,s in syst_groups.items() if not len(s)==0 }

    labels = hjson.load(open(config['labels'], 'r'))
    writeSystGroups(args, syst_groups, labels) # print NP groups to latex file

    dep_unfolding = []
    syst_groups_tmp = []

    for group, systs in syst_groups.items():
        outdir_tmp = '{}/{}/'.format(outdir, group)        
        # Create config file with NPs frozen to +1sigma or -1sigma
        syst_groups_tmp += systs
        makeConfig(outdir_tmp, args.config, syst_groups_tmp, NP_dict)

        # Run unfolding on batch
        outdir_output_tmp = '{}/results/'.format(outdir_tmp)
        os.makedirs(outdir_output_tmp, exist_ok=True)
        command = 'python scripts/runTRExFitter.py -i {} -o {} --config {}/config/config.json -v {} -c {} --campaign {} -d {} -m {}\n'.format(args.inpath, outdir_output_tmp, outdir_tmp, args.variable,
                                                                                                                                              ' '.join(args.channel), args.campaign, args.data, args.measurement)
        command += 'rm -rf {0}/{1}_{2}_{3}/{4}/Histograms/ {0}/{1}_{2}_{3}/{4}/Systematics/ {0}/{1}_{2}_{3}/{4}/RooStats/'.format(outdir_output_tmp, args.campaign, '_'.join(args.channel), args.variable, args.data)
        dep_unfolding.append( submit_job('unfold', group, outdir, command) )

    ## Launch job to compute breakdown
    command = 'python scripts/makeUncertaintyBreakdown.py -i {} -v {} -c {} --campaign {} -d {} -u {} -m {} -g {} --dobreakdown \n'.format( outdir, args.variable, ' '.join(args.channel), args.campaign,
                                                                                                                                            args.data, args.unconstrained_output, args.measurement, args.config )
    submit_job('breakdown', 'uncert', outdir, command, dependency=':'.join(dep_unfolding))

    ## Write down list of NP groups
    writejson( '{}/NP_groups.json'.format(outdir), list(syst_groups) )
    print('\nWriting systematic groups to {}/NP_groups.json'.format(outdir))
    print('\nThe summary table will be saved to {}/breakdown/breakdown_{}_{}_{}_{}.tex\n'.format(args.unconstrained_output, args.campaign, '_'.join(args.channel), args.variable, args.data))
        




if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description='Determine breakdown of uncertainties for measured spin parameter.', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--inpath', help='Path to input histograms', required=True)
    parser.add_argument('-o','--output', default='breakdown/', help='Path to output')
    parser.add_argument('-v', '--variable', type=str, help='Spin observable to run', required=True)
    parser.add_argument('-c', '--channel', nargs='*', help='Dileptonic channels to run', required=True)
    parser.add_argument('--campaign', type=str, help='MC campaign to run', required=True)
    parser.add_argument ('-d','--data', type=str, help='dataset to run on', default='mcdata')
    parser.add_argument('-g', '--config', type=str, default='config/config.json', help='Path to config file')
    parser.add_argument('-u', '--unconstrained-output', type=str, default='results/base/', help='Path to output of unconstrained fit')
    parser.add_argument('-b', '--breakdown-config', type=str, default='config/breakdown.json', help='Path to config file defining groups of systematics')
    parser.add_argument('-m', '--measurement', choices=['param', 'distrib'], default='param', help='o breakdown for the unfolded spin parameter of for the unfolded bin.')

    parser.add_argument('--dobreakdown', action='store_true', help='Perform uncertainty breakdown.')

    args = parser.parse_args()
    
    if not args.dobreakdown:
        runBreakdown(args)
    else:
        makeBreakdown(args)
