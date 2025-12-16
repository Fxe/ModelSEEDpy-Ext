import cobra
from cobra.core import Reaction


class BuildModel:

    @staticmethod
    def fix_template(template_neg):
        #template_neg = cobrakbase.io.load_kbase_zip_object('./GramNegModelTemplateV6_test.json.zip')
        from modelseedpy.core.mstemplate import MSTemplateReaction, MSTemplateSpecies
        template_neg.add_comp_compounds([
            MSTemplateSpecies('cpd25960_e', template_neg.compcompounds.cpd25960_c.charge, 'e', 'cpd25960')
        ])
        rxs43120_c0 = MSTemplateReaction('rxs43120_c', 'rxs43120', 'rxs43120', 'MePn Degradation', 0, 1000,
                                         gapfill_direction='>')
        rxs43120_c0.add_metabolites({
            template_neg.compcompounds.cpd00017_c: -1,
            template_neg.compcompounds.cpd26519_c: -1,
            template_neg.compcompounds.cpd01024_e: 1,
            template_neg.compcompounds.cpd03091_c: 1,
            template_neg.compcompounds.cpd00060_c: 1,
            template_neg.compcompounds.cpd26508_c: 1,

            template_neg.compcompounds.cpd00003_c: +1,
            template_neg.compcompounds.cpd00004_c: -1,
        })
        print(rxs43120_c0.check_mass_balance())
        rxs00001_c0 = MSTemplateReaction('rxs00001_c', 'rxs00001', 'MePn ABC Transporter', 'MePn Degradation', 0, 1000,
                                         gapfill_direction='>')
        rxs00001_c0.add_metabolites({
            template_neg.compcompounds.cpd00001_c: -1,
            template_neg.compcompounds.cpd00002_c: -1,
            template_neg.compcompounds.cpd00008_c: 1,
            template_neg.compcompounds.cpd00009_c: 1,
            template_neg.compcompounds.cpd00067_c: 1,

            template_neg.compcompounds.cpd25960_e: -1,
            template_neg.compcompounds.cpd25960_c: 1,
        })
        print(rxs00001_c0.check_mass_balance())
        template_neg.add_reactions([rxs43120_c0, rxs00001_c0])
        template_neg.reactions.rxn26447_c.add_metabolites({
            template_neg.compcompounds.cpd00067_c: -2
        })
        print(template_neg.reactions.rxn26447_c.check_mass_balance())
        disable = ['rxn23850_c', 'rxn09193_c', 'rxn05315_c', 'rxn05206_c']
        _remove = []
        for rxn_id in disable:
            if rxn_id in template_neg.reactions:
                _remove.append(rxn_id)
        template_neg.reactions -= _remove
        disable_lb = [
            'rxn09188_c', 'rxn00929_c',
            'rxn00145_c', 'rxn00499_c', 'rxn00146_c',
            'rxn08934_c', 'rxn08941_c', 'rxn08935_c', 'rxn08942_c',  # not sure about these
            'rxn11942_c', 'rxn11941_c',
            'rxn43329_c', 'rxn11940_c',
            'rxn40559_c', 'rxn04159_c', 'rxn04705_c',
            'rxn14415_c',
            'sul00003_c', 'rxn48579_c',
            'rxn39860_c'
        ]
        for rxn_id in disable_lb:
            template_reaction = template_neg.reactions.get_by_id(rxn_id)
            template_reaction.lower_bound = 0

    @staticmethod
    def _add_atpm(model):
        if 'ATPM_c0' not in model.reactions:
            atpm = Reaction(f'ATPM_c0', f'ATPM', 'ATPM', 0, 1000)
            atpm.add_metabolites({
                model.metabolites.cpd00001_c0: -1,
                model.metabolites.cpd00002_c0: -1,
                model.metabolites.cpd00008_c0: 1,
                model.metabolites.cpd00009_c0: 1,
                model.metabolites.cpd00067_c0: 1,
            })
            model.add_reactions([atpm])

    @staticmethod
    def conv_atp_results(atp_correction, tests):
        results = {'tests': {}, 'gapfill': {}}
        for x in tests:
            results['tests'][x['media'].id] = {
                'is_max_threshold': x['is_max_threshold'],
                'threshold': x['threshold'],
                'objective': x['objective'],
            }
        for x in atp_correction.media_gapfill_stats:
            gp_stats = atp_correction.media_gapfill_stats[x]
            # print(x.id, gp_stats)
            if gp_stats:
                results['gapfill'][x.id] = {
                    'reversed': gp_stats.get('reversed'),
                    'new': gp_stats.get('new'),
                    'target': gp_stats.get('target'),
                    'minobjective': gp_stats.get('minobjective'),
                    'binary_check': gp_stats.get('binary_check'),
                }
        return results

    @staticmethod
    def atp_test(model):
        disable = ['rxn23850_c0', 'rxn09193_c0', 'rxn05315_c0', 'rxn05206_c0']
        for rxn_id in disable:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = 0
                rxn.upper_bound = 0
        disable_lb = [
            'rxn09188_c0', 'rxn00929_c0',
            'rxn00145_c0', 'rxn00499_c0', 'rxn00146_c0',
            'rxn08934_c0', 'rxn08941_c0', 'rxn08935_c0', 'rxn08942_c0',  # not sure about these
            'rxn11942_c0', 'rxn11941_c0',
            'rxn43329_c0', 'rxn11940_c0',
            'rxn40559_c0', 'rxn04159_c0', 'rxn04705_c0',
            'rxn14415_c0',
            'sul00003_c0', 'rxn48579_c0',
            'rxn39860_c0'
        ]
        for rxn_id in disable_lb:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = 0
        rxn_balance_check = ['rxn15395_c0', 'rxn28338_c0']
        for rxn_id in rxn_balance_check:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                check = {'charge': -1, 'H': -1}
                if rxn.check_mass_balance() == check:
                    print(rxn)
                    rxn.add_metabolites({
                        model.metabolites.cpd00067_c0: +1
                    })
                    print(rxn.check_mass_balance())
        model.medium = {
            'EX_cpd00001_e0': 1000,
            'EX_cpd00067_e0': 1000,
            # 'EX_cpd00036_e0': 1,
            'EX_cpd00007_e0': 100,
        }
        model.objective = 'ATPM_c0'
        sol = cobra.flux_analysis.pfba(model)
        if sol.fluxes['ATPM_c0'] > 0.0:
            return False, sol
        return True, sol