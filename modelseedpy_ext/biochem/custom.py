from modelseedpy import MSMedia


def swap_template_compounds(t, swap_data):
    for _from_id in swap_data:
        _to_id = swap_data[_from_id]
        if _from_id in t.compcompounds and _to_id in t.compcompounds:
            _from = t.compcompounds.get_by_id(_from_id)
            _to = t.compcompounds.get_by_id(_to_id)
            for rxn in t.reactions:
                m = rxn.metabolites
                d = {}
                if _from in m:
                    v = m[_from]
                    d[_to] = v
                    d[_from] = -1 * v
                    rxn.add_metabolites(d)


def get_charge_transfer(r):
    cmp_charge = {}
    lhs = {}
    rhs = {}
    for m, v in r.metabolites.items():
        c = m.compartment
        if c not in cmp_charge:
            cmp_charge[c] = 0
        if v < 0:
            cmp_charge[c] += v * m.charge
            lhs[m] = v
        elif v > 0:
            cmp_charge[c] += v * m.charge
            rhs[m] = v
        else:
            raise ValueError('zero stoich: ' + str(r))
    for k in {k for k, v in cmp_charge.items() if v == 0}:
        del cmp_charge[k]
    return cmp_charge


def fix_template(t):
    swap = {
        'cpd27742_c': 'cpd00649_c',  # f420
        'cpd28067_c': 'cpd00792_c',  # f420
    }
    swap_template_compounds(t, swap)
    if 'rxn43076_c' in t.reactions:
        if len(t.reactions.rxn43076_c.check_mass_balance()) > 0:
            t.reactions.rxn43076_c.add_metabolites({
                t.compcompounds.cpd00067_c: 3
            })
            print('fix rxn43076_c', t.reactions.rxn43076_c.check_mass_balance())

    # fix charge balance
    if 'rxn05939_c' in t.reactions:
        mb = t.reactions.rxn05939_c.check_mass_balance()
        if len(mb) > 0:
            metabolites = t.reactions.rxn05939_c.metabolites
            t.reactions.rxn05939_c.add_metabolites({
                t.compcompounds.cpd11620_c: -1,
                t.compcompounds.cpd11621_c: 1
            })

    for rxn_id in REACTION_RENAME:
        if rxn_id in t.reactions:
            rxn = t.reactions.get_by_id(rxn_id)
            rxn.name = REACTION_RENAME[rxn_id]

    for cpd_id, name in COMPOUND_RENAME.items():
        t.compounds.get_by_id(cpd_id).name = name


def _hack_template(t):
    if 'rxn40355_c' in t.reactions:
        mb = t.reactions.rxn40355_c.check_mass_balance()
        if len(mb) > 0:
            t.reactions.rxn40355_c.add_metabolites({
                t.compcompounds.cpd00067_c: -1
            })
    if 'rxn00929_c' in t.reactions:
        t.reactions.rxn00929_c.upper_bound = 0
    if 'rxn00931_c' in t.reactions:
        t.reactions.rxn00931_c.upper_bound = 0
    if 'rxn09188_c' in t.reactions:
        t.reactions.rxn09188_c.upper_bound = 0
    if 'rxn00085_c' in t.reactions:
        t.reactions.rxn00085_c.upper_bound = 0
    if 'rxn00146_c' in t.reactions:
        t.reactions.rxn00146_c.lower_bound = 0
    if 'rxn45744_c' in t.reactions:
        t.reactions.rxn45744_c.upper_bound = 0
    if 'rxn39860_c' in t.reactions:
        t.reactions.rxn39860_c.lower_bound = 0
    if 'rxn05595_c' in t.reactions:
        t.reactions.rxn05595_c.lower_bound = 0
    if 'rxn05596_c' in t.reactions:
        t.reactions.rxn05596_c.upper_bound = 0
    if 'rxn43657_c' in t.reactions:
        t.reactions.rxn43657_c.upper_bound = 0
    if 'rxn09388_c' in t.reactions:
        t.reactions.rxn09388_c.upper_bound = 0


def add_tests(model):
    from cobra.core import Reaction
    tests = []
    metabolites_in_model = {m.id for m in model.metabolites}
    if 'test_atp_c0' not in model.reactions:
        rxn_test_atp = Reaction(f'test_atp_c0', f'Test ATP', 'TEST', 0, 0)
        rxn_test_atp.add_metabolites({
            model.metabolites.cpd00001_c0: -1,
            model.metabolites.cpd00002_c0: -1,
            model.metabolites.cpd00008_c0: 1,
            model.metabolites.cpd00009_c0: 1,
            model.metabolites.cpd00067_c0: 1,
        })
        tests.append(rxn_test_atp)
    if 'test_nadp_c0' not in model.reactions:
        need = {'cpd00006_c0', 'cpd00005_c0', 'cpd00067_e0'}
        if need & metabolites_in_model == need:
            test = Reaction(f'test_nadp_c0', f'Test NADP', 'TEST', 0, 0)
            test.add_metabolites({
                model.metabolites.cpd00006_c0: -1,
                model.metabolites.cpd00005_c0: 1,
                model.metabolites.cpd00067_e0: -1,
            })
            tests.append(test)
    if 'test_nad_c0' not in model.reactions:
        need = {'cpd00003_c0', 'cpd00004_c0', 'cpd00067_e0'}
        if need & metabolites_in_model == need:
            test = Reaction(f'test_nad_c0', f'Test NAD', 'TEST', 0, 0)
            test.add_metabolites({
                model.metabolites.cpd00003_c0: -1,
                model.metabolites.cpd00004_c0: 1,
                model.metabolites.cpd00067_e0: -1,
            })
            tests.append(test)
    model.add_reactions(tests)
    return tests


def add_anme_data_to_template(t):
    from cobrakbase.core.kbasefba.newmodeltemplate_metabolite import NewModelTemplateCompound, NewModelTemplateCompCompound
    from cobrakbase.core.kbasefba.newmodeltemplate_reaction import NewModelTemplateReaction
    template_anme_compounds = [
        NewModelTemplateCompound('cpdAQDS', 'C14H6O8S2', 'AQDS', 0),
        NewModelTemplateCompound('cpdAQH2DS', 'C14H8O8S2', 'AQH2DS', 0),
        NewModelTemplateCompound('cpdETCMe', 'R', 'ETC matrix (charged)', -1),
        NewModelTemplateCompound('cpdETCM', 'R', 'ETC matrix (uncharged)', 0)
    ]
    template_anme_ccompounds = [
        NewModelTemplateCompCompound('cpdAQDS_e', 0, 'e', 'cpdAQDS'),
        NewModelTemplateCompCompound('cpdAQH2DS_e', 0, 'e', 'cpdAQH2DS'),
        NewModelTemplateCompCompound('cpdETCMe_e', -1, 'e', 'cpdETCMe'),
        NewModelTemplateCompCompound('cpdETCM_e', 0, 'e', 'cpdETCM')
    ]
    t.add_compounds(template_anme_compounds)
    t.add_comp_compounds(template_anme_ccompounds)

    template_anme_reactions = {
        'rxnAQDS': NewModelTemplateReaction('rxnAQDS_c', 'rxn00000',
                                            'AQDS', 'ETC', 0, 1000, gapfill_direction='>'),
        'rxnMatMephen': NewModelTemplateReaction('rxnMatMephen_c', 'rxn00000',
                                                 'ANME Matrix', 'ETC', 0, 1000, gapfill_direction='>')
    }
    template_anme_reactions['rxnAQDS'].add_metabolites({
        t.compcompounds.cpd08702_c: -1,
        t.compcompounds.cpdAQDS_e: -1,
        t.compcompounds.cpd08701_c: 1,
        t.compcompounds.cpdAQH2DS_e: 1,
    })
    template_anme_reactions['rxnMatMephen'].add_metabolites({
        t.compcompounds.cpd08702_c: -1,
        t.compcompounds.cpd08701_c: 1,
        t.compcompounds.cpd00067_e: 2,
        t.compcompounds.cpdETCMe_e: 2,
        t.compcompounds.cpdETCM_e: -2,
    })
    # t.add_reactions([template_archaea.reactions.rxn40355_c])
    t.add_reactions(list(template_anme_reactions.values()))
    return t


media_co2_sulfate_matrix = MSMedia.from_dict({
    'cpd00149': (-100, 100),
    'cpd00099': (-100, 100),
    'cpd00067': (-100, 100),
    'cpd00063': (-100, 100),
    'cpd00058': (-100, 100),
    'cpd00048': (-10, 100),
    'cpd00034': (-100, 100),
    'cpd00030': (-100, 100),
    'cpd00013': (-100, 100),
    'cpd00009': (-100, 100),
    'cpd00001': (-100, 100),
    'cpd11640': (-100, 100),
    'cpd00205': (-100, 100),
    'cpd00254': (-100, 100),
    'cpd00971': (-100, 100),
    'cpd10515': (-100, 100),
    'cpd10516': (-100, 100),
    'cpd11574': (-100, 100),
    'cpd00244': (-100, 100),
    'cpd00528': (-100, 100),  # N2
    'cpd00011': (-100, 100),  # CO2
    # 'cpd00276': (-5, 100), # glucosamine
    #    'cpd00082': (-5, 100), # fructose
    #    'cpd00023': (-5, 100), # glutamate
    #    'cpd00794': (-5, 100), # trhalose
    #    'cpd00029': (-5, 100), # acetate

    'cpdETCMe': (-20, 10),  # matrix
})


REACTION_ANME_BOUNDS = {
    # 4 H+ [c] + 2 F420 [c] + CoM [c] + HTP [c] + 2 Reducedferredoxin [c] -->
    # 2 F420H2 [c] + CoM-S-S-CoB [c] + 2 Oxidizedferredoxin [c]
    'rxn43076_c': (0, 1000),
    # NAD [c] + CoA [c] + Pyruvate [c] <=> NADH [c] + CO2 [c] + Acetyl-CoA [c]
    'rxn00154_c': (-1000, 1000),
    # F420H2 [c] + MP [c] --> H+ [c] + F420 [c] + MPH2 [c]
    'rxn24614_c': (0, 1000),
}


MEDIUM_GENOME_SCALE_METHANE_AQDS = {
    'EX_cpd00067_e0': 1000.0,
    'EX_cpd00001_e0': 1000.0,
    'EX_cpd01024_e0': 10.0,
    'EX_cpdAQDS_e0': 1000.0,

    'EX_cpd00149_e0': 1000.0,
    'EX_cpd00099_e0': 1000.0,
    'EX_cpd00063_e0': 1000.0,
    'EX_cpd00058_e0': 1000.0,
    'EX_cpd00048_e0': 1000.0,

    'EX_cpd00034_e0': 1000.0,
    'EX_cpd00030_e0': 1000.0,

    'EX_cpd00009_e0': 1000.0,
    'EX_cpd11640_e0': 1000.0,
    'EX_cpd10515_e0': 1000.0,
    'EX_cpd10516_e0': 1000.0,

    'EX_cpd00205_e0': 1000.0,

    'EX_cpd00254_e0': 1000.0,
    'EX_cpd00971_e0': 1000.0,

    'EX_cpd11574_e0': 1000.0,
    'EX_cpd00244_e0': 1000.0,

    'EX_cpd00013_e0': 1000.0,  # nh4
    'EX_cpd00528_e0': 1000.0,  # n2
}

MEDIUM_ETC_METHANE_AQDS = {
    'EX_cpd01024_e0': 1,
    'EX_cpd00001_e0': 1000,
    'EX_cpd00067_e0': 1000,
    'EX_cpdAQDS_e0': 1000,
}

KO_ROLES = {
    'KO:K18209': 'fumarate reductase (CoM/CoB) subunit A [EC:1.3.4.1]',
    'KO:K18210': 'fumarate reductase (CoM/CoB) subunit B [EC:1.3.4.1]',
    'KO:K22171': 'F420H2:quinone oxidoreductase subunit A [EC:1.1.98.4]',
    'KO:K22172': 'F420H2:quinone oxidoreductase subunit B/C [EC:1.1.98.4]',
    'KO:K22173': 'F420H2:quinone oxidoreductase subunit D [EC:1.1.98.4]',
    'KO:K22174': 'F420H2:quinone oxidoreductase subunit F [EC:1.1.98.4]',
    'KO:K22175': 'F420H2:quinone oxidoreductase subunit H [EC:1.1.98.4]',
    'KO:K22176': 'F420H2:quinone oxidoreductase subunit I [EC:1.1.98.4]',
    'KO:K22177': 'F420H2:quinone oxidoreductase subunit J [EC:1.1.98.4]',
    'KO:K22178': 'F420H2:quinone oxidoreductase subunit K [EC:1.1.98.4]',
    'KO:K22179': 'F420H2:quinone oxidoreductase subunit L [EC:1.1.98.4]',
    'KO:K22180': 'F420H2:quinone oxidoreductase subunit M [EC:1.1.98.4]',
    'KO:K22181': 'F420H2:quinone oxidoreductase subunit N [EC:1.1.98.4]',
}

COMPOUND_RENAME = {
    'cpd00792': 'F420H2',
    'cpd00649': 'F420',
    'cpd08701': 'MP',
    'cpd08702': 'MPH2',
}

REACTION_RENAME = {
    'rxn03127_c': 'Mcr',
    'rxn03020_c': 'Mtr',
    'rxn03085_c': 'Mer',
    'rxn03079_c': 'Mtd',
    'rxn02480_c': 'Mch',
    'rxn03126_c': 'HdrDE',
    'rxn02431_c': 'Ftr',
    'rxn40505_c': 'cdh',
    'rxn40355_c': 'Fpo',
    'rxn05759_c': 'Frh',
    'rxn43076_c': 'HdrA',
    'rxn46184_c': 'Fwd/Fmd',

    'rxn15962_c': 'CODH/ACS',

    'rxa45615_c': 'Rnf',
    'rxn08173_c': 'ATPS',
    'rxn05469_c': 'Pyruvate Transport',
    'rxn05209_c': 'Na+/H+ Antiport',
    'rxn05467_c': 'CO2 Diffusion',
    'rxn05319_c': 'H2O Diffusion',

    'rxn13974_c': 'PFOR',

    'rxn00250_c': 'PC',
    'rxn00248_c': 'MDH NAD',
    'rxn00249_c': 'MDH NADP',
    'rxn00799_c': 'FH',
    'rxa09272_c': 'TrfAB',
    'rxn00285_c': 'SucCD',
    'rxn05939_c': 'KorABCD',
    'rxn14048_c': 'KorABCD',
    'rxn00974_c': 'acn',
    'rxn01388_c': 'acn',
    'rxn00199_c': 'icd',
    'rxn01387_c': 'icd',
    'rxn00256_c': 'ACLY',
}
