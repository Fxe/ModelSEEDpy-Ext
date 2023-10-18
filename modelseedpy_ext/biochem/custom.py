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
        'rxnAQDS': NewModelTemplateReaction('rxnAQDS_c', 'rxn00000', 'AQDS', 'ETC', 0, 1000, gapfill_direction='>'),
        'rxnMatMephen': NewModelTemplateReaction('rxnMatMephen_c', 'rxn00000', 'ANME Matrix', 'ETC', 0, 1000, gapfill_direction='>')
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