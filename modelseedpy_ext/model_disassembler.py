from modelseedpy.core.mstemplate import (
    MSTemplateReaction,
    MSTemplateMetabolite,
    MSTemplateSpecies,
    MSTemplateCompartment,
)


def to_template_compound(cpd, compartment):
    cpd_id = cpd.seed_id
    ccpd_id = f"{cpd.seed_id}_{compartment}"
    template_cpd = MSTemplateMetabolite(cpd_id, cpd.formula, cpd.name, cpd.charge)
    template_cpd.abbreviation = cpd.abbr
    template_ccpd = MSTemplateSpecies(ccpd_id, cpd.charge, compartment, template_cpd.id)
    return template_cpd, template_ccpd


class ModelDisassembler:
    def __init__(self, modelseed, template):
        self.modelseed = modelseed
        self.template = template

    def disassemble_metabolite(self, cpd, modelseed_id):
        print("disassemble_metabolite", cpd.id, modelseed_id)
        if modelseed_id in self.modelseed.compounds:
            seed_compound = self.modelseed.compounds.get_by_id(modelseed_id)
            return to_template_compound(seed_compound, cpd.compartment)
        else:
            print(modelseed_id, "not in modelseed database")
            return None, None

    def disassemble_reaction(self, rxn, modelseed_id):
        template_reaction = MSTemplateReaction(rxn.id, modelseed_id, rxn.name)
        template_reaction.lower_bound = rxn.lower_bound
        template_reaction.upper_bound = rxn.upper_bound

        metabolites = {}
        for m, v in rxn.metabolites.items():
            if m.id not in self.template.compcompounds:
                cpd, ccpd = self.disassemble_metabolite(m, m.id[:-2])
                # print(cpd.get_data())
                # print(ccpd.get_data())
                if cpd.id not in self.template.compounds:
                    self.template.add_compounds([cpd])
                if ccpd.id not in self.template.compcompounds:
                    self.template.add_comp_compounds([ccpd])

            template_compound = self.template.compcompounds.get_by_id(m.id)
            metabolites[template_compound] = v

        template_reaction.add_metabolites(metabolites)
        return template_reaction
