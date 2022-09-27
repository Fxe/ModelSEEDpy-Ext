import hashlib
import pandas as pd


def sha_hex(s: str):
    return hashlib.sha256(s.encode('utf-8')).hexdigest()


def test_core_template(template, media_file):
    medias = load_medias(media_file)
    model_core = build_full_model(template)
    test_core_model(model_core, medias)


def load_medias(filep):
    medias = pd.read_csv(filep, sep='\t', index_col=0).to_dict()
    return medias


def build_full_model(template):
    from modelseedpy import MSBuilder
    from cobra.core import Model, Reaction
    model = Model('modelo')
    model.add_reactions([r.to_reaction() for r in template.reactions if len(r.complexes) > 0 or not r.type == 'conditional'])
    MSBuilder.build_exchanges(model, 'e0')
    #model.add_reactions()
    bio_atp = Reaction('bio2', 'bio2', '', -10000, 10000)
    bio_atp.add_metabolites({
        model.metabolites.cpd00067_c0: 1,
        model.metabolites.cpd00002_c0: -1,
        model.metabolites.cpd00008_c0: 1,
        model.metabolites.cpd00001_c0: -1,
        model.metabolites.cpd00009_c0: 1,
    })
    model.add_reactions([bio_atp])
    model.objective = 'bio2'
    return model


def test_media(model, media, exs, comp):
    medium_ = {}
    for k in media:
        v = media[k]
        if v > 0:
            ex = k
            if ex in exs:
                medium_[ex] = v
            else:
                return None
    medium_.update(comp)
    model.medium = medium_
    sol = model.optimize()
    return sol


def test_core_model(model_core, medias, v=True):
    res = {}
    media_out = {}
    exs = set(map(lambda x: x.id, model_core.exchanges))
    comp = {'EX_cpd00001_e0': 1000, 'EX_cpd00067_e0': 1000}
    for media_id in medias:
        sol = test_media(model_core, medias[media_id], exs, comp)
        media_out[media_id] = None
        if sol and sol.status == 'optimal':
            
            ex_out = dict(filter(lambda x: x[0].startswith('EX_') and x[1] != 0, sol.fluxes.to_dict().items()))
            result = sol.objective_value
            for k in medias[media_id]:
                if medias[media_id][k] > 0:
                    if sol.fluxes[k] == 0:
                        result = 0
            res[media_id] = result
            if res[media_id] > 0:
                media_out[media_id] = ex_out
            if v:
                print(media_id, sol.objective_value)
        else:
            if v:
                print(media_id, 0)
            res[media_id] = 0
    return res, media_out


def xml_print(l_start, l_end, filename, line=True):
    l_read = 0
    with open(filename, 'rb') as fh:
        while l_read < l_end:
            l = fh.readline()
            if l_read >= l_start:
                if line:
                    print(f'[{l_read}]{l.decode("utf-8")[:-1]}')
                else:
                    print(f'{l.decode("utf-8")[:-1]}')
            l_read += 1
