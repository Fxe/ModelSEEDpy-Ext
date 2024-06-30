from cobra.util import format_long_string


class ReBiochemModelSEEDCompound:

    def __init__(self, doc, re, shallow=True):
        self.re = re
        self.doc = doc
        self._inchi = None
        self._smiles = None
        self.formula = None
        if not shallow:
            self._inchi = self.inchi
            self._smiles = self.smiles

        self.ontology = None

    @property
    def inchi(self):
        if self._inchi is None:
            aql = """
            FOR e IN re_biochem_inchi_key_has_metabolite
                FILTER e._to == @k
                RETURN e._from
            """
            res = self.re.db.AQLQuery(aql, bindVars={'k': self.doc['_id']}, rawResults=True)
            fetch_ids = {o for o in res}
            aql = """
            FOR o IN re_biochem_inchi_key
                FILTER o._id in @ids
                RETURN o
            """
            res = self.re.db.AQLQuery(aql, bindVars={'ids': list(fetch_ids)}, rawResults=True)
            records = [d for d in res]
            if len(records) == 0:
                self._inchi = (None, None)
            elif len(records) == 1:
                self._inchi = (records[0]['_key'], records[0]['inchi'])
            else:
                print('multiple inchis')

        return self._inchi

    @property
    def smiles(self):
        if self._smiles is None:
            aql = """
            LET edges = (
            FOR e IN re_biochem_smiles_has_metabolite
                FILTER e._to == @k
                RETURN e._from
            )
            FOR doc IN re_biochem_smiles
                FILTER doc._id IN edges
                RETURN doc
            """
            res = self.re.db.AQLQuery(aql, bindVars={'k': self.doc['_id']}, rawResults=True)
            records = [d for d in res]
            if len(records) == 0:
                self._smiles = (None, None)
            elif len(records) == 1:
                self._smiles = (records[0]['_key'], records[0]['_value'])
            else:
                print('multiple smiles')

        return self._smiles

    @property
    def id(self):
        return self.doc['metabolite_id']

    @property
    def label(self):
        return self.doc.collection.name

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Label</strong></td><td>{label}</td>
            </tr><tr>
                <td><strong>Compound identifier</strong></td><td>{id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Formula</strong></td><td>{formula}</td>
            </tr><tr>
                <td><strong>SMILES</strong></td><td>{smiles}</td>
            </tr><tr>
                <td><strong>Inchi</strong></td><td>{inchi}</td>
            </tr><tr>
                <td><strong>Inchi Key</strong></td><td>{inchi_key}</td>
            </tr><tr>
                <td><strong>In {n_species} species</strong></td><td>
                    {species}</td>
            </tr>
        </table>""".format(
            label=self.label,
            id=self.id,
            name=format_long_string(self.id),
            formula=self.id,
            smiles=self._smiles[1] if self._smiles else '?',
            inchi=self._inchi[1] if self._inchi else '?',
            inchi_key=self._inchi[0] if self._inchi else '?',
            address="0x0%x" % id(self),
            n_species=len(self.id),
            species=format_long_string(self.id, 200),
        )
