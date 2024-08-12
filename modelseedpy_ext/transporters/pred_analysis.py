import pandas as pd


def _count(d):
    res = {}
    for k in d:
        res[k] = len(d[k])
    return res


def _d_to_str(d):
    x = []
    for k, v in d.items():
        x.append(f'{k}:{v}')
    return ';'.join(x)


class PredAnalysis:

    def __init__(self, re_seq_dna_to_ko, dna_h_to_pred_feature, ko_info, ko_cat):
        """

        ko_cat ABC -> ['K10000', 'K20000']

        :param re_seq_dna_to_ko:
        :param dna_h_to_pred_feature:
        :param ko_info:
        :param ko_cat: dict[str] -> [str]
        """
        self.ko_info = ko_info
        self.ko_cat = ko_cat
        self.re_seq_dna_to_ko = re_seq_dna_to_ko
        self.dna_h_to_pred_feature = dna_h_to_pred_feature
        self.found_hashes = set(self.re_seq_dna_to_ko) & set(self.dna_h_to_pred_feature)
        print(len(self.found_hashes))

        self.ko_cat_in_genome_set = self.get_ko_cat_in_genome_set()

    def get_ko_cat_in_genome_set(self):
        ko_cat_in_genome_set = {}
        for cat in self.ko_cat:
            ko_cat_in_genome_set[cat] = {}
        for h in self.found_hashes:
            kos = {x.split(';')[-1] for x in self.re_seq_dna_to_ko[h]}
            if len(kos) == 1:
                _ko = list(kos)[0]
                for cat in self.ko_cat:
                    if _ko[3:] in self.ko_cat[cat]:
                        if _ko not in ko_cat_in_genome_set[cat]:
                            ko_cat_in_genome_set[cat][_ko] = set()
                        ko_cat_in_genome_set[cat][_ko].add(h)
        return ko_cat_in_genome_set

    def _get_family(self, cat, ko, fam_pred):
        kos_fam = {}
        for h in self.ko_cat_in_genome_set[cat][ko]:
            for feature_id in self.dna_h_to_pred_feature[h]:
                fam = fam_pred[feature_id]
                if fam not in kos_fam:
                    kos_fam[fam] = set()
                kos_fam[fam].add(feature_id)
        return kos_fam

    def get_df_bin_pred_scores(self):
        data = []
        for cat in self.ko_cat_in_genome_set:
            # print(cat, len(ko_cat_in_genome_set[cat]))
            for _ko in self.ko_cat_in_genome_set[cat]:
                scores = self._get_scores(cat, _ko)
                row = [cat, _ko, self.ko_info[_ko[3:]]['name'], len(scores[90]), len(scores[80]), len(scores[70]),
                       len(scores[60]), len(scores[50]), len(scores[0])]
                data.append(row)
        df = pd.DataFrame(data, columns=['cat', 'ko', 'name', '>=0.9', '>=0.8', '>=0.7', '>=0.6', '>=0.5', '>=0'])
        return df

    def _get_scores(self, cat, ko):
        kos_score = {
            90: set(),
            80: set(),
            70: set(),
            60: set(),
            50: set(),
            0: set()
        }
        for h in self.ko_cat_in_genome_set[cat][ko]:
            for feature_id in self.dna_h_to_pred_feature[h]:
                score = self.feature_id_to_bin_pred[feature_id]
                # print(feature_id, score, kos)
                if score >= 0.90:
                    kos_score[90].add(feature_id)
                elif score >= 0.80:
                    kos_score[80].add(feature_id)
                elif score >= 0.70:
                    kos_score[70].add(feature_id)
                elif score >= 0.60:
                    kos_score[60].add(feature_id)
                elif score >= 0.50:
                    kos_score[50].add(feature_id)
                else:
                    kos_score[0].add(feature_id)
        return kos_score

    def _get_ko_tc(self, ko_data):
        tc = None
        db_links = ko_data.get('db_links')
        if db_links:
            for s in db_links:
                if s.strip().startswith('TC:'):
                    tc = s
                    tc = tc.split(':')[1].strip()
                    break
        return tc

    def get_df_new_dataset_pred(self):
        data = []
        for cat in self.ko_cat_in_genome_set:
            # print(cat, len(ko_cat_in_genome_set[cat]))
            for _ko in self.ko_cat_in_genome_set[cat]:
                ko_data = self.ko_info[_ko[3:]]
                fam_pred_count = _count(self._get_family(cat, _ko, self.feature_id_to_pred_newdataset))
                row = [cat, _ko, ko_data['name'], _get_ko_tc(ko_data), _d_to_str(fam_pred_count)]
                data.append(row)
        df = pd.DataFrame(data, columns=['cat', 'ko', 'name', 'tc', 'newdataset'])
        return df

    def get_df_subfamily_pred(self):
        data = []
        for cat in self.ko_cat_in_genome_set:
            # print(cat, len(ko_cat_in_genome_set[cat]))
            for _ko in self.ko_cat_in_genome_set[cat]:
                ko_data = self.ko_info[_ko[3:]]
                fam_pred_count = []
                for k in [15, 20, 30, 50]:
                    fam_pred_count.append(_count(self._get_family(cat, _ko, self.feature_id_to_pred_subfamily[k])))
                row = [cat, _ko, ko_data['name'], _get_ko_tc(ko_data)]
                row += [_d_to_str(x) for x in fam_pred_count]
                data.append(row)
        df = pd.DataFrame(data, columns=['cat', 'ko', 'name', 'tc', 15, 20, 30, 50])
        return df

    def aaa(self):
        c = 0
        kos_score = {
            90: {},
            80: {},
            70: {},
            60: {},
            50: {},
            0: {}
        }
        for h in self.found_hashes:
            kos = {x.split(';')[-1] for x in self.re_seq_dna_to_ko[h]}
            if len(kos) == 1:
                _ko = list(kos)[0]
                for feature_id in self.dna_h_to_pred_feature[h]:
                    score = self.feature_id_to_bin_pred[feature_id]
                    #print(feature_id, score, kos)
                    if score > 0.90:
                        if _ko not in kos_score[90]:
                            kos_score[90][_ko] = set()
                        kos_score[90][_ko].add(feature_id)
                    elif score > 0.80:
                        if _ko not in kos_score[80]:
                            kos_score[80][_ko] = set()
                        kos_score[80][_ko].add(feature_id)
                    elif score > 0.70:
                        if _ko not in kos_score[70]:
                            kos_score[70][_ko] = set()
                        kos_score[70][_ko].add(feature_id)
                    elif score > 0.60:
                        if _ko not in kos_score[60]:
                            kos_score[60][_ko] = set()
                        kos_score[60][_ko].add(feature_id)
                    elif score > 0.50:
                        if _ko not in kos_score[50]:
                            kos_score[50][_ko] = set()
                        kos_score[50][_ko].add(feature_id)
                    else:
                        if _ko not in kos_score[0]:
                            kos_score[0][_ko] = set()
                        kos_score[0][_ko].add(feature_id)
            #print(dna_h_to_pred_feature[h], )
            c += 1
            #if c > 1000:
            #    break
        return kos_score