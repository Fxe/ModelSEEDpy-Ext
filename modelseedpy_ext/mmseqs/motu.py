import logging
from mOTUlizer.classes.MetaBin import MetaBin
from math import log10


logging.getLogger(__name__)

mean = lambda x : sum(x)/len(x)


def get_core_features(cluster, feature_to_genome, genomes, cov_cut=0.85):
    f_count = {}
    core_features = None
    for f_id in cluster:
        g_id = feature_to_genome[f_id]
        genome = genomes[g_id]

        feature = genome.features.get_by_id(f_id)
        terms = feature.ontology_terms.get('RAST', [])
        for t in terms:
            if t not in f_count:
                f_count[t] = 0
            f_count[t] += 1
    l = len(cluster)
    f_per = {k: (v / l) for k, v in f_count.items()}
    for f, cov in f_per.items():
        if cov > cov_cut:
            core_features = f
    return core_features


def load_mmseq_results(mmseqs_dat, precluster, mmseqs2, name):
    with open(mmseqs_dat) as handle:
        if precluster:
            recs = {g: l[:-1].split()[0] for l in handle for g in precluster[l[:-1].split()[1]]}
        else:
            recs = {l[:-1].split()[1]: l[:-1].split()[0] for l in handle}

    fill = len(str(len(set(recs.values()))))
    rep2clust = {k: name + str(i).zfill(fill) for i, k in enumerate(set(recs.values()))}
    gene_clusters2rep = {v: k for k, v in rep2clust.items()}

    logging.debug(f'For {len(recs)} CDSes we got {len(gene_clusters2rep)}  gene-clusters')

    recs = {k: rep2clust[v] for k, v in recs.items()}
    genome2gene_clusters = {k.id: set() for k in mmseqs2.genomes}

    for k, v in recs.items():
        g = mmseqs2.feature_genome[k]
        genome2gene_clusters[g.id].add(v)

    return genome2gene_clusters, recs, gene_clusters2rep


class FastMOTU:
    
    def __init__(self, tt):
        self.gene_clusters_dict = tt['genome2gene_clusterss']
        self.aa2gene_clusters = tt['aa2gene_clusters']
        self.quiet = False
        
    def __getitem__(self, i):
        return self.members[i]
        
    def __len__(self):
        return len(self.members)
    
    def run(self, genome_completion_dict = None, max_it=20):
        if genome_completion_dict is None:
            genome_completion_dict = {k: 100 for k in self.gene_clusters_dict.keys()}
        self.members = [MetaBin(bin_name, 
                                gene_clusterss = self.gene_clusters_dict[bin_name], 
                                faas = None, #self.faas.get(bin_name), 
                                fnas = None, 
                                complet = genome_completion_dict.get(bin_name)) for bin_name in self.gene_clusters_dict.keys()]
        
        self.gene_clustersCounts = {c : 0 for c in set.union(set([gene_clusters for mag in self.members for gene_clusters in mag.gene_clusterss]))}
        for mag in self.members:
            for gene_clusters in mag.gene_clusterss:
                self.gene_clustersCounts[gene_clusters] += 1
        self.mean_comp = mean(genome_completion_dict.values())
        """
        for gene_clusters, counts in self.gene_clustersCounts.items():
            if counts > 30:
                print(gene_clusters, counts)
                print(100*counts, len(self))
                print(100*counts/len(self) > self.mean_comp)
        """
        self.core = {gene_clusters for gene_clusters, counts in self.gene_clustersCounts.items() if (100*counts/len(self)) > self.mean_comp}
        self.likelies = self.__core_likelyhood(max_it = max_it)
        
    def __core_prob(self, gene_clusters, complet = "checkm"):
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        presence = [log10(comp(mag)) for mag in self if gene_clusters in mag.gene_clusterss]
        abscence = [log10(1 - comp(mag)) for mag in self if gene_clusters not in mag.gene_clusterss]
        return sum(presence + abscence)

    def __pange_prob(self, gene_clusters, core_size, complet = "checkm"):
#        pool_size = sum(self.gene_clustersCounts.values())
        pool_size = sum([c for k,c in  self.gene_clustersCounts.items()])
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        #presence = [1 - (1-self.gene_clustersCounts[gene_clusters]/pool_size)**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self if gene_clusters in mag.gene_clusterss]
        #abscence = [ (1-self.gene_clustersCounts[gene_clusters]/pool_size)**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self if gene_clusters not in mag.gene_clusterss]

#        mag_prob = {mag : ( 1-1/pool_size )**len(mag.gene_clusterss) for mag in self}
#        mag_prob = {mag : ( 1-1/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}
#        mag_prob = {mag : ( 1-self.gene_clustersCounts[gene_clusters]/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}

#        mag_prob = {mag : ( 1-self.gene_clustersCounts[gene_clusters]/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}
        mag_prob = {mag : ( 1-self.gene_clustersCounts[gene_clusters]/pool_size)**len(mag.gene_clusterss) for mag in self}

        presence = [ log10(1 -   mag_prob[mag]) if mag_prob[mag] < 1 else MIN_PROB                for mag in self if gene_clusters in mag.gene_clusterss]
        abscence = [ log10(      mag_prob[mag]) if mag_prob[mag] > 0 else log10(1-(10**MIN_PROB)) for mag in self if gene_clusters not in mag.gene_clusterss]

        #abscence = [ 1-self.gene_clustersCounts[gene_clusters]/len(self)*comp(mag) for mag in self if gene_clusters not in mag.gene_clusterss]
        #presence = [ self.gene_clustersCounts[gene_clusters]/len(self)*comp(mag) for mag in self if gene_clusters not in mag.gene_clusterss]

        return sum(presence + abscence)

    def __core_likely(self, gene_clusters, complet = "checkm", core_size = 0):
        pange_prob = self.__pange_prob(gene_clusters, core_size, complet)
        return self.__core_prob(gene_clusters, complet) - pange_prob
        
    def __core_likelyhood(self, max_it = 20, likeli_cutof = 0 ):
        likelies = {gene_clusters : self.__core_likely(gene_clusters) for gene_clusters in self.gene_clustersCounts}
        self.core = set([c for c, v in likelies.items() if v > likeli_cutof])
        core_len = len(self.core)
        i = 1
        if not self.quiet:
            print("iteration 1 : ", core_len, "sum_abs_LLHR:" , sum([l if l > 0 else -l for l in likelies.values()]), file = sys.stderr)
        for mag in self:
            if len(self.core) > 0:
                mag.new_completness = 100*len(mag.gene_clusterss.intersection(self.core))/len(self.core)
            else :
                mag.new_completness = 0
            mag.new_completness = mag.new_completness if mag.new_completness < 99.9 else 99.9
            mag.new_completness = mag.new_completness if mag.new_completness > 0 else 0.01
        for i in range(2,max_it):
            likelies = {gene_clusters : self.__core_likely(gene_clusters, complet = "new", core_size = core_len) for gene_clusters in self.gene_clustersCounts}
            old_core = self.core
            self.core = set([c for c, v in likelies.items() if v > likeli_cutof])
            new_core_len = len(self.core)
            for mag in self:
                if len(self.core) > 0:
                    mag.new_completness = 100*len(mag.gene_clusterss.intersection(self.core))/len(self.core)
                else :
                    mag.new_completness = 0
                mag.new_completness = mag.new_completness if mag.new_completness < 99.9 else 99.9
                mag.new_completness = mag.new_completness if mag.new_completness > 0 else 0.01
            if not self.quiet:
                print("iteration",i, ": ", new_core_len, "sum_abs_LLHR:" , sum([l if l > 0 else -l for l in likelies.values()]), file = sys.stderr)
            if self.core == old_core:
                break
            else :
                core_len = new_core_len

        if not self.quiet:
            pp =  "\nYour {name}-run for {nb_mags} genomes (with mean initial completeness {mean_start:.2f}) resulted\n"
            pp += "in a core of {core_len} traits with a total sum of loglikelihood-ratios {llhr:.2f} and a corrected \n"
            pp += "mean completness of {mean_new:.2f}, resulting to a estimated mean traits per genome count of {trait_count:.2f}\n"
            pp = pp.format(name = self.name, nb_mags = len(self), core_len = core_len, mean_start = mean([b.checkm_complet for b in self]),
                        mean_new =  mean([b.new_completness for b in self]), llhr =  sum([l if l > 0 else -l for l in likelies.values()]),
                        trait_count = mean([100*len(b.gene_clusterss)/b.new_completness for b in self]))
            print(pp, file = sys.stderr)
        self.iterations = i -1
        return likelies

    def get_stats(self):
        output = {
            "nb_genomes": len(self),
            "core": list(self.core) if self.core else None,
                     "aux_genome": [k for k, v in self.gene_clustersCounts.items() if
                                    k not in self.core] if self.core else None,
                     "singleton_gene_clusterss": [k for k, v in self.gene_clustersCounts.items() if k not in self.core
                                                  if v == 1] if self.core else None,
                     "gene_clusterss": None if self.gene_clusters_dict is None else {
                         'genome': {k: list(v) for k, v in self.gene_clusters_dict.items()},
                         'aa': self.aa2gene_clusters} if self.aa2gene_clusters else (
                         {k: list(v) for k, v in self.gene_clusters_dict.items()} if self.gene_clusters_dict else None),
                     "mean_ANI": self.get_mean_ani() if (
                                 hasattr(self, 'fastani_dict') or all([hasattr(g, "genome") for g in self])) else None,
                     "ANIs": [[k[0], k[1], v] for k, v in self.fastani_matrix().items()] if (
                                 hasattr(self, 'fastani_dict') or all([hasattr(g, "genome") for g in self])) else None,
                     "genomes": [v.get_data() for v in self],
                     "likelies": self.likelies
                     }
        return output

    def get_core_features(self, genomes, m_core, ontology='RAST', cov_cut=0.85):
        cluster_to_features = {}
        for f_id, c_id in self.aa2gene_clusters.items():
            if c_id not in cluster_to_features:
                cluster_to_features[c_id] = set()
            cluster_to_features[c_id].add(f_id)
        feature_id_to_genome = {}
        for g_id, genome in genomes.items():
            for f in genome.features:
                if f.id not in feature_id_to_genome:
                    feature_id_to_genome[f.id] = g_id
                else:
                    print('!!')

        core_features = {}
        for k in m_core:
            l = len(cluster_to_features[k])
            f_count = {}
            for f_id in cluster_to_features[k]:
                g_id = feature_id_to_genome[f_id]
                genome = genomes[g_id]
                feature = genome.features.get_by_id(f_id)
                terms = feature.ontology_terms.get(ontology, [])
                for t in terms:
                    if t not in f_count:
                        f_count[t] = 0
                    f_count[t] += 1
            f_per = {k:(v/l) for k,v in f_count.items()}
            for f, cov in f_per.items():
                if cov > cov_cut:
                    if f not in core_features:
                        core_features[f] = 0
                    core_features[f] += 1

        return core_features
