from mOTUlizer.classes.MetaBin import MetaBin
from math import log10

mean = lambda x : sum(x)/len(x)

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