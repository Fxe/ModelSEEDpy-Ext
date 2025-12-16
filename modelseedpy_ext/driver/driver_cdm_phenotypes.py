from scipy.stats import pearsonr, spearmanr
import pandas as pd
import matplotlib.pyplot as plt


def corr(df):
    x = df["d_ani"]
    y = df["d_phenotype"]
    # Pearson: linear correlation
    pearson_r, pearson_p = pearsonr(x, y)
    # Spearman: monotonic correlation (less sensitive to outliers)
    spearman_r, spearman_p = spearmanr(x, y)
    return {
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
    }


def build_g_pairs(dist_dict):
    import networkx as nx
    G = nx.Graph()
    for (g1, g2), d in dist_dict.items():
        G.add_edge(g1, g2, weight=d)
    matching = nx.algorithms.matching.min_weight_matching(G)
    pairs = list(matching)
    return pairs


def get_pheno_fit(_dist, _pairs, _od):
    corr_array = {}
    for phenotype in _od:
        data = []
        for g1, g2 in _pairs:
            p = (g1, g2)
            if p in _dist:
                v = _dist[p]
            else:
                v = _dist[(g2, g1)]
            d_pheno = abs(_od[phenotype][g1] - _od[phenotype][g2])
            data.append({
                'g1': g1,
                'g2': g2,
                'd_ani': v,
                'd_phenotype': d_pheno
            })
        df = pd.DataFrame(data)
        corr_array[phenotype] = corr(df)
    return corr_array


def plt_corr(df):
    # Extract Pearson p-values and r
    pearson_p = df.loc['pearson_p']
    pearson_r = df.loc['pearson_r']

    # Plot p-values
    plt.figure(figsize=(12, 10))
    bars = plt.bar(pearson_p.index, pearson_p.values, color='skyblue')
    plt.yscale('log')
    plt.ylabel('p-value')
    plt.title('Pearson Correlation')

    plt.xticks(rotation=90)

    # Annotate bars with correlation coefficient
    for bar, r in zip(bars, pearson_r):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.1,
                 f'r={r:.2f}', ha='center', va='bottom')

    # Draw horizontal lines for p-value thresholds
    plt.axhline(0.05, color='red', linestyle='--', label='p=0.05')

    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()
