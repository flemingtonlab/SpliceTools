import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns


path = sys.argv[1]

df = pd.read_table(path, index_col=0)

c = sns.clustermap(data=df, cmap='Blues', figsize=(30,30), z_score=None, \
    yticklabels=1, xticklabels=1, robust=True, method="ward", linewidth=.1)

# methods: ward, median, centroid, weighted, average, single

''' 
metric: ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’, ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘jensenshannon’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’
'''

length = len(df.index)
fontsize = 500 // length

plt.setp(c.ax_heatmap.xaxis.get_majorticklabels(),rotation=90)
plt.setp(c.ax_heatmap.yaxis.get_majorticklabels(),rotation=0)

c.ax_heatmap.set_xticklabels(c.ax_heatmap.get_xticklabels(), fontsize=fontsize, fontname='Arial')
c.ax_heatmap.set_yticklabels(c.ax_heatmap.get_yticklabels(), fontsize=fontsize, fontname='Arial')


row_index = c.dendrogram_row.reordered_ind
labels = [df.index[i] for i in row_index]
plt.tight_layout()
plt.savefig(path +'.svg')


with open(path + '.label_order', 'w') as outfile:
    for i in labels:
        outfile.write(f'{i}\n')   
