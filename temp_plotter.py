# coding: utf-8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

RAI_file = r'C:\Users\Samy\Dropbox\Samy_Dropbox\MSc\lncRNA_HOTAIRM1_Project_Dostie\jdostie-rnaseq\Qiagen\primary_analyses\QIAseqUltraplexRNA_90846\secondary_analysis\DGEA\DGEA_results_with_all_RAId_samples_siNC-siM1_fdr_0.1.tsv'
all_file = r'C:\Users\Samy\Dropbox\Samy_Dropbox\MSc\lncRNA_HOTAIRM1_Project_Dostie\jdostie-rnaseq\Qiagen\primary_analyses\QIAseqUltraplexRNA_90846\secondary_analysis\DGEA\comparisons-without-conditionwise-isolation\siNC_siM1\DGEA_results_with_all_samples_siNC-siM1_fdr_0.1.tsv'
RAI_df = pd.read_csv(RAI_file, sep='\t', header=0, index_col=None)

all_df = pd.read_csv(all_file, sep='\t', header=0, index_col=None)

def over_ts(df1, df2, trange):
    df1_counts, df2_counts, overlap_counts = [],[],[]
    for t in trange:
        df1_items = df1[ df1['edgeR.exactTest_FDR'] <= t ]['external_gene_name']
        df2_items = df2[ df2['edgeR.exactTest_FDR'] <= t ]['external_gene_name']
        df1_counts.append(len(df1_items))
        df2_counts.append(len(df2_items))
        overlap_counts.append( len(set(df1_items.values.tolist()).intersection(set(df2_items.values.tolist())) ))
    return df1_counts, df2_counts, overlap_counts
    
_1, _2, _o =  over_ts(all_df, RAI_df, np.arange(0.0, 0.3, 0.05))

fig,ax=plt.subplots(ncols=2, sharey=True)
ax[0].plot(np.arange(0.0,0.3,0.05), _1, 'k-', label="inc_NT", color='red')
ax[0].plot(np.arange(0.0,0.3,0.05), _2, 'k-', label="exc_NT", color='blue')
ax[1].plot(np.arange(0.0,0.3,0.05), _o, 'k-')
ax[0].legend()


# Show the major grid lines with dark grey lines
ax[0].grid(b=True, which='major', color='#666666', linestyle='-')
ax[1].grid(b=True, which='major', color='#666666', linestyle='-')

# Show the minor grid lines with very faint and almost transparent grey lines
ax[0].minorticks_on()
ax[0].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
ax[1].minorticks_on()
ax[1].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)



plt.show()
