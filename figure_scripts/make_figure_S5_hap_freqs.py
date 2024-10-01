'''
This program somewhat manually generates a figure S5 for the manuscript, using
outputs from graph_haplotype_freqs.py
'''
import pandas as pd
import plotly.express as px

amplicons=['heome-a rep1', 'heome-a rep2', 'heome-a rep1',
'heome-a rep2', 'dhps-540 rep1', 'dhps-540 rep2', 'dhps-540 rep1',
'dhps-540 rep2']

fracs=[0.994, 0.990, 0.006, 0.010, 0.948, 0.958, 0.052, 0.042]

point_names=['expected', 'expected', 'unexpected_hap1', 'unexpected_hap1',
'expected', 'expected', 'unexpected_hap1', 'unexpected_hap1']

data={'amplicons': amplicons, 'fracs':fracs, 'point_names':point_names}
df = pd.DataFrame(data)

fig1=px.bar(df, x='amplicons', y=fracs, color=point_names, text_auto=True)
fig1['layout']['xaxis']['title']='amplicon replicates (all at 1 parasite/ul'
fig1['layout']['yaxis']['title']='fraction of total reads per haplotype'

fig1.write_html('manuscript_figures/S5_hap_freqs.html')
fig1.write_image('manuscript_figures/S5_hap_freqs.svg')
fig1.write_image('manuscript_figures/S5_hap_freqs.png')