#Title: DAOx Illumina visualizations
#Author: Alexandre Schoepfer
#Version: 12.04.2021

import numpy as np
import pandas as pd
import altair as alt
import streamlit as st

aa_column  = 'aa_mutation'
n_aa_column = 'n_aa_substitutions'
bin_column = 'bin'

def up_f(types=['txt','csv','.xlsx']):
    up_file = st.file_uploader("Upload File",type=types)
    return up_file

@st.cache()
def build_df(df,heatmap, amplification_factor):
    
    global aa_column 
    global n_aa_column 
    global bin_column

    if heatmap == 'Expression':
        
        def min_max_norm(data):
            return (data - np.min(data)) / (np.max(data) - np.min(data))
        
        ceil_cells = np.array([7.6e7,5.2e7,4.8e7,3.7e7])
        count_cells = np.array([16709793,3754589,2834337,2249482])
        expr_scs  = min_max_norm(np.array([450,3500,9000,20000])) + 1
        take_cells = 4e7
        bins_filter = f"{bin_column} <= 4"

        af_sub = 1
        af_div = 1
    
    elif heatmap == 'PEG':
        
        ceil_cells = np.array([8.9e7,8.9e7,5.8e7,7.6e7,9.2e7])
        count_cells = np.array([187000,230000,340426,506977,969420])
        expr_scs  = np.array([.0,.05,.1,.15,.2]) + 1
        take_cells = 4e7
        bins_filter = f"{bin_column}%2 == 0 and {bin_column} > 4"

        af_sub = 6
        af_div = 2
    
    else:
        st.error(f"Subset \"{heatmap}\" not found.")
        st.stop()

    if amplification_factor == 'Size':
        amp_factors = count_cells
    
    elif amplification_factor == 'Size/OD':
        amp_factors = count_cells / ceil_cells

    elif amplification_factor == 'Size/Sample':
        amp_factors = count_cells / ( ( take_cells /ceil_cells ) * ceil_cells )

    else:
        st.error(f"Normalization \"{amplification_factor}\" not found.")
        st.stop()

    df = df.query(bins_filter)
    m_df = df.fillna('wt').groupby([aa_column,n_aa_column,bin_column], as_index=False).sum()

    wt = m_df.query(f"{aa_column} == 'wt'").sort_values(bin_column).copy()

    wt['i_size'] =  wt['size']
    wt['size'] = wt['size'] * amp_factors
    wt['i_bin'] = wt[bin_column]
    wt[bin_column] = expr_scs
    wt = wt.assign(bin_score=lambda x: x[bin_column] * x['size'])

    wt['bin_mean'] = wt['bin_score'].sum() / wt['size'].sum()

    wt_mean = wt['bin_mean'].mean()

    single_muts = m_df.query(f"{n_aa_column} == 1").sort_values(by=[aa_column,bin_column]).reset_index(drop=True)

    single_muts['i_size'] = single_muts['size']
    single_muts['i_bin'] = single_muts[bin_column]

    single_muts['size'] = single_muts['size'].astype('float')
    single_muts[bin_column] = single_muts[bin_column].astype('float')

    for i in single_muts.index:
            single_muts.at[i, 'size'] = single_muts.at[i, 'size'] * amp_factors[int( (single_muts.at[i, bin_column] - af_sub) / af_div )]
            single_muts.at[i, bin_column] = expr_scs[int( (single_muts.at[i, bin_column] - af_sub) / af_div )]

    single_muts = single_muts.assign(bin_score=lambda x: x[bin_column] * x['size'])      
            
    tm = single_muts[[aa_column,'i_size','size','bin_score']].groupby([aa_column],as_index=False).sum()
    single_muts = single_muts.merge(tm, on=aa_column)

    single_muts = single_muts.assign(bin_mean=lambda x: np.log2( ( x['bin_score_y'] / x['size_y']) / wt_mean ) )

    waas = single_muts[aa_column].str.extract(r'(^[A-Z*][0-9]+)')[0].unique()
    maas = single_muts[aa_column].str.extract(r'([A-Z*]$)')[0].unique()

    all_muts = np.empty([len(waas) * len(maas)], dtype=object)

    for i,waa in enumerate(waas):
        for j,maa in enumerate(maas):
            all_muts[i*len(maas)+j] = waa+maa
            
    single_muts = pd.DataFrame(data=all_muts, columns=[aa_column]).merge(single_muts, on=aa_column, how='left')

    single_muts['wt_aa'] = single_muts[aa_column].str.extract(r'(^[A-Z*])')
    single_muts['position'] = single_muts[aa_column].str.extract(r'([0-9]+)').astype(int)
    single_muts['mutated_aa'] = single_muts[aa_column].str.extract(r'([A-Z*]$)')

    single_muts['wt_t'] = single_muts.apply(lambda x: x['wt_aa'] if x['wt_aa'] == x['mutated_aa'] else '', axis=1)
    single_muts['bin_mean'] = single_muts.apply(lambda x: 0 if x['wt_aa'] == x['mutated_aa'] else x['bin_mean'], axis=1)

    return single_muts, wt, wt_mean

def wt_plot(wt):
    
    base = alt.Chart(wt).mark_bar().encode(
            x='i_bin:O'
    )

    wtd = alt.hconcat(
        base.encode(
        y='i_size:Q',
        tooltip=['i_size']
        ).properties(title="Post Amp."), 
        base.encode(
        y='size:Q',
        tooltip=['size']
        ).properties(title="Pre Amp."),  
    )

    return wtd


def sm_plot(single_muts, aa_order):

    global aa_column 
    global n_aa_column 
    global bin_column

    if aa_order == 'A':
        sort_order = ['*','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    elif aa_order == '1':
        sort_order = ['*','R','K','H','D','E','Q','N','S','T','Y','W','F','A','I','L','M','V','G','P','C']

    brush = alt.selection_interval(encodings=['x'])
    selector = alt.selection_single(empty='all', fields=[aa_column])

    bar = alt.Chart(single_muts).mark_rect(color='grey').encode(
    x = alt.X('position:O', axis=alt.Axis(labels = False, ticks=False))
    ).properties(width=600).add_selection(brush)

    base = alt.Chart(single_muts).encode(
        x=alt.X('position:O'),
        y=alt.Y('mutated_aa:O', sort=sort_order),
    )

    muts = base.mark_rect().encode(
        color=alt.Color('bin_mean:Q', scale=alt.Scale(scheme='redyellowblue', domainMid=0)),
        tooltip = ['i_size_y',aa_column,'bin_mean'],
        opacity=alt.condition(selector, alt.value(1), alt.value(0.2))
    ).transform_filter(brush).add_selection(selector)

    wt = base.mark_text().encode(
        text='wt_t:O',
        tooltip = [aa_column,'bin_mean'],
        color=alt.condition(selector, alt.ColorValue('#000000'), alt.ColorValue('#CCCCCC'))
    ).transform_filter(brush).add_selection(alt.selection_single())

    pts = base.mark_point().encode(
        size=alt.Size('i_size_y:Q'),
        color=alt.ColorValue('grey'),
        opacity=alt.condition(selector, alt.value(1), alt.value(0.1)),
        tooltip = ['i_size_y',aa_column,'bin_mean'],
    ).transform_filter(brush).add_selection(alt.selection_single())

    bar1 = alt.Chart(single_muts).mark_bar().encode(
        x='i_bin:O',
        y='i_size_x:Q',
        tooltip = ['i_size_x','size_x'],
    ).transform_filter(
        selector
    ).add_selection(alt.selection_single())

    bar1 = alt.Chart(single_muts, title="Pre Amp.").mark_bar().encode(
        x='i_bin:O',
        y='size_x:Q',
        tooltip = ['size_x','i_size_x'],
    ).transform_filter(
        selector
    ).add_selection(alt.selection_single())

    bar2 = alt.Chart(single_muts, title="Post Amp.").mark_bar().encode(
        x='i_bin:O',
        y='i_size_x:Q',
        tooltip = ['i_size_x','size_x'],
    ).transform_filter(
        selector
    ).add_selection(alt.selection_single())

    hm = alt.vconcat(

        bar,
        muts + wt + pts,
        alt.hconcat(
            bar1,
            bar2
        )

    )

    return hm

uploaded_file = up_f()

if uploaded_file:
    
    if uploaded_file.name.endswith('.xlsx'):
        df = pd.read_excel(uploaded_file)
    elif uploaded_file.name.endswith('.p'):
        df = pd.read.pickle(uploaded_file)
    else:
        df = pd.read_csv(uploaded_file)

    heatmap = st.radio('Subset', ['Expression','PEG'])
    amplification_factor = st.radio('Normalization', ['Size','Size/OD','Size/Sample'])
    aa_order = st.radio('Sort AA', ['1','A'])

    single_muts, wt, wt_mean = build_df(df,heatmap, amplification_factor)

    wtd = wt_plot(wt)
    st.header('Wild type')
    st.altair_chart(wtd)

    hm = sm_plot(single_muts, aa_order)
    st.header("Single mutations heatmap")
    st.altair_chart(hm)