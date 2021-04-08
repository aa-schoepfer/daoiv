#Title: DAOx Illumina visualizations
#Author: Alexandre Schoepfer
#Version: 08.04.2021

import numpy as np
import pandas as pd
import altair as alt
import streamlit as st

uploaded_file = st.file_uploader("Upload File",type=['txt','csv','.xlsx'])

if uploaded_file:
    
    if uploaded_file.name.endswith('.xlsx'):
        df = pd.read_excel(uploaded_file)
    elif uploaded_file.name.endswith('.p'):
        df = pd.read.pickle(uploaded_file)
    else:
        df = pd.read_csv(uploaded_file)

    heatmap = st.radio('Radio', ['Expression','PEG'])

    aa_column  = 'aa_mutation'
    n_aa_column = 'n_aa_substitutions'
    bin_column = 'bin'

    if heatmap == 'Expression':

        def min_max_norm(data):
            return (data - np.min(data)) / (np.max(data) - np.min(data))

        ceil_cells = np.array([7.6e7,5.2e7,4.8e7,3.7e7])
        count_cells = np.array([16709793,3754589,2834337,2249482])
        expr_scs  = min_max_norm(np.array([450,3500,9000,20000]))
        amp_factors = count_cells / ceil_cells
        bins_filter = "bin <= 4"

        df = df.query(bins_filter)
        m_df = df.fillna('wt').groupby([aa_column,n_aa_column,bin_column], as_index=False).sum()

        wt = m_df.query(f"{aa_column} == 'wt'").sort_values(bin_column).copy()

        wt['i_size'] =  wt['size']
        wt['size'] = wt['size'] * amp_factors
        wt['i_bin'] = wt['bin']
        wt['bin'] = expr_scs
        wt = wt.assign(bin_score=lambda x: x[bin_column] * x['size'])

        wt['bin_mean'] = wt['bin_score'].sum() / wt['size'].sum()

        wt_mean = wt['bin_mean'].mean()

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

        st.header('Wild type')
        st.altair_chart(wtd)

        single_muts = m_df.query(f"{n_aa_column} == 1").sort_values(by=[aa_column,bin_column]).reset_index(drop=True)

        single_muts['i_size'] = single_muts['size']
        single_muts['i_bin'] = single_muts[bin_column]

        single_muts['size'] = single_muts['size'].astype('float')
        single_muts[bin_column] = single_muts[bin_column].astype('float')


        for i in single_muts.index:
                single_muts.at[i, 'size'] = single_muts.at[i, 'size'] * amp_factors[int(single_muts.at[i, bin_column])-1]
                single_muts.at[i, bin_column] = expr_scs[int(single_muts.at[i, bin_column])-1]


        single_muts = single_muts.assign(bin_score=lambda x: x[bin_column] * x['size'])      
                
        tm = single_muts[[aa_column,'i_size','size','bin_score']].groupby([aa_column],as_index=False).sum()
        single_muts = single_muts.merge(tm, on=aa_column)

        single_muts = single_muts.assign(bin_mean=lambda x: x['bin_score_y'] / x['size_y'])

        single_muts['wt_aa'] = single_muts[aa_column].str.extract(r'(^[A-Z*])')
        single_muts['position'] = single_muts[aa_column].str.extract(r'([0-9]+)').astype(int)
        single_muts['mutated_aa'] = single_muts[aa_column].str.extract(r'([A-Z*]$)')

        brush = alt.selection_interval(encodings=['x'])
        selector = alt.selection_single(empty='all', fields=[aa_column])

        bar = alt.Chart(single_muts).mark_bar(size=2,color='grey').encode(
        alt.X('position:O', axis=alt.Axis(labels=False, ticks=False))
        ).properties(width=alt.Step(2)).add_selection(brush)

        muts = alt.Chart(single_muts).mark_rect().encode(
        x='position:O',
        y='mutated_aa:O',
        color=alt.Color('bin_mean:Q', scale=alt.Scale(domain=[0,1])),
        tooltip = ['i_size_y',aa_column,'bin_mean'],
        opacity=alt.condition(selector, alt.value(1), alt.value(0.2))
        ).transform_filter(brush).add_selection(selector)

        wtc = alt.Chart(single_muts).mark_text().encode(
        x='position:O', 
        y='wt_aa:O',
        text='wt_aa:O',
        tooltip = ['position'],
        color=alt.condition(selector, alt.ColorValue('#000000'), alt.ColorValue('#CCCCCC'))
        ).transform_filter(brush).add_selection(alt.selection_single())

        pts = alt.Chart(single_muts).mark_point().encode(
        x='position:O',
        y='mutated_aa:O',
        size=alt.Size('i_size_y:Q'),
        color=alt.condition(alt.FieldGTPredicate(field='bin_mean', gt=wt_mean), alt.ColorValue('#CC3333'), alt.ColorValue('#A0A0A0')),
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
            muts + wtc + pts,
            alt.hconcat(
                bar1,
                bar2
            )

        )
        st.header("Single mutations heatmap")
        st.markdown(f"Wild-type mean: **{wt_mean}**")
        st.altair_chart(hm)

    elif heatmap == 'PEG':

        ceil_cells = np.array([8.9e7,8.9e7,5.8e7,7.6e7,9.2e7])
        count_cells = np.array([187000,230000,340426,506977,969420])
        amp_factors = count_cells / ceil_cells
        bins_filter = "bin%2 == 0 and bin > 4"

        df = df.query(bins_filter)
        m_df = df.fillna('wt').groupby([aa_column,n_aa_column,bin_column], as_index=False).sum()

        wt = m_df.query(f"{aa_column} == 'wt'").reset_index(drop=True)
        wt['pre_amp'] = wt['size'] * amp_factors
        wt['peg'] = np.array([0,0.05,0.1,0.15,0.2])
        wt = wt.assign(bin_score=lambda x: x['peg'] * x['pre_amp'])
        wt['bin_mean'] = wt['bin_score'].sum() / wt['pre_amp'].sum()
        wt_mean = wt['bin_mean'].mean()

        base = alt.Chart(wt).mark_bar().encode(
            x='peg:O'
        )

        wtd = alt.hconcat(
            base.encode(
            y='size:Q',
            tooltip=['size']
            ).properties(title="Post Amp."), 
            base.encode(
            y='pre_amp:Q',
            tooltip=['pre_amp']
            ).properties(title="Pre Amp."),  
        )

        st.header('Wild type')
        st.altair_chart(wtd) 

        single_muts = m_df.query(f"{n_aa_column} == 1").copy()

        single_muts['pre_amp'] = single_muts['size'].astype('float')
        single_muts['peg'] = np.nan

        for i in single_muts.index:
                single_muts.at[i, 'pre_amp'] = single_muts.at[i, 'size'] * amp_factors[int((single_muts.at[i,'bin'] - 6) / 2)]
                single_muts.at[i, 'peg'] = round((0.025 * (single_muts.at[i,'bin'] - 6))*100)/100

        single_muts = single_muts.assign(bin_score=lambda x: x['peg'] * x['pre_amp'])
                
        tm = single_muts[[aa_column,'size','pre_amp','bin_score']].groupby([aa_column],as_index=False).sum()
        single_muts = single_muts.merge(tm, on=aa_column)
            
        single_muts = single_muts.assign(peg_mean=lambda x: x['bin_score_y'] / x['pre_amp_y'])
        
        single_muts['wt_aa'] = single_muts[aa_column].str.extract(r'(^[A-Z*])')
        single_muts['position'] = single_muts[aa_column].str.extract(r'([0-9]+)').astype(int)
        single_muts['mutated_aa'] = single_muts[aa_column].str.extract(r'([A-Z*]$)')

        brush = alt.selection_interval(encodings=['x'])
        selector = alt.selection_single(empty='all', fields=[aa_column])

        bar = alt.Chart(single_muts).mark_bar(size=2,color='grey').encode(
        alt.X('position:O', axis=alt.Axis(labels=False, ticks=False))
        ).properties(width=alt.Step(2)).add_selection(brush)

        muts = alt.Chart(single_muts).mark_rect().encode(
        x='position:O',
        y='mutated_aa:O',
        color=alt.Color('peg_mean:Q', scale=alt.Scale(domain=[0.2,0])),
        tooltip = ['size_y',aa_column,'peg_mean'],
        opacity=alt.condition(selector, alt.value(1), alt.value(0.2))
        ).transform_filter(brush).add_selection(selector)

        wtc = alt.Chart(single_muts).mark_text().encode(
        x='position:O', 
        y='wt_aa:O',
        text='wt_aa:O',
        tooltip = ['position'],
        color=alt.condition(selector, alt.ColorValue('#000000'), alt.ColorValue('#CCCCCC'))
        ).transform_filter(brush).add_selection(alt.selection_single())

        pts = alt.Chart(single_muts).mark_point().encode(
        #alt.ColorValue('grey'),
        x='position:O',
        y='mutated_aa:O',
        size=alt.Size('size_y:Q'),
        color=alt.condition(alt.FieldLTPredicate(field='peg_mean', lt=wt_mean), alt.ColorValue('#CC3333'), alt.ColorValue('#A0A0A0')),
        opacity=alt.condition(selector, alt.value(1), alt.value(0.1)),
        tooltip = ['size_y',aa_column,'peg_mean'],
        ).transform_filter(brush).add_selection(alt.selection_single())
        bar1 = alt.Chart(single_muts, title="Pre Amp.").mark_bar().encode(
            x='peg:O',
            y='pre_amp_x:Q',
            tooltip = ['pre_amp_x:Q','size_x'],
        ).transform_filter(
            selector
        ).add_selection(alt.selection_single())

        bar2 = alt.Chart(single_muts, title="Post Amp.").mark_bar().encode(
            x='peg:O',
            y='size_x:Q',
            tooltip = ['size_x','pre_amp_x:Q'],
        ).transform_filter(
            selector
        ).add_selection(alt.selection_single())

        hm = alt.vconcat(

            bar,
            muts + wtc + pts,
            alt.hconcat(
                bar1,
                bar2
            )

        )
        st.header("Single mutations heatmap")
        st.markdown(f"Wild-type mean: **{wt_mean}**")
        st.altair_chart(hm)
