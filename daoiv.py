#Title: DAOx Illumina visualizations
#Author: Alexandre Schoepfer
#Version: 07.05.2021

import numpy as np
import pandas as pd
import altair as alt
import streamlit as st

from moljs import mol_component

aa_column  = 'aa_mutation'
n_aa_column = 'n_aa_substitutions'
bin_column = 'bin'

def up_f(types=['txt','csv','.xlsx']):
    up_file = st.file_uploader("Upload File",type=types)
    return up_file

def wt_plot(df):
    
    base = alt.Chart(df).mark_bar().encode(
            x='bin:O'
    ).transform_filter(alt.FieldEqualPredicate(field=aa_column, equal='wt'))

    wtd = alt.hconcat(
        base.encode(
        y='size_x:Q',
        tooltip=['size_x']
        ).properties(title="Pre Amp."),
        base.encode(
        y='i_size:Q',
        tooltip=['i_size']
        ).properties(title="Post Amp."),   
    )

    return wtd

def sm_plot(df):

    global aa_column 
    global n_aa_column 
    global bin_column

    global aa_order
    global heatmap

    if aa_order == 'A':
        sort_order = ['*','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    elif aa_order == '1':
        sort_order = ['*','R','K','H','D','E','Q','N','S','T','Y','W','F','A','I','L','M','V','G','P','C']

    if heatmap == 'Expression':
        res = 'wt_bin_mean'
        a_res = 'wt_coef'
        dm = 0
    
    elif heatmap == 'PEG':
        res = 'wt_coef'
        a_res = 'wt_bin_mean'
        dm = 0
    
    elif heatmap == '0%':
        res = 'bin_score_y'
        a_res = 'wt_bin_mean'
        df = df.query("(bin >= 5 and bin <= 6 ) or mutated_aa == wt_t")
        dm = 1.101744 	

    elif heatmap == '5%':
        res = 'bin_score_y'
        a_res = 'wt_bin_mean'
        df = df.query("(bin >= 7 and bin <= 8 ) or mutated_aa == wt_t")
        dm = 1.129350

    elif heatmap == '10%':
        res = 'bin_score_y'
        a_res = 'wt_bin_mean'
        df = df.query("(bin >= 9 and bin <= 10 ) or mutated_aa == wt_t")
        dm = 1.219020

    elif heatmap == '15%':
        res = 'bin_score_y'
        a_res = 'wt_bin_mean'
        df = df.query("(bin >= 11 and bin <= 12 ) or mutated_aa == wt_t")
        dm = 1.290426

    elif heatmap == '20%':
        res = 'bin_score_y'
        a_res = 'wt_bin_mean'
        df = df.query("(bin >= 13 and bin <= 14 ) or mutated_aa == wt_t")
        dm = 1.528374

    brush = alt.selection_interval(encodings=['x'])
    selector = alt.selection_single(empty='all', fields=[aa_column])

    bar = alt.Chart(df).mark_rect(color='grey').encode(
    x = alt.X('position:O', axis=alt.Axis(labels = False, ticks=False))
    ).properties(width=714).add_selection(brush)

    base = alt.Chart(df).encode(
        x=alt.X('position:O'),
        y=alt.Y('mutated_aa:O', sort=sort_order),
    ).transform_filter(
        alt.FieldOneOfPredicate(field='mutated_aa',oneOf=sort_order)
    )

    muts = base.mark_rect().encode(
        color=alt.Color(f'{res}:Q', scale=alt.Scale(scheme='redyellowblue', domainMid=dm)),
        tooltip = [aa_column,res,a_res,'wt_inte','r2','sum(size_x)','tot_size_x','sum(i_size)','tot_i_size','e_conf','p_conf','t_conf'],
        opacity=alt.condition(selector, alt.value(1), alt.value(0.05)),
    ).transform_filter(brush).add_selection(selector)

    wtd = base.mark_text().encode(
        text='wt_t:O',
        tooltip = [aa_column],
        color=alt.condition(selector, alt.ColorValue('#000000'), alt.ColorValue('#CCCCCC'))
    ).transform_filter(brush).transform_filter(
        alt.FieldOneOfPredicate(field='wt_t',oneOf=sort_order)
    ).add_selection(alt.selection_single())

    bar1 = alt.Chart(df, title="Pre Amp.").mark_bar().encode(
        x='bin:O',
        y='size_x:Q',
        tooltip = [aa_column,'size_x','i_size'],
    ).transform_filter(
        selector
    ).add_selection(alt.selection_single())

    bar2 = alt.Chart(df, title="Post Amp.").mark_bar().encode(
        x='bin:O',
        y='i_size:Q',
        tooltip = [aa_column,'i_size','size_x'],
    ).transform_filter(
        selector
    ).add_selection(alt.selection_single())

    hm = alt.vconcat(

        bar,
        muts + wtd,
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

    st.sidebar.subheader("Heatmap")
    
    heatmap = st.sidebar.radio("Subset", ['Expression','PEG','0%','5%','10%','15%','20%'])
    
    aa_order = st.sidebar.radio("Sort AA", ['A','1'], index=1)

    wtd = wt_plot(df.query(f"{aa_column} == 'wt'"))
    st.header("Wild type")
    st.altair_chart(wtd)

    st.sidebar.subheader("Filters")
    ex_lm = st.sidebar.number_input(f"Expression lower limit (Default: {df['wt_bin_mean'].min()})", value=df['wt_bin_mean'].min()-0.01)
    ex_hm = st.sidebar.number_input(f"Expression higher limit (Default: {df['wt_bin_mean'].max()})", value=df['wt_bin_mean'].max()+0.01)
    co_lm = st.sidebar.number_input(f"PEG slope lower limit (Default: {df['wt_coef'].min()})", value=df['wt_coef'].min()-0.01)
    co_hm = st.sidebar.number_input(f"PEG slope higher limit (Default: {df['wt_coef'].max()})", value=df['wt_coef'].max()+0.01)
    in_lm = st.sidebar.number_input(f"PEG intercept lower limit (Default: {df['wt_inte'].min()})", value=df['wt_inte'].min()-0.01)
    in_hm = st.sidebar.number_input(f"PEG intercept higher limit (Default: {df['wt_inte'].max()})", value=df['wt_inte'].max()+0.01)
    
    as_lm = st.sidebar.number_input(f"Pre amplification, total frequency >= than", value=0)
    rs_lm = st.sidebar.number_input(f"Post amplification, total frequency >= than", value=0)

    ec_lm = st.sidebar.number_input(f"Expression confidence >= than", value=0)
    ep_lm = st.sidebar.number_input(f"PEG confidence >= than", value=0)
    et_lm = st.sidebar.number_input(f"Total confidence >= than", value=0)

    df = df.query(f"wt_bin_mean >= {ex_lm}")
    df = df.query(f"wt_bin_mean <= {ex_hm}")
    df = df.query(f"wt_coef >= {co_lm}")
    df = df.query(f"wt_coef <= {co_hm}")
    df = df.query(f"wt_inte >= {in_lm}")
    df = df.query(f"wt_inte <= {in_hm}")

    df = df.query(f"tot_size_x >= {as_lm}")
    df = df.query(f"tot_i_size >= {rs_lm}")

    df = df.query(f"e_conf >= {ec_lm}")
    df = df.query(f"p_conf >= {ep_lm}")
    df = df.query(f"t_conf >= {et_lm}")   

    hm = sm_plot(df)
    st.header(f"Single mutations heatmap for {heatmap}")
    st.altair_chart(hm)

    inpt = st.sidebar.text_input("Highlight residues, e.g. 50 or 50-60 or 50-60,300-310")
    st.header("DAOx PDB: 1c0p")
    mol_component(inpt)
