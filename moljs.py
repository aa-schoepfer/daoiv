import streamlit as st

import streamlit.components.v1 as components

def mol_component(resi, width=714, height=714):

    s_f = ""

    if resi != '':
        
        s_c = resi.split(',')
        for i,o in enumerate(s_c):
            if o != '':
                s_f += "'"
                s_m = o.split('-')
                for ii,oo in enumerate(s_m):
                    if len(s_m) == 1:
                        if int(oo):
                            s_f += str(int(oo)+1000)
                    elif len(s_m) == 2:
                        if ii == 0:
                            if int(oo):
                                s_f += str(int(oo)+1000) + '-'
                        elif ii == 1:
                            if int(oo):
                                s_f += str(int(oo)+1000)
                    else:
                        st.error("Only two numbers allowed. Make sure there is just one '-'")
                s_f += "',"
    

    components.html(f'''
<head>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>

    <script>
        var glviewer = null;
        var prot = null;

        $(document).ready(function () {{
            glviewer = $3Dmol.createViewer("gldiv", {{}});
            glviewer.setBackgroundColor(0xf4f4f4);
            prot = $3Dmol.download('pdb:1c0p', glviewer, {{doAssembly:true}}, function(){{
                glviewer.setStyle({{}},{{cartoon:{{color:'white'}}}});
                glviewer.setStyle({{resn:'FAD'}},{{stick:{{}}}});
                glviewer.setStyle({{resn:'DAL'}},{{stick:{{}}}});
                glviewer.setStyle({{resn:'PER'}},{{stick:{{}}}});
                glviewer.setStyle({{resi:[{s_f}]}},{{cartoon:{{color:'red'}}}});
                glviewer.render();
            }});
        }});
    </script>
</head>
<body>
    <div id="gldiv" style="width: 100%; height: 80vh;"></div>
</body>
'''
    , width=width, height=height)

