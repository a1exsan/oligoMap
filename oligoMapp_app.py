import pandas as pd
import streamlit as st

from bokeh.models.widgets import Button
from bokeh.models import CustomJS
from streamlit_bokeh_events import streamlit_bokeh_events

import oligoMap as map

#@st.cache
def load_map(fn):
    oMap = map.OligoMapReadWriter(fn)
    return oMap


st.set_page_config(layout="wide")

uploaded_xlsx = st.sidebar.file_uploader("Choose a map file (*.xlsx)")

if uploaded_xlsx is not None:
    with open(f'temp/{uploaded_xlsx.name}', 'wb') as f:
        f.write(uploaded_xlsx.getvalue())

    oMap = load_map(f'temp/{uploaded_xlsx.name}')
    df = oMap.synTab.astype(str)
    #print(df)
    st.dataframe(df, width=800, height=700)


    copy_button = Button(label="Copy to clipboard")
    copy_button.js_on_event("button_click", CustomJS(args=dict(df=oMap.copy_data().to_csv(sep='\t',
                                                                                #index=False,
                                                                                header=False)), code="""
    navigator.clipboard.writeText(df);
    """))

    no_event = streamlit_bokeh_events(
    copy_button,
    events="GET_TEXT",
    key="get_text",
    refresh_on_update=True,
    override_height=75,
    debounce_time=0)

    prog_name = st.text_area('Prog name: ', '', max_chars=100)

    if prog_name != '':
        oMap.convert_method2xml(out_name=prog_name)

        with open(f'temp/{prog_name}_main.pmp') as file:
            st.download_button(label='load main', data=file, file_name=f'temp/{prog_name}_main.pmp')

        with open(f'temp/{prog_name}_couple.pcp') as file:
            st.download_button(label='load couple', data=file, file_name=f'temp/{prog_name}_couple.pcp')

        with open(f'temp/{prog_name}_rmvDMT.prp') as file:
            st.download_button(label='load rmvDMT', data=file, file_name=f'temp/{prog_name}_rmvDMT.prp')

        with open(f'temp/{prog_name}_fin.pfp') as file:
            st.download_button(label='load finish', data=file, file_name=f'temp/{prog_name}_fin.pfp')

        with open(f'temp/{prog_name}_couple_dye.pcp') as file:
            st.download_button(label='load couple dye', data=file, file_name=f'temp/{prog_name}_couple_dye.pcp')


