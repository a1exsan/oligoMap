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

