import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, ElementTree
import pandas as pd

class polygenProgReader():
    def __init__(self, fn):
        self.step_list = None
        self.progXML = None
        self.name = None
        self.type = None
        self.descriprion = None
        self.CV1umol = 67 # colum volume ul 1 umol column slider
        self.CV5umol = 180 # colum volume ul 5 umol column slider
        self.pump_cal = {'RPM': 270, 'flow_ml_min': 2}
        self.pump_coef = self.pump_cal['flow_ml_min'] * 1000 / self.pump_cal['RPM'] # ul / rev

        with open(fn) as f:
            self.progXML = ET.parse(f).getroot()
            self.parse()

        self.__getDF()

    def __getDF(self):
        self.DF = pd.DataFrame(self.step_list)
        self.DF.fillna(value=0, inplace=True)
        self.DF['Time'] = self.DF['Time'].astype('float32')
        self.DF['RPM'] = self.DF['RPM'].astype('float32')
        self.DF['volume'] = self.DF['Time'] * self.DF['RPM'] * self.pump_coef / 60

    def parse(self):
        self.step_list = []
        for i in self.progXML.iter():
            if i.tag == 'Step':
                self.step_list.append(i.attrib)
            elif i.tag == 'MetaData':
                self.name = i.attrib['Name']
                self.type = i.attrib['Type']
                self.descriprion = i.attrib['Comment']

    def __call__(self, *args, **kwargs):
        self.DF['CV_1'] = self.DF['volume'] / self.CV1umol
        self.DF['CV_5'] = self.DF['volume'] / self.CV5umol
        return self.DF


class polygenProgConverter():
    def __init__(self, fn, main_n, couple_n, rDMT_n, finish_n,
                 file_names = [],
                 prog_types=['main', 'couple', 'removeDMT', 'finish'],
                 extentions = ['pmp', 'pcp', 'prp', 'pfp']):

        if fn.find('.csv') != -1:
            self.df = pd.read_csv(fn, sep='\t')
        elif fn.find('.xls') != -1:
            self.df = pd.read_excel(fn, sheet_name='synt_programs')

        print(self.df)
        if len(file_names) == 0:
            self.names = [main_n, couple_n, rDMT_n, finish_n]
        else:
            self.names = file_names
        self.extentions = extentions
        self.xmls = None
        self.prog_types = prog_types

    def get_prog_type_nums(self):
        ret = []
        d = {'pmp': '1', 'pcp': '2', 'prp': '4', 'pfp': '8'}
        for e in self.extentions:
            ret.append(d[e])
        return ret

    def __create_xml(self, df, meta_type, meta_name, meta_comment):
        # print(df.info())

        top = Element('ProgramSteps')

        for s, r, t, rpm, monitor in zip(df['Step'], df['Reagent'], df['Time'], df['RPM'], df['Monitor']):

            child = SubElement(top, 'Step')

            if r in ['WASH', 'DEBL']:
                child.attrib = {'Step': str(int(s)), 'Reagent': r, 'Time': str(t).replace(',', '.'),
                                'RPM': str(int(round(float(str(rpm).replace(',', '.'))))),
                                'Monitor': str(int(round(float(str(monitor).replace(',', '.')), 0)))}
            else:
                child.attrib = {'Step': str(int(s)), 'Reagent': r,
                                'Time': str(t).replace(',', '.'),
                                'RPM': str(int(round(float(str(rpm).replace(',', '.')))))}

        child = SubElement(top, 'MetaData')
        child.attrib = {'Name': meta_name, 'Type': meta_type, 'Comment': meta_comment}

        return top

    def convert_csv(self):

        self.xmls = []
        for prog, t, n in zip(self.prog_types, self.get_prog_type_nums(), self.names):
            df = self.df[self.df['type'] == prog]
            top = self.__create_xml(df, t, n, 'method')
            self.xmls.append(top)

    def write_xml(self):
        for top, name, ext in zip(self.xmls, self.names, self.extentions):
            # with open(f'{name}.{ext}', 'w') as f:
            ElementTree(top).write(f'{name}.{ext}', encoding='utf-8')


class ColumnSlider():
    def __init__(self, seq_list):
        self.CV1umol = 60  # colum volume ul 1 umol column slider
        self.CV5umol = 180  # colum volume ul 5 umol column slider
        self.seq_list = seq_list.copy()
        self.column_vec_count = None
        self.bases = pd.DataFrame([list(seq) for seq in self.seq_list])
        self.bases.fillna('null', inplace=True)
        self.bases = self.bases.T.values

    def get_column_count_vector(self):
        max_len = max([len(seq) for seq in self.seq_list])
        self.column_vec_count = [0 for i in range(max_len)]

        for seq in self.seq_list:
            for i, s in enumerate(seq):
                self.column_vec_count[i] += 1

        return self.column_vec_count

class Simulator():
    def __init__(self, fn, col_slider, df=None):
        if type(df) == type(pd.DataFrame()):
            self.df = df
        elif fn != '':
            self.df = pd.read_csv(fn, sep='\t')

        self.colSlider = col_slider
        if self.df['Time'].dtype == type(str):
            self.df['Time'] = self.df['Time'].str.replace(',', '.')
        self.df['Time'] = self.df['Time'].astype('float32')
        if self.df['RPM'].dtype == type(str):
            self.df['RPM'] = self.df['RPM'].str.replace(',', '.')
        self.df['RPM'] = self.df['RPM'].astype('float32')

        self.main = self.df[self.df['type'] == 'main']
        self.couple = self.df[self.df['type'] == 'couple']
        self.removeDMT = self.df[self.df['type'] == 'removeDMT']
        self.finish = self.df[self.df['type'] == 'finish']

        self.autoprime_sec = 3

        self.params = {}
        self.params['dbl_cycle_vol'] = 0.
        self.params['dbl_total_vol'] = 0.
        self.params['dbl_cycle_time'] = 0.
        self.params['total_time'] = 0.
        self.params['total_time, h'] = 0.
        self.params['acn_total_vol'] = 0.
        self.params['capA_cycle_vol'] = 0.
        self.params['capA_total_vol'] = 0.
        self.params['capA_cycle_time'] = 0.
        self.params['capB_cycle_vol'] = 0.
        self.params['capB_total_vol'] = 0.
        self.params['capB_cycle_time'] = 0.
        self.params['oxi_cycle_vol'] = 0.
        self.params['oxi_total_vol'] = 0.
        self.params['oxi_cycle_time'] = 0.
        self.params['act_cycle_vol'] = 0.
        self.params['act_total_vol'] = 0.
        self.params['act_cycle_time'] = 0.
        self.params['base_cycle_vol'] = 0.
        self.params['base_cycle_time'] = 0.
        self.params['base_A_vol'] = 0.
        self.params['base_C_vol'] = 0.
        self.params['base_G_vol'] = 0.
        self.params['base_T_vol'] = 0.
        self.params['total_waste'] = 0.

        self.repeats = {}

        self.debl_time = 0.
        self.capA_time = 0.
        self.capB_time = 0.
        self.oxid_time = 0.
        self.activ_time = 0.
        self.acn_time = 0.

    def get_VolTime(self, columns, time, rpm):
        #print(columns, time, rpm)
        vol = round(rpm * time * 2000 / (270 * 60), 1)
        sum_vol = columns * vol
        time = columns * time
        return vol, sum_vol, time

    def couple_cycle(self, bases, columns, cycle_count):
        for reag, time, rpm in zip(self.couple['Reagent'], self.couple['Time'], self.couple['RPM']):

            vol, sum_vol, t_time = self.get_VolTime(columns, time, rpm)

            if cycle_count == 1:
                self.repeats[reag] = 0
            elif cycle_count == 2:
                self.repeats[reag] += 1

            if reag == 'BASE':
                self.params['total_time'] += t_time + self.autoprime_sec

                if cycle_count == 1:
                    self.params['base_cycle_vol'] += vol
                    self.params['base_cycle_time'] += t_time

                for base in bases:
                    if base == 'A':
                        self.params['base_A_vol'] += vol
                    elif base == 'C':
                        self.params['base_C_vol'] += vol
                    elif base == 'G':
                        self.params['base_G_vol'] += vol
                    elif base == 'T':
                        self.params['base_T_vol'] += vol

            if reag == 'ACTIV':
                self.params['total_time'] += t_time + self.autoprime_sec
                self.params['act_total_vol'] += sum_vol

                if cycle_count == 1:
                    self.params['act_cycle_vol'] += vol
                    self.params['act_cycle_time'] += t_time

            if reag == 'ACNFLUSH':
                self.params['total_time'] += time
                self.params['acn_total_vol'] += vol

            if reag == 'WAIT':
                self.params['total_time'] += time

                if cycle_count == 1:
                    self.params['act_cycle_time'] += time
                    self.params['base_cycle_time'] += time

            if reag == 'CWAIT':
                self.params['total_time'] += time

                if cycle_count == 1:
                    self.params['act_cycle_time'] += time
                    self.params['base_cycle_time'] += time


    def main_cycle(self):
        cycle_count = 1
        step = 'INIT'
        for columns, bases in zip(self.colSlider.get_column_count_vector(),
                                  self.colSlider.bases):
            for reag, time, rpm in zip(self.main['Reagent'], self.main['Time'], self.main['RPM']):

                if cycle_count == 1:
                    self.repeats[reag] = 0
                elif cycle_count == 2:
                    self.repeats[reag] += 1

                if reag == 'COUPL':
                    self.couple_cycle(bases, columns, cycle_count)
                else:
                    vol, sum_vol, t_time = self.get_VolTime(columns, time, rpm)

                    if reag == 'DEBL':
                        step = 'DEBL'
                        self.debl_time = t_time
                        self.params['dbl_total_vol'] += sum_vol
                        self.params['total_time'] += t_time + self.autoprime_sec
                        if cycle_count == 1:
                            self.params['dbl_cycle_vol'] += vol

                    if reag == 'ACTIV':
                        step = 'ACTIV'
                        self.activ_time = t_time
                        self.params['act_total_vol'] += sum_vol
                        self.params['total_time'] += t_time + self.autoprime_sec
                        if cycle_count == 1:
                            self.params['act_cycle_vol'] += vol

                    if reag == 'CAPA':
                        step = 'CAP'
                        self.capA_time = t_time
                        self.params['capA_total_vol'] += sum_vol
                        self.params['total_time'] += t_time + self.autoprime_sec
                        if cycle_count == 1:
                            self.params['capA_cycle_vol'] += vol

                    if reag == 'CAPB':
                        step = 'CAP'
                        self.capB_time = t_time
                        self.params['capB_total_vol'] += sum_vol
                        self.params['total_time'] += t_time + self.autoprime_sec
                        if cycle_count == 1:
                            self.params['capB_cycle_vol'] += vol

                    if reag == 'OXID':
                        step = 'OXID'
                        self.oxid_time = t_time
                        self.params['oxi_total_vol'] += sum_vol
                        self.params['total_time'] += t_time + self.autoprime_sec
                        if cycle_count == 1:
                            self.params['oxi_cycle_vol'] += vol

                    if reag == 'WASH':
                        step = 'WASH'
                        self.activ_time = t_time
                        self.params['acn_total_vol'] += sum_vol
                        self.params['total_time'] += t_time + self.autoprime_sec

                    if reag == 'DWAIT':
                        if step == 'DEBL':
                            dt = time - self.debl_time
                            if dt > 0:
                                if cycle_count == 1:
                                    self.params['dbl_cycle_time'] += time
                                self.params['total_time'] += dt
                            else:
                                if cycle_count == 1:
                                    self.params['dbl_cycle_time'] += self.debl_time

                        if step == 'CAP':
                            dt = time - self.capA_time
                            if dt > 0:
                                if cycle_count == 1:
                                    self.params['capA_cycle_time'] += time
                            else:
                                if cycle_count == 1:
                                    self.params['capA_cycle_time'] += self.capA_time

                            dt = time - self.capB_time
                            if dt > 0:
                                if cycle_count == 1:
                                    self.params['capB_cycle_time'] += time
                                self.params['total_time'] += dt
                            else:
                                if cycle_count == 1:
                                    self.params['capB_cycle_time'] += self.capB_time

                        if step == 'OXID':
                            dt = time - self.oxid_time
                            if dt > 0:
                                if cycle_count == 1:
                                    self.params['oxi_cycle_time'] += time
                                self.params['total_time'] += dt
                            else:
                                if cycle_count == 1:
                                    self.params['oxi_cycle_time'] += self.oxid_time

                        if step == 'ACTIV':
                            dt = time - self.activ_time
                            if dt > 0:
                                if cycle_count == 1:
                                    self.params['act_cycle_time'] += time
                                self.params['total_time'] += dt
                            else:
                                if cycle_count == 1:
                                    self.params['act_cycle_time'] += self.oxid_time

            cycle_count += 1
        self.params['total_waste'] = self.params['dbl_total_vol'] + self.params['act_total_vol'] + \
                                     self.params['capA_total_vol'] + self.params['capB_total_vol'] +\
                                     self.params['oxi_total_vol'] + self.params['acn_total_vol'] + \
                                     self.params['base_A_vol'] + self.params['base_C_vol'] + \
                                     self.params['base_G_vol'] + self.params['base_T_vol']

        for key in ['DEBL', 'ACTIV', 'BASE', 'CAPA', 'CAPB', 'OXID']:
            if key in list(self.repeats.keys()):
                self.params[f'{key.lower()}_repeat_cycle'] = self.repeats[key]
            else:
                self.params[f'{key.lower()}_repeat_cycle'] = 0

    def rmvDMT_prog(self):
        columns = self.colSlider.get_column_count_vector()[0]
        for reag, time, rpm in zip(self.removeDMT['Reagent'], self.removeDMT['Time'], self.removeDMT['RPM']):

            vol, sum_vol, t_time = self.get_VolTime(columns, time, rpm)

            if reag == 'DEBL':
                self.params['dbl_total_vol'] += sum_vol
                self.params['total_time'] += t_time + self.autoprime_sec

            if reag == 'WASH':
                self.params['acn_total_vol'] += sum_vol
                self.params['total_time'] += t_time + self.autoprime_sec

            if reag == 'WAIT':
                self.params['total_time'] += time



    def finish_prog(self):
        columns = self.colSlider.get_column_count_vector()[0]
        for reag, time, rpm in zip(self.finish['Reagent'], self.finish['Time'], self.finish['RPM']):

            vol, sum_vol, t_time = self.get_VolTime(columns, time, rpm)

            if reag == 'WASH':
                self.params['acn_total_vol'] += sum_vol
                self.params['total_time'] += t_time + self.autoprime_sec

            if reag == 'WAIT':
                self.params['total_time'] += time

            if reag == 'GAS':
                self.params['total_time'] += t_time

        self.params['total_waste'] = self.params['dbl_total_vol'] + self.params['act_total_vol'] + \
                                     self.params['capA_total_vol'] + self.params['capB_total_vol'] + \
                                     self.params['oxi_total_vol'] + self.params['acn_total_vol'] + \
                                     self.params['base_A_vol'] + self.params['base_C_vol'] + \
                                     self.params['base_G_vol'] + self.params['base_T_vol']
        self.params['total_time, h'] = self.params['total_time'] / 3600


def calc_syn_params():

    def get_params(slider, df):
        sim = Simulator('', slider, df=df)
        sim.main_cycle()
        sim.rmvDMT_prog()
        sim.finish_prog()
        return sim.params

    slider1 = ColumnSlider(
        ['CCATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCCTGGTGCCGCGCGGCAGCCATATGGCTAGCATGACTGGTGGACAGC',
         'AAATGGGTCGCGGATCCATGCCGTCTGAACCGCCGTTCGGACGACATTTGATTTTTGCTTCTTTGACCTGCCTCATTGATGCGGTAT',
         'GCAAAAAAAGATACCATAACCAAAATGTTTATATATTATCTATTCTGCGTATGACTAGGAGTAAACCTGTAAATCGAACTGCCTTCT',
         'TTGCCGATGCACTAACCGCACCGCTCGACCATAAAGACAAAGGTTTGCAGTCTTTGACGCTGGATCAGTCCGTCAGGAAAAACGAGA',
         'cctcaGGTCTCtGAACGACAAGGTCAGCCGTTTCGACTTTATCCGCCAAATCGAAGTGGACGGGCAGCTCATTACCTTGGAGAGTGG',
         'AGAGTTCCAAGTATACAAACAAAGCCATTCCGCCTTAACCGCCTTTCAGACCGAGCAAATACAAGATTCGGAGCATTCCGGGAAGAT',
         'CCGAAGGCGGCAGGGCGACATATCGCGGGACGGCGTTCGGTTCAGACGATGCCGGCGGAAAACTGACCTACACCATAGATTTCGCCG',
         'CCAAGCAGGGAAACGGCAAAATCGAACATTTGAAATCGCTAGAACTCAATGTCGACCTGGCCGCCGCCGATATCAAGCCGGATGGAA',
         'CCCAGGAAGTTGCCGGCAGCGCGGAAGTGAAAACCGTAAACGGCATACGCCATATCGGCCTTGCCGCCAAGCAATAATGAACTCGAG',
         'GCTGCCTTTCTCTGACCACTGCCCTGATTCTGACCGCCTGCAGCAGCGGAGGGGGTGGTGTCGCCGCCGACATCGGTGCGGGGC',
         'AACGCCATGCCGTCATCAGCGGTTCCGTCCTTTACAACCAAGCCGAGAAAGGCAGTTACTCCCTCGGTATCTTTGGCGGAAAAG'])
    slider2 = ColumnSlider(
        ['GGTTGCGAAACGCCAGTTCAGAATCGGCGACATAGCGGGCGAACATACATCTTTTGACAAGCTTC',
         'AACTGAAGCTGGCGGCACAAGGTGCGGAAAAAACTTATGGAAACGGTGACAG'])
    slider3 = ColumnSlider(
        ['CCTCAATACGGGCAAATTGAAGAACtgagaccgaccgccg',
         'TATACTTGGAACTCTCCACTCTCCAAGGTAATGAGCTGC',
         'GGTATCTTTTTTTGCATACCGCATCAATGAGGCAGGTC',
         'GGTCAGAGAAAGGCAGCAGAAGGCAGTTCGAT'])
    slider4 = ColumnSlider(
        ['GCCGCCAGCTTCAGTTTCTCGTTTTTCCTGA',
         'ATCCGCGACCCATTTGCTGTCCACCAGTCA',
         'TTAGTGCATCGGCAAGCCCCGCACCGATGT',
         'TTGCCCGTATTGAGGCTGTCACCGTTTCCA',
         'TGGCGTTTCGCAACCATCTTCCCGGAATGC',
         'CCCTGCCGCCTTCGGGAAGCTTGTCAAAAG',
         'CGTTTCCCTGCTTGGCGGCGAAATCTATGG',
         'TGACGGCATGGCGTTTTCCATCCGGCTTGA',
         'CGGCAACTTCCTGGGCTTTTCCGCCAAAGA']
    )
    slider5 = ColumnSlider(
        ['CCATGGGCAGCAGCCATCAT',
         'CTCGAGTTCATTATTGCTTG']
    )
    # df = pd.read_excel('/home/alex/Documents/OligoSynt/DB/maps/synth_300322.xlsx', sheet_name='synt_programs')
    # df = pd.read_excel('/home/alex/Documents/OligoSynt/DB/maps/synth_110322_1.xlsx', sheet_name='synt_programs')
    df = pd.read_excel('/home/alex/Documents/OligoSynt/DB/maps/synth_220222_1.xlsx', sheet_name='synt_programs')

    syn_params = []

    syn_params.append(get_params(slider1, df))
    syn_params.append(get_params(slider2, df))
    syn_params.append(get_params(slider3, df))
    syn_params.append(get_params(slider4, df))
    syn_params.append(get_params(slider5, df))

    param_df = pd.DataFrame(syn_params)
    print(param_df.T)

    for key in param_df.keys():
        print(key, param_df[key].sum())

def append_to_excel(fpath, df, sheet_name):
    with pd.ExcelWriter(fpath, mode="a", if_sheet_exists='replace') as f:
        df.to_excel(f, sheet_name=sheet_name)

def order_analysis(fn):
    df_slider = pd.read_excel(fn, sheet_name='sliders')
    df_method = pd.read_excel(fn, sheet_name='method')
    df_reagents = pd.read_excel(fn, sheet_name='reagents')

    df_slider['length'] = df_slider['sequence'].str.len()
    print(f"total nucleotides: {df_slider['length'].sum()}")

    syn_params = []
    for sl in df_slider['slider'].unique():
        slider = ColumnSlider(df_slider[df_slider['slider'] == sl]['sequence'].values)
        #print(slider.bases)
        sim = Simulator('', slider, df=df_method)
        sim.main_cycle()
        sim.rmvDMT_prog()
        sim.finish_prog()

        print('@' * 50)
        print(sim.repeats)

        syn_params.append(sim.params)

    param_df = pd.DataFrame(syn_params)

    print(param_df.T)

    sum_list = []
    for key in param_df.keys():
        sum_list.append(param_df[key].sum())
    sum_df = pd.DataFrame([sum_list], columns=list(param_df.keys()))
    #print(sum_df.T)

    res = {}

    res['total_time, h'] = []
    res['ACN_vol, ml'] = []
    res['DCA, ml'] = []
    res['DCM, ml'] = []
    res['ACT, g'] = []
    res['Propionic anh, ml'] = []
    res['methylimidazole, ml'] = []
    res['I2, g'] = []
    res['Pyridine, ml'] = []
    res['DMT-dA(Bz)-CE phosphoramidite, g'] = []
    res['DMT-dC(Ac)-CE phosphoramidite, g'] = []
    res['DMT-dG(dmf)-CE phosphoramidite, g'] = []
    res['DMT-dT-CE phosphoramidite, g'] = []


    res['total_time, h'].append(sum_df['total_time'].values[0] / 3600)
    res['ACN_vol, ml'].append(sum_df['acn_total_vol'].values[0] / 1000)

    #print("res['ACN_vol, ml'].append(sum_df['acn_total_vol'].values[0] / 1000)", res['ACN_vol, ml'][0])

    dbl_conc = df_reagents[df_reagents['name']=='DEBL']['amount'].min() / df_reagents[df_reagents['name']=='DEBL']['amount'].sum()
    act_conc = df_reagents[df_reagents['name']=='ACT']['amount'].min() / df_reagents[df_reagents['name']=='ACT']['amount'].max()
    capA_conc = df_reagents[df_reagents['name']=='CAPA']['amount'].min() / df_reagents[df_reagents['name']=='CAPA']['amount'].sum()
    capB_conc = df_reagents[df_reagents['name']=='CAPB']['amount'].min() / df_reagents[df_reagents['name']=='CAPB']['amount'].sum()

    oxi_vol = df_reagents[(df_reagents['name']=='OXID')&(df_reagents['units']=='ml')]['amount'].sum()
    oxi_h2o = df_reagents[df_reagents['substance_name']=='Water']['amount'].min()
    I2_conc = df_reagents[df_reagents['substance_name']=='Iodine']['amount'].min() / oxi_vol
    pyr_conc = df_reagents[df_reagents['substance_name']=='Pyridine']['amount'].min() / oxi_vol

    dA_conc = df_reagents[df_reagents['name']=='BASE_A']['amount'].min() / df_reagents[df_reagents['name']=='BASE_A']['amount'].max()
    dG_conc = df_reagents[df_reagents['name']=='BASE_G']['amount'].min() / df_reagents[df_reagents['name']=='BASE_G']['amount'].max()
    dC_conc = df_reagents[df_reagents['name']=='BASE_C']['amount'].min() / df_reagents[df_reagents['name']=='BASE_C']['amount'].max()
    dT_conc = df_reagents[df_reagents['name']=='BASE_T']['amount'].min() / df_reagents[df_reagents['name']=='BASE_T']['amount'].max()

    res['DCA, ml'].append(sum_df['dbl_total_vol'].values[0] * dbl_conc / 1000)
    res['DCM, ml'].append(sum_df['dbl_total_vol'].values[0] / 1000 - res['DCA, ml'][0])
    res['ACT, g'].append(sum_df['act_total_vol'].values[0] * act_conc / 1000)
    res['ACN_vol, ml'][0] += sum_df['act_total_vol'].values[0] / 1000

    #print("res['ACN_vol, ml'][0] += sum_df['act_total_vol'].values[0] / 1000", res['ACN_vol, ml'][0])

    res['Propionic anh, ml'].append(sum_df['capA_total_vol'].values[0] * capA_conc / 1000)
    res['methylimidazole, ml'].append(sum_df['capB_total_vol'].values[0] * capB_conc / 1000)

    res['ACN_vol, ml'][0] += sum_df['capA_total_vol'].values[0] / 1000 - res['Propionic anh, ml'][0]

    #print("res['ACN_vol, ml'][0] += sum_df['capA_total_vol'].values[0] / 1000 - res['Propionic anh, ml'][0]", res['ACN_vol, ml'][0])

    res['ACN_vol, ml'][0] += sum_df['capB_total_vol'].values[0] / 1000 - res['methylimidazole, ml'][0]

    #print("res['ACN_vol, ml'][0] += sum_df['capB_total_vol'].values[0] / 1000 - res['methylimidazole, ml'][0]",
    #      res['ACN_vol, ml'][0])

    res['I2, g'].append(sum_df['oxi_total_vol'].values[0] * I2_conc / 1000)
    res['Pyridine, ml'].append(sum_df['oxi_total_vol'].values[0] * pyr_conc / 1000)
    res['ACN_vol, ml'][0] += sum_df['oxi_total_vol'].values[0] / 1000 - res['Pyridine, ml'][0] - oxi_h2o

    #print("res['ACN_vol, ml'][0] += oxi_vol - res['Pyridine, ml'][0] - oxi_h2o",
    #      res['ACN_vol, ml'][0], oxi_vol)

    res['DMT-dA(Bz)-CE phosphoramidite, g'].append(sum_df['base_A_vol'].values[0] * dA_conc / 1000)
    res['DMT-dC(Ac)-CE phosphoramidite, g'].append(sum_df['base_C_vol'].values[0] * dC_conc / 1000)
    res['DMT-dG(dmf)-CE phosphoramidite, g'].append(sum_df['base_G_vol'].values[0] * dG_conc / 1000)
    res['DMT-dT-CE phosphoramidite, g'].append(sum_df['base_T_vol'].values[0] * dT_conc / 1000)
    res['ACN_vol, ml'][0] += sum_df['base_A_vol'].values[0] / 1000
    res['ACN_vol, ml'][0] += sum_df['base_G_vol'].values[0] / 1000
    res['ACN_vol, ml'][0] += sum_df['base_C_vol'].values[0] / 1000
    res['ACN_vol, ml'][0] += sum_df['base_T_vol'].values[0] / 1000

    res['total_time, h'].append(0)
    res['ACN_vol, ml'].append(df_reagents[df_reagents['substance_name']=='ACN']['price'].min() * res['ACN_vol, ml'][0])
    res['ACT, g'].append(df_reagents[df_reagents['name']=='ACT']['price'].min() * res['ACT, g'][0])
    res['DCA, ml'].append(df_reagents[df_reagents['substance_name'].str.contains('DCA')]['price'].min() * res['DCA, ml'][0])
    res['DCM, ml'].append(df_reagents[df_reagents['substance_name'].str.contains('DCM')]['price'].min() * res['DCM, ml'][0])
    res['methylimidazole, ml'].append(df_reagents[df_reagents['substance_name'].str.contains('methylimidazole')]['price'].min() * res['methylimidazole, ml'][0])
    res['Propionic anh, ml'].append(df_reagents[df_reagents['substance_name'].str.contains('propionic')]['price'].min() * res['Propionic anh, ml'][0])
    res['I2, g'].append(df_reagents[df_reagents['substance_name'].str.contains('Iodine')]['price'].min() * res['I2, g'][0])
    res['Pyridine, ml'].append(df_reagents[df_reagents['substance_name'].str.contains('Pyridine')]['price'].min() * res['Pyridine, ml'][0])
    res['DMT-dA(Bz)-CE phosphoramidite, g'].append(df_reagents[df_reagents['substance_name'].str.contains('DMT-dA')]['price'].min() * res['DMT-dA(Bz)-CE phosphoramidite, g'][0])
    res['DMT-dC(Ac)-CE phosphoramidite, g'].append(df_reagents[df_reagents['substance_name'].str.contains('DMT-dC')]['price'].min() * res['DMT-dC(Ac)-CE phosphoramidite, g'][0])
    res['DMT-dG(dmf)-CE phosphoramidite, g'].append(df_reagents[df_reagents['substance_name'].str.contains('DMT-dG')]['price'].min() * res['DMT-dG(dmf)-CE phosphoramidite, g'][0])
    res['DMT-dT-CE phosphoramidite, g'].append(df_reagents[df_reagents['substance_name'].str.contains('DMT-dT')]['price'].min() * res['DMT-dT-CE phosphoramidite, g'][0])

    res_df = pd.DataFrame(res).T
    sum_cost = res_df[1].sum()
    res_df['cost%'] = res_df[1] * 100 / sum_cost

    print(res_df)

    out_keys, cycle_keys = [], []
    for key in sum_df.keys():
        if key.find('total') != -1:
            out_keys.append(key)
        if key.find('cycle') != -1:
            cycle_keys.append(key)
    out_keys.extend(['base_A_vol', 'base_C_vol', 'base_G_vol', 'base_T_vol'])
    print(cycle_keys)
    reag_vol_df = sum_df[out_keys] / 1000

    m_cycle_df = sum_df[cycle_keys].T
    slider = ColumnSlider(['AAA'])
    m_cycle_df['volume, cv'] = m_cycle_df[0] / slider.CV1umol
    m_cycle_df.reset_index(inplace=True)
    m_cycle_df.loc[m_cycle_df['index'].str.contains('time'), 'volume, cv'] = \
        m_cycle_df[m_cycle_df['index'].str.contains('time')][0]
    m_cycle_df.loc[m_cycle_df['index'].str.contains('repeat'), 'volume, cv'] = \
        m_cycle_df[m_cycle_df['index'].str.contains('repeat')][0]
    m_cycle_df['volume, cv'] = round(m_cycle_df['volume, cv'], 2)
    print(reag_vol_df.T)
    print(m_cycle_df)

    #try:
    append_to_excel(fn, res_df, 'reag_cost')
    #except:
    #    pass
    #try:
    append_to_excel(fn, reag_vol_df.T, 'reag_vol')
    #except:
    #    pass
    #try:
    append_to_excel(fn, m_cycle_df, 'cycle_params')
    #except:
    #    pass



def test_simul():

    def get_params(slider, df):
        sim = Simulator('', slider, df=df)
        sim.main_cycle()
        sim.rmvDMT_prog()
        sim.finish_prog()
        return sim.params

    #slider = ColumnSlider(
        #['CGAAGGTGTGACTTCCATG', 'TAATCAGACAAGGAACTGATTA', 'GAGCGGCTGTCTCCACAAGT', 'AGATTTGGACCTGCGAGCG'])
    #slider = ColumnSlider(
    #     ['TTTTTTTTTTTTTTTTTT', 'TTTTTTTTTTTTTTTTTT'])
    slider = ColumnSlider(
         ['CCCCCC', 'AAAAAA', 'GGGGGG', 'TTTTTT'])
    #print(slider.get_column_count_vector())

    #df = pd.read_excel('/home/alex/Documents/OligoSynt/DB/maps/synth_300322.xlsx', sheet_name='synt_programs')
    #df = pd.read_excel('/home/alex/Documents/OligoSynt/DB/maps/synth_110322_1.xlsx', sheet_name='synt_programs')
    df = pd.read_excel('/home/alex/Documents/OligoSynt/DB/maps/synth_220222_1.xlsx', sheet_name='synt_programs')

    sim = Simulator('', slider, df=df)
    sim.main_cycle()
    sim.rmvDMT_prog()
    sim.finish_prog()

    print('total time, h', sim.params['total_time'] / 3600)

    print('dbl_cycle_time', sim.params['dbl_cycle_time'])
    print('dbl_cycle_vol', sim.params['dbl_cycle_vol'])
    print('dbl_total_vol', sim.params['dbl_total_vol'])

    print('capA_cycle_time', sim.params['capA_cycle_time'])
    print('capA_cycle_vol', sim.params['capA_cycle_vol'])
    print('capA_total_vol', sim.params['capA_total_vol'])

    print('capB_cycle_time', sim.params['capB_cycle_time'])
    print('capB_cycle_vol', sim.params['capB_cycle_vol'])
    print('capB_total_vol', sim.params['capB_total_vol'])

    print('oxi_cycle_time', sim.params['oxi_cycle_time'])
    print('oxi_cycle_vol', sim.params['oxi_cycle_vol'])
    print('oxi_total_vol', sim.params['oxi_total_vol'])

    print('act_cycle_time', sim.params['act_cycle_time'])
    print('act_cycle_vol', sim.params['act_cycle_vol'])
    print('act_total_vol', sim.params['act_total_vol'])

    print('base_cycle_time', sim.params['base_cycle_time'])
    print('base_cycle_vol', sim.params['base_cycle_vol'])
    print('base_A_vol', sim.params['base_A_vol'])
    print('base_C_vol', sim.params['base_C_vol'])
    print('base_G_vol', sim.params['base_G_vol'])
    print('base_T_vol', sim.params['base_T_vol'])

    print('acn_total_vol', sim.params['acn_total_vol'])

    print('total_waste', sim.params['total_waste'])

def synprogs2params(fn):
    df = pd.read_csv(fn)
    out = []
    slider = ColumnSlider(
        ['CGAAGGTGTGACTTCCATG', 'TAATCAGACAAGGAACTGATTA', 'GAGCGGCTGTCTCCACAAGT', 'AGATTTGGACCTGCGAGCG'])
    for uid in df['uid'].unique():
        sim = Simulator('', slider, df=df[df['uid'] == uid])
        sim.main_cycle()
        sim.rmvDMT_prog()
        sim.finish_prog()

        out.append(sim.params)

    out = pd.DataFrame(out)
    out['uid'] = list(df['uid'].unique())

    out.to_csv(fn + '_params.csv', index=False)

def create_progs():
    #fn = '/home/alex/Documents/OligoSynt/DB/maps/nmap/synth_090622_gp_SIMA_mod.xlsx'
    fn = '/home/alex/Downloads/synth_050722_test_5umol.xlsx'
    main_fn = '/home/alex/Documents/OligoSynt/DB/maps/syn_programs/ID5L_1_main'
    coup_fn = '/home/alex/Documents/OligoSynt/DB/maps/syn_programs/ID5L_1_couple'
    remv_fn = '/home/alex/Documents/OligoSynt/DB/maps/syn_programs/ID5L_1_rmvDMT'
    fini_fn = '/home/alex/Documents/OligoSynt/DB/maps/syn_programs/ID5L_1_fin'
    copu_dye = '/home/alex/Documents/OligoSynt/DB/maps/syn_programs/ID5L_1_couple_dye'
    #copu_sima = '/home/alex/Documents/OligoSynt/DB/maps/syn_programs/ID31_couple_sima'
    #files = [main_fn, coup_fn, remv_fn, fini_fn, copu_sima, copu_fum6]
    files = [main_fn, coup_fn, remv_fn, fini_fn, copu_dye]
    prog_types = ['main', 'couple', 'removeDMT', 'finish', 'couple [6FAM]']
    extentions = ['pmp', 'pcp', 'prp', 'pfp', 'pcp']
    conv = polygenProgConverter(fn, main_fn, coup_fn, remv_fn, fini_fn, file_names=files,
                                prog_types=prog_types, extentions=extentions)
    conv.convert_csv()
    conv.write_xml()

def polygene_read():
    fn_main = '/home/alex/Documents/OligoSynt/S_P/1 umol Main.pmp'
    fn_couple = '/home/alex/Documents/OligoSynt/S_P/1 umol DNA dbl cpl.pcp'
    fn_rmvDMT = '/home/alex/Documents/OligoSynt/S_P/1 umol 2col Remove DMT.prp'
    fn_fin = '/home/alex/Documents/OligoSynt/S_P/1 umol Finish .pfp'
    main = polygenProgReader(fn_main)
    couple = polygenProgReader(fn_couple)
    rmv = polygenProgReader(fn_rmvDMT)
    fin = polygenProgReader(fn_fin)
    main.DF['type'] = 'main'
    couple.DF['type'] = 'couple'
    rmv.DF['type'] = 'removeDMT'
    fin.DF['type'] = 'finish'
    df = pd.concat([main.DF, couple.DF, rmv.DF, fin.DF])
    print(df)
    df.to_csv('/home/alex/Documents/OligoSynt/S_P/start_bundle.csv', index=False)


if __name__ == '__main__':
    #test_simul()
    #calc_syn_params()
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/order _analysis.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/order _analysis_large_scale.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/long_oligo_1col.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/order_Nikolay.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/order_timofey.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/opt_PS_2_RP_short.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/NR_HAPC_opt_1.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/FAM_SIMA_130422.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/HAPC_130422_opt3.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/HAPC_opt_4.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/HAPC_opt_5.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/HAPC_opt_6.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/HAPC_opt_7.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/HAPC_opt_6_dyes.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/HAPC_ID25_calc.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/HAPC_ID25_calc_10seq.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/test_orders/service_calc_1seq.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/NF_NR_2.xlsx')
    #order_analysis('/home/alex/Documents/OligoSynt/OligoOrders/test_orders/primer6.xlsx')
    create_progs()
    #synprogs2params('/home/alex/Documents/OligoSynt/DB/data/synprog.csv')

    #polygene_read()


