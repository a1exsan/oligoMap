import pandas as pd
import oligoMass.molmassOligo as mmo
import numpy as np

class OligoMapReadWriter():
    def __init__(self, fn):
        self.filePath = fn

        self.copy_list = [
                        'Molecular Mass, Da',
                       'Molar extinction, oe*L/mol',
                       'Total amount, OE',
                       'Total amount, nmol',
                       'Max yield, nmol',
                       'Yield%'
                          ]

        self.readSyntTab()


    def add_molProperties(self):
        mass, ext_cf = [], []
        self.synTab = self.synTab.T

        if self.synTab['Sequence (5-3)'].max() != 'null':
            for seq in self.synTab['Sequence (5-3)']:
                oligos = mmo.oligoNASequence(seq)
                mass.append(oligos.getAvgMass())
                ext_cf.append(oligos.getExtinction())

            self.synTab['Molecular Mass, Da'] = mass
            self.synTab['Molar extinction, oe*L/mol'] = ext_cf
        else:
            self.synTab['Molecular Mass, Da'] = self.getNullListDF(self.synTab)
            self.synTab['Molar extinction, oe*L/mol'] = self.getNullListDF(self.synTab)

        self.synTab = self.synTab.T


    def getNullListDF(self, df):
        return ['null' for i in range(df.shape[0])]


    def add_yieldParams(self):
        self.synTab = self.synTab.T
        if self.synTab['product_volume, ml'].max() != 'null' and self.synTab['product_conc, oe/ml'].max() != 'null':
            self.synTab['Total amount, OE'] = self.synTab['product_volume, ml'] * self.synTab['product_conc, oe/ml']
            self.synTab['Total amount, nmol'] = self.synTab['Total amount, OE'] * 1e6 \
                                            / self.synTab['Molar extinction, oe*L/mol']
            self.synTab['Max yield, nmol'] = self.synTab['support_amount, mg'] * self.synTab['support_capacity, nM/mg']
            self.synTab['Yield%'] = self.synTab['Total amount, nmol'] * 100 / self.synTab['Max yield, nmol']
        else:
            self.synTab['Total amount, OE'] = self.getNullListDF(self.synTab)
            self.synTab['Total amount, nmol'] = self.getNullListDF(self.synTab)
            self.synTab['Total amount, nmol'] = self.getNullListDF(self.synTab)
            self.synTab['Max yield, nmol'] = self.getNullListDF(self.synTab)
            self.synTab['Yield%'] = self.getNullListDF(self.synTab)

        self.synTab = self.synTab.T

    def add_reagentParams(self):

        self.reagTab = pd.read_excel(self.filePath, sheet_name='Reagents').fillna('null').reset_index()
        #print(self.reagTab)

        self.synTab = self.synTab.T

        synDate = self.synTab['Date'].max()

        if str(synDate.date()) != 'null':
            conc_d = {}
            conc_d['DEBL'] = self.reagTab[self.reagTab['name'] == 'DEBL']['amount'].min() / \
                       self.reagTab[self.reagTab['name'] == 'DEBL']['amount'].sum()
            conc_d['ACT'] = self.reagTab[self.reagTab['name'] == 'ACT']['amount'].min() / \
                       self.reagTab[self.reagTab['name'] == 'ACT']['amount'].max()
            conc_d['CAPA'] = self.reagTab[self.reagTab['name'] == 'CAPA']['amount'].min() / \
                        self.reagTab[self.reagTab['name'] == 'CAPA']['amount'].sum()
            conc_d['CAPB'] = self.reagTab[self.reagTab['name'] == 'CAPB']['amount'].min() / \
                        self.reagTab[self.reagTab['name'] == 'CAPB']['amount'].sum()

            oxi_vol = self.reagTab[(self.reagTab['name'] == 'OXID') & (self.reagTab['units'] == 'ml')]['amount'].sum()
            oxi_h2o = self.reagTab[self.reagTab['substance_name'] == 'Water']['amount'].min()
            conc_d['OXID'] = self.reagTab[self.reagTab['substance_name'] == 'Iodine']['amount'].min() / oxi_vol
            pyr_conc = self.reagTab[self.reagTab['substance_name'] == 'Pyridine']['amount'].min() / oxi_vol

            conc_d['BASE_A'] = self.reagTab[self.reagTab['name'] == 'BASE_A']['amount'].min() / \
                      self.reagTab[self.reagTab['name'] == 'BASE_A']['amount'].max()
            conc_d['BASE_G'] = self.reagTab[self.reagTab['name'] == 'BASE_G']['amount'].min() / \
                      self.reagTab[self.reagTab['name'] == 'BASE_G']['amount'].max()
            conc_d['BASE_C'] = self.reagTab[self.reagTab['name'] == 'BASE_C']['amount'].min() / \
                      self.reagTab[self.reagTab['name'] == 'BASE_C']['amount'].max()
            conc_d['BASE_T'] = self.reagTab[self.reagTab['name'] == 'BASE_T']['amount'].min() / \
                      self.reagTab[self.reagTab['name'] == 'BASE_T']['amount'].max()

            reag_list = ['ACT', 'DEBL', 'CAPA', 'CAPB', 'OXID', 'BASE_A', 'BASE_C', 'BASE_G', 'BASE_T']
            for reag in reag_list:
                df = self.reagTab[self.reagTab['name'] == reag]
                if str(df['prep_date'].max()) != 'null':
                    self.synTab[f'Reagent {reag} old, days'] = str(abs((self.reagTab[self.reagTab['name'] == reag]['prep_date']
                                                            - synDate).dt.days.values[0]))
                    self.synTab[f'Concentration {reag} units/ml'] = conc_d[reag]
                else:
                    self.synTab[f'Reagent {reag} old, days'] = 'null'
                    self.synTab[f'Concentration {reag} units/ml'] = 'null'
                self.copy_list.append(f'Reagent {reag} old, days')
                self.copy_list.append(f'Concentration {reag} units/ml')

        #print('@'*50)
        #print(synDate)

        self.synTab = self.synTab.T

    def add_synProgParams(self):
        self.synTab = self.synTab.T

        

        self.synTab = self.synTab.T

    def readSyntTab(self):
        self.synTab = pd.read_excel(self.filePath, sheet_name='synthesis').fillna('null').reset_index()
        self.synTab.set_index('param', inplace=True)
        self.synTab.drop('index', axis=1, inplace=True)

        drop_list = []
        for key in self.synTab.keys():
            if not self.synTab[key].loc['Sequence (5-3)'] != 'null':
                drop_list.append(key)

        self.synTab.drop(drop_list, axis=1, inplace=True)
        self.add_molProperties()
        self.add_yieldParams()
        self.add_reagentParams()
        self.add_synProgParams()

        #print(self.synTab)
        #self.synTab.reset_index(inplace=True)
        #self.add_sheet_to_excel(self.synTab, 'test')

    def add_sheet_to_excel(self, df, sheetname):
        with pd.ExcelWriter(self.filePath, mode="a", if_sheet_exists='overlay') as f:
            df.to_excel(f, sheet_name=sheetname, index=False, header=False)

    def copy_data(self):
        self.synTab = self.synTab.T
        column_list = self.copy_list
        real_list = []
        for key in self.synTab.keys():
            if key in column_list:
                real_list.append(key)

        out_df = self.synTab[real_list]
        self.synTab = self.synTab.T

        return out_df.T




def test1():
    map = OligoMapReadWriter('/home/alex/Documents/OligoSynt/DB/maps/nmap/synth_180422_6opt.xlsx')


if __name__ == '__main__':
    test1()
