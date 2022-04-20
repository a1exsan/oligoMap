import pandas as pd
import oligoMass.molmassOligo as mmo

class OligoMapReadWriter():
    def __init__(self, fn):
        self.filePath = fn
        self.readSyntTab()


    def get_molProperties(self):
        mass, ext_cf = [], []
        self.synTab = self.synTab.T
        for seq in self.synTab['Sequence (5-3)']:
            oligos = mmo.oligoNASequence(seq)
            mass.append(oligos.getAvgMass())
            ext_cf.append(oligos.getExtinction())

        self.synTab['Molecular Mass, Da'] = mass
        self.synTab['Molar extinction, oe*L/mol'] = ext_cf
        self.synTab = self.synTab.T


    def get_yieldParams(self):
        self.synTab = self.synTab.T

        self.synTab['Total amount, OE'] = self.synTab['product_volume, ml'] * self.synTab['product_conc, oe/ml']
        self.synTab['Total amount, nmol'] = self.synTab['Total amount, OE'] * 1e6 \
                                            / self.synTab['Molar extinction, oe*L/mol']
        self.synTab['Max yield, nmol'] = self.synTab['support_amount, mg'] * self.synTab['support_capacity, nM/mg']
        self.synTab['Yield%'] = self.synTab['Total amount, nmol'] * 100 / self.synTab['Max yield, nmol']
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
        self.get_molProperties()
        self.get_yieldParams()

        print(self.synTab)
        self.synTab.reset_index(inplace=True)
        self.add_sheet_to_excel(self.synTab, 'test')

    def add_sheet_to_excel(self, df, sheetname):
        with pd.ExcelWriter(self.filePath, mode="a", if_sheet_exists='overlay') as f:
            df.to_excel(f, sheet_name=sheetname, index=False, header=False)




def test1():
    map = OligoMapReadWriter('/home/alex/Documents/OligoSynt/DB/maps/nmap/synth_180422_6opt.xlsx')


if __name__ == '__main__':
    test1()
