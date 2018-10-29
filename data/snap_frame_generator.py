import numpy as np
import pandas as pd
from os import listdir
import scipy.io as sio


def main():
    folders = ['acrab_dnaq_soxs_lauv5_fis_rpst_marrab_purA_hdea_sodA_tolc_inaA_10ms','megaframe.csv']
    name_list=['acrAB','dnaQ','SoxS','lacuv5','Fis','rpsT','marRAB','purA','hdeA','sodA','tolC','InaA']

    second_list=['ompF', 'gadX',  'crp', 'sigma70', 'rpsT', 'rrnbp1', 'rob',  'tonB', 'CodB']

    for j,folder in enumerate(folders):
        if j==0:
            names = listdir(folder)

            if '.DS_Store' in names:
                names.remove('.DS_Store')

            names = [folder + '/' + name for name in names]

            for u,name in enumerate(names):
                mat_contents = sio.loadmat(name)
                df2 = pd.DataFrame(mat_contents['data'], columns=[x[0] for x in mat_contents['def'][0]])
                df2['Strain'] = name_list[int(np.floor(u/5))]

                # if len(df2.columns) >102:
                #     temp=df2['Fluor3 mean']
                #     df2 = df2.drop(columns=['Fluor3 mean death', 'Fluor3 sum death', 'Fluor3 mean', 'Fluor3 sum'])
                #     df2['Fluor3 mean']=temp
                # else:
                #     df2['Fluor3 mean'] = np.nan
                if 'df1' in locals():

                    print(df2['Strain'][0])
                    df1=pd.concat([df1,df2], ignore_index=True)
                    #df1 = df1.append(df2, ignore_index=True)
                else:
                    df1 = df2
        else:
            df2 = pd.read_csv('megaframe.csv')
            df2=df2[df2['Strain'].isin(second_list)]
            df1 = pd.concat([df1.ix[:, ['Strain','Fluor1 mean']], df2.ix[:, ['Strain','Fluor1 mean']]], ignore_index=True)

    df1.to_csv('snapframe.csv')


if __name__ == "__main__":
    main()
