import numpy as np
import pandas as pd
from os import listdir
import scipy.io as sio


def main():
    folders = ['sgadxk_cipro_2mg-L','srobk_cipro_2mg-L','smk_cipro_2mg-L']


    for folder in folders:

        names = listdir(folder)

        if '.DS_Store' in names:
            names.remove('.DS_Store')

        names = [folder + '/' + name for name in names]

        for name in names:
            mat_contents = sio.loadmat(name)
            df2 = pd.DataFrame(mat_contents['data'], columns=[x[0] for x in mat_contents['def'][0]])
            df2['Strain'] = mat_contents['Experiment'][0][0][0]

            if len(df2.columns) >102:
                temp=df2['Fluor3 mean']
                df2 = df2.drop(columns=['Fluor3 mean death', 'Fluor3 sum death', 'Fluor3 mean', 'Fluor3 sum'])
                df2['Fluor3 mean']=temp
            else:
                df2['Fluor3 mean'] = np.nan
            if 'df1' in locals():

                print(df2['Strain'][0])
                df1=pd.concat([df1,df2], ignore_index=True)
                #df1 = df1.append(df2, ignore_index=True)
            else:
                df1 = df2

    df1.to_csv('ciproframe.csv')
if __name__ == "__main__":
    main()
