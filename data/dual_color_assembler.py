import numpy as np
import pandas as pd
from os import listdir
import scipy.io as sio


def main():
    folders = ['ms_collected']




    for folder in folders:

        names = listdir(folder)

        if '.DS_Store' in names:
            names.remove('.DS_Store')

        names = [folder + '/' + name for name in names]

        for name in names:
            mat_contents = sio.loadmat(name)
            df2 = pd.DataFrame(mat_contents['data'], columns=[x[0] for x in mat_contents['def'][0]])
            df2['Strain'] = mat_contents['Experiment'][0][0][0]


            if 'df1' in locals():

                print(df2['Strain'][0])
                df1=pd.concat([df1,df2], ignore_index=True)
                #df1 = df1.append(df2, ignore_index=True)
            else:
                df1 = df2

    df1.to_csv('dual_megaframe.csv')
if __name__ == "__main__":
    main()
