import numpy as np
import pandas as pd
from os import listdir
import scipy.io as sio


def main():
    folders=['marrab_collected','gadx_collected','crp_collected','ompf_collected','sigma_collected','soxs_collected','purA_collected','hdea_collected','rob_collected']
    for folder in folders:

        names = listdir(folder)

        if '.DS_Store' in names:
            names.remove('.DS_Store')

        names = [folder + '/' + name for name in names]


        for name in names:
            mat_contents = sio.loadmat(name)
            df2 = pd.DataFrame(mat_contents['data'],columns=[x[0] for x in mat_contents['def'][0]])
            df2['Strain']=mat_contents['Experiment'][0][0][0]
            if 'df1' in locals():
                df1 = pd.concat([df1,df2])
            else:
                df1=df2

    df1.to_csv('megaframe.csv')
if __name__ == "__main__":
    main()
