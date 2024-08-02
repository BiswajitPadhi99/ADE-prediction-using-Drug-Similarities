import os
import numpy as np
from tqdm import tqdm
import pandas as pd

from mapping import global_variables 


class FAERSdata:
    def __init__(self, directory, method, year, quarter):

        Files = sorted(os.listdir(directory))

        file_name = f"{method}_Signal_Score_{int(year):02d}Q{quarter}.csv"
        
        if(method=="PRR"):
            score_col = "prr_score"
        elif(method=="ROR"):
            score_col = "ror_score"
        elif(method=="MGPS"):
            score_col = "gps_score"
        elif(method=="BCPNN"):
            score_col = "bcpnn_score"
        
        Files = [file_name]
        print(Files)
        X = {}
        Y = {}
        for i in tqdm(range(len(Files))):
            f = Files[i]
            x = np.zeros(shape=(len(global_variables.drug_list), len(global_variables.adr_list)))
            data = pd.read_csv(f'{directory}/{f}', delimiter=' ', quotechar='"')
            data.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
            if(method=="BCPNN"):
                data[score_col] = data[score_col].clip(upper=data[score_col].quantile(0.95))
                data[score_col] = np.exp(data[score_col])
            for _, row in data.iterrows():
                drug, adr, score = row["drug_code"], row["ade_code"], float(row[score_col])
                if drug in global_variables.drug_list and adr in global_variables.adr_list:
                    drug_id, adr_id = global_variables.drug2id.get(drug), global_variables.adr2id.get(adr)
                    x[drug_id, adr_id] = score

            

            y = np.zeros(shape=(len(global_variables.drug_list), len(global_variables.adr_list)))
            for drug, adr in global_variables.sider_eval_pairs:
                drug_id, adr_id = global_variables.drug2id.get(drug), global_variables.adr2id.get(adr)
                y[drug_id, adr_id] = 1


            y = np.asarray(y)


            X[i] = x
            Y[i] = y
 

        self.X = X
        self.Y = Y
        self.Files=Files














