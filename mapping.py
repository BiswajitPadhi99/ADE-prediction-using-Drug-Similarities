import pickle
import pandas as pd


class globalVaribales:
    def __init__(self):
        self.chem_sim_mat = None
        self.atc_sim_mat = None
        self.sw_sim_mat = None
        self.go_bp_sim_mat = None
        self.go_cc_sim_mat = None
        self.go_mf_sim_mat = None
        self.random_sim_mat_1 = None
        self.random_sim_mat_2 = None
        self.random_sim_mat_3 = None
        self.random_sim_mat_4 = None
        self.random_sim_mat_5 = None
        self.include_dbid_list = None
        self.similarity_matrix = None
        self.adrlist_freq = None
        self.sider_eval_pairs = None
        self.drug_list = None
        self.adr_list = None
        self.sider_adr_soc_pt_map = None
        self.id2drug = None
        self.drug2id = None
        self.id2adr = None
        self.adr2id = None
        
        
    def initialize(self, similarity):
    
        self.chem_sim_mat = pd.read_csv('./Similarity_Matrices/chem_similarity.csv', index_col=0)
        self.atc_sim_mat = pd.read_csv('./Similarity_Matrices/ATC_similarity.csv', index_col=0)
        self.sw_sim_mat = pd.read_csv('./Similarity_Matrices/SW_similarity.csv', index_col=0)
        self.go_bp_sim_mat = pd.read_csv('./Similarity_Matrices/Go_BP_similarity.csv', index_col=0)
        self.go_cc_sim_mat = pd.read_csv('./Similarity_Matrices/Go_CC_similarity.csv', index_col=0)
        self.go_mf_sim_mat = pd.read_csv('./Similarity_Matrices/Go_MF_similarity.csv', index_col=0)
        
        
        #chem_sim_sider = pd.read_csv('chem_sim_sider_1to1.csv', index_col=0)
        #exclude_drug = [col for col in chem_sim_sider.columns if chem_sim_sider[col].sum() == 0]


        if(similarity == 'chem'):
            self.similarity_matrix = self.chem_sim_mat
        elif(similarity == 'atc'):
            self.similarity_matrix = self.atc_sim_mat
        elif(similarity == 'sw'):
            self.similarity_matrix = self.sw_sim_mat
        elif(similarity == 'go_bp'):
            self.similarity_matrix = self.go_bp_sim_mat
        elif(similarity == 'go_cc'):
            self.similarity_matrix = self.go_cc_sim_mat
        elif(similarity == 'go_mf'):
            self.similarity_matrix = self.go_mf_sim_mat

        
        
        self.include_dbid_list = list(self.similarity_matrix.columns)
        
        self.sider_eval_pairs = pickle.load(open('./files/sider_sim_pkl_1to1.pkl', 'rb'))
        
        self.common_dbid_list = pickle.load(open('./files/common_dbid_list.pkl', 'rb'))
        
        
        
        self.adrlist_freq = pickle.load(open('./files/adrlist_freq_22q2.pkl', 'rb'))[:1000]
        self.sider_eval_pairs = [(drug, adr.lower()) for (drug, adr) in self.sider_eval_pairs if adr in self.adrlist_freq and drug in self.include_dbid_list] #and drug in self.include_dbid_list] #and drug in self.common_dbid_list]
        
        self.drug_list = sorted(list(set(drug for (drug, adr) in self.sider_eval_pairs)))
        self.adr_list = sorted(list(set(adr for (drug, adr) in self.sider_eval_pairs)))
        
        self.sider_adr_soc_pt_map = pickle.load(open('./files/sider_adr_pt_to_soc.pkl','rb'))
        
        self.similarity_matrix = self.similarity_matrix.loc[self.drug_list, self.drug_list]
        
        print('length of common dbid list: ', len(self.common_dbid_list))
        print('length of sider dbid list: ', len(self.drug_list))
        print('length of sider adr list: ', len(self.adr_list))
        
        self.id2drug = {i: drug for i, drug in enumerate(self.drug_list)}
        self.drug2id = {drug: i for i, drug in enumerate(self.drug_list)}
        
        self.id2adr = {i: adr for i, adr in enumerate(self.adr_list)}
        self.adr2id = {adr: i    for i, adr in enumerate(self.adr_list)}
        
        #chem_sim_sider = chem_sim_sider.loc[drug_list, drug_list]
        #print(chem_sim_sider.shape)


global_variables = globalVaribales()


