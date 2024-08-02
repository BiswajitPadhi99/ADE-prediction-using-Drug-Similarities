import numpy as np
import pickle
from collections import defaultdict

from Eval import Eval
from mapping import global_variables
from utils import split_data
from similarity import get_Jaccard_Similarity, matrix_normalize


class Model:
    def __init__(self, metrics, similarity):
        self.ALPHA = 0.1
        self.metrics = metrics
        self.similarity = similarity
        
    def get_dbid_similarity_matrix(self, similarity_matrix):
        if similarity_matrix.shape[0] == similarity_matrix.shape[1]:
            np.fill_diagonal(similarity_matrix.values, 0)
            similarity_matrix = matrix_normalize(np.matrix(similarity_matrix))
            print("Similarity matrix shape: ", similarity_matrix.shape)   
            
        return similarity_matrix

    def get_similarity_matrix(self, X):
        features_matrix = []
        for idx in range(X.shape[0]):
            drug = global_variables.id2drug.get(idx)
            rxnorm = drugid2rxnorm[drug]
            features = rxnorm2features[rxnorm]
            features_matrix.append(features)
        features_matrix = np.asarray(features_matrix)
        return get_Jaccard_Similarity(features_matrix)

    def label_propogation(self, X, alpha):
        # similarity_matrix = self.get_similarity_matrix(X)
        similarity_matrix = self.get_dbid_similarity_matrix(global_variables.similarity_matrix)
        #np.set_printoptions(threshold=np.inf)
        #print(self.similarity)
        #print("\n\nSimilarity Matrix: ")
        #print(similarity_matrix)
        #print("\n\nX: ")
        #print(X)
        score_matrix_drug = (1 - alpha) * np.matmul(np.linalg.pinv(
            np.eye(np.shape(X)[0]) - alpha * similarity_matrix), X)
        #print("\n\nScore matrix: ")
        #print(score_matrix_drug)
        return score_matrix_drug

    def validate(self, X, Y, idx):
        AUC = []
        for i in range(1, 10):
            alpha = i * 0.1
            Y_pred = self.predict(X, alpha)
            metrics = self.eval(Y_pred, Y, idx)
            auc = metrics[0]
            AUC.append(auc)
        print(AUC)
        max_auc = max(AUC)
        max_idx = AUC.index(max_auc)
        max_alpha = (max_idx + 1) * 0.1
        self.ALPHA = max_alpha

    def predict(self, X, alpha):
        Y_pred = self.label_propogation(X, alpha)
        return Y_pred

    def eval(self, Y_pred, Y, idx):
        y_pred, y_gold = [], []
        for r, c in zip(idx[0], idx[1]):
            y_pred.append(Y_pred[r, c])
            y_gold.append(Y[r, c])
        ev = Eval(y_pred, y_gold)
        return ev.Metrics(self.metrics)


    def eval_DME(self, Y_pred, Y, idx, DME):
        y_pred, y_gold = defaultdict(list), defaultdict(list)
        for r, c in zip(idx[0], idx[1]):
            adrid = global_variables.id2adr.get(c)
            if adrid in DME:
                y_pred[adrid].append(Y_pred[r, c])
                y_gold[adrid].append(Y[r, c])
        EV = {}
        for k in y_pred.keys():
            y_p, y_g = y_pred.get(k), y_gold.get(k)
            ev = Eval(y_p, y_g)
            EV[k] = ev.Metrics(self.metrics)
        return EV

    def eval_SOC_class(self, Y_pred, Y, idx):
        soc_class_eval_metrics = {}
        for soc_class, soc_adr_list in global_variables.sider_adr_soc_pt_map.items():
            y_pred, y_gold = [], []
            for r, c in zip(idx[0], idx[1]):
                if global_variables.id2adr.get(c) in soc_adr_list:
                    y_pred.append(Y_pred[r, c])
                    y_gold.append(Y[r, c])
            try:
                ev = Eval(y_pred, y_gold)
                soc_class_eval_metrics[soc_class] =  {'status': 'Success', 'metrics': ev.Metrics(self.metrics)}
            except:
                soc_class_eval_metrics[soc_class] =  {'status': 'Fail', 'metrics': ev.Metrics(self.metrics)}
                continue
        return soc_class_eval_metrics


