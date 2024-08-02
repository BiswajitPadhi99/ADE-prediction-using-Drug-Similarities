from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import numpy as np
import pandas as pd

from FAERSdata import FAERSdata
from Model import Model
from mapping import global_variables
from utils import split_data, sample_zeros
import random


def parse_args():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, conflict_handler='resolve')
    parser.add_argument('--input', default='signal_0513/', help='Input original signal scores file.')
    parser.add_argument('--method', default='prr', help='Signal detection algorithm')
    parser.add_argument('--year', default=4, help='Years of data used for model')
    parser.add_argument('--quarter', default=1, help='Years of data used for model')
    parser.add_argument('--eval_metrics', default=all, choices=['all', 'specificity-sensitivity'],
                        help='Evaluation metrics')
    parser.add_argument('--similarity', choices=['chem', 'atc', 'sw', 'go_bp', 'go_cc', 'go_mf', 'random1', 'random2', 'random3', 'random4', 'random5'], help='Type of similarity to use')
    parser.add_argument('--split', type=bool, default=True)
    parser.add_argument('--output', default='results/old_new_0513_prr_val5test95_sider_eval_pairs_final.csv')
    parser.add_argument('--soc_output')
    
    args = parser.parse_args()
    return args


def pretty_print_eval(res, metrics):
    if metrics == 'all':
        print('All metrics: ' + ','.join(np.round(res,3).astype(str)))
    else:
        print('fixed_sensitivity: ' + ','.join(np.round(res[1],3).astype(str)))
        print('fixed_specificity: ' + ','.join(np.round(res[2],3).astype(str)))
        
def save_drug_adr_scores(scores,path):
    idx=scores.shape
    drug_list = []
    ade_list = []
    for id in range(idx[0]):
        drug_list.append(global_variables.id2drug.get(id))
    for id in range(idx[1]):
        ade_list.append(global_variables.id2adr.get(id))
    drug_adr_df = pd.DataFrame(scores, index=drug_list, columns=ade_list)
    drug_adr_df.to_csv(path)

def set_seed(seed):
    np.random.seed(seed)
    random.seed(seed)


def main(args):

    set_seed(42)
    
    global_variables.initialize(args.similarity)

    print('#' * 50)
    print('Signal Detection Algorithm: {}'.format(args.method))
    print('#' * 50)


    data = FAERSdata(args.input, args.method, args.year, args.quarter)
    
    out = open(args.output, 'w')
    out.write('LP_auc,baseline_auc,LP_aupr,baseline_aupr,LP_f1,baseline_f1\n')
    
    sider_adr_soc_pt_map = global_variables.sider_adr_soc_pt_map
    out_soc = open(args.soc_output, 'w')
    for soc_class in sider_adr_soc_pt_map.keys():
        out_soc.write(f'LP_{soc_class}_AUC,baseline_{soc_class}_AUC,')
    out_soc.write('\n')
    

    for i in range(len(data.X.keys())):
        X, Y = data.X.get(i), data.Y.get(i)
        # all_idx = np.where(Y > -1)
        all_idx = sample_zeros(Y)
        if args.split:
            test, valid = split_data(Y)
            model = Model(args.eval_metrics,args.similarity)
            model.validate(X, Y, valid)
            Y_pred = model.predict(X, model.ALPHA)
            
            # save_drug_adr_scores(X,f'enhanced_signal_scores/{args.method}/original/{args.similarity}/{data.Files[i]}')
            # save_drug_adr_scores(Y_pred,f'enhanced_signal_scores/{args.method}/lp/{args.similarity}/{data.Files[i]}')
            
            valid_res = model.eval(Y_pred, Y, valid)
            test_res = model.eval(Y_pred, Y, test)
            lp_test_soc_class_res = model.eval_SOC_class(Y_pred, Y, test)
            print('LP-{}:'.format(args.method))
            print('alpha: {}'.format(model.ALPHA))
            print('valid:')
            pretty_print_eval(valid_res, args.eval_metrics)
            print('test:')
            pretty_print_eval(test_res, args.eval_metrics)
            print('soc_test:')
            print(lp_test_soc_class_res)
            lp_test_auc = test_res[0]
            lp_test_aupr = test_res[1]
            lp_test_f1 = test_res[-1]

            valid_res = model.eval(X, Y, valid)
            test_res = model.eval(X, Y, test)
            baseline_test_soc_class_res = model.eval_SOC_class(X, Y, test)
            print('baseline-{}:'.format(args.method))
            print('valid:')
            pretty_print_eval(valid_res, args.eval_metrics)
            print('test:')
            pretty_print_eval(test_res, args.eval_metrics)
            print('soc_test:')
            print(baseline_test_soc_class_res)
            baseline_test_auc = test_res[0]
            baseline_test_aupr = test_res[1]
            baseline_test_f1 = test_res[-1]
            out.write(f'{lp_test_auc},{baseline_test_auc},{lp_test_aupr},{baseline_test_aupr},{lp_test_f1},{baseline_test_f1} \n')
            
            for soc_class in sider_adr_soc_pt_map.keys():
                if lp_test_soc_class_res[soc_class]['status'] == 'Success':
                    soc_lp_res = lp_test_soc_class_res[soc_class]['metrics'][0]
                    out_soc.write(f'{soc_lp_res},')
                else:
                    out_soc.write(f'AUC calc. error,')
                if baseline_test_soc_class_res[soc_class]['status'] == 'Success' :
                    soc_baseline_res = baseline_test_soc_class_res[soc_class]['metrics'][0]
                    out_soc.write(f'{soc_baseline_res},')
                else:
                    out_soc.write(f'AUC calc. error,')
            out_soc.write('\n')
            
        else:
            model = Model(args.eval_metrics)
            model.validate(X, Y, all_idx)
            Y_pred = model.predict(X, model.ALPHA)
            res = model.eval(Y_pred, Y, all_idx)
            lp_test_soc_class_res = model.eval_SOC_class(Y_pred, Y, all_idx)
            print('LP-{}:'.format(args.method))
            pretty_print_eval(res, args.eval_metrics)
            lp_res_auc = res[0]
            lp_res_aupr = res[1]
            lp_res_f1 = res[-1]

            print('baseline-{}:'.format(args.method))
            res = model.eval(X, Y, all_idx)
            baseline_test_soc_class_res = model.eval_SOC_class(X, Y, all_idx)
            pretty_print_eval(res, args.eval_metrics)
            baseline_res_auc = res[0]
            baseline_res_aupr = res[1]
            baseline_res_f1 = res[-1]
            
            out.write(f'{lp_res_auc},{baseline_res_auc},{lp_res_aupr},{baseline_res_aupr},{lp_res_f1},{baseline_res_f1}\n')
            
            for soc_class in sider_adr_soc_pt_map.keys():
                if lp_test_soc_class_res[soc_class]['status'] == 'Success':
                    soc_lp_res = lp_test_soc_class_res[soc_class]['metrics'][0]
                    out_soc.write(f'{soc_lp_res},')
                else:
                    out_soc.write(f'AUC calc. error,')
                if baseline_test_soc_class_res[soc_class]['status'] == 'Success' :
                    soc_baseline_res = baseline_test_soc_class_res[soc_class]['metrics'][0]
                    out_soc.write(f'{soc_baseline_res},')
                else:
                    out_soc.write(f'AUC calc. error,')
            out_soc.write('\n')

def more_main():
    args = parse_args()
    main(args)

if __name__ == '__main__':
    more_main()
