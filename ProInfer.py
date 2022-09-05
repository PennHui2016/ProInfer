import pandas as pd
import numpy as np
import random
import os
import copy
from scipy import stats
import csv
import sys

def fdrsToQvals(fdrs):
  qvals = [0] * len(fdrs)
  qvals[len(fdrs)-1] = fdrs[-1]
  for i in range(len(fdrs)-2, -1, -1):
    qvals[i] = min(qvals[i+1], fdrs[i])
  return qvals

def fdr_2qvalue(score, label):
    ## input score is higher better
    ## fdr = (1+D)/(T+1)
    ##  call fdrsToQvals(fdrs)
    fdrs = [(len(np.where((label[np.where((score>=score[i]))[0]]==-1))[0])+1)/(len(np.where((score>=score[i]))[0])-len(np.where((label[np.where((score>=score[i]))[0]]==-1))[0])+1) for i in range(len(score))]
    fdrs = np.array(fdrs)
    idx = [i for i in range(len(score))]
    idx = np.array([idx]).T
    fdrs = np.array([fdrs]).T
    score = np.array([score]).T
    ifcb = np.column_stack((score, fdrs, idx))
    ifcb = ifcb[np.argsort(-ifcb[:, 0]), :]
    q_value = fdrsToQvals(ifcb[:, 1])
    ifcb[:, 1] = q_value
    ifcb = ifcb[np.argsort(ifcb[:, 2]), :]
    report_qvalue = ifcb[:, 1]
    return report_qvalue

def obtain_fasta(path, name):
    filepath = path # 'F:/NTU/proteomics/Human_database_including_decoys_(cRAP_added).fasta'
    fasta = {}
    with open(filepath) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = ''
                continue
            sequence = line
            fasta[active_sequence_name] += sequence

    protein_names = []
    protein_seqs = []
    protein_len = []
    for key in fasta.keys():
        strs = key.split(' ')
        protein_names.append(strs[0])
        protein_len.append(len(fasta[key]))
        protein_seqs.append(fasta[key])

    out = np.column_stack((np.array([protein_names]).T, np.array([protein_len]).T))
    return fasta, out

def complex_process(cpx_path, species):
    complex = pd.read_csv(cpx_path, sep='\t', header=0).values
    hu_complex = complex[np.where((complex[:, 2] == species))[0], 5]
    cpx_dict = {}
    pro_cpx_dict = {}
    for i in range(len(hu_complex)):
        pros_i = hu_complex[i].split(';')
        cpx_dict.update({str(i): pros_i})
        for j in pros_i:
            if j in pro_cpx_dict.keys():
                new_value = np.union1d(pro_cpx_dict[j], [i])
            else:
                new_value = [i]
            pro_cpx_dict.update({j: new_value})

    return cpx_dict, pro_cpx_dict

def read_store_dict_peptide(file_name, threshold):
    score_idx = 0
    #score_type_idx = 0
    sequence_idx = 0
    pep_pro_idx = 0
    dict_pep = {}
    dict_pro = {}
    pep_scores = {}
    f = open(file_name, 'r')
    lines = f.readlines()
    for line in lines:
        if line[0:8] == '#PEPTIDE':
            ss = line.replace('\n', '').replace('\r', '').split('\t')
            titles = ss
            score_idx = titles.index('score')
            score_type_idx = titles.index('score_type')
            sequence_idx = titles.index('sequence')
            pep_pro_idx = titles.index('accessions')
        elif line[0:7] == 'PEPTIDE':
            ss = line.replace('\n', '').replace('\r', '').split('\t')
            peptide = ss[sequence_idx].replace('(Oxidation)', '').replace('(Carbamidomethyl)', '').replace(
                '.(Glu->pyro-Glu)', '').replace('.(Gln->pyro-Glu)', '').replace('.(Acetyl)', '').replace('.(Ammonia-loss)', '')
            score = float(ss[score_idx])
            if float(score) > threshold:
                continue
            pep_pro = ss[pep_pro_idx].split(';')
            if peptide in dict_pep.keys(): # already exist
                new_value = np.union1d(pep_pro, dict_pep[peptide])
            else:
                new_value = pep_pro

            dict_pep.update({peptide: new_value})
            pep_scores.update({peptide: score})
            for j in pep_pro:
                if j in dict_pro.keys():
                    value = dict_pro[j]
                    if peptide in value.keys():
                        new_value = min(value[peptide], score)
                    else:
                        new_value = score
                    dict_pro[j].update({peptide: new_value})
                else:
                    dict_pro.update({j: {peptide: score}})

    return dict_pro, dict_pep, pep_scores

def accPEP_v1(dict_pro, dict_pep, proteins_info, decoy):
    pros = []
    scores = []
    lg_s = []
    label_pro = []
    gene_dict = {}
    for pro in dict_pro.keys():
        pros.append(pro)
        peps = dict_pro[pro]
        score = 1
        for pep in peps.keys():
            if len(dict_pep[pep])==1: # unique
                score *= dict_pro[pro][pep]
            elif len(dict_pep[pep]) > 1: # ambigious
                pro_j = dict_pep[pep]
                pro_j_len = [proteins_info[:, 1][np.where((proteins_info[:, 0] == pr))[0][0]] for pr in
                             np.unique(pro_j)]
                pro_j_len = np.array(pro_j_len, dtype='float')
                nor_pro_j_p = np.sum(pro_j_len) / np.array(pro_j_len)
                nor_pro_j_p = nor_pro_j_p / np.sum(nor_pro_j_p)
                p_from_pro = nor_pro_j_p[np.where((np.unique(pro_j) == pro))[0][0]]
                score *= 1 - (1 - dict_pro[pro][pep]) * p_from_pro

        scores.append(score)
        lg_s.append(-10 * np.log10(score + 0.00000000000001))
        gene_dict.update({pro: score})
        if pro[0:3] == decoy:
            label_pro.append(-1)
        else:
            label_pro.append(1)

    label_pro = np.array(label_pro)
    qvalue = fdr_2qvalue(np.array(lg_s), label_pro)

    return pros, scores, lg_s, qvalue, label_pro, gene_dict

def ProInfer_v1(per_pep_path, threshold, protein_database, save_path, decoy):
    dict_pro, dict_pep, pep_scores = read_store_dict_peptide(per_pep_path, threshold)
    fasta, proteins_info = obtain_fasta(protein_database, '')
    pros, scores, lg_s, qvalue, label_pro, gene_dict = accPEP_v1(dict_pro, dict_pep, proteins_info, decoy)
    out_all = np.column_stack((np.array([pros]).T, np.array([scores]).T, np.array([lg_s]).T, np.array([qvalue]).T,
                               np.array([label_pro]).T))
    with open(save_path, 'w', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(out_all)

    return out_all, gene_dict

def ProInfer_cpx_v2(per_pep_path, save_path_proinfer, save_path_cpx, psm_threshold, pro_qvalue_td, protein_database, complex_path, species, decoy):
    # get ProInfer results [id, accPEP score, lg10(accPEP), q-value, label]
    # gene_dict in the format of {protein: accPEP score}
    res_proinfer, gene_dict = ProInfer_v1(per_pep_path, psm_threshold, protein_database, save_path_proinfer, decoy)

    label_pro = np.array(res_proinfer[:, 4], dtype='int') # target label 1, decoy label -1
    pros = res_proinfer[:, 0] # candidate proteins
    qvalue = np.array(res_proinfer[:, 3], dtype='float') # ProInfer calculated q-values for candidate proteins

    # get complex networks
    # cpx_dict in format of {complex_id: [protein list]}
    # pro_cpx_dict in format of {protein_id: [complex ids]}
    cpx_dict, pro_cpx_dict = complex_process(complex_path, species)

    # get target proteins "sp|xxxx|xxxx_HUMAN"
    target_pros = np.array(pros)[np.where((np.array(label_pro) == 1))[0]]
    # get uniport ids of target proteins for searching complex
    target_pros_id = [m.replace('KKA1_ECOLX_CONTAMINANT', 'KKA1_ECOLX_CONTAMINANT||').split('|')[1] for m in target_pros]
    decoy_pros = np.array(pros)[np.where((np.array(label_pro) == -1))[0]]
    decoy_pros_id = [n.replace('KKA1_ECOLX_CONTAMINANT', 'KKA1_ECOLX_CONTAMINANT||').split('|')[1] for n in decoy_pros]

    res_out = res_proinfer # initialize proinfer_cpx outputs
    num0 = len(np.where((np.array(qvalue)[np.where((np.array(label_pro) == 1))[0]] < pro_qvalue_td))[0]) # initialize the number of reported target proteins
    flag = 0
    update = 0
    while flag != 1:
        # get confidences of complex presenting in sample
        p_complex_target = {}
        p_complex_decoy = {}
        for cpx in cpx_dict.keys():
            cpx_pros = cpx_dict[cpx]
            cpx_pros_target, idx_t1, idx_t2 = np.intersect1d(cpx_pros, target_pros_id, return_indices=True)
            cpx_pros_decoy, idx_d1, idx_d2 = np.intersect1d(cpx_pros, decoy_pros_id, return_indices=True)

            # probability of a complex present is calculated by the mean accPEP scores of its proteins
            # smaller better (in posterior error probability)
            if len(idx_t2) > 0:
                cpx_pros_targets = target_pros[idx_t2]
                complex_pros_ps_target = [gene_dict[t] for t in cpx_pros_targets]
                #p_complex_target.update({cpx: np.mean(complex_pros_ps_target)})
                p_complex_target.update({cpx: np.min(complex_pros_ps_target)})
            else:
                p_complex_target.update({cpx: 1})

            if len(idx_d2) > 0:
                cpx_pros_decoys = decoy_pros[idx_d2]
                complex_pros_ps_decoy = [gene_dict[t] for t in cpx_pros_decoys]
                #p_complex_decoy.update({cpx: np.mean(complex_pros_ps_decoy)})
                p_complex_decoy.update({cpx: np.min(complex_pros_ps_decoy)})
            else:
                p_complex_decoy.update({cpx: 1})

        # update probability of a protein present with min(ProInfer_accPEP, min(complexes' PEP)), smaller better
        res_complex_p = []
        res_complex_s = []

        for i in range(len(pros)):
            pro = pros[i]
            pro_id = pro.replace('KKA1_ECOLX_CONTAMINANT', 'KKA1_ECOLX_CONTAMINANT||').split('|')[1]
            label = label_pro[i]
            proinfer_score = gene_dict[pro]
            if pro_id in pro_cpx_dict.keys():
                pro_cpx = pro_cpx_dict[pro_id]
                if label == 1:
                    cpx_score = np.min([p_complex_target[str(x)] for x in pro_cpx])  # smaller the better
                elif label == -1:
                    cpx_score = np.min([p_complex_decoy[str(x)] for x in pro_cpx])  # smaller the better

                if proinfer_score >= cpx_score: # update gene_dict
                    gene_dict.update({pro: cpx_score})
                    res_complex_p.append(cpx_score)
                    res_complex_s.append(-10 * np.log10(cpx_score + 0.00000000000001))
                else:
                    res_complex_p.append(proinfer_score)
                    res_complex_s.append(-10 * np.log10(proinfer_score + 0.00000000000001))
            else:
                res_complex_p.append(proinfer_score)
                res_complex_s.append(-10 * np.log10(proinfer_score + 0.00000000000001))
                continue


        qvalue = fdr_2qvalue(np.array(res_complex_s), label_pro)  # psm_lfdr(np.array(lg_s), label_pro)

        num = len(np.where((np.array(qvalue)[np.where((np.array(label_pro) == 1))[0]] < pro_qvalue_td))[0]) # see whether more target proteins can be inferred
        # if num < (num0 + 10):
        if num <= num0:
            flag = 1
        else:
            num0 = num
            res_out = np.column_stack((res_out, np.array([res_complex_p]).T,
                                       np.array([res_complex_s]).T, np.array(qvalue)))
            update = update + 1

    with open(save_path_cpx, 'w', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(res_out)

    return res_out, update, res_proinfer

def main():
    # Training settings
    parser = argparse.ArgumentParser(description='run ProInfer and ProInfer_cpx')
    parser.add_argument('-i', '--input', default='./DDA1.tsv',
                        help='input peptide list from percolator, can obtain with the attached KNIME workflow in [preparing_peptides_workflow.knwf]')
    parser.add_argument('-t', '--run_type', type=int, default=1,
                        help='1 (default)--runs the ProInfer; 2--runs the ProInfer_cpx')
    parser.add_argument('-pt', '--psm_threshold', type=float, default=0.999,
                        help='Selection of posterior error probability (PEP) threshold ((0, 1], default 0.999) to filter low confidence peptides out.')
    parser.add_argument('-qt', '--pro_qvalue_td', type=float, default=0.01,
                        help='float (0, 1], default 0.01. Selection of False Discovery Rate (FDR) for reporting target proteins.')
    parser.add_argument('-sp','--save_path_proinfer', default='./res/proinfer_out.csv',
                        help='Indicate where to save the results from ProInfer')
    parser.add_argument('-spc','--save_path_cpx',default='./res/proinfer_cpx_out',
                        help='Indicate where to save the results from ProInfer')
    parser.add_argument('-db', '--protein_database', default='./2022-06-23-decoys-contam-uniprot-proteome_UP000005640_2022_5_5.fasta',
                        help='Indicate the protein database for database search. Should be the same as what you used in the KNIME workflow')
    parser.add_argument('-cp', '--complex_path', default='./allComplexes.txt',
                        help='The protein complexes used by ProInfer_cpx. The default file is downloaded from CORUM 3.0 (http://mips.helmholtz-muenchen.de/corum/).')
    args = parser.parse_args()

    per_pep_path = args.input # all_params[1]
    run_type = args.run_type #int(all_params[2])
    psm_threshold = args.psm_threshold # float(all_params[3])
    pro_qvalue_td = args.pro_qvalue_td # float(all_params[4])
    save_path_proinfer = args.save_path_proinfer # all_params[5]
    save_path_cpx = args.save_path_cpx # all_params[6]
    protein_database = args.protein_database # all_params[7]
    complex_path = args.complex_path # all_params[8]

    if run_type == 1:  # only run ProInfer
        ProInfer_v1(per_pep_path, psm_threshold, protein_database, save_path_proinfer)
    elif run_type == 2:  # run both ProInfer and ProInfer_cpx
        ProInfer_cpx_v2(per_pep_path, save_path_proinfer, save_path_cpx, psm_threshold, pro_qvalue_td, protein_database,
                        complex_path)

if __name__ == '__main__':
    # toy data Hela DDA1 with PSM pp>0.001
    # per_pep_path --- peptide list from percolator in .tsv format
    # save_path_proinfer --- where to save proinfer's outputs
    # save_path_cpx --- where to save proinfer_cpx's outputs
    # psm_threshold --- threshold filtering PSMs (peptides)
    # pro_qvalue_td --- FDR threshold
    # protein_database --- protein database used for database searching
    # complex_path --- complex list
    main()