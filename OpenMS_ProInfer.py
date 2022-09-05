import pandas as pd
import numpy as np
import random
import os
import copy
import ProInfer_core_v3
import re
import glob
import csv
import sys

def get_para_string(params):
    paraString = ''
    for i in range(len(params.keys())):#para in params.keys():
        para = list(params.keys())[i]
        if i!= len(params.keys())-1:
            if str(params[para]) == '':
                paraString += para + ' '
            else:
                paraString += para + ' ' + str(params[para]) + ' '
        else:
            paraString += para + ' ' + str(params[para])
    return paraString

def databaseSearch(tool, openms_path, msfragger_path, file, save_file, database, paras_tool):
    '''
    :param tool: database search tool, currently test MSFragger
    :param openms_path: OpenMs install folder
    :param msfragger_path: MSFragger path, .jar
    :param file: MS file, .mzML
    :param database: database, .fasta including decoy and contaminates
    :param paras: MSFragger parameters
    :paras_tool paras: save_file: output file path, .idXML
    :return: output file path, .idXML
    '''

    if paras_tool != '':
        paraString = get_para_string(paras_tool)

    if tool == 'MSFragger':
        os.system(openms_path + 'MSFraggerAdapter -license yes -executable ' + msfragger_path + ' -in ' + file +
                  ' -out ' + save_file + ' -database ' + database + ' ' + paraString)

    if os.path.exists(save_file):
        return save_file
    else:
        print('database search has been stopped')
        return ''

def peptideindexer(openms_path, in_file, out_file, database, para_pepidx):
    if para_pepidx != '':
        paraString = get_para_string(para_pepidx)

    print(openms_path + 'PeptideIndexer -in ' + in_file + ' -fasta ' + database + paraString + ' -out ' +
              out_file + ' ' + paraString)

    os.system(openms_path + 'PeptideIndexer -in ' + in_file + ' -fasta ' + database + paraString + ' -out ' +
              out_file + ' ' + paraString)

    if os.path.exists(out_file):
        return out_file
    else:
        print('peptideindexer search has been stopped')
        return ''

def PSMFeatureExtractor(openms_path, in_file, out_file, para_psmefea):
    if para_psmefea != '':
        paraString = get_para_string(para_psmefea)

    os.system(openms_path + 'PSMFeatureExtractor -in ' + in_file + ' -out ' +
              out_file + ' ' + paraString)

    if os.path.exists(out_file):
        #out_pin = out_file.replace('.idXML', '.tsv')
        #idXML2tsv(out_file, out_pin)
        return out_file
    else:
        print('PSMFeatureExtractor search has been stopped')
        return ''

def idXML2tsv(openms_path, in_file, out_file, para_i2t):
    # para_i2t = {'-id:proteins_only': 'false', ' -id:peptides_only': 'false', '-id:protein_groups': 'true',
    #             '-id:add_metavalues': 100, '-id:add_hit_metavalues': 100, '-id:add_protein_hit_metavalues': 100}
    if para_i2t != '':
        paraString = get_para_string(para_i2t)

    print(openms_path + 'TextExporter -in ' + in_file + ' -out ' +
              out_file + ' ' + paraString)
    os.system(openms_path + 'TextExporter -in ' + in_file + ' -out ' +
              out_file + ' ' + paraString)

    if os.path.exists(out_file):
        return out_file
    else:
        print('file format transfer has been stopped')
        return ''

def Percolator(openms_path, in_file, out_file, para_percolator):
    if para_percolator != '':
        paraString = get_para_string(para_percolator)

    os.system(openms_path + 'PercolatorAdapter -in ' + in_file + ' -out ' + out_file + ' ' + paraString)

    if os.path.exists(out_file):
        return out_file
    else:
        print('Percolator search has been stopped')
        return ''

def faslediscoveryrate(openms_path, in_file, out_file, para_fdr):
    if para_fdr != '':
        paraString = get_para_string(para_fdr)

    os.system(openms_path + 'FalseDiscoveryRate -in ' + in_file + ' -out ' + out_file + ' ' + paraString)

    if os.path.exists(out_file):
        return out_file
    else:
        print('FalseDiscoveryRate has been stopped')
        return ''

def idscoreswitch(openms_path, in_file, out_file, para_scsw):
    if para_scsw != '':
        paraString = get_para_string(para_scsw)

    os.system(openms_path + 'IDScoreSwitcher -in ' + in_file + ' -out ' + out_file + ' ' + paraString)

    if os.path.exists(out_file):
        return out_file
    else:
        print('IDScoreSwitcher has been stopped')
        return ''

def prepare_psms(openms_path, in_file, out_file, para_fdr, para_scsw):
    ### using pep to conduct protein identification but filtering psm with q-value
    out_fdr = in_file.replace('_percolator.idXML', '_fdr.idXML')
    faslediscoveryrate(openms_path, in_file, out_fdr, para_fdr)
    in_scsw= out_fdr
    out_scsw = in_scsw.replace('_fdr.idXML', '_scsw.idXML')
    idscoreswitch(openms_path, in_scsw, out_file, para_scsw)

def IDfilter(openms_path, in_file, out_file, threshold):
    para_idf = {
        '-score:pep': threshold
    }
    if para_idf != '':
        paraString = get_para_string(para_idf)
    os.system(openms_path + 'IDFilter -in ' + in_file + ' -out ' + out_file + ' ' + paraString)

    if os.path.exists(out_file):
        return out_file
    else:
        print('IDFilter search has been stopped')
        return ''

def report_proteins(pros_all, labels, qvalues, fdr):
    ind5 = np.where((np.array(qvalues, dtype='float') < fdr))[0]
    s = []
    if len(ind5) > 0:
        ind15 = np.where((np.array(labels)[ind5] == 1))[0]
        pros = np.array(pros_all)[ind5][ind15]
        tys = np.array([pro.split('|')[2] for pro in pros])
        ind_no_contanminant = np.where((tys != '_CONTAMINANT'))[0]
        # num_epi.append(len(ind_no_contanminant))
        s = np.array(pros_all)[ind5][ind15][ind_no_contanminant]
    return s

def Epifany(openms_path, in_file, psm_threshold, pro_qvalue_td, ggr, reg, decoy, species, save_path, cell_line, real, name):
    if reg == 'true':
        para_epifany = {
            '-protein_fdr': 'true',
            '-greedy_group_resolution':ggr, # 'none', 'remove_associations_only','remove_proteins_wo_evidence'
            '-threads': 8,
            '-algorithm:psm_probability_cutoff':1-psm_threshold, #0.001
            '-algorithm:model_parameters:regularize': '' #'true'
        }
    else:
        para_epifany = {
        '-protein_fdr': 'true',
        '-greedy_group_resolution': ggr,  # 'none', 'remove_associations_only','remove_proteins_wo_evidence'
        '-threads': 8,
        '-algorithm:psm_probability_cutoff': 1-psm_threshold}

    if para_epifany != '':
        paraString = get_para_string(para_epifany)

    out_file = save_path + str(psm_threshold) + '_' + ggr + '_' + reg + '_epifany.idXML'
    #print(openms_path + 'Epifany -in ' + in_file + ' -out ' + out_file + ' ' + paraString)
    os.system(openms_path + 'Epifany -in ' + in_file + ' -out ' + out_file + ' ' + paraString)
    in_file_tsv = out_file
    out_tsv = in_file_tsv.replace('.idXML', '.tsv')  # 'F:/NTU/proteomics/ProInfer_codes/test_res/CAP1_proinfer.tsv'
    psms_para = {
        '-id:add_protein_hit_metavalues': 100
    }
    idXML2tsv(openms_path, in_file_tsv, out_tsv, psms_para)
    epi_out = ProInfer.read_openMS_tsv_protein(out_tsv, 'epi_fdr')

    return epi_out

def Fido(openms_path, in_file, psm_threshold, pro_qvalue_td, ggr, decoy, species, save_path, cell_line, real, name):
    if ggr == 'false':
        para_fido = {
            '-threads': 4
        }
    else:
        para_fido = {
            '-greedy_group_resolution': '',
            '-threads': 4
        }

    in_filter = in_file
    out_filter = save_path + str(psm_threshold) + '_' + ggr + '_fido.idXML'
    IDfilter(openms_path, in_filter, out_filter, psm_threshold)

    if para_fido != '':
        paraString = get_para_string(para_fido)

    out_file = out_filter.replace(ggr, '_filtered_' + ggr)

    os.system(openms_path + 'FidoAdapter -in ' + out_filter + ' -out ' + out_file + ' ' + paraString)
    in_file_tsv = out_file
    out_tsv = in_file_tsv.replace('.idXML', '.tsv')  # 'F:/NTU/proteomics/ProInfer_codes/test_res/CAP1_proinfer.tsv'
    psms_para = {
        '-id:add_protein_hit_metavalues': 100
    }
    idXML2tsv(openms_path, in_file_tsv, out_tsv, psms_para)
    fido_out = ProInfer.read_openMS_tsv_protein(out_tsv, 'epi_fdr')

    return fido_out

def percolator_protein(openms_path, in_file, psm_threshold, pro_qvalue_td, decoy, species, save_path, cell_line, real, name):
    in_filter = in_file
    out_filter = save_path + 'filter_' + str(psm_threshold) + '.idXML'
    IDfilter(openms_path, in_filter, out_filter, psm_threshold)

    in_per = out_filter
    out_per = in_per.replace('filter_', '_per_')
    para_percolator = {
        '-score_type': 'q-value',
        '-peptide_level_fdrs': '',
        '-protein_level_fdrs': '',
        '-threads': 4,
        '-decoy_pattern':decoy
    }

    Percolator(openms_path, in_per, out_per, para_percolator)

    in_file_tsv = out_per
    out_tsv = in_file_tsv.replace('.idXML', '.tsv')  # 'F:/NTU/proteomics/ProInfer_codes/test_res/CAP1_proinfer.tsv'
    psms_para = {
        '-id:add_protein_hit_metavalues': 100
    }
    idXML2tsv(openms_path, in_file_tsv, out_tsv, psms_para)
    per_out = ProInfer.read_openMS_tsv_protein(out_tsv, 'epi_fdr')

    return per_out

def dbs_pi(openms_path, msfragger_path, in_file, database):
    params_frag = {
        '-tolerance:precursor_mass_tolerance_lower': 10,
        '-tolerance:precursor_mass_tolerance_upper': 10,
        '-tolerance:precursor_mass_unit': 'ppm',
        '-tolerance:precursor_true_tolerance': 10,
        '-tolerance:precursor_true_unit': 'ppm',
        '-tolerance:fragment_mass_tolerance': 0.05,
        '-tolerance:fragment_mass_unit': 'Da',
        '-digest:search_enzyme_name': 'Trypsin',
        '-digest:min_length': 7,
        '-digest:max_length': 50,
        '-varmod:masses': 15.994900,
        '-varmod:syntaxes': 'M',
        '-varmod:max_variable_mods_per_mod': 3,
        '-spectrum:minimum_peaks': 15,
        '-spectrum:use_topn_peaks': 150,
        '-spectrum:minimum_ratio': 0.01,
        '-search:min_fragments_modeling': 2,
        '-search:min_matched_fragments': 4,
        '-java_heapmemory': 8000,
        '-threads': 6
    }

    paras_tool = params_frag
    tool = 'MSFragger'

    in_file_dbs = in_file
    out_file_dbs = in_file.replace('.mzML', '.idXML')
    databaseSearch(tool, openms_path, msfragger_path, in_file_dbs, out_file_dbs, database, paras_tool)

    in_file_pepindx = out_file_dbs  # 'F:/NTU/quantification/PXD032186_true_test/DS/CAP1.idXML'
    out_file_pepindx = in_file_pepindx.replace('.idXML',
                                               '_idx.idXML')  # 'F:/NTU/proteomics/ProInfer_codes/test_res/CAP1_idx.idXML'
    # database = 'F:/NTU/quantification/PXD032186_true_test/mus_musculus_with_decoy_contam.fasta'
    para_pepidx = {}
    peptideindexer(openms_path, in_file_pepindx, out_file_pepindx, database, para_pepidx)

    in_file_PSMFea = out_file_pepindx
    out_file_PSMFea = in_file_PSMFea.replace('_idx',
                                             '_PSMFea')  # 'F:/NTU/proteomics/ProInfer_codes/test_res/CAP1_PSMFea.idXML'
    para_psmefea = {}
    PSMFeatureExtractor(openms_path, in_file_PSMFea, out_file_PSMFea, para_psmefea)

    in_file_percolator = out_file_PSMFea
    out_file_percolator = in_file_percolator.replace('_PSMFea',
                                                     '_percolator')  # 'F:/NTU/proteomics/ProInfer_codes/test_res/CAP1_percolator.idXML'
    para_percolator = {
        '-score_type': peroolator_score_type
    }
    Percolator(openms_path, in_file_percolator, out_file_percolator, para_percolator)

    psms_proinfer_para = {
        '-id:add_protein_hit_metavalues': 100
    }

    in_file_tsv = out_file_percolator
    out_percolator_tsv_proinfer = in_file_tsv.replace('_percolator.idXML',
                                                      '_proinfer.tsv')  # 'F:/NTU/proteomics/ProInfer_codes/test_res/CAP1_proinfer.tsv'

    idXML2tsv(openms_path, in_file_tsv, out_percolator_tsv_proinfer, psms_proinfer_para)
    for file_name in listdir('./'):
        if file_name.endswith('.idXML'):
            os.remove(folder_path + file_name)
    return out_percolator_tsv_proinfer

def main():
    # Training settings
    parser = argparse.ArgumentParser(description='run ProInfer and ProInfer_cpx')
    parser.add_argument('-e', '--openms', default='./openMS/OpenMS-2.7.0/bin/',
                        help='provide the installation path of OpenMS, e.g., C:/Users/xxx/Documents/openMS/OpenMS-2.7.0/bin/')
    parser.add_argument('-f', '--msfragger', default='./FragPipe-17.1/fragpipe/tools/MSFragger-3.4/MSFragger-3.4.jar',
                        help='provide the installation path of MSFragger for database search, e.g., C:/Users/xxx/Documents/FragPipe-17.1/fragpipe/tools/MSFragger-3.4/MSFragger-3.4.jar')
    parser.add_argument('-i', '--input', default='./HeLa_TechReps_Exp1_DDA_1.mzML',
                        help='input .mzML file for database search and peptide identification, raw file can be coverted to .mzML with MSConvert')
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
                        help='Indicate the protein database for database search.')
    parser.add_argument('-cp', '--complex_path', default='./allComplexes.txt',
                        help='The protein complexes used by ProInfer_cpx. The default file is downloaded from CORUM 3.0 (http://mips.helmholtz-muenchen.de/corum/).')
    parser.add_argument('-s', '--species', default='Human',
                        help='species of the sample')
    parser.add_argument('-d', '--decoy', default='rev',
                        help='prefix of the decoy proteins')
    parser.add_argument('-y', '--score_type', default='pep',
                        help='percolator score type')

    args = parser.parse_args()

    openms_path = args.openms
    input_file = args.input
    msfragger_path = args.msfragger
    run_type = args.run_type
    psm_threshold = args.psm_threshold
    pro_qvalue_td = args.pro_qvalue_td
    save_path_proinfer = args.save_path_proinfer
    save_path_cpx = args.save_path_cpx
    protein_database = args.protein_database
    complex_path = args.complex_path

    per_pep_path = dbs_pi(openms_path, msfragger_path, input_file, protein_database)


    if run_type == 1:  # only run ProInfer v1 (without complex info)
        ProInfer.ProInfer_v1(per_pep_path, psm_threshold, protein_database, save_path_proinfer)
    elif run_type == 2:  # run both ProInfer v1 and v2 (with complex info)
        ProInfer.ProInfer_cpx_v2(per_pep_path, save_path_proinfer, save_path_cpx, psm_threshold, pro_qvalue_td, protein_database,
                        complex_path)
if __name__ == '__main__':
    main()