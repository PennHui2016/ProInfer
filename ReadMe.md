# ProInfer
ProInfer is a protein inference tool.

If peptide identification results from percolator is already avaliable, then using ProInfer, otherwise, using OpenMS_ProInfer 

The ProInfer codes programmed with python 3.8.8.Required packages include:
1. pandas
2. numpy
3. random
4. os
5. copy
6. csv
7. sys
8. argparse
9. scipy
10. tkinter
11. PIL

# Part 1 KNIME workflow for ProInfer
ProInfer.py requires the iddentified peptides as input. To prepare the inputs, please
use the attached KNIME workflow (file name: preparing_peptides_workflow.knwf).

There are 10 required inputs:

1. per_pep_path: String, required. 
                 The path of the file containing input peptide list, can obtain with 
				 the workflow: preparing_peptides_workflow.knwf
				 we attached a toy data './DDA1.tsv'

2. run_type: int, 1 (default)--runs the ProInfer; 2--runs the ProInfer_cpx.
             Indicate which tool will be run. 

3. psm_threshold: float (0, 1], default 0.999. 
                  Selection of posterior error probability (PEP)
                  threshold to filter low confidence peptides out.

4. pro_qvalue_td: float (0, 1], default 0.01.
                  Selection of False Discovery Rate (FDR) for reporting target proteins.
				  
5. save_path_proinfer: String, default './res/proinfer_out.csv'.
                       Indicate where to save the results from ProInfer

6. save_path_cpx: String, default './res/proinfer_cpx_out'.
                       Indicate where to save the results from ProInfer

7. protein_database: String, default './Human_database_including_decoys_2022-06-23-decoys-contam-uniprot-proteome_UP000005640_2022_5_5.fasta'.
                     Indicate the protein database for database search. Should be the same as
					 what you used in the KNIME workflow.

8. complex_path: String, default './allComplexes.txt'.
                 The protein complexes used by ProInfer_cpx. The default file is downloaded
				 from CORUM 3.0 (http://mips.helmholtz-muenchen.de/corum/).
				 
9. species: String, default 'Human'.
	    Indicate which species does the sample come from

10. decoy: String, default 'rev'
	   Indicate the prefix of protein ID for decoy proteins
				 
The outputs from ProInfer are stored in .csv file. There are 5 columns, including 
Protein IDs (column 1), accPEP scores (column 2), confidence score (column 3),
q-values (column 4), proteins labels, 1 -- target; -1 -- decoy (column 5).

The outputs from ProInfer_cpx are stored in .csv file. There are at least 5 columns, including 
Protein IDs (column 1), accPEP scores (column 2), confidence score (column 3),
q-values (column 4), proteins labels, 1 -- target; -1 -- decoy (column 5), if more than
5 folumns in this file, then the following columns should be in the order of accPEP scores (column 6),
confidence score (column 7), q-values (column 8), accPEP scores (column 9), confidence score (column 10),
q-values (column 11), ..., the last group of accPEP scores, confidence score, q-values are the final outputs.
please use them for reporting proteins.

checking parameters:
    
    python ProInfer.py -h

The tool can be run in following command:
    
    python ProInfer.py -i [per_pep_path] -t [run_type] -pt [psm_threshold] -qt [pro_qvalue_td] -sp [save_path_proinfer] -spc [save_path_cpx] -db [protein_database] -cp [complex_path] -s [species] -d [decoy]
    
Example runnings with the toy data './DDA1.tsv':

Example 1: running the ProInfer:
     
     python ProInfer.py -i ./DDA1.tsv

which equals to:
    
    python ProInfer.py -i ./DDA1.tsv -t 1 -pt 0.999 -qt 0.01 -sp ./res/proinfer_out.csv -spc "" -db ./2022-06-23-decoys-contam-uniprot-proteome_UP000005640_2022_5_5.fasta -cp "" -s Human -d rev

Example 2: running the ProInfer_cpx:
     
     python ProInfer.py -i ./DDA1.tsv -t 2

which equals to:
    
    python ProInfer.py -i ./DDA1.tsv -t 2 -pt 0.999 -qt 0.01 -sp ./res/proinfer_out.csv -spc ./res/proinfer_cpx_out -db ./2022-06-23-decoys-contam-uniprot-proteome_UP000005640_2022_5_5.fasta -cp ./allComplexes.txt -s Human -d rev 

# Part 2 OpenMS for ProInfer
OpenMS_ProInfer.py accepts MS data in .mzML format as input. OpenMS (https://www.openms.de/downloads/) and MSFragger (https://github.com/Nesvilab/MSFragger) is required
to be installed. Parameters need to be specified including:
1. openms: string, OpenMS installation path
2. msfragger: string, MSFragger.jar package path
3. input: string, .mzML file, example data can be downloaded from https://www.ebi.ac.uk/pride/archive/projects/PXD022448 (.raw file can be converted to .mzML with MSConvert)
4. run_type: int, 1 (default)--runs the ProInfer; 2--runs the ProInfer_cpx.
             Indicate which tool will be run. 

5. psm_threshold: float (0, 1], default 0.999. 
                  Selection of posterior error probability (PEP)
                  threshold to filter low confidence peptides out.

6. pro_qvalue_td: float (0, 1], default 0.01.
                  Selection of False Discovery Rate (FDR) for reporting target proteins.
				  
7. save_path_proinfer: String, default './res/proinfer_out.csv'.
                       Indicate where to save the results from ProInfer

8. save_path_cpx: String, default './res/proinfer_cpx_out'.
                       Indicate where to save the results from ProInfer

9. protein_database: String, default './2022-06-23-decoys-contam-uniprot-proteome_UP000005640_2022_5_5.fasta'.
                     Indicate the protein database for database search.

10. complex_path: String, default './allComplexes.txt'.
                 The protein complexes used by ProInfer_cpx. The default file is downloaded
				 from CORUM 3.0 (http://mips.helmholtz-muenchen.de/corum/).
				 
11. species: String, default 'Human'. 
				Indicate which speies does the sample from.
				
12. decoy: String, default 'rev'. 
				Indicate what's the prefix of the decoy protein.	
13. score_type: String, default 'pep'. 
				Indicate which score is output by percolator.
				
checking parameters:
    python OpenMS_ProInfer.py -h
			
# Part 3 Stand-alone version ProInfer
We built a stand-alone version application ProInfer_tk, One can download the ProInfer_tk.exe file from https://drive.google.com/file/d/1GL3rXXADYKhjprHYwhGZ41_J-ThhsDOX/view?usp=sharing (ProInfer.exe is too big to upload to GitHub) or run ProInfer_tk.py with Python via following simple command:

	python ProInfer_tk.py
	
If OpenMS_ProInfer is used, please install OpenMS and MSFragger at first (see above part 2)

Any problems or requesting source codes for reproducing results in our paper please contact 
    Hui Peng: hui.peng@ntu.edu.sg
	Wilson Wen Bin Goh: wilsongoh@ntu.edu.sg
                        
