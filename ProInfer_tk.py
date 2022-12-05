import tkinter as tk
from tkinter.filedialog import (askopenfilename,
                                    askopenfilenames,
                                    askdirectory,
                               asksaveasfilename)

import ProInfer
import OpenMS_ProInfer
from PIL import Image, ImageTk


img = Image.open('./tk_photo.png')


class basedesk():
    def __init__(self, master):
        self.root = master
        self.root.config()
        self.root.title('Base page')
        self.root.geometry('1100x700')

        initface(self.root)


class initface():
    def __init__(self, master):
        self.master = master
        self.master.config(bg='white')
        # initface
        self.menubar = tk.Menu(self.master,)
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label='tools', menu=self.filemenu)
        self.filemenu.add_command(label='ProInfer', command=self.master.title('ProInfer'))
        self.filemenu.add_command(label='OpenMS_ProInfer', command=self.change)
        self.helpmenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label='Help', menu=self.helpmenu)
        self.helpmenu.add_command(label='Help', command=self.show_help)
        self.master.config(menu=self.menubar)
        self.initface = tk.Frame(self.master, bd=1,width=1100,height=700)
        self.initface.pack()

        self.path_var = tk.StringVar()
        self.path_var.set("./DDA1.tsv")
        #self.run_type_var = tk.StringVar()
        self.psm_threshold = tk.StringVar()
        self.psm_threshold.set('0.999')
        self.pro_fdr = tk.StringVar()
        self.pro_fdr.set('0.05')
        self.db_var = tk.StringVar()
        self.db_var.set('./2022-06-23-decoys-contam-uniprot-proteome_UP000005640_2022_5_5.fasta')
        self.cpx_var = tk.StringVar()
        self.cpx_var.set('./allComplexes.txt')
        self.dec_prefix = tk.StringVar()
        self.dec_prefix.set('rev')
        self.species_var = tk.StringVar()
        self.species_var.set('Human')
        self.save_var = tk.StringVar()
        self.save_var.set('./res')

        tk.Entry(self.initface, textvariable=self.path_var, width=80).place(x=60, y=45, anchor='nw')

        tk.Label(self.initface, text='please select the percolator result file:', font=('Arial', 12)).place(x=0, y=10, anchor='nw')
        tk.Button(self.initface, text='Select', command=self.select_pth).place(x=0, y=40, anchor='nw')
        #tk.Radiobutton(self.initface, text='without complex', variable=self.run_type_var, value='1', command=self.set_run_type).place(x=5, y=60, anchor='nw')

        #tk.Radiobutton(self.initface, text='with complex', variable=self.run_type_var, value='2', command=self.set_run_type).place(x=55, y=60, anchor='nw')
        ### PSM threshold
        tk.Label(self.initface, text='please set the PSM filtering threshold (0,1):', font=('Arial', 12)).place(x=0, y=80, anchor='nw')
        tk.Entry(self.initface, textvariable=self.psm_threshold, width=5).place(x=35, y=105, anchor='nw')

        tk.Button(self.initface, text='Set', command=self.set_psm_threshold).place(x=80, y=100, anchor='nw')

        ### protein report FDR
        tk.Label(self.initface, text='please set the protein reproting FDR (0,1):', font=('Arial', 12)).place(x=0,
                                                                                                                y=140,
                                                                                                                anchor='nw')
        tk.Entry(self.initface, textvariable=self.pro_fdr, width=5).place(x=35, y=165, anchor='nw')

        tk.Button(self.initface, text='Set', command=self.set_pro_fdr).place(x=80, y=160, anchor='nw')

        ### database path
        tk.Label(self.initface, text='please select the protein database file (contain decoy and contaminants):', font=('Arial', 12)).place(x=0,
                                                                                                              y=200,
                                                                                                              anchor='nw')
        tk.Entry(self.initface, textvariable=self.db_var, width=80).place(x=60, y=235, anchor='nw')

        tk.Button(self.initface, text='Select', command=self.select_db).place(x=0, y=230, anchor='nw')

        ### decy prefix string path
        tk.Label(self.initface, text='please specify the prefix string of decoy proteins:',
                 font=('Arial', 12)).place(x=0,
                                           y=280,
                                           anchor='nw')
        tk.Entry(self.initface, textvariable=self.dec_prefix, width=5).place(x=35, y=315, anchor='nw')

        tk.Button(self.initface, text='Set', command=self.set_dec_prefix).place(x=80, y=310, anchor='nw')

        ### complex path
        tk.Label(self.initface, text='please select the complex file:', font=('Arial', 12)).place(x=0,
                                                                                                           y=350,
                                                                                                           anchor='nw')
        tk.Entry(self.initface, textvariable=self.cpx_var, width=80).place(x=60, y=385, anchor='nw')

        tk.Button(self.initface, text='Select', command=self.select_cpx).place(x=0, y=380, anchor='nw')

        ### species
        tk.Label(self.initface, text='please specify the species, e.g., Human or Mouse:', font=('Arial', 12)).place(x=0,
                                                                                                  y=420,
                                                                                                  anchor='nw')
        tk.Entry(self.initface, textvariable=self.species_var, width=10).place(x=35, y=455, anchor='nw')

        tk.Button(self.initface, text='Set', command=self.set_species).place(x=90, y=450, anchor='nw')

        ### save path
        tk.Label(self.initface, text='please specify result saving folder:', font=('Arial', 12)).place(x=0,
                                                                                                  y=490,
                                                                                                  anchor='nw')
        tk.Entry(self.initface, textvariable=self.save_var, width=80).place(x=60, y=525, anchor='nw')

        tk.Button(self.initface, text='Select', command=self.save_pth).place(x=0, y=520, anchor='nw')

        ### run
        tk.Button(self.initface, text='Run ProInfer', command=self.run_Proinfer, width=20, bg = 'lightgray', font=('Arial', 12)).place(x=500, y=590, anchor='nw')

    def change(self, ):
        self.initface.destroy()
        face1(self.master)
    def show_help(self, ):
        self.initface.destroy()
        help_info(self.master)
        self.master.title('Help')
    def select_pth(self):
        file_pth = askopenfilename(title='Please choose a file',
                                   initialdir='/', filetypes=[('Percoltor result file', '*.tsv')])

        self.path_var.set(file_pth)
    def set_run_type(self):
        self.run_type_var.get()
        print(self.run_type_var.get())
    def set_psm_threshold(self):
        self.psm_threshold.get()
    def set_pro_fdr(self):
        self.pro_fdr.get()
    def select_db(self):
        file_pth = askopenfilename(title='Please choose a file',
                                   initialdir='/', filetypes=[('protein database file', '*.fasta')])
        self.db_var.set(file_pth)
    def select_cpx(self):
        file_pth = askopenfilename(title='Please choose a file',
                                   initialdir='/', filetypes=[('protein complex file', '*.*')])
        self.cpx_var.set(file_pth)
    def set_dec_prefix(self):
        self.dec_prefix.get()
        print(self.dec_prefix.get())

    def  set_species(self):
        self.species_var.get()
        print(self.species_var.get())
    def save_pth(self):
        file_pth = askdirectory(title='Please choose a directory',
                                   initialdir='/')
        self.save_var.set(file_pth)

    def run_Proinfer(self):
        print("running")

        per_pep_path = self.path_var.get()  # all_params[1]
        run_type = 2  # int(all_params[2])
        psm_threshold = float(self.psm_threshold.get())  # float(all_params[3])
        pro_qvalue_td = float(self.pro_fdr.get()) # float(all_params[4])
        save_path_proinfer = self.save_var.get() + '/proinfer_out.csv'  # all_params[5]
        save_path_cpx = self.save_var.get() + '/proinfer_cpx_out.csv'  # all_params[6]
        protein_database = self.db_var.get()  # all_params[7]
        complex_path = self.cpx_var.get()  # all_params[8]
        decoy = self.dec_prefix.get()
        species = self.species_var.get()

        if run_type == 1:  # only run ProInfer
            ProInfer.ProInfer_v1(per_pep_path, psm_threshold, protein_database, save_path_proinfer, decoy)
        elif run_type == 2:  # run both ProInfer and ProInfer_cpx
            ProInfer.ProInfer_cpx_v2(per_pep_path, save_path_proinfer, save_path_cpx, psm_threshold, pro_qvalue_td,
                            protein_database,
                            complex_path, species, decoy)

class face1():
    def __init__(self, master):
        self.master = master
        self.master.config(bg='white')
        self.face1 = tk.Frame(self.master, bd=1,width=1100,height=700)

        self.face1.pack()
        self.menubar = tk.Menu(self.master, )
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label='tools', menu=self.filemenu)
        self.filemenu.add_command(label='ProInfer', command=self.back)
        self.filemenu.add_command(label='OpenMS_ProInfer', command=self.master.title('OpenMS_ProInfer'))
        self.helpmenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label='Help', menu=self.helpmenu)
        self.helpmenu.add_command(label='Help', command=self.show_help)
        self.master.config(menu=self.menubar)

        self.path_var = tk.StringVar()
        self.path_var.set("F:/NTU/quantification/Hela/test/HeLa_TechReps_Exp1_DDA_1.mzML")
        # self.run_type_var = tk.StringVar()
        self.psm_threshold = tk.StringVar()
        self.psm_threshold.set('0.999')
        self.pro_fdr = tk.StringVar()
        self.pro_fdr.set('0.05')
        self.db_var = tk.StringVar()
        self.db_var.set('./2022-06-23-decoys-contam-uniprot-proteome_UP000005640_2022_5_5.fasta')
        self.cpx_var = tk.StringVar()
        self.cpx_var.set('./allComplexes.txt')
        self.dec_prefix = tk.StringVar()
        self.dec_prefix.set('rev')
        self.species_var = tk.StringVar()
        self.species_var.set('Human')
        self.save_var = tk.StringVar()
        self.save_var.set('./res')
        self.openms_pth = tk.StringVar()
        self.openms_pth.set('D:/openMS/OpenMS-2.7.0/bin/')
        self.msfragger_pth = tk.StringVar()
        self.msfragger_pth.set('D:/proteomics/2021.11.10/MSFragger-3.4/MSFragger-3.4.jar')
        self.score_type = tk.StringVar()
        self.score_type.set('pep')

        ### set openMS directory
        tk.Entry(self.face1, textvariable=self.openms_pth, width=80).place(x=65, y=45, anchor='nw')

        tk.Label(self.face1, text='please select the OpenMS installation directory:', font=('Arial', 12)).place(x=0, y=10,
                                                                                               anchor='nw')
        tk.Button(self.face1, text='Select', command=self.select_openms_pth).place(x=5, y=40, anchor='nw')

        ### set msfragger directory
        tk.Entry(self.face1, textvariable=self.msfragger_pth, width=80).place(x=65, y=115, anchor='nw')

        tk.Label(self.face1, text='please select msfragger .jar file path:', font=('Arial', 12)).place(x=0, y=80,
                                                                                               anchor='nw')
        tk.Button(self.face1, text='Select', command=self.select_msfragger_pth).place(x=5, y=110, anchor='nw')

        ### select mzML file
        tk.Entry(self.face1, textvariable=self.path_var, width=80).place(x=65, y=175, anchor='nw')

        tk.Label(self.face1, text='please select the mzML file:', font=('Arial', 12)).place(x=0, y=150,
                                                                                                            anchor='nw')
        tk.Button(self.face1, text='Select', command=self.select_pth).place(x=5, y=170, anchor='nw')
        # tk.Radiobutton(self.initface, text='without complex', variable=self.run_type_var, value='1', command=self.set_run_type).place(x=5, y=60, anchor='nw')

        # tk.Radiobutton(self.initface, text='with complex', variable=self.run_type_var, value='2', command=self.set_run_type).place(x=55, y=60, anchor='nw')
        ### PSM threshold
        tk.Label(self.face1, text='please set the PSM filtering threshold (0,1):', font=('Arial', 12)).place(x=0,
                                                                                                                y=220,
                                                                                                                anchor='nw')
        tk.Entry(self.face1, textvariable=self.psm_threshold, width=5).place(x=35, y=255, anchor='nw')

        tk.Button(self.face1, text='Set', command=self.set_psm_threshold).place(x=80, y=250, anchor='nw')

        ### protein report FDR
        tk.Label(self.face1, text='please set the protein reproting FDR (0,1):', font=('Arial', 12)).place(x=350,
                                                                                                              y=220,
                                                                                                              anchor='nw')
        tk.Entry(self.face1, textvariable=self.pro_fdr, width=5).place(x=385, y=255, anchor='nw')

        tk.Button(self.face1, text='Set', command=self.set_pro_fdr).place(x=460, y=250, anchor='nw')

        ### database path
        tk.Label(self.face1, text='please select the protein database file (contain decoy and contaminants):',
                 font=('Arial', 12)).place(x=0,
                                           y=290,
                                           anchor='nw')
        tk.Entry(self.face1, textvariable=self.db_var, width=80).place(x=60, y=325, anchor='nw')

        tk.Button(self.face1, text='Select', command=self.select_db).place(x=0, y=320, anchor='nw')

        ### decy prefix string path
        tk.Label(self.face1, text='please specify the prefix string of decoy proteins:',
                 font=('Arial', 12)).place(x=0,
                                           y=360,
                                           anchor='nw')
        tk.Entry(self.face1, textvariable=self.dec_prefix, width=5).place(x=35, y=395, anchor='nw')

        tk.Button(self.face1, text='Set', command=self.set_dec_prefix).place(x=80, y=390, anchor='nw')

        ### complex path
        tk.Label(self.face1, text='please select the complex file:', font=('Arial', 12)).place(x=0,
                                                                                                  y=430,
                                                                                                  anchor='nw')
        tk.Entry(self.face1, textvariable=self.cpx_var, width=80).place(x=60, y=465, anchor='nw')

        tk.Button(self.face1, text='Select', command=self.select_cpx).place(x=0, y=460, anchor='nw')

        ### species
        tk.Label(self.face1, text='please specify the species, e.g., Human or Mouse:', font=('Arial', 12)).place(x=0,
                                                                                                                    y=500,
                                                                                                                    anchor='nw')
        tk.Entry(self.face1, textvariable=self.species_var, width=10).place(x=35, y=535, anchor='nw')

        tk.Button(self.face1, text='Set', command=self.set_species).place(x=90, y=530, anchor='nw')

        ### percolator score type
        tk.Label(self.face1, text='please specify the score type from percolator, e.g., pep or fdr:', font=('Arial', 12)).place(x=500,
                                                                                                                 y=500,
                                                                                                                 anchor='nw')
        tk.Entry(self.face1, textvariable=self.score_type, width=10).place(x=535, y=535, anchor='nw')

        tk.Button(self.face1, text='Set', command=self.set_score_type).place(x=590, y=530, anchor='nw')

        ### save path
        tk.Label(self.face1, text='please specify result saving folder:', font=('Arial', 12)).place(x=0,
                                                                                                       y=570,
                                                                                                       anchor='nw')
        tk.Entry(self.face1, textvariable=self.save_var, width=80).place(x=60, y=605, anchor='nw')

        tk.Button(self.face1, text='Select', command=self.save_pth).place(x=0, y=600, anchor='nw')

        ### run
        tk.Button(self.face1, text='Run OpenMS_ProInfer', command=self.run_Proinfer, width=20, bg='lightgray',
                  font=('Arial', 12)).place(x=500, y=650, anchor='nw')
        #btn_back = tk.Button(self.face1, text='face1 back', command=self.back)
        #btn_back.pack()

    def back(self):
        self.face1.destroy()
        initface(self.master)
    def show_help(self, ):
        self.face1.destroy()
        help_info(self.master)
        self.master.title('Help')
    def select_openms_pth(self):
        file_pth = askdirectory(title='Please choose a directory',
                                   initialdir='/')

        self.openms_pth.set(file_pth)
    def select_msfragger_pth(self):
        file_pth = askopenfilename(title='Please choose the msfragger jar file',
                                   initialdir='/', filetypes=[('.jar file', '*.jar')])

        self.msfragger_pth.set(file_pth)
    def select_pth(self):
        file_pth = askopenfilename(title='Please choose a file',
                                   initialdir='/', filetypes=[('mzML file', '*.mzML')])

        self.path_var.set(file_pth)
    def set_run_type(self):
        self.run_type_var.get()
        print(self.run_type_var.get())
    def set_psm_threshold(self):
        self.psm_threshold.get()
    def set_score_type(self):
        self.score_type.get()
    def set_pro_fdr(self):
        self.pro_fdr.get()
    def select_db(self):
        file_pth = askopenfilename(title='Please choose a file',
                                   initialdir='/', filetypes=[('protein database file', '*.fasta')])
        self.db_var.set(file_pth)
    def select_cpx(self):
        file_pth = askopenfilename(title='Please choose a file',
                                   initialdir='/', filetypes=[('protein complex file', '*.*')])
        self.cpx_var.set(file_pth)
    def set_dec_prefix(self):
        self.dec_prefix.get()
        print(self.dec_prefix.get())

    def  set_species(self):
        self.species_var.get()
        print(self.species_var.get())
    def save_pth(self):
        file_pth = askdirectory(title='Please choose a directory',
                                   initialdir='/')
        self.save_var.set(file_pth)

    def run_Proinfer(self):
        print("running")
        openms_path = self.openms_pth.get()
        msfragger_path = self.msfragger_pth.get()
        input_file = self.path_var.get()  # all_params[1]
        run_type = 2  # int(all_params[2])
        psm_threshold = float(self.psm_threshold.get())  # float(all_params[3])
        pro_qvalue_td = float(self.pro_fdr.get()) # float(all_params[4])
        save_path_proinfer = self.save_var.get() + '/proinfer_out.csv'  # all_params[5]
        save_path_cpx = self.save_var.get() + '/proinfer_cpx_out.csv'  # all_params[6]
        protein_database = self.db_var.get()  # all_params[7]
        complex_path = self.cpx_var.get()  # all_params[8]
        decoy = self.dec_prefix.get()
        species = self.species_var.get()
        score_type = self.score_type.get()

        per_pep_path = OpenMS_ProInfer.dbs_pi(openms_path, msfragger_path, input_file, protein_database, score_type)

        if run_type == 1:  # only run ProInfer
            ProInfer.ProInfer_v1(per_pep_path, psm_threshold, protein_database, save_path_proinfer, decoy)
        elif run_type == 2:  # run both ProInfer and ProInfer_cpx
            ProInfer.ProInfer_cpx_v2(per_pep_path, save_path_proinfer, save_path_cpx, psm_threshold, pro_qvalue_td,
                            protein_database,
                            complex_path, species, decoy)

def ctrlEvent(event):
    if(12==event.state and event.keysym=='c' ):
        return
    else:
        return "break"
class help_info():
    def __init__(self, master):
        self.master = master
        self.master.config(bg='lightgray')
        self.help = tk.Frame(self.master, bd=2,width=1100, height=700)

        self.help.pack()
        self.menubar = tk.Menu(self.master, )
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label='tools', menu=self.filemenu)
        self.filemenu.add_command(label='ProInfer', command=self.back_proinfer)
        self.filemenu.add_command(label='OpenMS_ProInfer', command=self.back_openmsproinfer)
        self.helpmenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label='Help', menu=self.helpmenu)
        self.helpmenu.add_command(label='Help', command=None)
        self.master.config(menu=self.menubar)

        tk.Message(self.help, text='Tool Description                                              ', width=1100, justify='left',
                   font=('Arial', 15, 'bold')).pack(anchor='nw')
        tk.Message(self.help, text = 'ProInfer requires the iddentified peptides as input. '
                                     'To prepare the inputs, please use the attached KNIME workflow: '
                                     'https://github.com/PennHui2016/ProInfer/preparing_peptides_workflow.knwf.', width=1100, justify='left',
                   font=('Arial', 12)).pack(anchor='nw')
        tk.Message(self.help, text='OpenMS_ProInfer accepts MS data in .mzML format as input. OpenMS'
                                   ' (https://www.openms.de/downloads/) and MSFragger (https://github.com/Nesvilab/MSFragger) '
                                   'is required to be installed and their path should be provided.', width=1100, justify='left',
                   font=('Arial', 12)).pack(anchor='nw')
        tk.Message(self.help, text='More details and python source codes can be found in https://github.com/PennHui2016/ProInfer', width=1100,
                   justify='left',
                   font=('Arial', 12)).pack(anchor='nw')
        tk.Message(self.help, text='Contact', font=('Arial', 15, 'bold'), width=1100, justify='left').pack(anchor='nw')
        tk.Message(self.help,
                   text='Any problems or requesting source codes for reproducing results in our paper please contact Hui Peng: hui.peng@ntu.edu.sg Wilson Wen Bin Goh: wilsongoh@ntu.edu.sg\n\n\n',
                   width=1100,
                   justify='left',
                   font=('Arial', 12)).pack(anchor='nw')
        global photo
        photo= ImageTk.PhotoImage(img)
        tk.Label(self.help,image=photo).pack(anchor='s')

    def back_proinfer(self):
        self.help.destroy()
        self.master.title('ProInfer')
        initface(self.master)
    def back_openmsproinfer(self):
        self.help.destroy()
        self.master.title('OpenMS_ProInfer')
        face1(self.master)

if __name__ == '__main__':
    root = tk.Tk()
    basedesk(root)
    root.mainloop()
