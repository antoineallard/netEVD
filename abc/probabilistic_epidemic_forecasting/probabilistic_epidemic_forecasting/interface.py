'''
File that creates a GUI that simplify the usage of the whole package (except results refinement)
'''

import tkinter as tk
from tkinter import Tk, ttk, filedialog


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Simulation d\'épidémie')
        self.geometry("500x700")

        self.degree_dist_var = tk.StringVar()

        self.param_A = tk.StringVar()
        self.param_A.set('  ')
        self.param_B = tk.StringVar()
        self.param_B.set('  ')
        self.param_C = tk.StringVar()
        self.param_C.set('  ')
        self.sigma = tk.StringVar()
        self.sigma.set('  ')
        self.gamma_shape = tk.StringVar()
        self.gamma_shape.set('  ')

        self.path_chosen_var = tk.StringVar()
        self.path_chosen_var.set('No data chosen yet')

        self.param_A_type_var = tk.StringVar()
        self.param_A_type_var.set('uniform')
        self.param_A_low_label_var = tk.StringVar()
        self.param_A_low_label_var.set('low')
        self.param_A_high_label_var = tk.StringVar()
        self.param_A_high_label_var.set('high')
        self.param_A_low_entry_var = tk.StringVar()
        self.param_A_low_entry_var.set(1)
        self.param_A_high_entry_var = tk.StringVar()
        self.param_A_high_entry_var.set(2)

        self.param_B_type_var = tk.StringVar()
        self.param_B_type_var.set('uniform')
        self.param_B_low_label_var = tk.StringVar()
        self.param_B_low_label_var.set('low')
        self.param_B_high_label_var = tk.StringVar()
        self.param_B_high_label_var.set('high')
        self.param_B_low_entry_var = tk.StringVar()
        self.param_B_high_entry_var = tk.StringVar()
        self.param_B_low_entry_var.set(0.1)
        self.param_B_high_entry_var.set(0.9)

        self.param_C_type_var = tk.StringVar()
        self.param_C_type_var.set('uniform')
        self.param_C_low_label_var = tk.StringVar()
        self.param_C_low_label_var.set('low')
        self.param_C_high_label_var = tk.StringVar()
        self.param_C_high_label_var.set('high')
        self.param_C_low_entry_var = tk.StringVar()
        self.param_C_high_entry_var = tk.StringVar()
        self.param_C_low_entry_var.set(0.1)
        self.param_C_high_entry_var.set(0.9)

        self.sigma_type_var = tk.StringVar()
        self.sigma_type_var.set('uniform')
        self.sigma_low_label_var = tk.StringVar()
        self.sigma_low_label_var.set('low')
        self.sigma_high_label_var = tk.StringVar()
        self.sigma_high_label_var.set('high')
        self.sigma_low_entry_var = tk.StringVar()
        self.sigma_high_entry_var = tk.StringVar()
        self.sigma_low_entry_var.set(10)
        self.sigma_high_entry_var.set(30)

        self.gamma_shape_type_var = tk.StringVar()
        self.gamma_shape_type_var.set('uniform')
        self.gamma_shape_low_label_var = tk.StringVar()
        self.gamma_shape_low_label_var.set('low')
        self.gamma_shape_high_label_var = tk.StringVar()
        self.gamma_shape_high_label_var.set('high')
        self.gamma_shape_low_entry_var = tk.StringVar()
        self.gamma_shape_high_entry_var = tk.StringVar()
        self.gamma_shape_low_entry_var.set(0.2)
        self.gamma_shape_high_entry_var.set(2)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)

        self.create_widgets()
        self.param_init()



    def create_widgets(self):

        padding = {'padx': 5, 'pady': 5}

        # Combobox for degree distribution
        ttk.Label(self, text='Excess Degree Distribution (EDD):', font=('TkDefaultFont', 10, 'bold')).grid(column=0, row=0, **padding)
        self.degree_distribution = ttk.Combobox(self, value =('NBinom','Binom','ShiftedPowerLaw','ShiftedPowExpcut', 'Poisson','Exponential'), textvariable=self.degree_dist_var)
        self.degree_distribution.grid(column=1, row=0, **padding)
        self.degree_distribution.current(0)
        self.degree_distribution.bind("<<ComboboxSelected>>", self.fix_params)
        
        
        # Entry for param_A
        self.param_A_label = ttk.Label(self, textvariable=self.param_A, font=('TkDefaultFont', 9, 'bold'))
        self.param_A_label.grid(column=0, row=2, **padding)
        self.param_A_type =  ttk.Combobox(self, value =('uniform','gamma','normal','loguniform'), textvariable=self.param_A_type_var)
        self.param_A_type.bind("<<ComboboxSelected>>", self.fix_A_subparams)
        self.param_A_type.grid(column=1, row=2, **padding)
        self.param_A_low_label = ttk.Label(self, textvariable=self.param_A_low_label_var)
        self.param_A_high_label = ttk.Label(self, textvariable=self.param_A_high_label_var)
        self.param_A_low_label.grid(column=0, row=3, **padding)
        self.param_A_high_label.grid(column=0, row=4, **padding)
        self.param_A_low_entry = ttk.Entry(self, textvariable=self.param_A_low_entry_var)
        self.param_A_high_entry = ttk.Entry(self, textvariable=self.param_A_high_entry_var)
        self.param_A_low_entry.grid(column=1, row=3, **padding)
        self.param_A_high_entry.grid(column=1, row=4, **padding)

        # Entry for param_B
        self.param_B_label = ttk.Label(self, textvariable=self.param_B,font=('TkDefaultFont', 9, 'bold'))
        self.param_B_label.grid(column=0, row=5, **padding)
        self.param_B_type =  ttk.Combobox(self, value =('uniform','gamma','normal','loguniform'), textvariable=self.param_B_type_var)
        self.param_B_type.bind("<<ComboboxSelected>>", self.fix_B_subparams)
        self.param_B_type.grid(column=1, row=5, **padding)
        self.param_B_low_label = ttk.Label(self, textvariable=self.param_B_low_label_var)
        self.param_B_high_label = ttk.Label(self, textvariable=self.param_B_high_label_var)
        self.param_B_low_label.grid(column=0, row=6, **padding)
        self.param_B_high_label.grid(column=0, row=7, **padding)
        self.param_B_low_entry = ttk.Entry(self, textvariable=self.param_B_low_entry_var)
        self.param_B_high_entry = ttk.Entry(self, textvariable=self.param_B_high_entry_var)
        self.param_B_low_entry.grid(column=1, row=6, **padding)
        self.param_B_high_entry.grid(column=1, row=7, **padding)
        
        # Entry for param_C
        self.param_C_label = ttk.Label(self, textvariable=self.param_C,font=('TkDefaultFont', 9, 'bold'))
        self.param_C_label.grid(column=0, row=8, **padding)
        self.param_C_type =  ttk.Combobox(self, value =('uniform','gamma','normal','loguniform'), textvariable=self.param_C_type_var)
        self.param_C_type.bind("<<ComboboxSelected>>", self.fix_C_subparams)
        self.param_C_type.grid(column=1, row=8, **padding)
        self.param_C_low_label = ttk.Label(self, textvariable=self.param_C_low_label_var)
        self.param_C_high_label = ttk.Label(self, textvariable=self.param_C_high_label_var)
        self.param_C_low_label.grid(column=0, row=9, **padding)
        self.param_C_high_label.grid(column=0, row=10, **padding)
        self.param_C_low_entry = ttk.Entry(self, textvariable=self.param_C_low_entry_var)
        self.param_C_high_entry = ttk.Entry(self, textvariable=self.param_C_high_entry_var)
        self.param_C_low_entry.grid(column=1, row=9, **padding)
        self.param_C_high_entry.grid(column=1, row=10, **padding)

        # Entry for sigma
        self.sigma_label = ttk.Label(self, textvariable=self.sigma,font=('TkDefaultFont', 9, 'bold'))
        self.sigma_label.grid(column=0, row=11, **padding)
        self.sigma_type =  ttk.Combobox(self, value =('uniform','gamma','normal','loguniform'), textvariable=self.sigma_type_var)
        self.sigma_type.bind("<<ComboboxSelected>>", self.fix_sigma_subparams)
        self.sigma_type.grid(column=1, row=11, **padding)
        self.sigma_low_label = ttk.Label(self, textvariable=self.sigma_low_label_var)
        self.sigma_high_label = ttk.Label(self, textvariable=self.sigma_high_label_var)
        self.sigma_low_label.grid(column=0, row=12, **padding)
        self.sigma_high_label.grid(column=0, row=13, **padding)
        self.sigma_low_entry = ttk.Entry(self, textvariable=self.sigma_low_entry_var)
        self.sigma_high_entry = ttk.Entry(self, textvariable=self.sigma_high_entry_var)
        self.sigma_low_entry.grid(column=1, row=12, **padding)
        self.sigma_high_entry.grid(column=1, row=13, **padding)

        # Entry for gamma_shape
        self.gamma_shape_label = ttk.Label(self, textvariable=self.gamma_shape,font=('TkDefaultFont', 9, 'bold'))
        self.gamma_shape_label.grid(column=0, row=14, **padding)
        self.gamma_shape_type =  ttk.Combobox(self, value =('uniform','gamma','normal','loguniform'), textvariable=self.gamma_shape_type_var)
        self.gamma_shape_type.bind("<<ComboboxSelected>>", self.fix_gamma_subparams)
        self.gamma_shape_type.grid(column=1, row=14, **padding)
        self.gamma_shape_low_label = ttk.Label(self, textvariable=self.gamma_shape_low_label_var)
        self.gamma_shape_high_label = ttk.Label(self, textvariable=self.gamma_shape_high_label_var)
        self.gamma_shape_low_label.grid(column=0, row=15, **padding)
        self.gamma_shape_high_label.grid(column=0, row=16, **padding)
        self.gamma_shape_low_entry = ttk.Entry(self, textvariable=self.gamma_shape_low_entry_var)
        self.gamma_shape_high_entry = ttk.Entry(self, textvariable=self.gamma_shape_high_entry_var)
        self.gamma_shape_low_entry.grid(column=1, row=15, **padding)
        self.gamma_shape_high_entry.grid(column=1, row=16, **padding)

        # Entry for initial cases
        self.seed_label = ttk.Label(self, text='Initial cases',font=('TkDefaultFont', 9, 'bold'))
        self.seed_label.grid(column=0, row=17, **padding)
        self.seed_entry_var = tk.StringVar()
        self.seed_entry_var.set(1)
        self.seed_entry = ttk.Entry(self, textvariable=self.seed_entry_var)
        self.seed_entry.grid(column=1, row=17, **padding)

        # Entry for max cases
        self.max_cases_label = ttk.Label(self, text='Max cases',font=('TkDefaultFont', 9, 'bold'))
        self.max_cases_label.grid(column=0, row=18, **padding)
        self.max_cases_entry_var = tk.StringVar()
        self.max_cases_entry_var.set(10000)
        self.max_cases_entry = ttk.Entry(self, textvariable=self.max_cases_entry_var)
        self.max_cases_entry.grid(column=1, row=18, **padding)

        # Entry for max time
        self.max_time_label = ttk.Label(self, text='Max time',font=('TkDefaultFont', 9, 'bold'))
        self.max_time_label.grid(column=0, row=19, **padding)
        self.max_time_entry_var = tk.StringVar()
        self.max_time_entry_var.set(90)
        self.max_time_entry = ttk.Entry(self, textvariable=self.max_time_entry_var)
        self.max_time_entry.grid(column=1, row=19, **padding)

        # Entry for number of simulatons
        self.N_simulation_label = ttk.Label(self, text='Number of simul.',font=('TkDefaultFont', 9, 'bold'))
        self.N_simulation_label.grid(column=0, row=20, **padding)
        self.N_simulation_entry_var = tk.StringVar()
        self.N_simulation_entry_var.set(1000)
        self.N_simulation_entry = ttk.Entry(self, textvariable=self.N_simulation_entry_var)
        self.N_simulation_entry.grid(column=1, row=20, **padding)

        # Entry for days of tolerance check
        self.windows_label = ttk.Label(self, text='Days of tolerance check',font=('TkDefaultFont', 9, 'bold'))
        self.windows_label.grid(column=0, row=21, **padding)
        self.windows_entry_var = tk.StringVar()
        self.windows_entry_var.set('45,90')
        self.windows_entry = ttk.Entry(self, textvariable=self.windows_entry_var)
        self.windows_entry.grid(column=1, row=21, **padding)

        # Entry for tolerances
        self.tolerances_label = ttk.Label(self, text='Tolerances',font=('TkDefaultFont', 9, 'bold'))
        self.tolerances_label.grid(column=0, row=22, **padding)
        self.tolerances_entry_var = tk.StringVar()
        self.tolerances_entry_var.set('30,30')
        self.tolerances_entry = ttk.Entry(self, textvariable=self.tolerances_entry_var)
        self.tolerances_entry.grid(column=1, row=22, **padding)
        
        # File search button
        self.path_label = ttk.Label(self, text='Raw data path',font=('TkDefaultFont', 9, 'bold'))
        self.path_label.grid(column=0, row=23, **padding)
        self.explore_button = ttk.Button(self, text='Choose raw data', command=self.explore)
        self.explore_button.grid(column=1, row=23, **padding)
        self.path_chosen = ttk.Label(self, textvariable=self.path_chosen_var,font=('TkDefaultFont', 9))
        self.path_chosen.grid(column=0, row=24, **padding, columnspan=2)

        # Submit Button
        submit_button = ttk.Button(self, text='Submit', command=self.submit)
        submit_button.grid(column=1, row=40, **padding)


    def submit(self):
        self.destroy()
        self.config_dict = {'param_A': {'distribution': self.param_A_type_var.get(),  
                                self.param_A_type_var.get(): {self.param_A_low_label_var.get(): float(self.param_A_low_entry_var.get()),
                                      self.param_A_high_label_var.get(): float(self.param_A_high_entry_var.get())}},
                          

                   'param_B': {'distribution': self.param_B_type_var.get(),  
                                self.param_B_type_var.get(): {self.param_B_low_label_var.get(): float(self.param_B_low_entry_var.get()),
                                     self.param_B_high_label_var.get(): float(self.param_B_high_entry_var.get())}},

                   'param_C': {'distribution': self.param_C_type_var.get(),  
                                self.param_C_type_var.get(): {self.param_C_low_label_var.get(): float(self.param_C_low_entry_var.get()),
                                     self.param_C_high_label_var.get(): float(self.param_C_high_entry_var.get())}},

                   'sigma': {'distribution': self.sigma_type_var.get(),  
                                self.sigma_type_var.get(): {self.sigma_low_label_var.get(): float(self.sigma_low_entry_var.get()),
                                     self.sigma_high_label_var.get(): float(self.sigma_high_entry_var.get())}},

                   'gamma_shape': {'distribution': self.gamma_shape_type_var.get(),  
                                self.gamma_shape_type_var.get(): {self.gamma_shape_low_label_var.get(): float(self.gamma_shape_low_entry_var.get()),
                                     self.gamma_shape_high_label_var.get(): float(self.gamma_shape_high_entry_var.get())}},

                   'seed': int(self.seed_entry_var.get()),                            
                   'max_cases': int(self.max_cases_entry_var.get()), 
                   'max_time': int(self.max_time_entry_var.get()),                   
                   'N_simulations': int(self.N_simulation_entry_var.get()),                
                   'excess_outward_degree_distributions': self.degree_dist_var.get(),  
                   'manual_p0':False,
                   'tolerances' : [[int(i) for i in self.windows_entry_var.get().split(',')],
                                    [float(i) for i in self.tolerances_entry_var.get().split(',')]],
                   'data_path': self.path,
                   'folder_path': None,
                   'folder_name': None,
                   'custom_p0': None} 

    def explore(self):
        self.path = filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = (("Text files",
                                                            "*.txt*"),
                                                        ("all files",
                                                            "*.*")))
        self.path_chosen_var.set(self.path)
    
    def remove_C(self):
        self.param_C_label.grid_remove()
        self.param_C_type.grid_remove()
        self.param_C_low_entry.grid_remove()
        self.param_C_high_entry.grid_remove()
        self.param_C_low_label.grid_remove()
        self.param_C_high_label.grid_remove()

    def remove_B(self):
        self.param_B_label.grid_remove()
        self.param_B_type.grid_remove()
        self.param_B_low_entry.grid_remove()
        self.param_B_high_entry.grid_remove()
        self.param_B_low_label.grid_remove()
        self.param_B_high_label.grid_remove()

    def show_B(self):
        self.param_B_label.grid()
        self.param_B_type.grid()
        self.param_B_low_entry.grid()
        self.param_B_high_entry.grid()
        self.param_B_low_label.grid()
        self.param_B_high_label.grid()

    def param_init(self):
        self.param_A.set('R_0 prior distribution')
        self.param_B.set('k prior distribution') 
        self.sigma.set('sigma prior distribution')
        self.gamma_shape.set('gamma_shape prior distribution')
        self.remove_C()

    def fix_params(self, event):
        distr = event.widget.get()
        if distr == 'NBinom':
            self.param_A.set('R_0 prior distribution')
            self.param_B.set('k prior distribution')
            self.show_B()
            self.remove_C()
        if distr == 'Binom':
            self.param_A.set('R_0 prior distribution')
            self.param_B.set('k prior distribution')
            self.show_B()
            self.remove_C()
        if distr == 'ShiftedPowExpcut':
            self.param_A.set('Kappa prior distribution')
            self.param_B.set('tau prior distribution')
            self.show_B()
            self.remove_C()
        if distr == 'ShiftedPowerLaw':
            self.param_A.set('Tau prior distribution')
            self.remove_B()
            self.remove_C()
        if distr == 'Poisson':
            self.param_A.set('z')
            self.remove_B()
            self.remove_C()
        if distr == 'Exponential':
            self.param_A.set('kappa')
            self.remove_B()
            self.remove_C()
    
    def fix_A_subparams(self, event):
        distr = event.widget.get()
        if distr == 'uniform':
            self.param_A_low_label_var.set('low')
            self.param_A_high_label_var.set('high')
        if distr == 'loguniform':
            self.param_A_low_label_var.set('low')
            self.param_A_high_label_var.set('high')
        if distr == 'gamma':
            self.param_A_low_label_var.set('a')
            self.param_A_high_label_var.set('b')
        if distr == 'normal':
            self.param_A_low_label_var.set('mean')
            self.param_A_high_label_var.set('std')
    
    def fix_B_subparams(self, event):
        distr = event.widget.get()
        if distr == 'uniform':
            self.param_B_low_label_var.set('low')
            self.param_B_high_label_var.set('high')
        if distr == 'loguniform':
            self.param_B_low_label_var.set('low')
            self.param_B_high_label_var.set('high')
        if distr == 'gamma':
            self.param_B_low_label_var.set('a')
            self.param_B_high_label_var.set('b')
        if distr == 'normal':
            self.param_B_low_label_var.set('mean')
            self.param_B_high_label_var.set('std')
    
    def fix_C_subparams(self, event):
        distr = event.widget.get()
        if distr == 'uniform':
            self.param_C_low_label_var.set('low')
            self.param_C_high_label_var.set('high')
        if distr == 'uniform':
            self.param_C_low_label_var.set('low')
            self.param_C_high_label_var.set('high')
        if distr == 'gamma':
            self.param_C_low_label_var.set('a')
            self.param_C_high_label_var.set('b')
        if distr == 'normal':
            self.param_C_low_label_var.set('mean')
            self.param_C_high_label_var.set('std')
    
    def fix_sigma_subparams(self, event):
        distr = event.widget.get()
        if distr == 'uniform':
            self.sigma_low_label_var.set('low')
            self.sigma_high_label_var.set('high')
        if distr == 'loguniform':
            self.sigma_low_label_var.set('low')
            self.sigma_high_label_var.set('high')
        if distr == 'gamma':
            self.sigma_low_label_var.set('a')
            self.sigma_high_label_var.set('b')
        if distr == 'normal':
            self.sigma_low_label_var.set('mean')
            self.sigma_high_label_var.set('std')
    
    def fix_gamma_subparams(self, event):
        distr = event.widget.get()
        if distr == 'uniform':
            self.gamma_shape_low_label_var.set('low')
            self.gamma_shape_high_label_var.set('high')
        if distr == 'loguniform':
            self.gamma_shape_low_label_var.set('low')
            self.gamma_shape_high_label_var.set('high')
        if distr == 'gamma':
            self.gamma_shape_low_label_var.set('a')
            self.gamma_shape_high_label_var.set('b')
        if distr == 'normal':
            self.gamma_shape_low_label_var.set('mean')
            self.gamma_shape_high_label_var.set('std')



if __name__ == "__main__":
    print('\n')