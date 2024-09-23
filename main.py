from misc import *
from simulator.ode import ODE
from simulator.graph import Simulator
from experiment import Experiment
import sys, getopt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# file locations
mRNA_file = './test_data/Ecoli-1_wildtype.tsv'
protein_file = './test_data/Ecoli-1_proteins_wildtype.tsv'
perturbs_file = './test_data/Ecoli-1_dream4_timeseries_perturbations.tsv'
graph_file = './test_data/Ecoli-1_goldstandard_signed.tsv'
grn_file = './test_data/Ecoli-1.xml'

# data frames
mRNA_df = pd.read_csv(mRNA_file, sep='\t', decimal=',')
protein_df = pd.read_csv(protein_file, sep='\t', decimal=',')
perturbs_df = pd.read_csv(perturbs_file, sep='\t', decimal='.')

# time series parameters
N = 51       # number of samples
max_t = 1000 # max time

t = np.linspace(0, max_t, N)


def run_kinetic():
    # simulate kinetic model and plot timeseries
    genes_dict = sbml_to_dict(grn_file)
    genes = list(genes_dict.keys())
    x = list(mRNA_df.iloc[0,:])
    y = list(protein_df.iloc[0,:])
    perturbs = np.array(perturbs_df.iloc[0])

    ode = ODE(genes_dict, x, y)
    ode.solve_ode(t, perturbs)

    plt.figure()
    for i in range(10):
        plt.plot(np.array(ode.y_vecs)[:,i], label=genes[i])
        plt.legend()
        plt.title('Timeseries plot of normalized mRNA expression levels')
        plt.xlabel('time [s]')
        plt.ylabel('mRNA expression level')
        plt.savefig('./figures/example-kinetic.pdf') 
    plt.show()

def run_boolean():
    # simulate boolean model and display GRN evolution
    graph_df = pd.read_csv(graph_file, header=None, sep='\t', decimal=',')
    graph_config = [tuple(item) for item in graph_df.to_numpy()]
    genes_dict = sbml_to_dict(grn_file)
    perturbs = np.array(perturbs_df.iloc[0])

    init_states = mRNA_df.iloc[0,:] >= 0.5

    sim = Simulator(genes_dict, graph_config, init_states)
    out = sim.run_simulation(perturbs,visualize=True)

def main(argv):
    model = ''
    try:
        opts, args = getopt.getopt(argv,"m:he",["model=", "help", "experiment"])
    except getopt.GetoptError:
        print('main.py --model [boolean|kinetic]')
        sys.exit(2)

    if opts == []:
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('main.py --model [boolean|kinetic]')
            model = arg
        elif opt in ("-m", "--model"):
            model = arg
        elif opt in ("-e", "--experiment"):
            model = None
            experiment = Experiment()
            data = experiment.run()
            print(data)
            experiment.plot_data(data)
    
    if model == 'boolean':
        run_boolean()
    elif model == 'kinetic':
        run_kinetic()


if __name__ == "__main__":
   main(sys.argv[1:])