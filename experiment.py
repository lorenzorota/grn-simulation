from misc import *
from simulator.ode import ODE
from simulator.graph import Simulator
from scipy.stats import f_oneway
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Experiment:

    def __init__(self):
        pass

    def run_kinetic(self, genes_dict, mRNA_df, protein_df, perturbs_df, t, exp):
        # simulate kinetic model and plot timeseries
        genes = list(genes_dict.keys())
        x = list(mRNA_df.iloc[0,:])
        y = list(protein_df.iloc[0,:])
        perturbs = np.array(perturbs_df.iloc[exp])

        ode = ODE(genes_dict, x, y)
        ode.solve_ode(t, perturbs)

        # return the discretized state vector
        return np.array(ode.y_vecs)[-1] >= 0.5

    def run_boolean(self, genes_dict, mRNA_df, graph_df, perturbs_df, exp):
        # simulate boolean model and display GRN evolution
        graph_config = [tuple(item) for item in graph_df.to_numpy()]
        init_states = mRNA_df.iloc[0,:] >= 0.5
        perturbs = np.array(perturbs_df.iloc[exp])

        sim = Simulator(genes_dict, graph_config, init_states)
        out = sim.run_simulation(perturbs,visualize=False)
        
        # return last state vector (assuming it is a fixed point state)
        return out[:,-1].astype('bool')

    def run(self):
        # experiment parameters (should not be changed unless data is generated for it)
        node_counts = [5, 10, 20, 50, 100]
        subnets = [i for i in range(1,11)]
        experiments = [i for i in range(0,10)]

        # time series parameters
        M = 10          # number of different experiments
        N = 21          # number of samples
        max_t = 1000    # max time
        t = np.linspace(0, max_t, N)

        node_count_dependent_errors = []
        for node_count in node_counts:
            print('----------------')
            print('node count: {}'.format(node_count))
            print('----------------')

            subnet_errors = []
            for subnet in subnets:
                print('----------------')
                print('subnet: {}'.format(subnet))
                print('----------------')

                experiment_errors = []
                for experiment in experiments:

                    # location for files of subnet
                    mRNA_file = './experiment_data/{}_node/{}/Ecoli-{}_wildtype.tsv'.format(str(node_count), str(subnet), str(subnet))
                    protein_file = './experiment_data/{}_node/{}/Ecoli-{}_proteins_wildtype.tsv'.format(str(node_count), str(subnet), str(subnet))
                    perturbs_file = './experiment_data/{}_node/{}/Ecoli-{}_dream4_timeseries_perturbations.tsv'.format(str(node_count), str(subnet), str(subnet))
                    graph_file = './experiment_data/{}_node/{}/Ecoli-{}_goldstandard_signed.tsv'.format(str(node_count), str(subnet), str(subnet))
                    grn_file = './experiment_data/{}_node/{}/Ecoli-{}.xml'.format(str(node_count), str(subnet), str(subnet))

                    # extracting info
                    mRNA_df = pd.read_csv(mRNA_file, sep='\t', decimal=',')
                    protein_df = pd.read_csv(protein_file, sep='\t', decimal=',')
                    perturbs_df = pd.read_csv(perturbs_file, sep='\t', decimal='.')
                    graph_df = pd.read_csv(graph_file, header=None, sep='\t', decimal=',')
                    genes_dict = sbml_to_dict(grn_file)


                    out1 = self.run_kinetic(genes_dict, mRNA_df, protein_df, perturbs_df, t, experiment)
                    out2 = self.run_boolean(genes_dict, mRNA_df, graph_df, perturbs_df, experiment)
                    
                    # compute error as hamming distance between the two outcomes
                    errors = np.count_nonzero(out1!=out2)
                    experiment_errors.append(errors)
                print('min:\t', np.min(experiment_errors))
                print('max:\t', np.max(experiment_errors))
                print('mean:\t', np.mean(experiment_errors))
                print('median:\t', np.median(experiment_errors))
                print('std:\t', np.std(experiment_errors))
                subnet_errors.append(experiment_errors)
            F, p = f_oneway(*subnet_errors)
            print('================')
            print('f-statistic:\t', F)
            print('p-value:\t', p)
            print('================')
            node_count_dependent_errors.append(subnet_errors)

        return node_count_dependent_errors

    def plot_data(self, data):
        # experiment parameters (should not be changed unless data is generated for it)
        node_counts = [5, 10, 20, 50, 100]
        subnets = [i for i in range(1,11)]
        experiments = [i for i in range(0,10)]

        plt.figure(figsize=(6, 4))
        for j in range(len(subnets)):
            subnet_errors = []
            for i in range(len(node_counts)):
                mean_error = np.mean(data[i][j])
                subnet_errors.append(mean_error)
            plt.plot(node_counts, subnet_errors, 'x', label='subnet {}'.format(j+1))
            plt.legend()
            plt.title('Mean error of Boolean network steady state')
            plt.xlabel('Size of subnet GRN [nodes]')
            plt.ylabel('Error [hamming distance]')
            plt.savefig('./figures/error-analysis.pdf')  
        plt.show()

if __name__ == '__main__':
    experiment = Experiment()
    data = experiment.run()
    print(data)
    experiment.plot_data(data)