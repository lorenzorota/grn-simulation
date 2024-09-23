from model.boolean import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


class Graph:
    
    def __init__(self, config, init_state):
        self.config = config
        self.init_state = init_state
        self.graph, self.layout = self.init_network(config, init_state)
    
    def init_network(self, config, init_state):
        """Initializes a Gene Regulatory Network from the list of triples"""

        # separate edges and weights
        edges = [elem[0:2] for elem in config]
        weights = [elem[2] for elem in config]

        G = nx.DiGraph()
        G.add_edges_from(edges)

        # assign weights to edges with values:
        # True: activation
        # False: repression
        for tup in config:
            if tup[2] == '+':
                G[tup[0]][tup[1]]['weight'] = True
            else:
                G[tup[0]][tup[1]]['weight'] = False

        # True: means active
        # False: means inactive
        node_states = {k : init_state[k] for k in G.nodes()}
        nx.set_node_attributes(G, node_states, 'state')
        # the shell layout is arguably the nicest one they offer
        nodes_pos = nx.shell_layout(G)
        
        return G, nodes_pos
        
    def draw(self):
        """Generates the visual representation of the network"""
        
        G = self.graph
        pos = self.layout
        
        for e in G.edges():
            if int(G[e[0]][e[1]]['weight']):
                nx.draw_networkx_edges(G, pos=pos, arrows=True, edgelist=[e], arrowstyle='-|>')
            else:
                nx.draw_networkx_edges(G, pos=pos, arrows=True, edgelist=[e], arrowstyle='-[', edge_color='r')

        # define the color scheme for active and inactive nodes
        color_table = {0: "grey", 1: "green"}
        node_colors = [color_table[G.nodes[i]['state']] for i in G.nodes()]

        # draw the graph and labels
        _ = nx.draw_networkx_nodes(G, pos=pos, node_color=node_colors)
        _ = nx.draw_networkx_labels(G, pos=pos)
    
    def reset(self):
        self.graph, self.layout = self.init_network(self.config, self.init_state)


class Simulator:
    
    def __init__(self, genes_dict, network_config, init_states):
        self.init_states = init_states
        self.graph = Graph(network_config, init_states)
        self.grn = BooleanGRN(genes_dict)
    
    def apply_function(self):
        """ Input function """
        
        G = self.graph.graph
        
        updated_states = []
        for node in G.nodes():
            update = self.grn.genes[node].input_function(G)
            updated_states.append(update)

        # apply the state transition synchronously (namely at the end)
        for node in G.nodes():
            G.nodes[node]['state'] = updated_states.pop(0)
    
    def run_simulation(self, perturbs=None, visualize=False, max_epoch=100):

        # initialize simulation
        matched = False
        epoch_sleep = 1
        self.graph.reset()
        G = self.graph.graph
        
        # create initial matrix
        matr = np.array([[G.nodes[node]['state'] for node in G.nodes()]]).T

        # apply perturbations (if any)
        if perturbs is not None:
            for idx, gene in enumerate(self.grn.genes):
                self.grn.genes[gene].perturb_activation(perturbs[idx])
        
        # simulate the GRN
        for epoch in range(max_epoch):
            if visualize:
                self.graph.draw()
                plt.title('epoch {}'.format(epoch))
                plt.savefig('./figures/example-boolean({}).pdf'.format(str(epoch))) 
                plt.pause(epoch_sleep)

            self.apply_function()
            new_states = [G.nodes[node]['state'] for node in G.nodes()]
            new_vec = np.array([new_states])
            matr = np.concatenate((matr, new_vec.T), axis=1)

            # check for elementwise for each column if vector corresponds
            equiv_matr = np.array([np.equal(col, new_vec[0]) for col in matr.T]).T

            # check for oscillation or stable convergence
            last_idx = equiv_matr.shape[1] - 1
            for idx, col in enumerate(equiv_matr.T):
                # check if new_vec is equivalent to penultimate column
                if np.all(col) and idx == last_idx - 1:
                    # we have stable convergence
                    if not visualize:
                        matched = True
                        break
                    else:
                        self.graph.draw()
                        plt.title('epoch {}, detected fixed point attractor'.format(epoch+1))
                        plt.savefig('./figures/example-boolean({}).pdf'.format(str(epoch+1))) 
                        plt.pause(epoch_sleep)
                        matched = True
                        break

                # check if new_vec appears anywhere except for the last column
                elif np.all(col) and idx != last_idx:
                    # we have an oscillation
                    if not visualize:
                        matched = True
                        break
                    else:
                        self.graph.draw()
                        plt.title('epoch {}, detected limit cycle attractor'.format(epoch+1))
                        plt.savefig('./figures/example-boolean({}).pdf'.format(str(epoch+1)))
                        plt.pause(epoch_sleep)
                        matched = True
                        break

            if matched == True:
                break

        plt.show()
        matr = matr.astype('uint8')
        return matr