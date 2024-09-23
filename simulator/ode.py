from model.kinetic import *
import numpy as np


class ODE:
    
    def __init__(self, grn_config, x0, y0):
        """
        
        x0: initial protein concentration
        y0: initial mRNA concentration
        """
        
        self.grn = KineticGRN(grn_config)
        self.x_vecs = [np.array(x0)]
        self.y_vecs = [np.array(y0)]
    
    def reset(self):
        self.x_vecs = [self.x_vecs[0]]
        self.y_vecs = [self.y_vecs[0]]

    def deriv(self, z, t, i):
        """
        Compute the derivative at time t for gene i
        
        z: vector of the form array([x[i], y[i]])
        t: time at which the derivative is evaluated
        i: the index for the i-th gene
        """
        
        grn = self.grn
        
        xi = z[0]
        yi = z[1]
        mi = grn.genes[i].max_transcription
        ri = grn.genes[i].max_translation
        lambda_mrna = grn.genes[i].delta_mRNA_deg
        lambda_prot = grn.genes[i].delta_protein_deg
        fi = grn.genes[i].input_func(self.y_vecs[-1])

        return np.array([mi*fi - lambda_mrna*xi, ri*xi - lambda_prot*yi])
    
    def solve_ode(self, t, perturbs=None, perturb_stop=False):
        """
        Use the standard Runge-Kutta 4 method
        """
        
        grn = self.grn
        
        n = len(t) # number of timesteps
        # reset timeseries in ODE object
        self.reset()
        x_vecs = [self.x_vecs[0]]
        y_vecs = [self.y_vecs[0]]
        
        if perturbs is not None:
            for i in range(len(grn.genes)):
                grn.genes[i].perturb_activation(perturbs[i])

        for j in range(n - 1):
            if j == n//2 and perturb_stop:
                for k in range(len(grn.genes)):
                    grn.genes[k].restore_activation()
                
            x_vec = []
            y_vec = []
            z = np.array([x_vecs[j], y_vecs[j]]).T

            # for each gene evaluate solution at t_j
            for i in range(len(grn.genes)):
                # stepsize
                h = t[j+1] - t[j]
                # rk4 parameters
                k1 = self.deriv(z[i], t[j], i)
                k2 = self.deriv(z[i] + k1 * h / 2, t[j] + h / 2, i)
                k3 = self.deriv(z[i] + k2 * h / 2, t[j] + h / 2, i)
                k4 = self.deriv(z[i] + k3 * h, t[j] + h, i)
                # evaluation of z at timestep j+1
                out = z[i] + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)
                x_vec.append(out[0])
                y_vec.append(out[1])
            # store solution of all mRNA and protein concentrations in the timeseries list
            x_vecs.append(np.array(x_vec))
            y_vecs.append(np.array(y_vec))
        self.x_vecs = x_vecs
        self.y_vecs = y_vecs