import numpy as np
import math 
import itertools
import sys
import os as os
import cPickle as pickle
import pandas as pd







class SignalNet:

    def __init__(self, parameters = {'proteins':4, 'input_grid': 10, 'n_samples': 10**5, 'Emin': -8.0, 'Emax': 2.0} ):

        self.n = parameters['proteins']
        self.pin_grid = parameters['input_grid']
        self.ns = parameters['n_samples']
        self.emin = parameters['Emin']
        self.emax = parameters['Emax']


        self.deltaEbar = 10**-4
        self.delta = 10**-6
        self.y = np.array(list(itertools.product([[0],[1]], repeat = self.n)))[:,:,0]
        
        self.thinit = (self.emax - self.emin) * np.random.random_sample((self.n, 2)) + self.emin
        self.etinit = (self.emax - self.emin) * np.random.random_sample((2, self.n)) + self.emin
        self.pi, self.cor = self.Pin_n_cor()
    

    
    def Pin_n_cor(self):

            
        q1 = np.linspace(0.01, 0.99, num = self.pin_grid)
        q2 = np.linspace(0.01, 0.99, num = self.pin_grid)
        q3 = np.linspace(0.01, 0.99, num = self.pin_grid)

        Q1, Q2, Q3 = np.meshgrid(q1, q2, q3)
        Q1r = np.reshape(np.moveaxis(Q1, 0, 1), self.pin_grid**3)
        Q2r = np.reshape(np.moveaxis(Q2, 0, 1), self.pin_grid**3)
        Q3r = np.reshape(np.moveaxis(Q3, 0, 1), self.pin_grid**3)

        to_delete = np.where(1.0-Q1r-Q2r-Q3r < 0.0)[0]
        Q1r = np.append(np.delete(Q1r, to_delete), 0.25)
        Q2r = np.append(np.delete(Q2r, to_delete), 0.25)
        Q3r = np.append(np.delete(Q3r, to_delete), 0.25)

        C = Q1*(1-Q1-Q2-Q3)-Q2*Q3
        C = np.reshape(np.moveaxis(C, 0, 1), self.pin_grid**3)
        Cr = np.append(np.delete(C, to_delete), 0.0) 


        return [ np.array([ Q1r, Q2r, Q3r, 1.0-Q1r-Q2r-Q3r ]), Cr]
    
    def MIcalculate(self, th, et):
            
        f1  = self.func_f(th, 0)
        f2  = self.func_f(th, 1)
        f12 = self.func_f(th, None)
        g1  = self.func_g(et, 0)
        g2  = self.func_g(et, 1)

        P = np.zeros([2**self.n, 4], dtype = float)
        Q = np.zeros([4, 2**self.n], dtype = float)
        P[0, 0] = 1.0
        Q[0, 0] = 1.0
        Q[3, 2**self.n -1] = 1.0
        

        
        u = list(np.argwhere(self.y[i] == 1) for i in xrange(len(self.y))) 
        v = list(np.argwhere(self.y[i] == 0) for i in xrange(len(self.y)))

     

        for mu in range(1, 2**self.n -1): # loop over all config. of y, except boundary cases
            
            P[mu, 1] =  np.prod(f1[u[mu]][:, 0]) *  np.prod(1.0 - f1[v[mu]][:, 0])
            P[mu, 2] =  np.prod(f2[u[mu]][:, 0]) *  np.prod(1.0 - f2[v[mu]][:, 0])
            P[mu, 3] =  np.prod(f12[u[mu]][:, 0]) *  np.prod(1.0 -f12[v[mu]][:, 0])

        
            Q[0, mu] =  (1 - g1[mu]) * (1 - g2[mu])
            Q[1, mu] =  g1[mu] * (1 - g2[mu])
            Q[2, mu] =  (1 - g1[mu]) * g2[mu]
            Q[3, mu] =  g1[mu] *  g2[mu]

        P[0, 1] =  np.prod(1.0 - f1[v[0]][:, 0])
        P[0, 2] =  np.prod(1.0 - f2[v[0]][:, 0])
        P[0, 3] =  np.prod(1.0 - f12[v[0]][:, 0])
        P[2**self.n -1, 1] = np.prod(f1[v[0]][:, 0])
        P[2**self.n -1, 2] = np.prod(f2[v[0]][:, 0])
        P[2**self.n -1, 3] = np.prod(f12[v[0]][:, 0])

        Pio = np.dot(Q,P) + self.delta
        
        MI = []
        
        for i in range(0, self.pi.shape[1]):
			Po =  np.dot(Pio, self.pi[:,i]) + self.delta
			tmp_val = np.multiply(Pio, np.log2(Pio / Po [:, None]))
			tmp_val = np.nan_to_num(tmp_val)
			MItmp = np.sum(tmp_val, axis = 0)
			MI = np.append(MI,  np.dot(MItmp, self.pi[:,i]))

		
        return MI


    
    def MI_anneal(self):

        th = self.thinit
        et = self.etinit
        MI = self.MIcalculate(th, et)


        istep = 1

    

        while istep < self.ns and (th > self.emin).all() and (th < self.emax).all() and (et > self.emin).all() and (et < self.emax).all():
            
            delta_th = np.random.choice([-1.0, 0.0, 1.0], size = (self.n, 2), p=[1.0/3.0, 1.0/3.0, 1.0/3.0]) * np.random.random_sample((self.n, 2)) * self.deltaEbar
            delta_et = np.random.choice([-1.0, 0.0, 1.0], size = (2, self.n), p=[1.0/3.0, 1.0/3.0, 1.0/3.0]) * np.random.random_sample((2, self.n)) * self.deltaEbar

            MItmp = self.MIcalculate(th + delta_th, et + delta_et)
            if (np.random.rand(1) < np.exp(np.log(istep) * ( MItmp - MI))).all():
                MI = MItmp
                th = th + delta_th
                et = et + delta_et
                istep = istep + 1
            
            else:
                istep = istep + 1
            
        

        return [MI, th, et]

	

    def func_f(self, theta, c):

        if c is None:
            return np.sum(np.exp(-theta), axis = 1)/ ( 1 +  np.sum(np.exp(-theta), axis = 1))
        else:
            return np.exp(-theta)[:,c]/ ( 1 +  np.exp(-theta)[:,c])



    def func_g(self, eta, c):

        return np.sum(np.exp(-eta[c,:])* self.y.astype(float), axis = 1) / ( 1 + np.sum(np.exp(-eta[c,:])* self.y.astype(float), axis =1))







	
#*****************
# batch tools
#*****************


def create_name_short(n_samples, n_pro, n_prog):

	str_r = 'Ns' + str(n_samples) + '_np' + str(n_pro) + '_progindex_' + str(n_prog)

	return str_r


#*****************
# Main program
#*****************

# Directory 
dir_name='./'

# Network parameters
n_proteins =	int(sys.argv[1])   
prog_index =    int(sys.argv[2])
n_samples = 10**6
affinity_lb = -8.0
affinity_up = 2.0
q_grid = 20




parameters = {'proteins': n_proteins,'input_grid': q_grid , 'n_samples': n_samples, 'Emin': affinity_lb , 'Emax': affinity_up}
tn = SignalNet(parameters = parameters)
input_cor = tn.cor
input_pdf = tn.pi
MI_opt, theta_opt, eta_opt = tn.MI_anneal()


df1 = pd.DataFrame({'MI_opt': list(MI_opt), 'input_corr': list(input_cor)})
df2 = pd.DataFrame({'I1': list(theta_opt[:,0]), 'I2': list(theta_opt[:,1]), 'O1': list(eta_opt[0,:]), 'O2': list(eta_opt[1,:])})



df1_name = create_name_short(n_samples, n_proteins, prog_index)  + 'q_n_I.pkl'
df2_name = create_name_short(n_samples, n_proteins, prog_index)  + 'theta_n_eta.pkl'

df1_name = os.path.join(dir_name, df1_name)
df2_name = os.path.join(dir_name, df2_name)

df1.to_pickle(df1_name)
df2.to_pickle(df2_name)
