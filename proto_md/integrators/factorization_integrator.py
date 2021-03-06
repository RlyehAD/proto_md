'''
Created on Jan 30, 2013

@author: andy
'''

import integrator
import proto_md.dynamics as dynamics
import numpy as np

class FactorizationIntegrator(integrator.Integrator):
    def cg_step(self):
        """
        perform a forward euler step (the most idiotic and unstable of all possible integrators)
        to integrate the state variables. 
        
        reads self.cg_momenta to advance the cg state variables. 
        
        X[n+1] = X[n] + dt*dX/dt[n], and dX/dt is cg_moment/cg_mass.
        """
        
        # forward euler 
        # cg_velocities: nensemble x n subsystem x n_step x n_cg
        
        # inner mean: average over ensembles - results in a  n subsystem x n_step x n_cg array
        # next mean: average over steps - results in a n subsystem x n_cg array
        # flatten(): makes a 1D array 
        # [:,np.newaxis]: makes a 1xN column vector.
	# this used to be: avg_velocities = np.mean(np.mean(self.system.cg_velocities, axis = 0), axis = 1).flatten()[:,np.newaxis]

	# this is hackish ~ but more reliable than using the very chaotic velocities
	avg_velocities = (self.system.cg_positions[0,:,-1,:] - self.system.cg_positions[0,:,0,:]) / \
			(self.system.config["md_steps"]	- self.system.config["nstxout"])

	avg_velocities = avg_velocities.flatten()[:,np.newaxis]

	# this is hackish. The 1000 factor should be replaced with 'dt' from the mdp file in case the user
	# takes more than 1fs timestep

        cg_translate = self.system.dt * avg_velocities * 1000.0 # (ps to fs)
	
        self.system.translate(cg_translate)
