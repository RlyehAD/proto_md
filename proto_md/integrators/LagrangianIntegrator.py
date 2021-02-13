'''

@author： Glacier 

'''

import integrator
import proto_md.subsystems.rigid_subsystem as rigid_subsystem
import proto_md.dynamics as dynamics
import numpy as np 

class LagrangianIntegrator(integrator.Integrator):
	def cg_step(self):

	#Perform a forward euler step to extrapolate the cg variable (center of mass)
	#The code here is copied from the code in factorization integrator (18~22)
	#Flatten is a function used to make a 1d vector
	avg_velocities = (self.system.cg_positions[0,:,-1,:] - self.system.cg_positions[0,:,0,:]) / \
			(self.system.config["md_steps"]	- self.system.config["nstxout"])
	avg_velocities = avg_velocities.flatten()[:,np.newaxis]
	dt = self.system.config["dt"]
	cg_translate = dt * avg_velocities * 1000.

	#Get the masses and the positions of the atoms as two column vectors
	masses = self.atoms.masses()[:, np.newaxis]
	pos = self.atoms.positions
	#The outer loop is for subsystems
	#for nss in range(len(self.system.cg_positions[0,:,-1,0]))
		
	updating_cm = self.system.cg_positions[0,:,-1,:]
	updating_cm = updating_cm.flatten()[:, np.newaxis] + cg_translate


	pos = pos + np.transpose(cg_translate)
	#Turns cg_translate to a row vector otherwise cannot broadcast
	#When there's only 1 cg var, this is useless

	spos = pos / self.system.box * 2.0 * pi - pi
	updating_cm = (arctan2(sum(masses * sin(spos), axis=0), sum(masses * cos(spos), axis=0)) + pi) * \
	self.system.box / 2.0 / pi#row vec

	 

		







