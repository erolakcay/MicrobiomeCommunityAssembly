#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python script that simulates microbial communities.
Eeman Abbasi and Erol Akçay, 2022. "Host control and species interactions jointly determine
microbiome community structure"
Extention of code implemented by Qian and Akçay 2020.
"""

import os
os.environ["OMP_NUM_THREADS"]="1"
os.environ["MKL_NUM_THREADS"]="1"
os.environ["NUMEXPR_NUM_THREADS"]="1"
import numpy as np
from numpy import random
from scipy.stats import halfnorm
from scipy import integrate
from statsmodels.tsa.stattools import adfuller
#import statsmodels.tsa.stattools as ts
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import itertools
import sys
import warnings
import math
from random import sample



class Community(object):
	"""
	An ecological community whose assembly is simulated over time. The class
	is initialized by passing in eleven mandatory parameters and some optional parameters.
	"""

	def __init__(self, C, Pc, Pm, h, kappa,immune_microbial_load, immuneResponse,change_kappa,new_kappa,repeat_num,path_to_save,plotbool=False, plot_failed_eq=False,
				 study_selection=False, relative_eq_criteria=False,
				 IIM_bool = False, VUM_bool=False,
				 lamb = 0.5):
		"""
		Takes in parameters that are varied across communities.
		Initialize default parameters of the community as attributes.s
		Input:
			C (float): desired connectivity of the community variable_r
			Pc (float): proportion of competitive interactions
			Pm (float): proportion of mutualistic interactions
			h (float): half-saturation constant for type II functional response
			plotbool (boolean): whether or not to plot populations over time
			plot_failed_eq (boolean): whether to plot dynamics at equilibria that do not converge before the time limit
			study_selection (boolean): set True to run the analysis of selection as described for Figure 5, S14-17
			relative_eq_criteria (boolean): set True to use equilibration criteria that is relative to population sizes
			IIM_bool (boolean): set True to implement the interchangeable interactions model, otherwise use unique interactions model (make sure VUM_bool set to False)
			variable_r (boolean): set True to vary r between -1 and 1 based on the proportion of predatory interactions
			VUM_bool (boolean): set True to implement the variable uniqueness model (make sure IIM_bool is set to False)
			lamb (float): uniqueness coefficient (lambda) in the VUM
			immuneResponse (boolean): indicates if the community is exposed to host immune control
			kappa (int): the host immune control strength, or constraint on microbial species growth
			repeat_num (int): the number repeat simulation runs for a given microbial community 
			change_kappa (boolean): determines whether the kappa value will be updated once the microbial community has reached the steady state 
			new_kappa (int): new kappa value imposed on the microbial community 
			original_kappa (int): the kappa value that the microbial community is exposed to at the start of the simualtion. 
			path_to_save (string): reference folder where all the simulation outputs will be saved
		"""
		# Simulation parameters that are varied across communities. Taken as inputs
		self.C = C
		self.Pc = Pc
		self.Pm = Pm
		self.halfsat = h
		self.plotbool = plotbool
		self.plot_failed_eq = plot_failed_eq
		self.study_selection = study_selection
		self.relative_eq_criteria = relative_eq_criteria
		self.IIM_bool = IIM_bool
		self.immune_microbial_load = immune_microbial_load
		self.VUM_bool = VUM_bool
		self.lamb = lamb
		self.immuneResponse = immuneResponse
		self.kappa = kappa
		self.repeat_num = repeat_num
		self.change_kappa = change_kappa
		self.new_kappa = new_kappa
		self.original_kappa = self.kappa
		self.path_to_save = path_to_save
		 

		# Additional simulation parameters and variables
		self.S = 10 # number of initial species
		self.numspecies = self.S # number of total species, this changes over time
		self.numspecies_eq = [self.numspecies] # list of species richness at each equilibrium
		self.selfinteraction = -1.0 # self-regulation coefficient (diagonals of the interaction matrix)
		sigma = 0.5 # half-normal distribution scale parameter
		self.K = 100 # carrying capacity of each species
		sigma_c = sigma / self.K # half-normal distribution scale parameter for competitive interaction
		self.strengthdist = halfnorm(scale=sigma) # half-normal distribution for mutualistic and exploitative interactions
		self.strengthdist_comp = halfnorm(scale=sigma_c) # half-normal distribution for competitive interactions
		self.r = 1 # intrinsic growth rate of all species
		self.eqcriteria = 0.001 # equilibration criteria for changes in population
		self.eqcriteriamatrix = np.full((1,self.S),self.eqcriteria) # 1 by S matrix of equilibrium criteria
		self.extinctthres = 0.01 # extinction threshold for population sizes
		self.timelim = 5000 # time limit to wait for equilibrium
		self.failedinvasionlim = 1000 # failed invasion threshold
		self.numspecieslim = 750 # limit on species richness
		self.steadytime = 0 # equilibrium at which species richness converges
		self.pvalthres = 0.01 # p-value threshold to reject null hypothesis
		self.window = 1000 # number of equilibria analyzed in augmented Dickey-Fuller test
		self.poststeadywindow = 1500 # number of equilibria simulated after steady state is reached
		self.t_start = 1. # simulation start time
		self.t_end = 10**9 # time limit of simulation; this is never reached
		self.currenttime=self.t_start
		self.integratedtime = self.t_start
		self.extinctiontimes = [] # list of extinction times
		self.numextinct = 0 # number/counter of total extinctions
		self.eqs = [0] # list of times that the system is at equilibrium
		self.failedinvasiontimes = [] # list of times a failed invasion occurs
		self.failedinvasionsteady=0 # number of failed invasions during the steady state
		self.extinctionsteady=0 # number of extinctions during the steady state
		self.numspeciesatsteady=0 # species richness when steady state is reached
		introtime=np.full((1,self.S), 0) # 1 by S matrix describing introduction times of species
		self.introtime = introtime.flatten()
		self.ages=[] # list of the species persistence for every species that existed
		self.ages_steady=[] # list of species persistence at the steady state 
		self.wait=[] # times the time limit to wait for equilibrium is reached
		self.num_invasions_entire = 0 # num invasions for the entire community history 
		self.num_invasions_steady = 0 # num invasions that take place at the steady state 
		self.competitive_invasion_condition = []  #list of competitive invasion coefficent scaled by competitive species abundance in the community 
		self.mutualistic_invasion_condition = []  #list of mutualistic invasion coefficent scaled by mutualistic species abundance in the community 
		self.exploitative_invasion_condition_positive = [] #list of exploitative invasion coefficent scaled by exploitator species abundance in the community 
		self.exploitative_invasion_condition_negative = [] #list of exploitative invasion coefficent scaled by exploitated species abundance in the community 
		self.popSteady = [] #list of community population at steady state 
		self.hostControl = 0 #global host immune control imposed on the community 
		


		#self.mut_int = 0
		#self.comp_int = 0
		#self.exp_int = 0
		#self.sample_mut_int = 0 
		#self.sample_comp_int = 0
		#self.sample_exp_int = 0

		# variables to be saved once an invasion event occurs 
		#self.extinct_inv = [] 
		#self.host_control_inv = [] 
		#self.comm_pop = [] 
		#self.A_comp_beforeInv =[]  
		#self.A_mut_beforeInv = []
		#self.A_comp_afterInv = []
		#self.A_mut_afterInv = []


		# initializing the community property variables for the microbial community/we use these for plotting 
		self.invasion_steady = [] # list of invasions at steady state 
		self.inv_entire = [] #list of invasions taking place during entire community history 
		self.growth_rate = [] #list of intrinsic growth rate of each species
		self.host_control = [] #list of host immune control imposed on the microbial community
		self.species_list = [] #list of species richness at steady state 
		self.comm_abun = [] #list of community abundance at entire community history 
		self.mean_species_persistence_steady = []  #list of mean species that persisted during the steady state 
		self.mean_species_persistence = [] #list of mean species that persisted during entire community history 
		self.median_species_persistence = [] #list of median species that persisted during the steady state 
		self.median_species_persistence_steady = [] #list of median species that persisted during entire community history 
		self.extinction_steady = [] #list of num extinctions taking place at the steady state 
		self.extinction_entire = [] #list of num extinctions taking place during the entire community history
		self.competitive = [] #list of mean of competitive_invasion_condition
		self.mut_invasion_condition = [] #list of mean of mutualistic_invasion_condition
		self.exp_invasion_condition_positive = [] #list of mean of exploitative_invasion_condition_positive
		self.exp_invasion_condition_negative = [] #list of mean of exploitative_invasion_condition_negative
		self.exp_invasion_condition = [] #list of sum of exp_invasion_condition_positive and exp_invasion_condition_negative
		self.failed_invasions = [] # list of mean of failedinvasiontimes
		self.pop_size_steady = [] # list of mean of popSteady

		self.start_eq = 0  
		self.kappa_eq = True #set to false when equlibria == 1500 and change_kappa == True to allow community to reach equilbrium. 
	


		if plotbool:
			self.time = [self.t_start]
		if plot_failed_eq:
			self.failed_eq_counter = 0
		if study_selection:
			self.sampledmatrices = [] # list of sampled interaction matrices
			self.sampledinvaderrows=[] # list of interaction matrix rows of sampled successful invaders
			self.sampledinvadercols=[] # list of interaction matrix cols of sampled successful invaders
			self.sampledextinctionrows=[] # list of interaction matrix rows of sampled species that go extinct
			self.sampledextinctioncols=[] # list of interaction matrix cols of sampled species that go extinct

			# sampling scheme
			self.inv_ext_sampleinterval = 10 # sampling interval as described in Methods
			self.numcommunitysamples=100 # how many matrices we sample after steady state
			self.ext_samplecounter = 0 # counter to sample species that go exctinct with given interval
		if relative_eq_criteria:
			self.rel_eq_criteria = 1/10**6

		# create integration object
		ode = integrate.ode(self.eq_system)
		ode.set_integrator('dopri5',nsteps=1)
		warnings.filterwarnings("ignore", category=UserWarning) # ignore warnings of integrating step by step
		self.ode = ode

		# Boolean attributes reporting whether various simulation limits are reached
		self.steady = False # boolean indicating whether community has come to the steady state
		self.broken = False # boolean indicating whether community has come to un-invadable equilibrium or exceeds limit on species richness

	def start(self):
		"""
		Construct the initial community, create uninitialized interaction
		matrix, initialize interaction matrix, and impose desired connectivity
		onto interaction matrix. Must call this method after initialization.
		"""
		S = self.S
		# construct the initial community
		X0 = np.random.rand(S) # matrix describing initial population size of each species
		X0 = X0 * 10 # initialize populations of initial species with uniform random samples from 0 to 10
		#self.R = np.full((S,1),self.r) # S by 1 matrix describing intrinsic growth rate of each species
	   
		self.population = X0 # current population
		self.prevpop = X0 # previous population
		if self.plotbool:
			self.populations = X0 # save full population histories
		self.R  = np.full((self.S,1),self.r) # so no change once we start with the community
		#self.R = self.popDenGrowthRate()

		# create uninitialized interaction matrix
		A = np.zeros((S,S)) # uninitialized S by S interaction matrix, a_ij is j's effect on i
		A_c = np.zeros((S,S)) # uninitialized interaction matrix with competitive terms only
		A_m = np.zeros((S,S))  # uninitialized interaction matrix with mutualistic terms only
		A_e_pos = np.zeros((S,S))  # uninitialized interaction matrix with positive exploitative terms only
		A_e_neg = np.zeros((S,S))  # uninitialized interaction matrix with negative exploitative terms only

		# initialize interaction matrix, impose desired connectivity onto interaction matrix
		for i in range(0,S):

			num_interaction = 0
			num_mutualistic = 0
			num_predatory = 0

			for j in range(i+1,S): # shouldn't be range(i,S) since j and i shouldn't be equal
				if np.random.uniform(0,1) <= self.C: # then there is an interaction between i and j
					randnum = np.random.uniform(0,1)
				   
					num_interaction = num_interaction + 1
					if randnum >= (1-self.Pm): # then mutualistic interaction
						# inverse transform sampling: use percent point fxn of the half normal
						A[i,j] = self.strengthdist.ppf(np.random.uniform(0,1)) # interaction strength
						A[j,i] = self.strengthdist.ppf(np.random.uniform(0,1)) # not necessarily symmetric

						A_m[i,j] = A[i,j]
						A_m[j,i] = A[j,i]

						# update the mut tracker for both species i and j
						num_mutualistic = num_mutualistic + 1
					   
					elif randnum <= self.Pc: # then competitive interaction
						A[i,j] = -1 * self.strengthdist_comp.ppf(np.random.uniform(0,1))
						A[j,i] = -1 * self.strengthdist_comp.ppf(np.random.uniform(0,1))

						A_c[i,j] = A[i,j]
						A_c[j,i] = A[j,i]

					else: # exploitative interaction
						if np.random.uniform(0,1) <= 0.5: # i benefits, j suffers
							A[i,j] = self.strengthdist.ppf(np.random.uniform(0,1))
							A[j,i] = -1 * self.strengthdist.ppf(np.random.uniform(0,1))

							A_e_pos[i,j] = A[i,j]
							A_e_neg[j,i] = A[j,i]

							# update the pred tracker for both species i and j
							num_predatory = num_predatory + 1


						else: # i suffers, j benefits
							A[i,j] = -1 * self.strengthdist.ppf(np.random.uniform(0,1))
							A[j,i] = self.strengthdist.ppf(np.random.uniform(0,1))

							A_e_neg[i,j] = A[i,j]
							A_e_pos[j,i] = A[j,i]



		np.fill_diagonal(A, self.selfinteraction ) # all self-interactions are set to the same value
		

		self.A = A
		self.A_c = A_c
		self.A_m = A_m
		self.A_e_pos = A_e_pos
		self.A_e_neg = A_e_neg
		self.mut_species_track = np.count_nonzero(self.A_m,axis = 1)
		self.pred_species_track = np.count_nonzero(self.A_e_pos, axis = 1)

		# populate the meu based on the method
		

	def immuneDenDependence(self, population_list, kappa, immune_microbial_load):

		"""
		Determine the growth of each microbial species in the community as a function of abundance of 
		the entire microbial community. 
		Returns the 1 by S matrix of population growth, and the host immune control imposed on the microbial 
		community 
		"""
		pop = population_list
		sum_pop = np.sum(pop)
		ratio = sum_pop/immune_microbial_load # host control imposed on all microbial species 
		ratio_kappa = (kappa * ratio)
		growth_rate = 1 - ratio_kappa

		return growth_rate,ratio_kappa


	def eq_system(self,t,X):
		"""
		Solve the differential equation over t for initial condition X.
		X is a 1 by S array.
		Returns the 1 by S matrix of population growth.
		"""
		Xnew = X.reshape((-1,1)) # make it into 2D S by 1 matrix
		X_m = Xnew / (self.halfsat + Xnew)
		output =  Xnew * (self.R  + self.selfinteraction * Xnew / self.K + \
						 np.dot(self.A_c, Xnew)+ np.dot(self.A_m, X_m) + \
						 np.dot(self.A_e_pos, X_m) + \
						 np.dot(self.A_e_neg, Xnew)/(self.halfsat + Xnew))
		output = output.T # transpose back to 1 by S
		return output.flatten()

	def introduce(self, currenteq):
		"""
		Add new species once equilibrium is reached and iterate until
		next equilibrium. Input is integer number of the current equilibrium.
		"""

		if self.numspecies>self.numspecieslim:
			self.broken = True

		if self.broken:
			return

		self.currenttime = self.integratedtime

		if self.currenttime >= self.t_end: # do not iterate if t_end has already been reached
			return
		if self.currenttime > self.t_start: # don't add species/check for extinction if first iteration

			# check for extinction, and remove the species immediately if extinct
			boolmask = self.population > self.extinctthres

			# extract the ages of extinct species
			extinctbirth = self.introtime[np.logical_not(boolmask)]

			# remove extinct species
			self.population = self.population[boolmask]
			self.prevpop = self.prevpop[boolmask]
			self.introtime = self.introtime[boolmask]
			self.R = self.R[boolmask]

			if self.plotbool:
				self.populations = self.populations[boolmask]

			# update running total of number of extinctions
			newextinct = self.numspecies - self.population.shape[0]
			self.numextinct = self.numextinct + newextinct

			if self.steady == True:
				self.extinctionsteady = self.extinctionsteady + newextinct

				if self.study_selection:
					self.ext_samplecounter = self.ext_samplecounter + newextinct

					# sample the species that go extinct
					if self.ext_samplecounter >= self.inv_ext_sampleinterval:
						A_extinct_rows = self.A[np.logical_not(boolmask)]
						A_extinct_cols = self.A.T[np.logical_not(boolmask)]
						self.sampledextinctionrows.append(A_extinct_rows[0])
						self.sampledextinctioncols.append(A_extinct_cols[0])
						self.ext_samplecounter = self.ext_samplecounter % self.inv_ext_sampleinterval

			self.numspecies = self.population.size
			self.eqcriteriamatrix = np.full((1,self.numspecies),self.eqcriteria)

			# remove the extinct species rows from the interaction matrix
			self.A = self.A[boolmask]
			self.A_c = self.A_c[boolmask]
			self.A_m = self.A_m[boolmask]
			self.A_e_pos = self.A_e_pos[boolmask]
			self.A_e_neg = self.A_e_neg[boolmask]
			self.mut_species_track = self.mut_species_track[boolmask]
			self.pred_species_track = self.pred_species_track[boolmask]

			# remove the extinct species cols from interaction matrix
			self.A = self.A.T[boolmask]
			self.A_c = self.A_c.T[boolmask]
			self.A_m = self.A_m.T[boolmask]
			self.A_e_pos = self.A_e_pos.T[boolmask]
			self.A_e_neg = self.A_e_neg.T[boolmask]

			# invert back after indexing the target col as a row in transpose
			self.A = self.A.T
			self.A_c = self.A_c.T
			self.A_m = self.A_m.T
			self.A_e_pos = self.A_e_pos.T
			self.A_e_neg = self.A_e_neg.T


			if newextinct > 0:
				self.extinctiontimes.append(self.ode.t)
				self.ages.extend(currenteq - extinctbirth)# save the persistence of extinct species
				if self.steady == True:
					self.ages_steady.extend(currenteq - extinctbirth)

			

			# Changing the community growth rate, adding the invader population size beforehand and determining the growth rate
			new_population = []
	   
			invade_newpop = 10 * np.random.uniform(0,1) # initial population of new species (rare in abundance)
			if self.immuneResponse == True:
				new_population = np.append(self.population, invade_newpop)
				output = self.immuneDenDependence(new_population, self.kappa, self.immune_microbial_load)# vary for  10^4, 10^15
				self.hostControl = output[1]
				new_R = output[0]
				if currenteq == 1501:
					print("updating R with pop 1501: " , np.sum(self.population))
				if currenteq == 1502:
					print("updating R with pop 1502 : " , np.sum(self.population))
				self.R = np.full((self.population.size,1),new_R)
			else:
				new_R = self.r

			if self.kappa_eq == True:

				invasionfailed = True
				failedcounter = 0 # number of failed invasions at this particular equilibrium
				num_invasions = 0 # resetting for this particular eq.
				failed_invasion_counter = False


				while invasionfailed: # create a new species, and check if it will invade

					# initialize
					newrow = np.zeros((1,self.numspecies)) # new row to be added into A
					newcol = np.zeros((self.numspecies+1,1)) # new column to be added into A, extra term because hstack after newrow
					newrow_c = np.zeros((1,self.numspecies)) # new row only with incoming competitive interactions
					newrow_m = np.zeros((1,self.numspecies)) # new row only with incoming mutualistic interactions
					newrow_e_pos = np.zeros((1,self.numspecies)) # new row only with incoming positive exploitative interactions
					newrow_e_neg = np.zeros((1,self.numspecies)) # new row only with incoming negative exploitative interactions
					newcol_c = np.zeros((self.numspecies+1,1))
					newcol_m = np.zeros((self.numspecies+1,1))
					newcol_e_pos = np.zeros((self.numspecies+1,1))
					newcol_e_neg = np.zeros((self.numspecies+1,1))
				   

					for j in range(0,self.numspecies):

							if np.random.uniform(0,1) <= self.C: # then there is an interaction between i and j
								
								randnum = np.random.uniform(0,1)
								

								if randnum >= (1-self.Pm): # then mutualistic interaction
										newrow[0,j]= self.strengthdist.ppf(np.random.uniform(0,1))
										newcol[j,0]= self.strengthdist.ppf(np.random.uniform(0,1))

										newrow_m[0,j]= newrow[0,j]
										newcol_m[j,0]= newcol[j,0]



								elif randnum <= self.Pc: # then competitive interaction
										newrow[0,j]= -1 * self.strengthdist_comp.ppf(np.random.uniform(0,1))
										newcol[j,0]= -1 * self.strengthdist_comp.ppf(np.random.uniform(0,1))

										newrow_c[0,j]= newrow[0,j]
										newcol_c[j,0]= newcol[j,0]


								else: # exploitative interaction
									if np.random.uniform(0,1) <= 0.5: # i benefits, j suffers
										newrow[0,j]= self.strengthdist.ppf(np.random.uniform(0,1))
										newcol[j,0]= -1 * self.strengthdist.ppf(np.random.uniform(0,1))

										newrow_e_pos[0,j]= newrow[0,j]
										newcol_e_neg[j,0]= newcol[j,0]


									else: # i suffers, j benefits
										newrow[0,j]= -1 * self.strengthdist.ppf(np.random.uniform(0,1))
										newcol[j,0]= self.strengthdist.ppf(np.random.uniform(0,1))

										newrow_e_neg[0,j]= newrow[0,j]
										newcol_e_pos[j,0]= newcol[j,0]




					# Check to see if invasion can occur:
					#creating a temporary interaction matrix G

					xjcol = self.population
					new_r =  new_R # self.r
					X_m = xjcol / (self.halfsat + xjcol)
					# invasion condition
					if (new_r + np.dot(newrow_c, xjcol) + np.dot(newrow_m, X_m) + \
						np.dot(newrow_e_pos, X_m) + \
						np.dot(newrow_e_neg, xjcol)/self.halfsat) > 0:
						invasionfailed = False
						#print("invasion_pass:", new_r)
						if self.steady == True:
							self.num_invasions_steady +=1
						self.num_invasions_entire += 1
					 
					# Keeping track of the interaction type for each invader 
						

					 # keeping track of the assmebly of the species in the community 

						self.competitive_invasion_condition.append(np.dot(newrow_c, xjcol))
						self.mutualistic_invasion_condition.append(np.dot(newrow_m, X_m))
						self.exploitative_invasion_condition_positive.append(np.dot(newrow_e_pos, X_m))
						self.exploitative_invasion_condition_negative.append( np.dot(newrow_e_neg, xjcol)/self.halfsat)

						self.comm_pop.append(np.sum(xjcol))
						self.host_control_inv.append(self.hostControl)
						self.extinct_inv.append(newextinct)
						self.A_comp_beforeInv.append(np.count_nonzero(self.A_c,axis = 1))
						self.A_mut_beforeInv.append(np.count_nonzero(self.A_m,axis = 1))
						
						

					else:

						if self.change_kappa == True:

							if currenteq > 1500 and currenteq != self.steadytime:
								print("intermission invasion")
								break 


							self.failedinvasiontimes.append(self.currenttime) # track the times of failed invasions
							failedcounter += 1

							if self.steady == True:
								self.failedinvasionsteady += 1

							#in the case of an uninvadable eq when change_kappa == True, reset the failed invasion counter 
							if failedcounter > self.failedinvasionlim:

								print("uninvadable equilibrium")
								print(self.kappa, currenteq, self.R[-1], np.sum(new_population),np.sum(X_m),self.numspecies)
								print("failed_invasion", new_r)
								
								if self.change_kappa == True and failed_invasion_counter == False:
									
									failedcounter = 0 
									failed_invasion_counter = True
									print("counter is being reset")

									if self.kappa == self.original_kappa:
										
										self.kappa = self.new_kappa
										self.steady = False 

										if self.kappa == 0:
											self.immuneResponse = False 
											self.R = np.full((self.population.size,1),self.r)

										else:
											self.immuneResponse = True
									     
								else:
									self.broken = True
									return
						else:

							self.failedinvasiontimes.append(self.currenttime)
							failedcounter = failedcounter + 1

							if self.steady == True:
								self.failedinvasionsteady +=1

							if failedcounter > self.failedinvasionlim:
								self.broken = True
								print("uninvadable eq")
								return

				# add new species into interaction matrix once invasion condition is met
				newcol[self.numspecies,0] = self.selfinteraction
				self.A = np.vstack((self.A,newrow))
				self.A = np.hstack((self.A,newcol))

				self.A_c = np.vstack((self.A_c,newrow_c))
				self.A_c = np.hstack((self.A_c,newcol_c))
				self.A_m = np.vstack((self.A_m,newrow_m))
				self.A_m = np.hstack((self.A_m,newcol_m))
				self.A_e_pos = np.vstack((self.A_e_pos,newrow_e_pos))
				self.A_e_pos = np.hstack((self.A_e_pos,newcol_e_pos))
				self.A_e_neg = np.vstack((self.A_e_neg,newrow_e_neg))
				self.A_e_neg = np.hstack((self.A_e_neg,newcol_e_neg))
				self.mut_species_track = np.count_nonzero(self.A_m,axis = 1)
				self.pred_species_track = np.count_nonzero(self.A_e_pos,axis = 1)
				self.A_mut_afterInv.append(self.mut_species_track)
				self.A_comp_afterInv.append(np.count_nonzero(self.A_c,axis = 1))

				if (self.study_selection and self.steady == True and
					((currenteq - self.steadytime) % self.inv_ext_sampleinterval == 0)):
					newrow = np.hstack((newrow.flatten(),self.selfinteraction))
					self.sampledinvaderrows.append(newrow)
					self.sampledinvadercols.append(newcol)

				#newpop = 10 * np.random.uniform(0,1) # initial population of new species
				newpop = invade_newpop
				self.prevpop = np.append(self.population,0)
				self.population = np.append(self.population, newpop)

				if self.plotbool:
					newrowpop = np.zeros((1,self.populations.shape[1]))
					newrowpop[0,self.populations.shape[1]-1]= newpop
					self.populations = np.vstack((self.populations, newrowpop))

				self.introtime = np.append(self.introtime, currenteq)
				self.eqs.append(self.currenttime)

				# extend matrices for eqcriteria and R and extinction threshold
				self.eqcriteriamatrix=np.append(self.eqcriteriamatrix, self.eqcriteria)
				self.R = np.vstack((self.R,np.array([new_r])))
				self.numspecies += 1
				self.numspecies_eq.append(self.numspecies)
				lastpop = self.population
			else:
				lastpop = self.population
				print("No invasions took place: ", currenteq, "self.kappa_eq: ", self.kappa_eq)


		else: # this is the first iteration of simulation
			lastpop = self.population

		if currenteq == 1502 and self.kappa_eq == False:
			self.kappa_eq = True 

		if self.plot_failed_eq:
			self.eq_population = self.population
			self.this_eq_time = [self.currenttime]

		# integrate population dynamics until the next equilibrium
		self.ode.set_initial_value(lastpop.flatten(),self.currenttime)
		while (self.ode.t<self.t_end and (self.ode.t - self.currenttime)<self.timelim):
			if self.ode.t > self.t_start: # do not check for equilibrium/extinction on first pass
				check = np.absolute(self.population - self.prevpop)
				if self.relative_eq_criteria:
					relative_criteria = self.population * self.rel_eq_criteria
					relative_criteria = np.clip(relative_criteria, a_min=self.eqcriteria, a_max=None)
					check = check < relative_criteria
				else:
					check = check < self.eqcriteriamatrix
				if np.all(check): # if all growth rates are less than equilibrium criteria
					if currenteq == 1500:
						print("for 1500:" , np.sum(self.population))
					if currenteq == 1501:
						print("for 1501 : ", np.sum(self.population))
					if currenteq == 1502:
						print("for 1502: ", np.sum(self.population))

					return

			self.ode.integrate(self.t_end)
			self.integratedtime = self.ode.t
			self.prevpop = self.population
			self.population = self.ode.y # update population
			if (self.ode.t - self.currenttime) > self.timelim: # did not converge within time limit
				self.wait.append(currenteq)

				if self.plot_failed_eq: # plot the dynamics
					self.plot_dynamics_single_eq()
					self.failed_eq_counter = self.failed_eq_counter+1

			if self.plotbool:
				self.time.append(self.ode.t)
				self.populations = np.column_stack((self.populations,self.ode.y))
			if self.plot_failed_eq:
				self.this_eq_time.append(self.ode.t)
				self.eq_population = np.column_stack((self.eq_population,self.ode.y))


	def iterate(self, limit=5000):
		"""
		Introduce species for the specified number of iterations.
		Input (optional): limit on the number of equilibria.
		Saves data as Numpy .npz file.
		"""

		self.iterationslimit=limit
		for i in range(0,limit):
			
			self.introduce(i)

			if self.broken: # has reached an unsteady state
				break
				

			if self.change_kappa == False:

				if self.steady == False:
					if (i >= self.window - 1): # start checking for steady state
						series=np.array(self.numspecies_eq)[(i+1-self.window):(i+1)]
						pval = adfuller(series)[1] # calculate p-value of augmented Dickey-Fuller test
						if pval< self.pvalthres:
							self.steadytime = i # equilibrium number in which steady state has been reached
							print(self.steadytime)
							self.steady = True
							self.numspeciesatsteady = self.numspecies
							self.popSteady = self.population

				if self.steady == True:

					if i >= self.steadytime+self.poststeadywindow: # stop the simulation
						break


			else: 

				if self.steady == False:
					if (i >= self.window - 1 + self.start_eq): # start checking for steady state at 2500 , so compare between 1501 - 2501
						series=np.array(self.numspecies_eq)[(i+1-self.window):(i+1)]
						pval = adfuller(series)[1] # calculate p-value of augmented Dickey-Fuller test
						if pval< self.pvalthres:
							self.steadytime = i # equilibrium number in which steady state has been reached
							print(self.steadytime)
							self.steady = True
							print('steady state is true')
							self.numspeciesatsteady = self.numspecies
							self.popSteady = self.population
						
				if self.steady == True:
					if(i >= 1500  and self.kappa == self.original_kappa):
						self.kappa = self.new_kappa
						self.steady = False
						self.start_eq = 1501 # updating the window to check the steady state 
						print('steady state has been set to false', i)
						if self.new_kappa == 0:
							self.immuneResponse = False 
							self.R = np.full((self.population.size,1),self.r)
						else:
							self.immuneResponse = True # this is imp if we start with no immune response, when kappa = 1

						if i == 1500:
							self.kappa_eq = False
					
					if(i>=3000 and i<=3500):

						self.commProperty() # append the last 500 simuulations results 

						if(i==3500):
							np.savez(self.path_to_save+'/last_500_simm/output_'+str(self.repeat_num) +'/output_'+str(self.kappa)+'_'+str(self.new_kappa)+"_"+str(self.Pc)+'_'+str(self.Pm)+'_'+str(self.halfsat),\
							species_list = np.mean(self.species_list), invasion_steady = np.mean(self.invasion_steady),inv_entire = np.mean(self.inv_entire),\
							extinction_steady = np.mean(self.extinction_steady), extinction_entire = np.mean(self.extinction_entire), mean_species_persistence = np.mean(self.mean_species_persistence),\
							median_species_persistence = np.mean(self.median_species_persistence), growth_rate = np.mean(self.growth_rate), pop_size_steady = np.mean(self.pop_size_steady),\
							failed_invasions = np.mean(self.failed_invasions), competitive = np.mean(self.competitive), mut_invasion_condition = np.mean(self.mut_invasion_condition),\
							exp_invasion_condition_positive = np.mean(self.exp_invasion_condition_positive), exp_invasion_condition_negative = np.mean(self.exp_invasion_condition_negative),\
							exp_invasion_condition = np.mean(self.exp_invasion_condition), comm_abun = np.mean(self.comm_abun),\
							host_control = np.mean(self.host_control), mean_species_persistence_steady = np.mean(self.mean_species_persistence_steady))
					
					if(i >= self.steadytime+self.poststeadywindow + 500): # stop the simulation, add 2000 more timesteps
						break
					
				
				
			#paramvals = np.array([self.C, self.Pc, self.Pm, self.halfsat, i])
			#aftersteady = np.array([self.extinctionsteady, self.failedinvasionsteady, self.num_invasions_steady, self.numspeciesatsteady])

		# save the persistence of the species in the final population
		self.ages.extend(len(self.eqs) - self.introtime)
		if self.steady == True:
			self.ages_steady.extend(len(self.eqs) - self.introtime)
		print(self.steady, self.steadytime)
		self.output(self.path_to_save+ '/final_community/output_') # check if this saving for the final population. 
		



	def output(self, location):
		# export the data from the simulation
		paramvals = np.array([self.C, self.Pc, self.Pm, self.halfsat])
		brokenarray = np.array([self.broken])
		steadyarray = np.array([self.steady])
		waitarray = np.array([self.wait])
		aftersteady = np.array([self.extinctionsteady, self.failedinvasionsteady, self.num_invasions_steady, self.numspeciesatsteady])
		np.savez(location +str(self.repeat_num) +'/output_'+str(self.kappa)+'_'+str(self.new_kappa)+"_"+str(self.Pc)+'_'+str(self.Pm)+'_'+str(self.halfsat),\
				 paramvals = paramvals, population = self.population,\
				 failedinvasiontimes=self.failedinvasiontimes, extinctiontimes = self.extinctiontimes,\
				 eqs = self.eqs, numspecies_eq = self.numspecies_eq, broken=brokenarray, \
				 steadytime = self.steadytime, ages=self.ages, ages_steady=self.ages_steady, steady=steadyarray,\
				 aftersteady = aftersteady, waitarray = waitarray, numExtinct = self.numextinct,\
				 growthRate = self.R, invasion_entire = self.num_invasions_entire, eq_criteria = self.eqcriteriamatrix,\
				 mutualistic_invasion_condition = self.mutualistic_invasion_condition,\
				 exploitative_invasion_condition_positive = self.exploitative_invasion_condition_positive,\
				 exploitative_invasion_condition_negative = self.exploitative_invasion_condition_negative,\
				 competitive_invasion_condition = self.competitive_invasion_condition, comm_pop = self.comm_pop, host_control_inv = self.host_control_inv,\
				 extinct_inv = self.extinct_inv, steady_pop = self.popSteady, final_pop = self.population, hostControl = self.hostControl,\
				 A_comp_beforeInv = self.A_comp_beforeInv, A_mut_beforeInv = self.A_mut_beforeInv, A_comp_afterInv = self.A_comp_afterInv,\
				 A_mut_afterInv = self.A_mut_afterInv, interaction_mat = self.A)



	def prev_output(self, location):
		# export the data from the simulation, for a specific equilibrium 
		paramvals = np.array([self.C, self.Pc, self.Pm, self.halfsat])
		brokenarray = np.array([self.broken])
		steadyarray = np.array([self.steady])
		waitarray = np.array([self.wait])
		aftersteady = np.array([self.extinctionsteady, self.failedinvasionsteady, self.num_invasions_steady, self.numspeciesatsteady])
		np.savez(location +str(self.repeat_num) +'/output_'+str(self.kappa)+'_'+str(self.new_kappa)+"_"+str(self.Pc)+'_'+str(self.Pm)+'_'+str(self.halfsat),\
				 paramvals = paramvals, population = self.population,\
				 failedinvasiontimes=self.failedinvasiontimes, extinctiontimes = self.extinctiontimes,\
				 eqs = self.eqs, numspecies_eq = self.numspecies_eq, broken=brokenarray, \
				 steadytime = self.steadytime, ages=self.ages, ages_steady=self.ages_steady, steady=steadyarray,\
				 aftersteady = aftersteady, waitarray = waitarray, numExtinct = self.numextinct,\
				 growthRate = self.R, invasion_entire = self.num_invasions_entire, eq_criteria = self.eqcriteriamatrix,\
				 comm_pop = self.comm_pop, host_control_inv = self.host_control_inv,\
				 steady_pop = self.popSteady, hostControl = self.hostControl,\
				 A_c = self.A_c, A_m = self.A_m, A_e_pos = self.A_e_pos, A_e_neg = self.A_e_neg,\
				 interaction_mat = self.A)


		

	def plot_richness(self):
		"""
		Plot the species richness over time over entire community history.
		"""
		# plot the number of species over time using eqs
		plt.figure()
		fig, ax = plt.subplots(nrows=1, ncols=1)
		eq_num = np.arange(0,len(self.eqs))
		ax.plot(eq_num, self.numspecies_eq, 'r')
		if self.steady == True:
			ax.axvline(x=self.steadytime, color = 'grey', ls='--') # plot vertical line when steady state is reached
		ax.set_ylabel('Species richness')
		ax.set_xlabel('Equilibria')
		ax.set_title(r'C = %.1f, $P_c = $ %.1f, $P_m = $%.1f' %(self.C, self.Pc, self.Pm))
		fig.savefig(self.path_to_save+'/Species_eqNum/output_'+str(self.repeat_num) +'/output_'+str(self.kappa)+'_'+str(self.new_kappa)+"_"+str(self.Pc)+'_'+str(self.Pm)+'_'+str(self.halfsat)+'.eps')

	def plot_dynamics(self):
		"""
		Plot the species population dynamics over the entire community history.
		Community must be initialized with plotbool = True.
		"""
		if self.plotbool == False:
			raise RuntimeError('Must have plotbool=True to plot dynamics.')

		plt.figure()
		fig, ax = plt.subplots(nrows=1, ncols=1)
		# plot extant species
		for i in range(0,self.numspecies):
			ax.plot(self.time, self.populations[i])
		for xc in self.eqs: # plot vertical line whenever new species is introduced
			plt.axvline(x=xc, color = 'grey', ls=':')
		ax.set_ylabel('Population size')
		ax.set_xlabel('Time')
		ax.set_title(r'C = %.1f, $P_c = $ %.1f, $P_m = $%.1f' %(self.C, self.Pc, self.Pm))
		fig.savefig(self.path_to_save+'/Population/output_'+str(self.repeat_num)+'/output_'+str(self.iterationslimit)+'_'+str(self.C)+'_'\
					+str(self.Pc)+'_'+str(self.Pm)+str(self.halfsat)+'.eps')

	def plot_dynamics_single_eq(self):
		"""
		Plot the species population dynamics for a single equilibrium.
		Community must be initialized with plot_failed_eq = True.
		"""
		if self.plot_failed_eq == False:
			raise RuntimeError('Must have plot_failed_eq=True to plot dynamics.')

		plt.figure()
		fig, ax = plt.subplots(nrows=1, ncols=1)

		# plot extant species
		for i in range(0,self.eq_population.shape[0]):
			ax.plot(self.this_eq_time, self.eq_population[i],alpha=0.5)
		ax.set_ylabel('Population size')
		ax.set_xlabel('Time')
		ax.set_title(r'C = %.1f, $P_c = $ %.1f, $P_m = $%.1f, h = %.1f' %(self.C, self.Pc, self.Pm, self.halfsat))
		fig.savefig(self.path_to_save+'/Failed-eq/output_'+str(self.repeat_num)+'/output_'+str(self.failed_eq_counter)+'_'+str(self.iterationslimit)+'_'+str(self.C)+'_'\
					+str(self.Pc)+'_'+str(self.Pm)+'_'+str(self.halfsat)+'.eps')
		
	def num_interactions(self):

		"""
		Determines how much of each interaction exists within the resident microbial 
		community at the end of the community simulation. 
		"""

		# create an adjacency matrix with only ones and zero 
		#make copies of A and pop_size
		num_species = self.population.size
		A_mat = self.A.copy()
		#A_mat_len = len(A_mat[0])
		A_mat[A_mat > 0] = 1 # creating an adjacency matrix 
		A_mat[A_mat < 0] = 1 # creating an adjacency matrix 
		np.fill_diagonal(A_mat,0) # setting all diagonal values to be zero 
		
		A_interaction = A_mat.copy()

		# crearting a mat that stores the interaction type between each species 
		index_m_rows,index_m_cols  = np.where(self.A_m > 0)
		index_c_rows,index_c_cols = np.where(self.A_c < 0)
		index_e_pos_rows,index_e_pos_cols = np.where(self.A_e_pos > 0)
		index_e_neg_rows,index_e_neg_cols = np.where(self.A_e_neg < 0)

		A_interaction[index_m_rows,index_m_cols] = 99    # for M
		A_interaction[index_c_rows,index_c_cols] = -99   # for C
		A_interaction[index_e_pos_rows,index_e_pos_cols] = 45   # for e_pos
		A_interaction[index_e_neg_rows,index_e_neg_cols] = -45  # for e_neg
		
		
	   
		for i in range(0, len(A_interaction[0])):
			for j in range(0,len(A_interaction[0])):
				
			   
				if A_interaction[i,j] == 99.: 
					self.mut_int = self.mut_int + 1 
				elif A_interaction[i,j]== -45. or A_interaction[i,j] == 45.:
					
					self.exp_int = self.exp_int + 1
				elif A_interaction[i,j]== -99.: 
					self.comp_int = self.comp_int + 1
				else:
					"pass"
				
		#paramvals = np.array([self.C, self.Pc, self.Pm, self.halfsat])
		#np.savez(self.path_to_save+'/random_walk/output_'+str(self.repeat_num)+'/output_'+str(self.C)+'_'+str(self.Pc)+'_'+str(self.Pm)+'_'+str(self.halfsat),\
				 #paramvals=paramvals,mut=self.mut_int,comp=self.comp_int,exp=self.exp_int)
		return [self.comp_int,self.mut_int] 
	


	#def sample_interactions(self):
		"""
		Function to perfom a weighted walk and sample the type of microbial interactions among the 
		chosen species and its most weighted neighbor. 
		"""

		# create an adjacency matrix with only ones and zero 
		#make copies of A and pop_size
		#num_species = self.population.size
		#sample_species = sample(list(range(num_species)), num_species)
		#A_mat = self.A.copy()
	   
		#A_mat[A_mat > 0] = 1 # creating an adjacency matrix 
		#A_mat[A_mat < 0] = 1 # creating an adjacency matrix 
		#np.fill_diagonal(A_mat,0) # setting all diagonal values to be zero 
		
		#A_interaction = A_mat.copy()

		# crearting a mat that stores the interaction type between each species 
		#index_m_rows,index_m_cols  = np.where(self.A_m > 0)
		#index_c_rows,index_c_cols = np.where(self.A_c < 0)
		#index_e_pos_rows,index_e_pos_cols = np.where(self.A_e_pos > 0)
		#index_e_neg_rows,index_e_neg_cols = np.where(self.A_e_neg < 0)

		#A_interaction[index_m_rows,index_m_cols] = 99    # for M
		#A_interaction[index_c_rows,index_c_cols] = -99   # for C
		#A_interaction[index_e_pos_rows,index_e_pos_cols] = 45   # for e_pos
		#A_interaction[index_e_neg_rows,index_e_neg_cols] = -45  # for e_neg
		
		#chosen_neib = []
		#chosen_neib_weight=[]
		#chosen_partners = []
		#chosen_partnerW = []
		#interaction_type = []
		
		#for i in zip(sample_species):
			#C_mat = A_mat.copy()
			#neighbors = np.where(C_mat[:, i] == 1)[0].tolist()
			#if(len(neighbors)!=0):
				#chosenNeighbor = i
				#chosen_neib.append(chosenNeighbor)
				#chosenNeighborWeight = np.array([self.population[index] for index in neighbors])
				#chosen_neib_weight.append(np.sum(chosenNeighborWeight))
				#max_weight = max(chosenNeighborWeight)
				#chosen_partnerW.append(max_weight)
				#index_chosen_weight = (np.where(chosenNeighborWeight == max_weight)[0]).tolist()
				#chosen_partner = neighbors[index_chosen_weight[0]]
				#chosen_partners.append(chosen_partner)
				#interaction = A_interaction[chosen_partner,chosenNeighbor]
				#interaction_type.append(interaction)
				#if interaction == 99.: 
					#self.sample_mut_int = self.sample_mut_int + 1 
				#elif interaction == -45. or interaction == 45.:
					#self.sample_exp_int = self.sample_exp_int + 1
				#elif interaction == -99.: 
					#self.sample_comp_int = self.sample_comp_int + 1
				#else:
					#"error, there should be an interaction"
			#else:
				#print("no neighbor dedected for: ", i)
				
		#paramvals = np.array([self.C, self.Pc, self.Pm, self.halfsat])
		#np.savez(self.path_to_save+'/sample_walk/output_'+str(self.repeat_num)+'/output_'+str(self.C)+'_'+str(self.Pc)+'_'+str(self.Pm)+'_'+str(self.halfsat),\
				 #paramvals=paramvals,mut=self.sample_mut_int,comp=self.sample_comp_int,exp=self.sample_exp_int, chosen_neib = chosen_neib,\
				 #chosen_neib_weight = chosen_neib_weight, chosen_partner=chosen_partners,chosen_partnerW=chosen_partnerW,\
				 #interaction_type=interaction_type)
		
	

									
	def commProperty(self):
		"""
		This function determines all the community properties of the microbial community, that we later use for plotting purposes. 
		"""
	
		self.species_list.append(self.numspeciesatsteady)
		self.invasion_steady.append(self.num_invasions_steady/(self.num_invasions_steady+ self.failedinvasionsteady))
		self.inv_entire.append(int(self.num_invasions_entire)/(len(self.failedinvasiontimes)+int(self.num_invasions_entire)))
		self.extinction_steady.append(self.extinctionsteady/(self.extinctionsteady+self.numspeciesatsteady))
		self.mean_species_persistence_steady.append(np.mean(self.ages_steady))
		self.median_species_persistence_steady.append(np.median(self.ages_steady))
		self.extinction_entire.append(self.numextinct/(self.numextinct+ (self.numspecies_eq[-1])))
		self.mean_species_persistence.append(np.mean(self.ages))
		self.median_species_persistence.append(np.median(self.ages))
		self.growth_rate.append(float(self.R[-1]))
		self.pop_size_steady.append(np.sum(self.popSteady)/10000)
		self.failed_invasions.append(self.failedinvasionsteady)
		self.competitive.append(np.mean(self.competitive_invasion_condition))
		self.mut_invasion_condition.append(np.mean(self.mutualistic_invasion_condition))
		self.exp_invasion_condition_positive.append(np.mean(self.exploitative_invasion_condition_positive))
		self.exp_invasion_condition_negative.append(np.mean(self.exploitative_invasion_condition_negative))
		self.exp_invasion_condition.append(self.exp_invasion_condition_positive[-1] + self.exp_invasion_condition_negative[-1])
		self.comm_abun.append(np.sum(self.population)/10000)
		self.host_control.append(self.hostControl)

		   
	
def run_simulation(Pc,Pm,kappa,immune_microbial_load,immuneResponse,change_kappa,new_kappa,repeat_num,path_to_save):
	"""
	Create combinations of parameters and simulate communities.
	Takes an index from 0 to 593 as an argument. This index determines the
	parameter values that are used to create the community.
	"""
	# create lists all possible parameter values for Pc, Pm, C, h
	h = 100
	C = 0.5
	

	simulation = Community(C, Pc, Pm, h, kappa,immune_microbial_load, immuneResponse,change_kappa,new_kappa,repeat_num,path_to_save,study_selection=False, IIM_bool = False,
						   VUM_bool = False, lamb=0.5)
	simulation.start()
	simulation.iterate()
	simulation.plot_richness()
	#simulation.plot_dynamics()
	#simulation.num_interactions()
	#simulation.sample_interactions()

	return
