"""
Best fit to a target density given a set of eigenstates |Psi_n>
The output are the coefficients c_n of the state |Psi_b> = sum c_n|Psi_n> 
that better approximate the target density rho_t
"""

import random
import numpy as np
import math

class ga_fit:
   def __init__(self,bset,eigenstates, tdens, npop=100, ngen=10, 
                mort=0.01, nmin=10, selection_pow=3, pmutation = 0.05):
      """
      INPUT:
         bset = basis set object
	 eigenstates = numpy matrix where Psi_n = eigenstates[:,n]
	 tdens = target density. List of tuples consisting of (x,y,dens).
	         This density is intended normalized to 1.
	 
      """
      
      self.bset = bset
      self.eigs = eigenstates
      self.targ = tdens
      self.sel = selection_pow
      self.mort=mort

      self.ngenes = eigenstates.shape[0]
      self.npop=npop
      self.ngen=ngen
      self.nmin=nmin
      self.pmutation = pmutation

      self.gen_init_pop()
      
      print "Initial population loaded. Best fitness is: ",self.fitness[0]," Now running..."
      
      for i in xrange(self.ngen):
         cbest = 1*self.pop[0]
	 self.crosstrm()
	 self.mutation()
	 self.npop = self.npop - int((self.npop-self.nmin)*self.mort)   
         self.npop = self.npop - self.npop%2
	 if self.npop<self.nmin:
	    self.npop=self.nmin
	    
	 self.pop = self.pop[:self.npop]
	 self.pop[0] = cbest
	 self.ffitness()
         self.sortit()
         self.best = self.pop[0]
         print "Gen ", i,". Current population: ",self.npop," Best Fitness: ",self.fitness[0]
	 
   
      psi = np.zeros(self.ngenes)
      for i,cn in enumerate(self.best):
         psi = psi + float(cn)*self.eigs[:,i]   

      self.beststate = psi.tolist()	 

   def gen_init_pop(self):
      """
      Initialize random pops
      """      
      self.pop=[]
      
      for i in xrange(self.npop):
         
	 cpop=[]
	 nrm=0
	 for j in xrange(self.ngenes):
            cvl = 2*random.random() - 1.0
	    cpop.append(cvl)
	    nrm = nrm + cvl*cvl
	    
	 cpop = map(lambda x : float(x)/nrm, cpop)
	 self.pop.append(cpop)
	 
      self.ffitness()
      self.sortit()
         	 


   def crosstrm(self):
      """
      Create the new generation of pops
      """
      
      npops=[]
      for i in xrange(self.npop/2):
         ii = int( (random.random()**self.sel) * self.npop )
         jj = int( (random.random()**self.sel) * self.npop )
         while ii == jj:
	    jj = int( (random.random()**self.sel) * self.npop )
         
         pop1=[]
	 pop2=[]
	 npop1=0
	 npop2=0
	 for cni, cnj in zip(self.pop[ii], self.pop[jj]):
	    if random.random() < 0.5:
	       c1,c2 = cni,cnj
	    else:
	       c1,c2 = cnj,cni
         
	    pop1.append(c1)
	    npop1 = npop1 + c1*c1
	    pop2.append(c2)
	    npop2 = npop2 + c2*c2
	 
	 npops.append(map(lambda x : float(x)/npop1, pop1))
	 npops.append(map(lambda x : float(x)/npop2, pop2))
	 
      self.pop = npops
      
   
   def mutation(self):
      """
      Make a mutation over a pair of chars so that the normalization is preserved
      """  
      for i in xrange(self.npop):
         if random.random() > self.pmutation:
            continue
	    
	 vm = 2*random.random() - 1.0
	 ichr1 = int(random.random()*self.ngenes)
	 ichr2 = int(random.random()*self.ngenes)
	 while ichr1 == ichr2:
	    ichr2 = int(random.random()*self.ngenes)
	 
	 sq = self.pop[i][ichr1]*self.pop[i][ichr1] + self.pop[i][ichr2]*self.pop[i][ichr2]
	 self.pop[i][ichr1]= vm
	 self.pop[i][ichr2]= sq - vm*vm
	 

   def ffitness(self):
      """
      Compute fitness of pops
      """
      
      self.fitness=[]
      for cp in self.pop:
         psi = np.zeros(self.ngenes)
	 for i,ng in enumerate(cp):
	    psi = psi + float(ng)*self.eigs[:,i]
	    
	 ftn=0
	 for entry in self.targ:
	    x=entry[0]
	    y=entry[1]
	    val=entry[2]
	    model = self.bset.density((x,y),psi.tolist())
	    ftn = ftn + (val-model)**2
	    
	 ftn = 1.0 - math.exp(-ftn/len(self.targ))
	 self.fitness.append(ftn)  	 


   def sortit(self):
      """
      Sort the pop array according to the fitness value
      """
      
      ev=[]
      for cp,cf in zip(self.pop, self.fitness):
         ev.append((cp,cf))
	 
      ev.sort(key = lambda e : e[1], reverse=False)
      
      self.pop=[]
      self.fitness=[]
      for el in ev:
         self.pop.append(el[0])
	 self.fitness.append(el[1])
	 
      
