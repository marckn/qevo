"""
User-defined units.
From IS to userdef units [X]:

1m = cl [L]  (i.e.: for nm cl = 1E9)
1s = ct [t]
1Kg = cm [M]
"""

import math

class units:
   si_hbar = 1.0545718E-34
   si_enunit = 1.0
   si_kb = 1.38064852E-23
   si_e = 1.60217662E-19
   si_eps0 = 8.854187817E-12
   si_coulfac = si_e**2/(4*math.pi*si_eps0)

   def __init__(self,pcl,pct,pcm):
      self.cl=pcl
      self.ct=pct
      self.cm=pcm
      self.hbar = self.si_hbar*pcm*pcl*pcl/pct
      self.enunit = self.si_enunit*pcm*pcl*pcl/(pct*pct)
      self.kb = self.si_kb*self.enunit
      self.coulfac = self.si_coulfac*pcm*pcl**3/pct**2
      self.toev = 1.0/(self.si_e*self.enunit)
      self.tomev = 1000.0/(self.si_e*self.enunit)
      print(
      """
      #######################
      #######################
      ######## UNITS ########
      #######################
      #######################

      # [L] -> {0} m
      # [T] -> K
      # [t] -> {1} s
      # [M] -> {2} Kg

      # So that:
      # [E] = {3} J = {4} eV
      # hbar = {5} [E][t]
      # Kb = {6} [E]/K
      # e^2/(4pi*eps_o[L]) = {7}
      #
      ########################
      ########################
      """.format(1.0/self.cl,1.0/self.ct,1.0/self.cm,1.0/self.enunit,self.toev,self.hbar,self.kb,self.coulfac))


