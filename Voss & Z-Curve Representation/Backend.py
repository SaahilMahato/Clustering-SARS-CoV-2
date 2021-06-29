import numpy as np 
from scipy.fft import fft

def read_dna_seq(file_name):
  fil = open(file_name,'r')
  fil_list = fil.readlines()
  fil.close

  genome = {}
  acession = ''
  protien_name = ''
  gene_seq = ''
  for i in fil_list:
    if i[0] == '>':
      if list(genome.keys()) != []:
        gene_seq = gene_seq.replace('\n','')
        genome[acession].append(gene_seq)
      gene_seq = ''
      acession_st = i.find('M')
      acession_end = i.find('|')
      date_st = i.find('2020-') 
      date_end = i.find('T00')

      if acession_st > 0 and acession_end > 0:
        acession = i[acession_st:acession_end]
        date = i[date_st:date_end]
        genome[acession] = []
        genome[acession].append(date)
    else:
      gene_seq += i
  gene_seq = gene_seq.replace('\n','')
  genome[acession].append(gene_seq)    
  return genome


class dna:
  def __init__(self, dna_seq):
    self.date = dna_seq[0]
    self.dna_seq = dna_seq[1].upper()
  
  def _vectorization(self):
    v1 = []
    v2 = []
    v3 = []
    v4 = []
    for G in self.dna_seq:
      if G == 'A':
        v1.append(1)
        v2.append(0)
        v3.append(0)
        v4.append(0)
      if G == 'C':
        v1.append(0)
        v2.append(1)
        v3.append(0)
        v4.append(0)
      if G == 'G':
        v1.append(0)
        v2.append(0)
        v3.append(1)
        v4.append(0)
      if G == 'T':
        v1.append(0)
        v2.append(0)
        v3.append(0)
        v4.append(1)
    return [v1,v2,v3,v4]
  
  def _ZCurve(self):
    x = [0]
    y = [0]
    z = [0]
    for i,G in enumerate(self.dna_seq):
      if G == 'A':
        x.append(x[i]+1)
        y.append(y[i]+1)
        z.append(z[i]+1)
      elif G == 'G':
        x.append(x[i]+1)
        y.append(y[i]-1)
        z.append(z[i]-1)
      elif G == 'C':
        x.append(x[i]-1)
        y.append(y[i]+1)
        z.append(z[i]-1)
      elif G == 'T':
        x.append(x[i]-1)
        y.append(y[i]-1)
        z.append(z[i]+1)
      else:
        x.append(x[i]-1)
        y.append(y[i]-1)
        z.append(z[i]-1)   
    return [x[1:],y[1:],z[1:]]

  def _DFT(self,V):
    Trans_V = []
    for v in V:
      Trans = fft(v)
      Trans_V.append(Trans)
    return Trans_V

  def PowerSpectrum(self, V):
    Trans_V = self._DFT(V)
    Ps = []
    for i in range(len(Trans_V[0])):
      value = 0
      for j in range(len(Trans_V)):
        value += abs(Trans_V[j][i])**2
      Ps.append(value)
    return Ps

  def PAPR(self, V):
    Ps = self.PowerSpectrum(V)
    peak_2_average_pwr_ratio = max(Ps)/(sum(Ps)/len(Ps)) 
    return peak_2_average_pwr_ratio
  
  def extract(self):
    V_vector = self._vectorization()
    V_zcurve = self._ZCurve()
    Feature_V = self.PAPR(V_vector)
    Feature_Z = self.PAPR(V_zcurve)
    return Feature_V, Feature_Z
