#Argument 1: training data file (separated by tab) [checked.bed]
#Argument 2: file containing parameter values for hyperbola curve (separated by spaces), if any, or file where parameter values have to be saved [best1.txt]

from random import choice
from numpy import array, dot
import time
import sys, os, subprocess, linecache
import numpy.random as rand
import re
import random

unit_step = lambda x: 0 if x < 0 else 1

training = open(sys.argv[1], 'r')
training_data=[]

#getting training data
for line in training.readlines():
  line_list=line.strip().split('\t')
  print(line_list)
  coord=array([line_list[3],line_list[4]])
  tupl=(coord,line_list[2])
  training_data.extend([tupl])
  
w = [0.0, 0.0] 
b = 80000.0
errors = []
eta = 100 #training speed
n = 500 #number of times random training values are considered
rp=100 #number of times the whole code repeats

if(os.path.isfile(sys.argv[2])):
  #taking best parameters if saved earlier
  best=open('best1.txt', 'r+')
  l=linecache.getline('best1.txt', 1).strip().split()
  minct=int(float(l[0]))
  bw=[float(l[1]),float(l[2])]
  w=bw
  bb=float(l[3])
  b=bb
else:
  minct=10000
  bw=[0.0,0.0]
  bb=100000.0

print('The starting parameters are:')
print(minct)  
print(bw)
print(bb)

for num in range(rp):
  for i in xrange(n):
    x, expected = choice(training_data)
    #print(x)
    #print(expected)
    result = float((float(x[0])-float(w[0]))*(float(x[1])-float(w[1]))-float(b))
    r=-1
    if(result>0):
      r=1
    error = float(expected) - r
    #print(error)
    if(error!=0):
      #print('here')
      b=b-float(eta*error)
      #print(b)
  
  count=0
  for (x, exp) in training_data:
    #print(x,exp)
    if(exp==0):
      exp=-1
    result = float((float(x[0])-float(w[0]))*(float(x[1])-float(w[1]))-float(b))
    r=-1
    if(result>0):
      r=1
    if (r!=int(exp)):
      count=count+1
    print("{}: {} -> {}".format(x[:2], exp, r))
    #prints expected and actual value of all training points
  
  if(count<=minct):
    bw=w
    bb=b
    minct=count

training.close()
#best.truncate(0)
if(os.path.isfile(sys.argv[2])):
  best.close()
best=open(sys.argv[2], 'w')
best.close()
best=open(sys.argv[2], 'w')
best.write(repr(minct)+'\t'+repr(bw[0])+'\t'+repr(bw[1])+'\t'+repr(bb))
best.write('\n')
best.close()

print(minct)  
print(bb)
