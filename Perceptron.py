from random import choice
from numpy import array, dot
import time
import sys, os, subprocess, linecache
import numpy.random as rand
import re
import random

unit_step = lambda x: 0 if x < 0 else 1

training = open(sys.argv[1], 'r')
best=open('best.txt', 'r+')

training_data=[]

for line in training.readlines():
  line_list=line.strip().split('\t')
  coord=array([line_list[1],line_list[2]])
  tupl=(coord,line_list[3])
  training_data.extend([tupl])

#print(training_data)

#training_data = [ #insert from file
#  (array([0,0,1]), 0), 
#  (array([0,1,1]), 1), 
#  (array([1,0,1]), 1), 
#  (array([1,1,1]), 1), 
#]

w = [0.8, 0.7] 
b = 0.06
errors = []
eta = 0.01
n = 10000
eta2= 0.02
rp=10

print('49 -11751844.038644727 2491654.468133597 -2.1359449978471168'.split())

l=linecache.getline('best.txt', 1).strip().split()
#print(l_list)
#l_list=re.split(r'\s+', l)

#print(l[0])
#print(l[1],l[2])
#print(l[3])

minct=int(float(l[0]))
bw=[float(l[1]),float(l[2])]
bb=float(l[3])

print(minct)  
print(bw)
print(bb)

for num in range(rp):
  for i in xrange(n):
    x, expected = choice(training_data)
    #print(x)
    #print(expected)
    result =  float(w[0])*float(x[1]) + float(w[1])*float(x[0]) + float(b)
    #print(float(b))
    #r=-1
    #if(result>0):
    #  r=1
    error = float(expected) - result
    #print(error)
    #errors.append(error)
    if(error!=0):
      #print('here')
      #print(float(eta) * float(error)*float(x[0]))
      w[0] = float(float(w[0]) + float(float(eta) * float(error) * float(x[1])))
      #print(float(w[0]))
      w[1] = float(w[1]+ float(float(eta) * float(error) * float(x[0])))
      #print(w[1])
      b=b+float(eta2*error)
      #print(b)
  
  count=0
  for (x, exp) in training_data:
    #print(x,exp)
    result = float(w[0])*float(x[1]) + float(w[1])*float(x[0]) + float(b)
    r=-1
    if(result>0):
      r=1
    if (r!=int(exp)):
      count=count+1
    #print("{}: {} -> {}".format(x[:2], exp, r))
  
  print(num)
  print(count)  
  print(w)
  print(b)
  
  if(count<minct):
    bw=w
    bb=b
    minct=count

training.close()
#best.truncate(0)
best.close()
best=open('best.txt', 'w')
best.close()
best=open('best.txt', 'w')
best.write(repr(minct)+'\t'+repr(bw[0])+'\t'+repr(bw[1])+'\t'+repr(bb))
best.close()

#print(minct)  
#print(bw)
#print(bb)

#print('over')

#process=subprocess.Popen(['Rscript','Plots.R',repr(float(bb/(-bw[1]))),repr(float(bw[0]/(-bw[1])))])
#process=subprocess.Popen(['Rscript','Plots.R','best.txt']
#process.wait()
