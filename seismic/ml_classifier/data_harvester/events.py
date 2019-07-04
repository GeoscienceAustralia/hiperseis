from mat4py import *

#load GA's pick database

GA=loadmat("GA.mat")['GA']

print GA['picks'].keys()
