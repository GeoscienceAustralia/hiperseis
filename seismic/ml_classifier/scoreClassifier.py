from model import *
from data import *

sGen=SGenerator(16)
nGen=NGenerator(16)

model=shakenet(pretrained_weights='shakenet-model.hdf5')

sRes= model.predict_generator(sGen,steps=len(sGen))
nRes= model.predict_generator(nGen,steps=len(nGen))

tS=0
fS=0
tN=0
fN=0

for result in sRes:
    if np.argmax(result)==0:
        tS+=1
    else:
        fN+=1

for result in nRes:
    if np.argmax(result)==0:
        fS+=1
    else:
        tN+=1

print("S recall = "+str(float(tS)/len(sRes))+", S precision = "+str(float(tS)/(tS+fS)))
print("N recall = "+str(float(tN)/len(nRes))+", N precision = "+str(float(tN)/(tN+fN)))

tGen=testGenerator(16)

results = model.evaluate_generator(tGen,steps=len(tGen))

print("Test set loss: "+str(results[0]))
print("Test set accuracy: "+str(results[1]))
