from model import *
from data import *
from time import time
from tensorflow.python.keras.callbacks import TensorBoard


#os.environ["CUDA_VISIBLE_DEVICES"] = "0"



myGene = trainGenerator(128)
valGene = valGenerator(128)
testGene = testGenerator(128)

model = shakenet()
model_checkpoint = ModelCheckpoint('shakenet-model.hdf5', monitor='val_acc',verbose=1, save_best_only=True, mode='max')
tensorboard = TensorBoard(log_dir="logs/{}".format(time()))

#model.fit_generator(myGene,steps_per_epoch=len(myGene),epochs=60,callbacks=[model_checkpoint],validation_data=valGene,nb_val_samples=len(valGene))
model.fit_generator(myGene,steps_per_epoch=len(myGene),epochs=30,callbacks=[model_checkpoint,tensorboard],validation_data=valGene,nb_val_samples=len(valGene))



results = model.evaluate_generator(testGene,steps=len(testGene))

print("Test set loss: "+str(results[0]))
print("Test set accuracy: "+str(results[1]))
