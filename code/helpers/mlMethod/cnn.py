import numpy as np

samplesNumber = 20000
windowSize = 14*14
trainPer = 90

X, y = None, None

with open("../mapSequence/contigData.txt", "r") as f:
    lines = f.readlines()
    X = [[float(j) for j in lines[i].strip("[]/\n").split(", ")] for i in range(1, len(lines), 2)]
    y = [lines[i].strip() for i in range(0, len(lines), 2)]

X = np.array([np.array(i) for i in X])

classes = set()
for i in y:
    if i not in classes:
        classes.add(i)
classes = list(classes)

classToInt, intToClass = {}, {}
for i in range(len(classes)):
    classToInt[classes[i]] = i
    intToClass[i] = classes[i]

y = np.array([classToInt[i] for i in y], dtype=np.int)

newX, newY = [], []

for i in range(len(X)):
    for j in range(len(X[i])//windowSize):
        newX.append(X[i][j*windowSize:j*windowSize+windowSize])
        newY.append(y[i])

X = np.array(newX)
y = np.array(newY)

X = X[:samplesNumber]
y = y[:samplesNumber]

######################
from sklearn.utils import shuffle
from keras.utils import to_categorical

y = np.array([0 if intToClass[i] == "mtDNA" else 1 for i in y])

print("SHape is")
print(X.shape)
print(y.shape)

posIndex = (y==0)
negIndex = (y==1)

X_pos = X[posIndex]
X_neg = X[negIndex]

y_pos = y[posIndex]
y_neg = y[negIndex]
X = np.append(X_pos[:900], X_neg[:900]).reshape(1800, windowSize)
y = np.append(y_pos[:900], y_neg[:900])

print("SHape is" + str(y_pos.shape) + " " + str(y_neg.shape))
print(X.shape)
print(y.shape)

X, y = shuffle(X, y, random_state=0)

num_of_classes = len(set(y))

y = to_categorical(y)

X = (X.T[:] - np.mean(X, axis=1, dtype='float64')).T
X[:] -= X.mean(axis=0, dtype = 'float64')
X[:] /= X.std(axis=0, dtype = 'float64')

X = X.reshape((X.shape[0], 14, 14, 1))

trainPart = (trainPer*X.shape[0])//100

X_train = X[:trainPart]
X_valid = X[trainPart:]

y_train = y[:trainPart]
y_valid = y[trainPart:]

######################

from keras.models import Sequential
from keras.layers import Dense, Conv2D, Flatten, Dropout, MaxPooling2D
from keras.layers.normalization import BatchNormalization
from keras.layers.core import Activation
from keras.optimizers import Adam

import matplotlib.pyplot as plt

model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3), input_shape=X.shape[1:], padding = "valid"))
model.add(Activation("relu"))
model.add(BatchNormalization())

model.add(Dropout(0.2))

model.add(Conv2D(32, (3, 3), padding="same"))
model.add(Activation("relu"))
model.add(BatchNormalization())

model.add(MaxPooling2D(pool_size=(2, 2), padding = "valid"))

model.add(Conv2D(64, (3, 3), padding="same"))
model.add(Activation("relu"))
model.add(BatchNormalization())

model.add(Dropout(0.4))

model.add(Flatten())

model.add(Dense(64))
model.add(Activation("relu"))
model.add(BatchNormalization())

model.add(Dense(num_of_classes, activation='softmax'))

model.compile(optimizer='adam', loss='categorical_crossentropy',  metrics=['accuracy'])

model.summary()

history = model.fit(X_train, y_train, validation_data=(X_valid, y_valid), batch_size=128, epochs=20)

plt.figure()
plt.plot(history.history['accuracy'], label='train accuracy')
plt.plot(history.history['val_accuracy'], label='validation accuracy')
plt.legend(loc='best')
plt.show()

