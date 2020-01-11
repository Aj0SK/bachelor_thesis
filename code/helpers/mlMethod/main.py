import numpy as np

windowSize = 250
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

print(len(X))
X = X[:10000]
y = y[:10000]

from sklearn import svm
from sklearn.model_selection import cross_val_score

svc = svm.SVC(C=1, kernel='rbf')
scores = cross_val_score(svc, X, y, cv=5)
print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
