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

X = X[:10000]
y = y[:10000]

y = [0 if intToClass[i] == "mtDNA" else 1 for i in y]

from sklearn import svm
from sklearn.model_selection import cross_val_score
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.utils import shuffle

X, y = shuffle(X, y, random_state=0)

print(len(classes))
print(classes)

X = (X.T[:] - np.mean(X, axis=1, dtype='float64')).T
X[:] -= X.mean(axis=0, dtype = 'float64')
X[:] /= X.std(axis=0, dtype = 'float64')

svc = svm.SVC(C=1.0, kernel='rbf', gamma='auto')
scores = cross_val_score(svc, X, y, cv=5)
print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

plt.plot(X[0])
plt.show()

'''
pca = PCA(n_components = 3)
prinCom = pca.fit_transform(X)

x_axis = prinCom[:, 1]
y_axis = prinCom[:, 2]

t = np.zeros((X.shape[0], 3))
t[y[:]==0] = (1, 0, 0)
t[y[:]==1] = (0, 1, 0)

plt.scatter(x_axis, y_axis, c=t, marker = 'o')
plt.xlabel("second_component")
plt.ylabel("third_component")
plt.show()
'''
