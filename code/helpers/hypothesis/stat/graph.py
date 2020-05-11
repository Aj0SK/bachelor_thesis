import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# sphinx_gallery_thumbnail_number = 2

klens = []
levels = []

d = {}

with open("out2.txt", "r") as f:
    for line in f.readlines():
        line = line.split()
        k = int(line[2])
        level = int(line[4])
        hits = int(line[6])
        
        klens.append(k)
        levels.append(level)
        
        if line[0] == "+":
            d[(k, level)] = d.get((k, level), 1) * hits
        if line[0] == "-":
            d[(k, level)] = d.get((k, level), 1) / hits

data = []

klens = list(set(klens))
levels = list(set(levels))

for i in range(1, 100):
    x = []
    for j in range(1, 100):
        if (i, j) in d:
            x.append(d[(i, j)])
    if len(x) != 0:
        data.append(x)

data = np.array(data)
data = np.log(data)

fig, ax = plt.subplots()
im = ax.imshow(np.transpose(data), cmap = "hot", interpolation = "nearest")
fig.colorbar(im, ax=ax, orientation='vertical')

ax.set_xticks(np.arange(len(klens)))
ax.set_yticks(np.arange(len(levels)))

ax.set_xticklabels(klens)
ax.set_yticklabels(levels)

plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
#for i in range(len(klens)):
#    for j in range(len(levels)):
#        text = ax.text(i, j, harvest[i, j],
#                       ha="center", va="center", color="w")

ax.set_title("Log of positive read hits / negative read hits (> 0 is better)")
fig.tight_layout()
plt.show()
