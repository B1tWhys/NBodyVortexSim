import matplotlib.pyplot as plt
import re

f = open("timeData", 'r')

totT = 0
step = 0

reg = re.compile("""^[^\d]+(\d+)[^\d]+([\d.]+).*""")

data = {int(x):float(y) for x, y in [reg.search(line).groups([1, 2]) for line in f]}
print("Total compute time: %fs"%round(sum(data.values()), 2))


plt.plot(data.keys(), data.values())
plt.show()