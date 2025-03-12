import matplotlib.pyplot as plt

f, ax = plt.subplots()
y1 = 4.178996099
y2 = 1
y3 = 2.132797971
L = 1

ax.plot([0,0],[-y1/2, y1/2], color='black', linestyle='--')
ax.plot([L,L],[-y2/2, y2/2], color='black', linestyle='--')
ax.plot([2*L,2*L],[-y3/2, y3/2], color='black', linestyle='--')
plt.show()