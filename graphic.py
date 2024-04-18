import matplotlib.pyplot as plt

data_x = []
data_y = []
with open("./out.txt") as file:
    for i in file.readlines():
        data_x.append(float(i.split()[0]))
        data_y.append(float(i.split()[1]))

plt.plot(data_x, data_y, "ro")
plt.grid(True)
plt.show()
