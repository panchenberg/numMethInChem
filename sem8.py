# %%
import numpy as np
import matplotlib.pyplot as plt


#здесь мы смотрели на то как графики соотносятся друг с другом
# %%
f_data = open("NO2_Cp.txt")
T = []
Cp_exp = []
for line in f_data:
    t = line.split() # " "
    T.append(float(t[0]))
    Cp_exp.append(float(t[1]))
    

# %%
plt.scatter(T, Cp_exp)

# %%
factors = [lambda T: 1,
          lambda T: T,
          lambda T: (T**2)]


# %%
n = len(Cp_exp)
k = len(factors)
F = []
y = []
for i in range(k):
    row = []
    for j in range(k):
        sum = 0
        for l in range(n):
            sum += factors[i](T[l]) * factors[j](T[l]) 
        row.append(sum)
    F.append(row)

    sum = 0
    for l in range(n):
        sum += factors[i](T[l]) * Cp_exp[l]
    y.append(sum)

# %%
F = np.matrix(F)
y = np.array(y)
F_inv = np.linalg.inv(F)
b = np.asarray(np.dot(F_inv, y))[0]
print(b)

# %%
T_fine = np.arange(300, 1000, 1)
Cp_calc = []
for Ti in T_fine:
    Cp = 0
    for i in range(k):
        Cp += b[i] * factors[i](Ti)
    Cp_calc.append(Cp)

# %%
plt.scatter(T, Cp_exp)
plt.plot(T_fine, Cp_calc, color="red")

# %%
SS_res = 0
SS_tot = 0
y_mean = np.mean(Cp_exp)
for i in range(n):
    Cp = 0 
    for j in range(k):
        Cp += b[j] * factors[j](T[i])
    SS_res += (Cp - Cp_exp[i])**2
    SS_tot += (Cp - y_mean)**2
R2 = 1 - SS_res/SS_tot
print(f"R2 = {R2}")
print(R2)

# %%
