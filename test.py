from funcAndClas import *

O2_test = Compound("O2", "g")
print(O2_test.s_formula)


print(O2_test.dfH0_298)


print(O2_test.cp_powers)


print(O2_test.cp0_298)


print(O2_test)

#%%
r_test = Reaction(["C4H10", "O2", "CO2", "H2O"], ["g", "g", "g", "g"], [-1, -13/2 , 4, 5])

print(r_test.drH0_298)
print(r_test)

print(r_test.calc_gibbs())

xt = np.arange(0, 10, 1e-5)
yt = xt**2
yt_int = integrate_trapezoid(xt, yt)
print(yt_int)

print(r_test.calc_heat_effect(400))
print(r_test.calc_react_entropy(400))
print(r_test.calc_gibbs(400))

print(r_test)

r_test2 = Reaction(["SO2", "O2", "SO3"],["g","g","g"],[-1,-1/2,1])

print(r_test2)

print(r_test2.calc_gibbs(298.15))

T = 900
drG0_T = r_test2.calc_gibbs(T)
print(drG0_T)

#%%
xi = np.arange(0.01, 0.99, 0.01)
drG = [calculate_drG(xi_i, r_test2, drG0_T, T, 1, [1, 0.5, 0]) for xi_i in xi]
#%%
#Построим график
import matplotlib.pyplot as plt
#%%
plt.plot(xi, drG)
plt.xlim(0, 1)
plt.plot([0,1], [0,0], color="black", linestyle="--")
