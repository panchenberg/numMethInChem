# %%
from funcAndClas import *


# %%
def compare_compounds(c1, c2):
    if c1.s_state != c2.s_state:
        return False
    if set(c1.composition.keys()) != set(c2.composition.keys()):
        return False
    for element in c1.composition:
        if c1.composition[element] != c2.composition[element]:
            return False
    return True


# %%
class PT_equilibrium:
    def __init__(self, reactions, n0_s, P, T):
        self.reactions = reactions
        self.n0_s = n0_s
        self.P = P
        self.T = T

        self.comm_parts = []
        for i in range(len(self.reactions)):
            react_vec = []
            for j in range(len(self.reactions[i].parts)):
                comp_vec = []
                for k in range(len(self.reactions)):
                    for l in range(len(self.reactions[k].parts)):
                        if compare_compounds(self.reactions[i].parts[j], self.reactions[k].parts[l]):
                            comp_vec.append(k)

                react_vec.append(comp_vec)
            self.comm_parts.append(react_vec)

        print(self.comm_parts)


# %%
text1 = "1*SO2(g) + 1*O2(g) = 1*SO3(g)"
text2 = "1*SO3(g) + 1*N2O(g) = 1*SO2(g) + 2*NO(g)"
eq_1 = PT_equilibrium([getReaction(text1),
                       getReaction(text2)],
                      [[1, 1, 1], [1, 10, 1, 2]],
                      1,
                      900)

# %%
