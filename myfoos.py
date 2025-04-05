def composition_to_mass(composition):
    weights = load_weights("elements.txt")
    mass = 0.0
    for element in composition:
        index = composition[element]
        Al = weights[element]
        mass += index * Al
    return mass

def decompose_formula(s_formula):
    composition = dict()
    s_buffer = ""
    s_formula += "\r"
    for ch in s_formula:
        if (ch.isupper() or ch == "\r") and s_buffer != "":
            s_element = ""
            s_index = ""
    
            for ch_2 in s_buffer:
                if ch_2.isdigit():
                    s_index += ch_2
                else:
                    s_element += ch_2
    
            if s_index == "":
                s_index = "1"
    
            if not s_element in composition:
                composition[s_element] = int(s_index)
            else:
                composition[s_element] += int(s_index)
        
            s_buffer = ""
        
        if ch == "\r":
            break
    
        s_buffer += ch

    return composition

def load_weights(s_filename):
    weights = dict()
    f = open(s_filename)
    for line in f:
        tmp = line.split("\t")
        weights[tmp[1]] = float(tmp[3])
    return weights

def calc_M(s_formula):
    return composition_to_mass(decompose_formula(s_formula))

def integrate_trapezoid(x, y):
    sum = 0
    for i in range(len(x)-1):
        sum += (y[i+1] + y[i])/2 * (x[i+1] - x[i])
    return sum