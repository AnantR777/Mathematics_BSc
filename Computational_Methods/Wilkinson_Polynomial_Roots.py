import numpy as np
import matplotlib.pyplot as plt

#  coefficients of perturbed polynomial stored in list
epsilon = 1e-10
P_perturbed_coeffs = [
1 + epsilon,
-210,
20615,
-1256850,
53327946,
-1672280820,
40171771630,
-756111184500,
11310276995381,
-135585182899530,
1307535010540395,
-10142299865511450,
63030812099294896,
-311333643161390640,
1206647803780373360,
-3599979517947607200,
8037811822645051776,
-12870931245150988800,
13803759753640704000,
-8752948036761600000,
2432902008176640000,
]

#array containing zeros of P, given in question that roots are from 1 to 20
P_zeros = np.arange(1, 21)

#array containing zeros of P_perturbed
P_perturbed_zeros = np.roots(P_perturbed_coeffs)

# extract real part of roots of P
p_real = [proot.real for proot in P_zeros]

# extract imaginary part of roots of P
p_imag = [proot.imag for proot in P_zeros]

# extract real part for roots of perturbed P
pert_real = [pertroot.real for pertroot in P_perturbed_zeros]

# extract imaginary part for roots of perturbed P
pert_imag = [pertroot.imag for pertroot in P_perturbed_zeros]


plt.scatter(p_real, p_imag, c = 'red', marker= "o")
plt.scatter(pert_real, pert_imag, c = 'black', marker= "x")
plt.xlabel('Real')
plt.ylabel('Imaginary')
plt.show()
