import numpy as np
from scipy.optimize import fsolve
import time
import matplotlib.pyplot as plt
import cantera as ct
print(f"Running Cantera version: {ct.__version__}")

input_file = "sofc.yaml"


F = 96485e3 #C/Kmol

# Temperature and pressure
T = 300  # K
P = ct.one_atm

gas = ct.Solution(input_file, "gas")
metal = ct.Solution(input_file, "metal")                    # Ni bulk (anode)
oxide = ct.Solution(input_file, "oxide_bulk")               # YSZ bulk (cathode)
metal_surface = ct.Solution(input_file, "metal_surface")    # Ni surface
oxide_surface = ct.Solution(input_file, "oxide_surface")    # YSZ surface
tpb = ct.Solution(input_file, "tbp")                        # triple phase boundary

phases = [gas, metal, oxide, metal_surface, oxide_surface, tpb]
for ph in phases:
    ph.TP = T,P

# --- Constants ---
R = ct.gas_constant  # J/kmol-K
n = 2                # electrons transferred
beta = 0.5           # from yaml

# --- Get species indices for TPB reaction ---
tpb_species = tpb.species_names

i_Hm = tpb_species.index('H(m)')
i_Oox = tpb_species.index("O''(ox)")
i_OHm = tpb_species.index('OH(m)')
i_m = tpb_species.index('(m)')

# --- Get rate constants ---
k_f = tpb.forward_rate_constants[0]
k_r = tpb.reverse_rate_constants[0]

# --- Get concentrations ---
C = tpb.concentrations

C_Hm = C[i_Hm]
C_Oox = C[i_Oox]
C_OHm = C[i_OHm]

# --- Compute equilibrium potential difference ---
# Using mass-action equilibrium condition
delta_phi_eq = (R * T) / (n * F) * np.log((k_f * C_Hm * C_Oox) / (k_r * C_OHm))

print(f"Delta phi_eq = {delta_phi_eq:.4e} V")

# --- Sweep overpotential ---
eta = np.linspace(0, 0.3, 100)
current = []

for e in eta:
    delta_phi = delta_phi_eq + e

    # Set potentials (only difference matters)
    metal.electric_potential = delta_phi
    oxide.electric_potential = 0.0

    # Update rates
    rop_f = tpb.forward_rates_of_progress[0]
    rop_r = tpb.reverse_rates_of_progress[0]

    rop_net = rop_f - rop_r

    # Current density
    i_val = n * F * rop_net
    current.append(i_val)

current = np.array(current)

# --- Plot ---
plt.figure()
plt.plot(eta, current, label='Cantera')
plt.xlabel('Overpotential η (V)')
plt.ylabel('Current density (A/m^2)')
plt.legend()
plt.grid()
plt.show()