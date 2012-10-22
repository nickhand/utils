# assumes everything in cgs
from numpy import pi

c_light = 2.9979e10
h_planck = 6.6261e-27
h_bar = h_planck/(2.*pi)
G = 6.6726e-8
q_e = 4.8032e-10
m_e = 9.1094e-28
m_p = 1.6726e-24
N_a = 6.0221e23
k_b = 1.3807e-16
a_rad = 7.5646e-15
sigma_sb = 5.6705e-5
sigma_T = 6.6524e-25  # Thomson cross section
T_CMB = 2.725
H_0 = 100.0 # in units of h km/s/Mpc

a_0 = 5.292e-9

# unit conversions
cmPerAU = 1.496e13
ergsPerCm = 1.6022e-12
cmPerParsec = 3.086e18
cmPerLightYear = 9.463e17
MpcPerKm = 3.241e-20
eVPerErg = 6.242e11
radiansToArcmin = 180.0/pi*60.0
cgsToJansky = 1e23

Ry = 13.6 # in eV

# solar parameters
M_sun = 1.99e33
L_sun = 3.839e33
T_sun = 5.78e3
R_sun = 6.96e10

M_earth = 5.97219e27 # in grams

# AB magnitude zero point
AB_flux_zero = 3631 # in Janskys

