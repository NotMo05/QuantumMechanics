import math
import matplotlib.pyplot as plt
from utils import *
import os

clearFolder("PotentialBarrier", "RealParts")
clearFolder("PotentialBarrier", "ImParts")
clearFolder("PotentialBarrier", "ProbDistributions")

PI = math.pi
e = math.e
i = 1j

m = 9.11e-31 #mass of quantum particle

L = replace_default("Enter the width of the barrier (in tenths of nanometers): ", 1) * 1e-10
V0 = replace_default("Enter the potential of the barrier (in electronVolts): ", 2) * 1.6e-19
E = replace_default("Enter the energy of the quantum particle (in electronVolts): ", 1) * 1.6e-19
mode = input("Enter e to simulate an electron or p to simulate a proton:")
if mode == "p":
    m = 1.67e-27  # Proton mass in kg

NO_OF_STEPS = 500 #Number of plots inside the barrier
rel_width = 1/18 #Relative width of the barrier on the plot
STEP = L / NO_OF_STEPS


A = 1 #Set A as 1 as this is arbitrary

#CALCULATE CONSTANTS
k = calculateWV(m, E)
q = calculateWV(m, V0-E)
F = A * 4*i*k*q*e**(q*L - i*k*L)
F /= (q + i*k)**2 - (e**(2*q*L))*((q-i*k)**2)
C = 0.5 * (1 + i*k/q) * e**(i*k*L-q*L) * F
D = 0.5 * (1 - i*k/q) * e**(i*k*L+q*L) * F
B = C + D - A


print("Plotting Region 1")
#REGION 1
minx = -(((1/rel_width) - 1) / 2) * L
x = minx
while x < 0:
  phi = A*e**(i*k*x) + B*e**(-i*k*x)
  real_phi = phi.real
  imag_phi = phi.imag
  phi_sq = real_phi**2 + imag_phi**2
  plt.figure(1)
  plt.plot(x, real_phi, marker = ',', color = 'r')
  plt.figure(2)
  plt.plot(x, imag_phi, marker = ',', color = 'r')
  plt.figure(3)
  plt.plot(x, phi_sq, marker = ',', color = 'orange')
  x += STEP

print("Plotting Region 2")
#REGION 2
while x < L:
  phi = C*e**(q*x) + D*e**(-q*x)
  real_phi = phi.real
  imag_phi = phi.imag
  phi_sq = real_phi**2 + imag_phi**2
  plt.figure(1)
  plt.plot(x, real_phi, marker = ',', color = 'b')
  plt.figure(2)
  plt.plot(x, imag_phi, marker = ',', color = 'b')
  plt.figure(3)
  plt.plot(x, phi_sq, marker = ',', color = 'green')
  x += STEP

print("Plotting Region 3")
#REGION 3
while x < L - minx:
  phi = F*e**(i*k*x)
  real_phi = phi.real
  imag_phi = phi.imag
  phi_sq = real_phi**2 + imag_phi**2
  plt.figure(1)
  plt.plot(x, real_phi, marker = ',', color = 'r')
  plt.figure(2)
  plt.plot(x, imag_phi, marker = ',', color = 'r')
  plt.figure(3)
  plt.plot(x, phi_sq, marker = ',', color = 'orange')
  x += STEP

#PLOT DETAILS
plt.figure(1)
plt.axvline(x=0, color = 'black')
plt.axvline(x=L, color = 'black')
plt.xlabel('Position / m')
plt.ylabel('Real part of wave function')
file_path = os.path.join("PotentialBarrier", "RealParts", "Real_Part_Of_WF.png")
plt.savefig(file_path)
plt.figure(2)
plt.axvline(x=0, color = 'black')
plt.axvline(x=L, color = 'black')
plt.xlabel('Position / m')
plt.ylabel('Imaginary part of wave function')
file_path = os.path.join("PotentialBarrier", "ImParts", "Im_Part_Of_WF.png")
plt.savefig(file_path)
plt.figure(3)
plt.axvline(x=0, color = 'black')
plt.axvline(x=L, color = 'black')
plt.xlabel('Position / m')
plt.ylabel('Magnitude of wave function squared')
file_path = os.path.join("PotentialBarrier", "ProbDistributions", "Probability_Distribution")
plt.savefig(file_path)