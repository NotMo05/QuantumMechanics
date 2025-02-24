import math
import matplotlib.pyplot as plt
from utils import *
import os

clearFolder("PotentialWell", "ProbDistributions")
clearFolder("PotentialWell", "WaveFuncs")
clearFolder("PotentialWell", "EnergyLevels")

PI = math.pi
e = math.e
TRIES = 1000
SF_PRECISION = 10
STABILISE_REPEATS = 15
NO_OF_STEPS = 1000 #Number of plots inside the well


m = 1.67e-27

L = replace_default("Enter the width of the well (in tenths of nanometers): ", 1) * 1e-10  # Convert to meters
V1 = replace_default("Enter the LHS potential (in electronVolts): ", 0.1) * 1.6e-19  # Convert to Joules
V2 = replace_default("Enter the RHS potential (in electronVolts): ", 1) * 1.6e-19  # Convert to Joules
mode = input("Enter e to simulate an electron or p to simulate a proton:")
if mode == "e":
    m = 9.11e-31  # Proton mass in kg

a = L*calculateWV(m, 1)

STEP = L / NO_OF_STEPS

if V1 < V2:
  smaller = 'V1'
else:
  smaller = 'V2'

max_z = math.sqrt(min([V1, V2]))
Asymptote = math.sqrt((V1 * V2) / (V1 + V2))
Period = PI / a

def save_plot(n, thing, ylabel):
  if ylabel == 'Wave function':
    plt.figure(2*n + 1)
  else:
     plt.figure(2*n + 2)
  plt.axvline(x=0, color = 'black')
  plt.axvline(x=L, color = 'black')
  plt.xlabel('Position / m')
  plt.ylabel(ylabel)
  plt.text(-0.75*L, thing, f'V1 = {V1}', color = 'black')
  plt.text(1.25*L, thing, f'V2 = {V2}', color = 'black')
  if ylabel == 'Wave function':
    file_path = os.path.join("PotentialWell", "WaveFuncs", f"Wavefunction_{n+1}.png")
  else:
    file_path = os.path.join("PotentialWell", "ProbDistributions", f"Probability_Distribution_{n+1}.png")
  plt.savefig(file_path)

def sf(x, SF):
  i = 0
  while str(x)[i] == '0':
    i += 1
  return str(x)[:SF + i + 1]

def round_solution(solution, SF_PRECISION):
  significant_digits = sf(solution, SF_PRECISION)
  mult = 0
  for i in range(len(str(solution))):
    if str(solution)[i] == 'e':
      mult = int(str(solution)[i+1:])
  return float(significant_digits) * (10**mult)

#FUNCTION TO NUMERICALLY SOLVE
def f(z, check = 0):
  zs = z**2
  if check == 'V1':
    zs = V1
  elif check == 'V2':
    zs = V2
  LHS = math.tan(a*z)
  RHS = math.sqrt(V1 - zs) + math.sqrt(V2 - zs)
  RHS /= z - math.sqrt(z**2 - (V1+V2) + ((V1*V2)/zs))
  return LHS - RHS

#DERIVATIVE OF FUNCTION
def f_prime(z):
  S = math.sqrt(z**2 - (V1+V2) + ((V1*V2)/(z**2)))
  sec = 1/math.cos(a*z)
  term1 = (a*sec**2)
  term2 = z/math.sqrt(V1 - z**2) + z/math.sqrt(V2 - z**2)
  term2 /= z - S
  term3 = math.sqrt(V1 - z**2) + math.sqrt(V2 - z**2)
  term3 *= 1 - ((z - (V1*V2)/(z**3)) / S)
  term4 = (z - S) ** 2
  return term1 + term2 + (term3 / term4)

def Newton_Raphson(start):
  count = 0
  OK = False
  solution = -1
  prev_x = -1
  x_n = start
  rounded_solution = -1
  for cycle in range(TRIES):
    while x_n > max_z:
      x_n = 0.5 * (x_n + prev_x)
    f_x = f(x_n)
    fp_x = f_prime(x_n)
    x_n = x_n - (f_x / fp_x)
    if sf(x_n, SF_PRECISION) == sf(prev_x, SF_PRECISION):
      count += 1
    else:
      count = 0
      prev_x = x_n
    if count > STABILISE_REPEATS:
      solution = x_n
      rounded_solution = round_solution(solution, SF_PRECISION)
      OK = True
      break
  return OK, rounded_solution

def normalise_to_find_B(k, l, u):
  term1 = 1/(2*k)
  term2 = (l**2 - k**2) * math.sin(2*l*L)
  term2 -= 4*l*L * (math.cos(l*L))**2
  term2 += 2*L * l**3
  term2 += 2*l*L * k**2
  term2 += 4*k*L
  term2 /= 4 * l**3
  term3 = (math.cos(l*L) + (k/l) * math.sin(l*L)) ** 2
  term3 /= 2*u
  B = 1 / math.sqrt(term1 + term2 + term3)
  return B

#FIND NUMBER OF SOLUTIONS
number_below_asymptote = round(Asymptote / Period)
number_above_asymptote = round(max_z / Period) - round(Asymptote / Period)
if f(max_z, check = smaller) > 0:
  number_above_asymptote += 1
number_of_solutions = number_below_asymptote + number_above_asymptote
solutions = ()

#BEFORE ASYMPTOTE
start = 0
for n in range(number_below_asymptote - 1):
  start += Period
  OK, solution = Newton_Raphson(start)
  while not OK:
    start *= 1.01
    OK, solution = Newton_Raphson(start)
  solutions += (solution,)

#JUST BEFORE ASYMPTOTE
prev_tan_root = start
start = 0.5 * (prev_tan_root + Asymptote)
OK = False
while not OK:
  OK, solution = Newton_Raphson(start)
  if solution in solutions:
    OK = False
  start = 0.5 * (start + Asymptote)
solutions += (solution,)

#JUST AFTER ASYMPTOTE
next_tan_root = prev_tan_root + Period
start = 0.5*(next_tan_root + Asymptote)
OK = False
while not OK:
  OK, solution = Newton_Raphson(start)
  if solution in solutions:
    OK = False
  start = 0.5 * (start + Asymptote)
solutions += (solution,)

#AFTER ASYMPTOTE
start = next_tan_root
for n in range(number_above_asymptote - 1):
  start += Period
  OK, solution = Newton_Raphson(start)
  while not OK:
    start *= 1.01
    OK, solution = Newton_Raphson(start)
  solutions += (solution,)
energy_levels = [x**2 for x in solutions]

#SET UP PLOTS
plt.figure(0)
for n in range(number_of_solutions):
  E = energy_levels[n]
  ax1 = plt.axes(frameon=False)
  ax1.axes.get_xaxis().set_visible(False)
  plt.axhline(y=E, xmin=0, xmax=1, color = 'b')
  plt.plot(1, E, marker = ',', color='b')
  plt.text(1.06, E, f'n = {n+1}')
  plt.axhline(y=V1, xmin=0, xmax=1, color = 'r')
  plt.plot(1, V1, marker = ',', color='r')
  plt.text(1.06, V1, 'V1')
  plt.axhline(y=V2, xmin=0, xmax=1, color = 'r')
  plt.plot(1, V2, marker = ',', color='r')
  plt.text(1.06, V2, 'V2')
  file_path = os.path.join("PotentialWell", "EnergyLevels", f"Energy_Levels_{n+1}.png")
  plt.savefig(file_path)

#SOLVE EACH ENERGY LEVEL
for n in range(number_of_solutions):
  max_phi = 0
  record_x = 0
  E = energy_levels[n]
  k = calculateWV(m, V1-E)
  l = calculateWV(m, E)
  u = calculateWV(m, V2-E)
  B = normalise_to_find_B(k, l, u)
  normalisation_sum_check = 0

  #REGION 1
  x = -L
  while x < 0:
    phi = B * (e ** (k * x))
    plt.figure(2*n + 1)
    plt.plot(x, phi, marker = ',', color = 'r')
    plt.figure(2*n + 2)
    plt.plot(x, phi**2, marker = ',', color = 'orange')
    normalisation_sum_check += STEP * phi**2
    x += STEP

  #REGION 2
  while x < L:
    phi = B * (math.cos(l*x) + (k/l) * math.sin(l*x))
    plt.figure(2*n + 1)
    plt.plot(x, phi, marker = ',', color = 'b')
    plt.figure(2*n + 2)
    plt.plot(x, phi**2, marker = ',', color = 'g')
    normalisation_sum_check += STEP * phi**2
    x += STEP
    if n == 0:
      if phi > max_phi:
        max_phi = phi
        record_x = x

  #REGION 3
  while x < 2*L:
    phi = B * (math.cos(l*L) + (k/l) * math.sin(l*L)) * e**(u*(L-x))
    plt.figure(2*n + 1)
    plt.plot(x, phi, marker = ',', color = 'r')
    plt.figure(2*n + 2)
    plt.plot(x, phi**2, marker = ',', color = 'orange')
    normalisation_sum_check += STEP * phi**2
    x += STEP

  #PLOT DETAILS
  print(f"Plotting solution {n+1} out of {number_of_solutions}")
  save_plot(n, 0.75*max_phi, 'Wave function')
  save_plot(n, 0.75*max_phi**2, 'Probility density')


