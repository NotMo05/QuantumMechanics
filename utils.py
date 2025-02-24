import math
from scipy.constants import hbar
import os

def clearFolder(mainFolder, subFolder):
  folder_path = os.path.join(mainFolder, subFolder)
  os.makedirs(folder_path, exist_ok=True)
  for filename in os.listdir(folder_path):
      file_path = os.path.join(folder_path, filename)
      if os.path.isfile(file_path):
          os.remove(file_path)

def replace_default(msg, default):
    user_input = input(msg)
    if user_input:
        return float(user_input)  # Convert the input to float if the user provides a value
    return default

def calculateWV(m, E):
  return (1 / hbar) * math.sqrt(2*m*E)
