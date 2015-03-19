"""Define simulation variables
"""

import itertools as it

# Define the electrodes
electrodes = {
    "positive": [1, 2, 3, 4],
    "negative": [5, 6]
}

# Define the pairings of electrodes
electrode_pairs = list(it.product(electrodes["positive"], electrodes["negative"]))
# QLJ: Why is there this final pair?
electrode_pairs.append((5, 6))

# For each relevant pair, define a voltage
voltages = [1300, 1500, 1300, 1900, 1300, 1300, 1300, 1900, 1300]
tissues = {"liver": {
    "indices": (11,),
    "sigma": (0.08, 0.32),
    "threshold reversible": 35000,
    "threshold irreversible": 70000,
    "relative permittivity": 1.,
}, "needles": {
    "indices": (1, 2, 3, 4, 5, 6),
    "relative permittivity": 1.,
    "threshold reversible": 35000,
    "threshold irreversible": 70000,
    "sigma": (1e3, 1e3),
}, "vessels": {
    "indices": (8, 9, 10,),
    "sigma": (0.7, 1.5),
    "threshold reversible": 10000,
    "threshold irreversible": 45000,
    "relative permittivity": 1.,
}, "tumour": {
    "indices": (7,),
    "sigma": (0.2, 0.7),
    "threshold reversible": 40000,
    "threshold irreversible": 80000,
    "relative permittivity": 1.,
}, "background": {
    "indices": (12,),
    "sigma": (0.02, 0.08),
    "threshold reversible": 10000,
    "threshold irreversible": 40000,
    "relative permittivity": 1.,
}}

# Define the domain and voxel size
dim_width = 512
dim_height = 352
dim_depth = 21
dim_bits = 16
delta_width = 0.68359 / 1000.
delta_height = 0.68359 / 1000.
delta_depth = 2.7500 / 1000.
offset_x = 0.147
offset_y = 0.12
offset_z = 0.021
voxel_size = delta_width * delta_height * delta_depth

max_restarts = 1

patient_data_location = "../../scratch/IRE/data/"
