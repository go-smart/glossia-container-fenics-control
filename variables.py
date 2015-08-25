"""Define simulation variables
"""

from gosmart.parameters import P, R

# Define the pairings of electrodes
# For each relevant pair, define a voltage
electrode_triples = P.CONSTANT_IRE_NEEDLEPAIR_VOLTAGE

tissues_utility = {
    "liver": ("organs", "TISSUE"),
    "needles": ("needles", "NEEDLES"),
    "vessels": ("vessels", "VESSELS"),
    "tumour": ("tumours", "TUMOURS"),
    "background": ("tissues", "BACKGROUND")
}

tissues = {}

for name, (group, suffix) in tissues_utility.items():
    tissues[name] = {
        "indices": [r.idx for r in R.group(group)],
        "relative permittivity": P["CONSTANT_IRE_RELATIVE_PERMITTIVITY_%s" % suffix],
        "sigma": [P["CONSTANT_IRE_ELECTRIC_CONDUCTIVITY_%s_%s" % (limit, suffix)] for limit in ("LOWER", "UPPER")],
        "threshold reversible": P["CONSTANT_IRE_ELECTRIC_CONDUCTIVITY_THRESHOLD_LOWER_%s" % suffix],
        "threshold irreversible": P["CONSTANT_IRE_ELECTRIC_CONDUCTIVITY_THRESHOLD_UPPER_%s" % suffix],
    }

# Define the domain and voxel size
# dim_width = 512
# dim_height = 352
# dim_depth = 21
# dim_bits = 16
# delta_width = 0.68359 / 1000.
# delta_height = 0.68359 / 1000.
# delta_depth = 2.7500 / 1000.
# offset_x = 0.147
# offset_y = 0.12
# offset_z = 0.021
# voxel_size = delta_width * delta_height * delta_depth

max_restarts = 10
