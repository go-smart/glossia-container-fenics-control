"""Define simulation variables
"""

from gosmart.parameters import P, R
import sys

# Define the pairings of electrodes
# For each relevant pair, define a voltage
electrode_triples = P.CONSTANT_IRE_NEEDLEPAIR_VOLTAGE

tissues_utility = {
    "liver": ("organs", "TISSUE"),
    "needles": ("needles", "NEEDLE"),
    "vessels": ("vessels", "VESSELS"),
    "tumour": ("tumours", "TUMOUR"),
    "background": ("tissues", "BACKGROUND")
}

tissues = {}

for name, (group, suffix) in tissues_utility.items():
    map = R.group(group)
    if map is None:
        print("Failed to find a required surface")
        sys.exit(1)
    tissues[name] = {
        "indices": [r.idx for r in R.group(group)],
        "relative permittivity": P["CONSTANT_IRE_RELATIVE_PERMITTIVITY_%s" % suffix]
    }
    try:
        tissues[name]["sigma"] = [P["CONSTANT_IRE_ELECTRIC_CONDUCTIVITY_%s_%s" % (limit, suffix)] for limit in ("LOWER", "UPPER")]
    except KeyError:
        tissues[name]["sigma"] = P["CONSTANT_IRE_ELECTRIC_CONDUCTIVITY_%s" % suffix]
    else:
        tissues[name].update({
            "threshold reversible": P["CONSTANT_IRE_ELECTRIC_CONDUCTIVITY_THRESHOLD_LOWER_%s" % suffix],
            "threshold irreversible": P["CONSTANT_IRE_ELECTRIC_CONDUCTIVITY_THRESHOLD_UPPER_%s" % suffix],
        })

max_restarts = 10
