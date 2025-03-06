import luminescent as lumi
import gdsfactory as gf
import numpy as np
import os

radius = 5
c = gf.components.bend_circular(radius=radius, allow_min_radius_violation=True)
# c.plot()

path = os.path.join("test", f"bend_R{radius}")
wavelengths = 1.55
keys = ["2,1"]  # same as keys=["o2@0,o1@0"]
nres = 8
gpu = None

lumi.make_pic_sim_problem(
    path, c, wavelengths=wavelengths, nres=nres, keys=keys, gpu=gpu)
