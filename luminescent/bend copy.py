import luminescent as lumi
import gdsfactory as gf
import numpy as np
import os

radius = 5
# c.plot()

wavelengths = 1.55
keys = ["2,1"]  # same as keys=["o2@0,o1@0"]
nres = 10
gpu = None

c = gf.components.bend_circular(radius=radius, allow_min_radius_violation=True)
path = os.path.join("test", f"bend_R{radius}")
lumi.make_prob(
    path, c, wavelengths=wavelengths, nres=nres, keys=keys, gpu=gpu)

c = gf.components.bend_euler(radius=radius, )
path = os.path.join("test", f"euler_bend_R{radius}")
lumi.make_prob(
    path, c, wavelengths=wavelengths, nres=nres, keys=keys, gpu=gpu)
