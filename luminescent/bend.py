import luminescent as lumi
import gdsfactory as gf
import numpy as np
import os

radius = 5
# wavelengths = np.linspace(1.4, 1.6, 5)
keys = ["2,1"]  # same as keys=["o2@0,o1@0"]
nres = 10
gpu = None
materials = lumi.MATERIALS
layer_stack = lumi.SOI

c = gf.components.bend_circular(radius)

path = os.path.join("test", f"bend_R{radius}")
wavelengths = np.linspace(1.45, 1.65, 7)
# wavelengths = 1.55
lumi.make_pic_sim_problem(path, c, wavelengths=wavelengths, nres=nres, keys=keys, gpu=gpu,
                          materials=materials, layer_stack=layer_stack)

# c = gf.components.straight(1)
# path = os.path.join("test", f"straight")
# nres = 8
# lumi.make_pic_sim_problem(
#     path, c, wavelengths=wavelengths, nres=nres, keys=keys, gpu=gpu)
