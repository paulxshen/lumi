import luminescent as lumi
import gdsfactory as gf
import numpy as np
import os

radius = 5
keys = ["2,1"]  # same as keys=["o2@0,o1@0"]
nres = 6
gpu = None
materials = lumi.MATERIALS
layer_stack = lumi.SOI
wavelengths = np.linspace(1.5, 1.6, 5)

c = gf.components.bend_circular(radius)
path = "build/precompile_execution/back_float16"

# lumi.make_prob(
#     path, c, wavelengths=wavelengths, nres=nres, keys=keys, gpu=gpu,
#     materials=materials, layer_stack=layer_stack)
lumi.load(path)
