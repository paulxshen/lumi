import itertools
import os
import luminescent as lumi
import gdsfactory as gf

# dir = "runs"
dir = os.path.join("build", "precompile_execution")
# lumi.solve(os.path.join(dir, "tiny_2_float32_None"))


c = gf.components.straight(length=0.5,)
wavelengths = [1, 1.5]

# for N, dtype, gpu in itertools.product([2, 3], ["float32"], [None, "CUDA"]):
for N, dtype, gpu in itertools.product([2, 3], ["float32"], [None]):
    path = os.path.join(dir, f"tiny_{N}_{dtype}_{gpu}")
    approx_2D_mode = "TE" if N == 2 else None
    lumi.make_prob(path, c, wavelengths=wavelengths, keys=[
        "2,1"], nres=4, approx_2D_mode=approx_2D_mode, gpu=gpu, dtype=dtype)


# raise NotImplementedError("This is a stub")
c = lumi.mimo(west=1, east=1, l=.5, w=.5,  wwg=.5)
targets = {"tparams": {
    1.5: {
        "2,1": 1.0
    }}}
for dtype in ['float16']:
    path = os.path.join(dir, f"back_{dtype}")
    lumi.make_design_prob(
        path,  c, targets,
        fill_material="Si", void_material="SiO2",
        lvoid=0.2, iters=2, nres=4,
        approx_2D_mode="TE", dtype=dtype)


# lumi.finetune(os.path.join("runs", "back"), iters=2)
0
