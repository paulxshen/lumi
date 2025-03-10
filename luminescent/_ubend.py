import luminescent as lumi
import os

path = os.path.join('runs', "ubend")
wwg = .5
gap = .5
c = lumi.mimo(west=[wwg/2, 1.5*wwg+gap], l=6,
              w=2*wwg+gap,  wwg=wwg, taper=0.0)
targets = {"tparams": {1.55: {"2,1": 1.0}}}
# c.show()

prob = lumi.make_pic_inv_prob(
    path,   c, targets=targets,
    symmetries=[1,],
    nres=30,
    lvoid=0.15, lsolid=.1,
    iters=50, stoploss=.03, Ttrans='2.5x',
    dtype='float16', approx_2D_mode="TE")
sol = lumi.solve(path)
lumi.load_res(path)
