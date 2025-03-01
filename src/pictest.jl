# using Luminescent
include("main.jl")
ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
Random.seed!(1234)
using CUDA

# picrun(joinpath("runs", "straight");)# array=cu)
# picrun(joinpath("runs", "bend_R5"), array=cu)
# picrun(joinpath("runs", "mode_converter"))

# picrun(joinpath("runs", "splitter"); array=cu)
# picrun(joinpath("runs", "splitter"))
# using GLMakie: volume
# picrun(joinpath("build", "precompile_execution", "tiny_2_float32_CUDA"))#; framerate=10)
picrun(joinpath("build", "precompile_execution", "tiny_3_float32_CUDA"), cu)
# picrun(joinpath("build", "precompile_execution", "tiny_3_float32_None"))
# picrun(joinpath("build", "precompile_execution", "back_float32"))
# picrun(joinpath("runs", "tiny3"))
# picrun(joinpath("runs", "back"))# array=cu)
# models[1]()

# picrun(joinpath("runs", "demux"))
# picrun(joinpath("runs", "straight"))

# for p = readdir("build/precompile_execution", join=true)
#     # if !contains(string(p), "16") && !contains(string(p), "back")
#     if !contains(string(p), "16") #&& contains(string(p), "back")
#         # if !contains(string(p), "16") && contains(string(p), "back")
#         picrun(p)
#     end
# end
# # 

# using Pkg
# pkg"
# dev C:\Users\pxshe\OneDrive\Desktop\beans\Porcupine.jl;
# dev C:\Users\pxshe\OneDrive\Desktop\beans\ArrayPadding.jl;
#  dev C:\Users\pxshe\OneDrive\Desktop\beans\Jello.jl;
# up"

using GLMakie: volume
# volume(prob.geometry)
# volume(prob.source_instances[1].sigmodes[1][2].Jx |> cpu .|> abs)
# volume(_gf2[1].Jx |> cpu .|> abs)
# volume(sols[1].u.H.Hz |> cpu .|> abs)
# prob.geometry.invϵ |> cpu |> first |> extrema
# volume(prob.geometry.invϵ |> cpu |> first)