# using Luminescent
include("main.jl")
ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
ENV["JULIA_DEBUG"] = "Main"
Random.seed!(1234)
# using CUDA
RUNS = joinpath("build", "precompile_execution")
# picrun(joinpath("test", "straight");)# array=cu)
# picrun(joinpath("test", "bend_R5"), cu)
# picrun(joinpath("test", "euler_bend_R5"), cu)
# picrun(joinpath("runs", "mode_converter"))

# picrun(joinpath(RUNS, "tiny_2_float32_CUDA"), cu)#; framerate=10)
# # picrun(joinpath(RUNS, "tiny_3_float32_CUDA"), cu)
# picrun(joinpath(RUNS, "tiny_3_float32_None"))
# picrun(joinpath(RUNS, "back_float32"))
# picrun(joinpath(RUNS, "back_float32"))# array=cu)

# prob.canvas_instances[1]._frame[:ϵ]|>extrema
# picrun(joinpath("runs", "demux"))

for p = readdir("build/precompile_execution", join=true)
    if contains(string(p), "16") #&& contains(string(p), "back")
        picrun(p)
    end
end
# # 

# using Pkg
# pkg"
# dev C:\Users\pxshe\OneDrive\Desktop\beans\Porcupine.jl;
# dev C:\Users\pxshe\OneDrive\Desktop\beans\ArrayPadding.jl;
#  dev C:\Users\pxshe\OneDrive\Desktop\beans\Jello.jl;
# up"

# using GLMakie: volume
# volume(prob.geometry)
# volume(prob.source_instances[1].sigmodes[1][2].Jx |> cpu .|> abs)
# volume(_gf2[1].Jx |> cpu .|> abs)
# volume(sols[1].u.H.Hz |> cpu .|> abs)
# prob.geometry.invϵ |> cpu |> first |> extrema
# volume(prob.geometry.invϵ |> cpu |> first)
