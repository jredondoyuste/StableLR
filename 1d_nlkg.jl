using ApproxFun 
using HDF5
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--N"
            help = "Number of modes"
            arg_type = Int64
            default = 128
        "--l"
            help = "First mode"
            arg_type = Int64
            default = 1
        "--A"
            help = "Amplitude of first mode"
            arg_type = Float64
            default = 0.5
        "--tf"
            help = "Final time"
            arg_type = Float64
            default = 30.0
    end
    return parse_args(s)
end

parsed_args = parse_commandline()

#############################

Nmodes  = parsed_args["N"]         
l0      = parsed_args["l"]              
Amp     = parsed_args["A"]         
t_end   = parsed_args["tf"]         

dt = 0.5/Nmodes     # Time step
A_tot = sqrt(Amp^2)
id = "/groups/astro/jredondo/2D_NLKG/data_1d/Convergence/1d_nlkg_N_$(Nmodes)_l_$(l0)_A_$(round(A_tot,sigdigits=2))"
#############################

S = Fourier(0..2π)

freqs = [l for l in 0:Nmodes]
is2 = 1/sqrt(2)
ics_phi = [is2*Amp*(l==2*l0)+is2*Amp*(l==2*l0+2) for l in 0:Nmodes]
# ics_phi = [Amp*(l==2*l0+1) for l in 0:Nmodes]
ics_psi = [0.0 for l in 0:Nmodes]
ics = [ics_phi, ics_psi]

function nl(phi::AbstractArray)
    fun = Fun(S, phi)
    nl_term = -fun^3
    coefs = 0*phi
    length_nl = length(nl_term.coefficients)
    for i in 1:min(Nmodes+1,length_nl)
        coefs[i] = nl_term.coefficients[i]
    end
    if length_nl < Nmodes+1
        coefs[length_nl+1:end] .= 0.0
    end
    coefs[div(Nmodes, 2)+1:end] .= 0.0          # Dealiasing
    return coefs
end

function rhs!(du, u)
    phi = u[1]
    du[1] = u[2]
    nlterm = nl(phi)
    du[2] = @. - freqs^2 * phi + nlterm
end 

function rk4!(u_new, u, h)
    k1 = 0*u
    k2 = 0*u
    k3 = 0*u
    k4 = 0*u
    rhs!(k1, u)
    rhs!(k2, u + h/2*k1)
    rhs!(k3, u + h/2*k2)
    rhs!(k4, u + h*k3)
    u_new .= u + h/6*(k1 + 2*k2 + 2*k3 + k4)
end

#############################

function create_file(Nsaves)
    fid=h5open("$(id).h5", "w")
    fid["l"]=[l for l in 0:div(Nmodes,2)]
    create_dataset(fid, "t", datatype(Float64), dataspace((Nsaves,)))
    len2=div(Nmodes,2)+1
    create_dataset(fid, "phi", datatype(Float64), dataspace((Nsaves,len2)))
    close(fid)
end

function write_chk(u, t, idx)
    fid=h5open("$(id).h5", "r+")
    fid["t"][idx] = t 
    fid["phi"][idx,:]= u[1][1:2:end]
    close(fid)
end

#############################

u_new = ics
t = 0.0
t_chk = 0.1
create_file(10*Int(t_end)+1)
write_chk(u_new, t, 1)
counter = 1

while t < t_end
    global u = u_new
    rk4!(u_new, u, dt)
    global t += dt
    if mod(t,t_chk) < dt
        max_phi = maximum(abs.(u_new[1]))
        if max_phi > 100 * A_tot
            println("Amplitude grew too much at t = $t")
            break
        end
        println("t = $t, max(phi) = $(max_phi)")
        global counter += 1
        write_chk(u_new, t, counter)
    end
end