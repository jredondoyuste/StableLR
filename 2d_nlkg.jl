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
        "--l1"
            help = "First angular mode"
            arg_type = Int64
            default = 2
        "--l2"
            help = "Second angular mode"
            arg_type = Int64
            default = 3
        "--A1"
            help = "Amplitude of first mode"
            arg_type = Float64
            default = 1.0
        "--A2"
            help = "Amplitude of second mode"
            arg_type = Float64
            default = 0.0
        "--tf"
            help = "Final time"
            arg_type = Float64
            default = 50.0
    end
    return parse_args(s)
end

parsed_args = parse_commandline()

#############################

Nmodes  = parsed_args["N"]         
l0_1    = parsed_args["l1"]              
l0_2    = parsed_args["l2"]              
Amp_1   = parsed_args["A1"]         
Amp_2   = parsed_args["A2"] 
t_end   = parsed_args["tf"]         

dt = 0.5/Nmodes     # Time step
A_tot = sqrt(Amp_1^2 + Amp_2^2)
id = "/groups/astro/jredondo/2D_NLKG/data/Norm_growth/2d_nlkg_N_$(Nmodes)_l_$(l0_1)_$(l0_2)_A_$(round(A_tot,sigdigits=2))"
#############################

S = Legendre(-1..1)

freqs = [sqrt(l*(l+1)) for l in 0:Nmodes]

ics_phi = [Amp_1*(l==l0_1) + Amp_2*(l==l0_2) for l in 0:Nmodes]
ics_psi = [0.0 for l in 0:Nmodes]
ics = [ics_phi, ics_psi]

function nl(phi::AbstractArray)
    fun = Fun(S, phi)
    nl_term = -fun^3
    coefs = nl_term.coefficients[1:Nmodes+1]
    coefs[div(Nmodes, 2)+1:end] .= 0.0          # Dealiasing
    return coefs
end

function rhs!(du, u)
    phi = u[1]
    psi = u[2]
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
    fid["l"]=[l for l in 0:Nmodes]
    create_dataset(fid, "t", datatype(Float64), dataspace((Nsaves,)))
    create_dataset(fid, "phi", datatype(Float64), dataspace((Nsaves,Nmodes+1)))
    close(fid)
end

function write_chk(u, t, idx)
    fid=h5open("$(id).h5", "r+")
    fid["t"][idx] = t 
    fid["phi"][idx,:] = u[1]
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