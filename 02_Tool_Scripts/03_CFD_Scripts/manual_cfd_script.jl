using Plots
using XCALibre
using TOML

grids_dir = joinpath(@__DIR__, "..", "..", "04_Meshes")
toml_dir  = joinpath(@__DIR__, "..", "..", "01_Input_Files")

mesh_file = joinpath(grids_dir, "EXAMPLE_3_el_wing.unv")
toml_file = joinpath(toml_dir, "EXAMPLE_mesh_3_el_wing.toml")

toml_config = TOML.parsefile(toml_file)
flow_config = toml_config["flow"]
total_chord = sum(toml_config["foils"]["CHORD"]) # sum of all element chords, m (reference length for Cl/Cd normalisation)
mesh = UNV2D_mesh(mesh_file, scale=1)

backend = CPU(); workgroup = 1024; activate_multithread(backend)
hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

rho = flow_config["DENSITY"]           # kg/m3
nu = flow_config["VISCOSITY"]          # kinematic viscosity, m^2/s
u_mag = flow_config["VELOCITY"]        # m/s
velocity = [u_mag, 0.0, 0.0]
Tu = 0.05
nuR = 100
k_inlet = 3/2*(Tu*u_mag)^2
ω_inlet = k_inlet/(nuR*nu)
νt_inlet = k_inlet/ω_inlet
Re = velocity[1]*total_chord/nu

model = Physics(
    time = Steady(),
    fluid = Fluid{Incompressible}(nu = nu),
    turbulence = RANS{KOmega}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )

BCs = assign(
    region = mesh_dev,
    (
        U = [
            Dirichlet(:inlet, velocity),
            Zerogradient(:outlet),
            Wall(:airfoil, [0.0, 0.0, 0.0]),
            Symmetry(:farfield)
        ],
        p = [
            Zerogradient(:inlet),
            Dirichlet(:outlet, 0.0),
            Wall(:airfoil),
            Symmetry(:farfield)
        ],
        k = [
            Dirichlet(:inlet, k_inlet),
            Zerogradient(:outlet),
            KWallFunction(:airfoil),
            Symmetry(:farfield)
        ],
        omega = [
            Dirichlet(:inlet, ω_inlet),
            Zerogradient(:outlet),
            OmegaWallFunction(:airfoil),
            Symmetry(:farfield)
        ],
        nut = [
            Dirichlet(:inlet, νt_inlet),
            Zerogradient(:outlet),
            NutWallFunction(:airfoil),
            Symmetry(:farfield)
        ]
    )
)

schemes = (
    U = Schemes(divergence=Upwind),
    p = Schemes(divergence=Upwind),
    k = Schemes(divergence=Upwind),
    omega = Schemes(divergence=Upwind)
)

solvers = (
    U = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(),
        convergence = 1e-7,
        relax       = 0.5,
        rtol = 1e-2,
        atol = 1e-10
    ),
    p = SolverSetup(
        solver      = Cg(), #Gmres(), #Cg(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(),
        convergence = 1e-7,
        relax       = 0.2,
        rtol = 1e-3,
        atol = 1e-10
    ),
    k = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(),
        convergence = 1e-7,
        relax       = 0.7,
        rtol = 1e-2,
        atol = 1e-10
    ),
    omega = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(),
        convergence = 1e-7,
        relax       = 0.7,
        rtol = 1e-2,
        atol = 1e-10
    )
)

runtime = Runtime(iterations=1000, write_interval=100, time_step=1)
# runtime = Runtime(iterations=2, write_interval=-1, time_step=1)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)


GC.gc()

initialise!(model.momentum.U, velocity)
initialise!(model.momentum.p, 0.0)
initialise!(model.turbulence.k, k_inlet)
initialise!(model.turbulence.omega, ω_inlet)
initialise!(model.turbulence.nut, νt_inlet)

residuals = run!(model, config)

# Lift and drag calculation
Fp = pressure_force(:airfoil, model.momentum.p, rho)
Fv = viscous_force(:airfoil, model.momentum.U, rho, nu, model.turbulence.nut, config)
Ft = Fp + Fv # total force per unit span, N/m (inflow is horizontal, so no rotation needed)

drag = Ft[1]
lift = Ft[2]
q = 0.5*rho*u_mag^2*total_chord

Cl = lift/q
Cd = drag/q

println("Lift: ", round(lift, sigdigits=4), " N/m   Drag: ", round(drag, sigdigits=4), " N/m")
println("Cl: ", round(Cl, sigdigits=4), "   Cd: ", round(Cd, sigdigits=4))

println("L/D ratio: ", round(Cl/Cd, sigdigits=4))

println("Reynolds number: ", round(Re, sigdigits=4))