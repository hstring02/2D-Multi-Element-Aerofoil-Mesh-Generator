# Persistent k-omega solve worker for the batch runner (02_Tool_Scripts/01_Batch_Run_Script)
# and RUN_CASE.py (01_Input_Files).
#
# Runs once and solves any number of meshes fed to it over stdin, so
# XCALibre's Physics/run! machinery only pays its JIT/precompile cost once
# for the whole sweep instead of once per mesh. Protocol (plain text, no
# extra Julia package dependency):
#
#   stdout, on startup : ##READY##
#   stdin,  per case    : <mesh_path>|<toml_path>, or ##QUIT## to exit
#   stdout, on success  : ##RESULT## <mesh_path>|<Cl>|<Cd>|<lift>|<drag>
#   stdout, on failure  : ##ERROR## <mesh_path>|<message>
#
# toml_path is the flow-config file for the case, supplied explicitly by the
# caller (the batch runner writes per-case TOMLs under 01_Input_Files/Batch,
# not next to a same-named mesh, so it can't be guessed from the mesh path).

using XCALibre
using TOML

const BACKEND   = CPU()
const HARDWARE  = Hardware(backend=BACKEND, workgroup=1024)
activate_multithread(BACKEND)

function run_case(mesh_file::AbstractString, toml_file::AbstractString)
    if !isfile(mesh_file)
        error("Mesh file not found: $mesh_file")
    end
    if !isfile(toml_file)
        error("Flow config file not found: $toml_file")
    end
    flow_config = TOML.parsefile(toml_file)["flow"]

    mesh = UNV2D_mesh(mesh_file, scale=1)
    mesh_dev = adapt(BACKEND, mesh)

    rho = flow_config["DENSITY"]
    nu = flow_config["VISCOSITY"]
    u_mag = flow_config["VELOCITY"]
    chord = flow_config["REYNOLDS_LENGTH"]
    velocity = [u_mag, 0.0, 0.0]
    Tu = 0.05
    nuR = 100
    k_inlet = 3/2*(Tu*u_mag)^2
    ω_inlet = k_inlet/(nuR*nu)
    νt_inlet = k_inlet / ω_inlet

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
            solver = Bicgstab(), preconditioner = Jacobi(),
            convergence = 1e-7, relax = 0.5, rtol = 1e-2, atol = 1e-10
        ),
        p = SolverSetup(
            solver = Cg(), preconditioner = Jacobi(),
            convergence = 1e-7, relax = 0.2, rtol = 1e-3, atol = 1e-10
        ),
        k = SolverSetup(
            solver = Bicgstab(), preconditioner = Jacobi(),
            convergence = 1e-7, relax = 0.7, rtol = 1e-2, atol = 1e-10
        ),
        omega = SolverSetup(
            solver = Bicgstab(), preconditioner = Jacobi(),
            convergence = 1e-7, relax = 0.7, rtol = 1e-2, atol = 1e-10
        )
    )

    runtime = Runtime(iterations=1000, write_interval=-1, time_step=1)

    config = Configuration(
        solvers=solvers, schemes=schemes, runtime=runtime, hardware=HARDWARE, boundaries=BCs)

    GC.gc()

    initialise!(model.momentum.U, velocity)
    initialise!(model.momentum.p, 0.0)
    initialise!(model.turbulence.k, k_inlet)
    initialise!(model.turbulence.omega, ω_inlet)
    initialise!(model.turbulence.nut, νt_inlet)

    run!(model, config)

    Fp = pressure_force(:airfoil, model.momentum.p, rho)
    Fv = viscous_force(:airfoil, model.momentum.U, rho, nu, model.turbulence.nut, config)
    Ft = Fp + Fv

    drag = Ft[1]
    lift = Ft[2]
    q = 0.5 * rho * u_mag^2 * chord

    return (Cl = lift / q, Cd = drag / q, lift = lift, drag = drag)
end

function main()
    println("##READY##")
    flush(stdout)

    for line in eachline(stdin)
        line = strip(line)
        isempty(line) && continue
        line == "##QUIT##" && break

        mesh_file, toml_file = split(line, "|"; limit=2)
        try
            r = run_case(mesh_file, toml_file)
            println("##RESULT## $(mesh_file)|$(r.Cl)|$(r.Cd)|$(r.lift)|$(r.drag)")
        catch e
            msg = replace(sprint(showerror, e), "\n" => " | ")
            println("##ERROR## $(mesh_file)|$(msg)")
        end
        flush(stdout)
    end
end

main()
