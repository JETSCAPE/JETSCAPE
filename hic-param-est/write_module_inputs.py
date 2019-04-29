#a function to write module input files given design params
import math

def write_module_inputs(
    outdir = 'test/',
    #design point
    design_pt = 0,

    #trento params
    projectile = 'Pb',
    target = 'Pb',
    sqrts = 2760,
    cross_section = 6.4,
    normalization = 30.0,
    reduced_thickness = 1.0,
    fluctuation = 0.3,
    nucleon_width = 2.0,
    nucleon_min_dist = 1.27,

    #freestreaming params
    tau_s = 1.0,

    #hydro params
    T_s = 0.120,

    #shear viscosity p'zation
    eta_over_s_min = 1.0,
    eta_over_s_slope = 0.0,
    eta_over_s_curv = 0.0,
    #bulk viscosity p'zation
    bulk_viscosity_normalisation = 0.05,
    bulk_viscosity_width_in_GeV = 0.02,
    bulk_viscosity_peak_in_GeV = 0.18
    ):

    #choose grid spacing based on nucleon width
    dx = 0.15 * nucleon_width #[fm]
    dy = 0.15 * nucleon_width #[fm]

    #choose a grid size large enough to capture central events
    #Does [-15fm, 15fm] work?
    Lx = 15.0 #[fm]
    Ly = 15.0 #[fm]

    nx = 1.0 + (2.0 * Lx)/dx
    ny = 1.0 + (2.0 * Ly)/dy

    #get nearest integer
    nx = math.ceil(nx)
    ny = math.ceil(ny)

    #recalculate grid size
    Lx = ( (nx-1) / 2.0) * dx
    Ly = ( (nx-1) / 2.0) * dx

    #TRENTo needs slightly larger for odd number of points
    max_x = Lx + 0.5*dx #max x [fm]
    max_y = Ly + 0.5*dy #may y [fm]

    #Parameters common to modules
    #nx = 301 #num grid points in x
    #ny = 301 #num grid points in y
    #dx = 0.2 #grid spacing x [fm]
    #dy = 0.2 #grid spacing y [fm]
    #L_x = (nx - 1)/2.0 * dx #size of grid in x[fm]
    #L_y = (ny - 1)/2.0 * dy #size of grid in y[fm]
    #set max_x = L_x + 0.5dx
    #max_x = L_x + 0.5*dx #max x [fm]
    #max_y = L_y + 0.5*dy #may y [fm]

    dtau = 0.02 #hydro time step [fm/c]

    #dummy parameter unused
    e_s = 1.7

    #write appropriate input files
    cent_low = 0
    cent_high = 100

    #freestream file
    fs_file = open(outdir + 'freestream_input_' + str(design_pt),'w')

    fs_file.write("OUTPUTFORMAT 2\n")
    fs_file.write("BARYON 0\n")
    fs_file.write("IC_ENERGY 5\n")
    fs_file.write("IC_BARYON 1\n")
    fs_file.write("ETA_WIDTH 0.5\n")
    fs_file.write("ETA_FLAT 0.5\n")
    fs_file.write("SIGMA 0.5\n")
    fs_file.write("SIGMA_B 0.5\n")
    fs_file.write("DIM_X " + str(nx) + "\n")
    fs_file.write("DIM_Y " + str(ny) + "\n")
    fs_file.write("DIM_ETA 1\n")
    fs_file.write("DIM_RAP 1\n")
    fs_file.write("DIM_PHIP 5000\n")
    fs_file.write("DX " + str(dx) + "\n")
    fs_file.write("DY " + str(dy) + "\n")
    fs_file.write("DETA 0.1\n")
    fs_file.write("DRAP 0.07\n")
    fs_file.write("DTAU " + str(tau_s) + "\n")
    fs_file.write("TAU0 0.0\n")
    fs_file.write("EOS_TYPE 1\n")
    fs_file.write("E_FREEZE " + str(e_s))
    fs_file.write("VISCOUS_MATCHING 1")

    fs_file.close()

    #MUSIC file
    music_file = open(outdir + 'music_input_' + str(design_pt),'w')

    music_file.write("echo_level  1\n")                  # control the mount of message output to screen
    music_file.write("mode 2\n")                         # MUSIC running mode
    music_file.write("Initial_profile 42\n")             # type of initial condition
    music_file.write("initialize_with_entropy 0\n")      # init with entropy density or energy density
    music_file.write("s_factor 1.00\n")                  # normalization factor read in
    music_file.write("boost_invariant  1\n")             # whether the simulation is boost-invariant
    music_file.write("Initial_time_tau_0 " +str(tau_s) + "\n")# starting time of the hydro
    music_file.write("Total_evolution_time_tau 30.\n")   # the maximum allowed running time
    music_file.write("Delta_Tau " + str(dtau) + "\n")    # time step to use in the evolution [fm/c]
    music_file.write("Eta_grid_size 1.0\n")              # spatial rapidity range
    music_file.write("Grid_size_in_eta 1\n")             # number of the grid points in spatial
    music_file.write("X_grid_size_in_fm " + str(max_x) + "\n")# spatial range along x direction in the
    music_file.write("Y_grid_size_in_fm " + str(max_y) + "\n")# spatial range along y direction in the
    music_file.write("Grid_size_in_y " + str(nx) + "\n")             # number of the grid points in y direction
    music_file.write("Grid_size_in_x " + str(ny) + "\n")             # number of the grid points in x direction
    music_file.write("EOS_to_use 9\n")                   # type of the equation of state
    music_file.write("reconst_type  1\n")                # 0: solve energy density for hydro eqns. 1: solve flow velocity for hydro eqns.
    music_file.write("Minmod_Theta 1.8\n")               # theta parameter in the min-mod like limiter
    music_file.write("Runge_Kutta_order 2\n")            # order of Runge_Kutta for temporal evolution
    music_file.write("Viscosity_Flag_Yes_1_No_0 1\n")    # turn on viscosity in the evolution
    music_file.write("Include_Shear_Visc_Yes_1_No_0 1\n")# include shear viscous effect

    music_file.write("T_dependent_Shear_to_S_ratio  2\n")# flag to use temperature dep. \eta/s(T)
    music_file.write("T_dependent_Bulk_to_S_ratio  2\n")# flag to use temperature dep. \zeta/s(T)
    music_file.write("eta_over_s_min " + str(eta_over_s_min) + "\n")
    music_file.write("eta_over_s_slope " + str(eta_over_s_slope) + "\n")
    music_file.write("eta_over_s_curv " + str(eta_over_s_curv) + "\n")
    music_file.write("bulk_viscosity_normalisation " + str(bulk_viscosity_normalisation) + "\n")
    music_file.write("bulk_viscosity_width_in_GeV " + str(bulk_viscosity_width_in_GeV) + "\n")
    music_file.write("bulk_viscosity_peak_in_GeV " + str(bulk_viscosity_peak_in_GeV) + "\n")

    music_file.write("Include_Bulk_Visc_Yes_1_No_0 1\n") # include bulk viscous effect
    music_file.write("Include_second_order_terms 1\n")   # include second order non-linear coupling terms
    music_file.write("store_hydro_info_in_memory 0\n")   # flag to store hydro info in memory
    music_file.write("output_evolution_data 0\n")        # flag to output evolution history to file
    music_file.write("Do_FreezeOut_Yes_1_No_0 1\n")      # flag to find freeze-out surface
    music_file.write("Do_FreezeOut_lowtemp 0\n")         # flag to include cold corona
    music_file.write("freeze_out_method 4\n")            # method for hyper-surface finder
    music_file.write("average_surface_over_this_many_time_steps 5\n")   # the step skipped in the tau
    music_file.write("epsilon_freeze " + str(e_s) + "\n")            # the freeze out energy density (GeV/fm^3)
    music_file.write("use_eps_for_freeze_out 0\n")       # 0: use temperature, 1: use energy density
    music_file.write("T_freeze " + str(T_s) + "\n")                 # freeze-out temperature (GeV)
    music_file.write("EndOfData\n")

    music_file.close()

    #the iS3D files, one for each delta_f
    #delta_f_mode = 4 # 1: 14 moment, 2: C.E., 3: McNelis feq_mod, 4: Bernhard feq_mod
    for df_mode in range(1,5):
        iS3D_file = open(outdir + 'iS3D_parameters_' + str(design_pt) + '_df_' + str(df_mode) + '.dat','w')

        iS3D_file.write("operation = 2\n")
        iS3D_file.write("mode      = 6\n")
        iS3D_file.write("hrg_eos   = 3\n")
        iS3D_file.write("dimension = 2\n")
        iS3D_file.write("df_mode   = " + str(df_mode) + "\n")
        iS3D_file.write("include_baryon            	= 0\n")
        iS3D_file.write("include_bulk_deltaf       	= 1\n")
        iS3D_file.write("include_shear_deltaf      	= 1\n")
        iS3D_file.write("include_baryondiff_deltaf 	= 0\n")
        iS3D_file.write("regulate_deltaf           	= 0\n")
        iS3D_file.write("outflow 			= 1\n")
        iS3D_file.write("deta_min 			= 1.e-2\n")
        iS3D_file.write("group_particles            = 0\n")
        iS3D_file.write("particle_diff_tolerance    = 0.01\n")
        iS3D_file.write("mass_pion0		    = 0.138\n")
        iS3D_file.write("do_resonance_decays = 0\n")
        iS3D_file.write("lightest_particle 	 = 111\n")
        iS3D_file.write("oversample		     = 0\n")
        iS3D_file.write("min_num_hadrons     = 1.e+6\n")
        iS3D_file.write("sampler_seed	     = -1\n")
        #these only used for testing, are dummys
        iS3D_file.write("test_sampler = 0\n")
        iS3D_file.write("pT_lower_cut = 0.0\n")
        iS3D_file.write("pT_upper_cut = 3.0\n")
        iS3D_file.write("pT_bins = 100\n")
        iS3D_file.write("y_cut = 5.0\n")
        iS3D_file.write("tau_min = 0.0\n")
        iS3D_file.write("tau_max = 12.0\n")
        iS3D_file.write("tau_bins = 120\n")
        iS3D_file.write("r_min = 0.0\n")
        iS3D_file.write("r_max = 10.0\n")
        iS3D_file.write("r_bins = 50\n")

        iS3D_file.close()


    #the jetscape init xml file
    #note that this file can potentially override parameters set in the MUSIC input file
    js_file = open(outdir + 'jetscape_init_' + str(design_pt) + '.xml', 'w')

    js_file.write("<?xml version=\"1.0\"?>\n")
    js_file.write(" <jetscape>\n")
    js_file.write("  <debug> on </debug>\n")
    js_file.write("  <remark> off </remark>\n")
    js_file.write("  <vlevel> 0 </vlevel>\n")
    js_file.write("  <Random>\n")
    js_file.write("    <seed>0</seed>\n")
    js_file.write("  </Random>\n")

    #parameters common to TRENTo and MUSIC
    js_file.write("  <IS>\n")
    js_file.write("    <grid_max_x> " + str(max_x) + " </grid_max_x>\n")
    js_file.write("    <grid_max_y> " + str(max_y) + " </grid_max_y>\n")
    js_file.write("    <grid_max_z> 0 </grid_max_z>\n")
    js_file.write("    <grid_step_x> " + str(dx) + " </grid_step_x>\n")
    js_file.write("    <grid_step_y> " + str(dy) + " </grid_step_y>\n")
    js_file.write("    <grid_step_z> 0.5 </grid_step_z>\n")

    #TRENTo parameters
    js_file.write("    <Trento>\n")
    js_file.write("		<PhysicsInputs	projectile=\' " + str(projectile) + " \'\n")
    js_file.write("						target=\' " + str(target) + "\'\n")
    js_file.write("						sqrts=\' " + str(sqrts) + "\'\n")
    js_file.write("						cross-section=\' " + str(cross_section) + "\'\n")
    js_file.write("						normalization=\'" + str(normalization) + "\'>\n")
    js_file.write("		</PhysicsInputs>\n")
    js_file.write("		<CutInputs	centrality-low=\'" + str(cent_low) + "\'\n")
    js_file.write("					centrality-high=\'" + str(cent_high) + "\'>\n")
    js_file.write("		</CutInputs>\n")
    js_file.write("		<TransInputs	reduced-thickness=\'" + str(reduced_thickness) + "\'\n")
    js_file.write("						fluctuation=\'" + str(fluctuation) + "\'\n")
    js_file.write("						nucleon-width=\'" + str(nucleon_width) + "\'\n")
    js_file.write("						nucleon-min-dist=\'" + str(nucleon_min_dist) + "\'>\n")
    js_file.write("		</TransInputs>\n")
    js_file.write("		<LongiInputs	mean-coeff=\'1.0\'\n")
    js_file.write("						std-coeff=\'3.0\'\n")
    js_file.write("						skew-coeff=\'0.0\'\n")
    js_file.write("						skew-type=\'1\'\n")
    js_file.write("						jacobian=\'0.8\'>\n")
    js_file.write("		</LongiInputs>\n")
    js_file.write("    </Trento>\n")
    js_file.write("    <initial_profile_path>../examples/test_hydro_files</initial_profile_path>\n")
    js_file.write("  </IS>\n")

    #these are dummies, not actually read by freestream-milne
    js_file.write("  <Preequilibrium>\n")
    js_file.write("    <tau0>0.0</tau0>\n")
    js_file.write("    <taus>0.5</taus>\n")
    js_file.write("    <FreestreamMilne>\n")
    js_file.write("      <name>FreestreamMilne</name>\n")
    js_file.write("      <freestream_input_file>freestream_input</freestream_input_file>\n")
    js_file.write("    </FreestreamMilne>\n")
    js_file.write("  </Preequilibrium>\n")

    #fixed params for MUSIC
    js_file.write("  <Hydro>\n")
    js_file.write("    <MUSIC>\n")
    js_file.write("      <name>MUSIC</name>\n")
    js_file.write("      <MUSIC_input_file>music_input</MUSIC_input_file>\n")
    js_file.write("      <Perform_CooperFrye_Feezeout>0</Perform_CooperFrye_Feezeout>\n")
    js_file.write("    </MUSIC>\n")
    js_file.write("  </Hydro>\n")

    #fixed params for SMASH
    js_file.write("  <Afterburner>\n")
    js_file.write("    <SMASH>\n")
    js_file.write("      <name>SMASH</name>\n")
    js_file.write("      <SMASH_config_file>smash_input/config.yaml</SMASH_config_file>\n")
    js_file.write("      <SMASH_particles_file>smash_input/particles.txt</SMASH_particles_file>\n")
    js_file.write("      <SMASH_decaymodes_file>smash_input/decaymodes.txt</SMASH_decaymodes_file>\n")
    js_file.write("      <end_time>300.0</end_time>\n")
    js_file.write("      <only_decays>0</only_decays>\n")
    js_file.write("    </SMASH>\n")
    js_file.write("  </Afterburner>\n")
    js_file.write("</jetscape>\n")

    js_file.close()
