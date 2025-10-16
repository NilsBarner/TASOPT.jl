using TASOPT
include(TASOPT.__TASOPTindices__)

#---------------------------------     
# hxsize!()
#---------------------------------
HXgas = TASOPT.engine.HX_gas()
HXgeom = TASOPT.engine.HX_tubular()

HXgas.fluid_p = "air"
HXgas.fluid_c = "h2"
HXgas.mdot_p = 1144/60
HXgas.mdot_c = 9.95/60
HXgas.ε = 0.8
HXgas.Tp_in = 778
HXgas.Tc_in = 264
HXgas.Mp_in  = 0.19
HXgas.Mc_in = 0.0285
HXgas.pp_in = 40e3
HXgas.pc_in = 1515e3

HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
HXgas.igas_c = 40

HXgeom.is_concentric = true
HXgeom.has_recirculation = false
HXgeom.has_shaft = false
HXgeom.D_i = 0.564
HXgeom.l = 0.6084530646014857 #tube length
HXgeom.n_stages = 4
HXgeom.xt_D = 6
HXgeom.xl_D = 1.25
HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
HXgeom.material = TASOPT.StructuralAlloy("SS-304")
HXgeom.Δpdes = 3e6

Q = 1.5e6  # pare_sl[iQheat]
println(HXgas.mdot_c)
HXgas.mdot_c = TASOPT.engine.radiator_coolant_mass(HXgas, Q)
println(HXgas.mdot_c)

TASOPT.engine.hxsize!(HXgas, HXgeom)

size_out = [HXgas.Tp_out, HXgas.Tc_out, HXgas.Δp_p, HXgeom.N_t, HXgeom.n_passes, HXgeom.tD_o, HXgeom.A_cs]
# println(size_out)

size_out_check = 
[731.5888605437423, 665.8848846504773, 1414.022586064259, 62.03322510460286, 8.040614281031441, 0.004760508726403918, 1.0189779296746375]
# println(size_out_check)

println("HXgeom.D_i = ", HXgeom.D_i)
println("HXgeom.D_o = ", HXgeom.D_o)
println("HXgeom.l = ", HXgeom.l)