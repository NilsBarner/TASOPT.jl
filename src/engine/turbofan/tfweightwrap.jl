"""
    tfweightwrap!(ac)

General function to estimate and store the weight of a turbofan engine. 
This function is basically a wrapper on tfweight, going from the
basic aircraft inputs to those required by the function and storing the outputs.
      
!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine weights.
"""
function tfweightwrap!(ac)
    parg = ac.parg
    wing = ac.wing
    neng = parg[igneng]
    
    Weng, Wnace, Webare, W_HXs, Snace1 = tfweight(ac)
        
    custom_weight_delta = ac.engine.model.custom_weight_delta  # line added by NILS
    parg[igWeng] = Weng + custom_weight_delta  # second term added by NILS
    # println("custom_weight_delta INTERNAL =", custom_weight_delta)
    parg[igWebare] = Webare
    parg[igWnace] = Wnace
    parg[igWHXs] = W_HXs

    # set new nacelle area / reference area  fraction fSnace
    S = wing.layout.S

    # # NILS: correct nacelle volume (and thus surface area) by FCS volume
    # dfan = ac.parg[igdfan]  # NILS
    # rSnace = parg[igrSnace]  # NILS
    # Vfcs = parg[igVfcsnac]  # NILS
    # lfcs = Vfcs / (pi * dfan^2)  # NILS

    # Snace = Snace1 * (rSnace + lfcs) / rSnace * neng  # second term added by NILS
    Snace = Snace1 * neng  # second term added by NILS
    fSnace = Snace / S
    parg[igfSnace] = fSnace

    # set new nacelle area / reference area  fraction fSnace
    # Snace = Snace1 * (rSnace + lfcs) / rSnace * neng  # second term added by NILS
    Snace = Snace1 * neng  # second term added by NILS
    fSnace = Snace / S
    parg[igfSnace] = fSnace
    lnace = parg[igdfan] * parg[igrSnace] * 0.15# + lfcs  # last term added by NILS
    parg[iglnace] = lnace

    # NILS: calculate and log total available nacelle volume for FCS
    dfan = ac.parg[igdfan]  # NILS
    HTR_f = parg[igHTRf]
    V_fcs_av_nace = 1/3 * pi * lnace * ((dfan/2 * HTR_f)^2 + (dfan/4 * HTR_f)^2 + (dfan/2 * HTR_f) * (dfan/4 * HTR_f))  # volume of truncated cone 
    # println("lnace, dfan, HTR_f, V_fcs_av_nace = ", lnace, ", ", dfan, ", ", HTR_f, ", ", V_fcs_av_nace)
    parg[igVfcsavnacetot] = neng * V_fcs_av_nace

end