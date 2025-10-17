"""
This script contains plotting methods written by Nils to visualise
the effect of varying FCS weights, locations, and distributions
on aircraft structural weight and drag.
"""

module PointLoadPlots

    export compare_stickfig,
        plot_wrel_vs_segment,
        plot_wrel_vs_sigma_fcs,
        plot_wmto_vs_sigma_fcs,
        plot_weight_drag_vs_weight_split,
        plot_wmto_vs_weight_position,
        plot_drag_vs_sigma_fcs,
        scatter_weight_drag_segment
    using Revise
    using Plots
    using TASOPT
    include(TASOPT.__TASOPTindices__)

    """
        compare_stickfig(aircrafts::Vector{Any}; colour_list::Vector, kwargs...)

    Plot multiple aircraft stick figures on the same figure by overlaying their 2D sketches.
    Any keyword arguments are passed down to `stickfig`.
    """
    function compare_stickfig(aircrafts::Vector{Any}; colour_list::Vector, kwargs...)

        # Start a new plot with the first aircraft
        plt = TASOPT.stickfig(
            aircrafts[1];
            alpha = 0.5, show_seats = false, fill_tank = false, annotate_sketch = false, annotate_length = false, annotate_group = false, airframe_colour = :black, engine_colour = :black,
            kwargs...
        )

        # Overlay the remaining aircraft on the same plot
        for (i, ac) in enumerate(aircrafts[2:end])
            colour = colour_list[i]
            TASOPT.stickfig(
                ac; plot_obj = plt,
                show_seats = false, fill_tank = false, annotate_sketch = false, annotate_length = false, annotate_group = false, airframe_colour = colour, engine_colour = colour,
                kwargs...
            )
        end

        # UNCOMMENT below two lines if want to plot only two aircraft at a time (one case)
        # xlims!(plt, -2.5, 62.5)
        # ylims!(plt, -27.5, 27.5)

        # UNCOMMENT below lines to add a custom legend if want to plot all cases

        # # Create a unique case list with one entry per case
        # case_list = ["Case 1", "Case 21", "Case 22", "Case 31", "Case 32", "Case 4", "Case 5"]

        # # read current axis limits
        # x_min, x_max = xlims(plt)
        # y_min, y_max = ylims(plt)
        # dx = x_max - x_min
        # dy = y_max - y_min

        # # enlarge x-range on the right to make room for text (28% of current width)
        # xr_new = x_max + 0.28 * dx
        # xlims!(plt, x_min, xr_new)

        # # legend anchor x coordinate (inside the newly added right area)
        # x_marker = x_max + 0.04 * dx
        # x_text   = x_max + 0.08 * dx

        # # top y coordinate for first label, and vertical spacing
        # y_top = y_max - 0.02 * dy
        # spacing = 0.06 * dy

        # # draw one line per case (marker + coloured text)
        # for (k, col) in enumerate(colour_list[1:2:end])  # Use every second color to match unique cases
        #     yk = y_top - (k-1) * spacing
        #     # small coloured square marker (no legend label)
        #     scatter!(plt, [x_marker], [yk]; marker=:rect, ms=8, color=col, label="")

        #     # coloured text label; you can change halign/valign/fontsize here
        #     annotate!(plt, (x_text, yk, text(case_list[k], 10, halign=:left, valign=:center, color=col)))
        # end

        return plt

    end

    ###

    """
        plot_wrel_vs_segment(
            ac_refs::Vector{Any}, result_dict::Dict, case_indices::Vector{Int}, N_points::Int,
            sigma_fc_min::Number, sigma_fc_max::Number,
        )
    
    Plot relative weights of fuselage, wing, and tail for regional vs narrowbody aircraft
    as a function of FC system specific power.
    """
    function plot_wrel_vs_segment(
        ac_refs::Vector{Any}, result_dict::Dict, case_indices::Vector{Int}, N_points::Int,
        sigma_fc_min::Number, sigma_fc_max::Number,
    )

        ##

        # Aircraft weights

        # Define colors and linestyles for cases and components
        colors = [:black, :blue, :red, :green, :purple]

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            fuselage_weight_ref = ac_refs[i].fuselage.weight
            wing_weight_ref = ac_refs[i].wing.weight
            htail_weight_ref = ac_refs[i].htail.weight
            vtail_weight_ref = ac_refs[i].vtail.weight
            x_values = range(sigma_fc_max, stop=sigma_fc_min, length=N_points) / 1e3

            ac_list = result_dict[case_idx]["ac"]
            Wpointload_fuselage_list = result_dict[case_idx]["Wpointload_fuselage"]

            plot!(x_values, [ac.fuselage.weight - Wpointload_fuselage_list[j] for (j, ac) in enumerate(ac_list)] / fuselage_weight_ref, color=color, linestyle=:solid, label="", yscale=:log10)
            plot!(x_values, [ac.wing.weight for ac in ac_list] / wing_weight_ref, color=color, linestyle=:dash, label="", yscale=:log10)
            # Plot horizontal and vertical tail weights separately
            # plot!(x_values, [ac.htail.weight for ac in ac_list] / htail_weight_ref, color=color, linestyle=:dot, label="", yscale=:log10)
            # plot!(x_values, [ac.vtail.weight for ac in ac_list] / vtail_weight_ref, color=color, linestyle=:dashdot, label="", yscale=:log10)
            # Plot horizontal and vertical tail weights combined
            plot!(x_values, [ac.htail.weight + ac.vtail.weight for ac in ac_list] / (htail_weight_ref + vtail_weight_ref), color=color, linestyle=:dashdot, label="", yscale=:log10)
        end

        # yticks = ([1, 2, 3, 5, 10], string.([1, 2, 3, 5, 10]))  # choose ticks on case-by-case basis
        yticks = ([1, 1.1, 1.2, 1.3, 1.5], string.([1, 1.1, 1.2, 1.3, 1.5]))  # choose ticks on case-by-case basis
        plot!(yticks=yticks)

        # Custom legend

        x_dummy = [NaN]

        # Cases (colors)
        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]
            plot!(x_dummy, x_dummy, color=color, linestyle=:solid, label="Case $case_idx")
        end

        # Components (linestyles)
        plot!(x_dummy, x_dummy, color=:black, linestyle=:solid, label="Fuselage")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dash, label="Wing")
        # Plot horizontal and vertical tail weights separately
        # plot!(x_dummy, x_dummy, color=:black, linestyle=:dot, label="Horizontal Tail")
        # plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Vertical Tail")
        # Plot horizontal and vertical tail weights combined
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Tail")

        # Final plot settings
        plot!(
            xlabel="FC system specific power (kW/kg)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            # title="Normalised component weights"
            dpi=400
        )
        display(p)
        # savefig("normalised_aircraft_weights_segment.png")

    end

    ###

    """
        plot_wrel_vs_sigma_fcs(
            ac_ref::aircraft, result_dict::Dict, case_indices::Vector{Int}, N_points::Int,
            sigma_fc_min::Number, sigma_fc_max::Number,
        )
    
    Plot relative weights of fuselage, wing, and tail air aircraft, wing, and fuselage level.
    """
    function plot_wrel_vs_sigma_fcs(
        ac_ref::aircraft, result_dict::Dict, case_indices::Vector{Int}, N_points::Int,
        sigma_fc_min::Number, sigma_fc_max::Number,
    )

        # Aircraft weights

        fuselage_weight_ref = ac_ref.fuselage.weight
        wing_weight_ref = ac_ref.wing.weight
        htail_weight_ref = ac_ref.htail.weight
        vtail_weight_ref = ac_ref.vtail.weight
        x_values = range(sigma_fc_max, stop=sigma_fc_min, length=N_points) / 1e3

        # Define colors and linestyles for cases and components
        colors = [:black, :blue, :red, :green, :purple]

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            ac_list = result_dict[case_idx]["ac"]
            Wpointload_fuselage_list = result_dict[case_idx]["Wpointload_fuselage"]

            plot!(x_values, [ac.fuselage.weight - Wpointload_fuselage_list[j] for (j, ac) in enumerate(ac_list)] / fuselage_weight_ref, color=color, linestyle=:solid, label="", yscale=:log10)
            plot!(x_values, [ac.wing.weight for ac in ac_list] / wing_weight_ref, color=color, linestyle=:dash, label="", yscale=:log10)
            # Plot horizontal and vertical tail weights separately
            # plot!(x_values, [ac.htail.weight for ac in ac_list] / htail_weight_ref, color=color, linestyle=:dot, label="", yscale=:log10)
            # plot!(x_values, [ac.vtail.weight for ac in ac_list] / vtail_weight_ref, color=color, linestyle=:dashdot, label="", yscale=:log10)
            # Plot horizontal and vertical tail weights combined
            plot!(x_values, [ac.htail.weight + ac.vtail.weight for ac in ac_list] / (htail_weight_ref + vtail_weight_ref), color=color, linestyle=:dashdot, label="", yscale=:log10)
        end

        # yticks = ([1, 2, 3, 5, 10], string.([1, 2, 3, 5, 10]))  # choose ticks on case-by-case basis
        yticks = ([1, 1.1, 1.2, 1.3, 1.5], string.([1, 1.1, 1.2, 1.3, 1.5]))  # choose ticks on case-by-case basis
        plot!(yticks=yticks)

        # Custom legend

        x_dummy = [NaN]

        # Cases (colors)
        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]
            plot!(x_dummy, x_dummy, color=color, linestyle=:solid, label="Case $case_idx")
        end

        # Components (linestyles)
        plot!(x_dummy, x_dummy, color=:black, linestyle=:solid, label="Fuselage")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dash, label="Wing")
        # Plot horizontal and vertical tail weights separately
        # plot!(x_dummy, x_dummy, color=:black, linestyle=:dot, label="Horizontal Tail")
        # plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Vertical Tail")
        # Plot horizontal and vertical tail weights combined
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Tail")

        # Final plot settings
        plot!(
            xlabel="FC system specific power (kW/kg)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            # title="Normalised component weights"
            dpi=400
        )
        display(p)
        # savefig("normalised_aircraft_weights_pspec.png")

        ##

        # Fuselage weights

        fuselage_shell_weight_ref = ac_ref.fuselage.shell.weight.W
        fuselage_bendingmaterial_h_weight_ref = ac_ref.fuselage.bendingmaterial_h.weight.W
        fuselage_bendingmaterial_v_weight_ref = ac_ref.fuselage.bendingmaterial_v.weight.W
        x_values = range(sigma_fc_max, stop=sigma_fc_min, length=N_points) / 1e3

        # Define colors for cases
        colors = [:black, :green, :blue, :red, :purple]

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            ac_list = result_dict[case_idx]["ac"]
            plot!(x_values, [ac.fuselage.shell.weight.W for ac in ac_list] / fuselage_shell_weight_ref, label="", color=color, linestyle=:solid)
            plot!(x_values, [ac.fuselage.bendingmaterial_h.weight.W for ac in ac_list] / fuselage_bendingmaterial_h_weight_ref, label="", color=color, linestyle=:dash)
            plot!(x_values, [ac.fuselage.bendingmaterial_v.weight.W for ac in ac_list] / fuselage_bendingmaterial_v_weight_ref, label="", color=color, linestyle=:dot)
        end

        # Custom legend

        x_dummy = [NaN]

        # Cases (linestyles)
        plot!(x_dummy, x_dummy, color="black", linestyle=:solid, label="Case 1")
        plot!(x_dummy, x_dummy, color="blue", linestyle=:solid, label="Case 21")
        plot!(x_dummy, x_dummy, color="red", linestyle=:solid, label="Case 22")

        # Components (colors)
        plot!(x_dummy, x_dummy, color="black", linestyle=:solid, label="Shell")
        plot!(x_dummy, x_dummy, color="black", linestyle=:dash, label="Vert. bending material")
        plot!(x_dummy, x_dummy, color="black", linestyle=:dot, label="Horz. bending material")

        # Final plot settings
        plot!(
            xlabel="FC system specific power (kW/kg)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            # title="Normalised Fuselage component weights"
        )
        display(p)
        # savefig("normalised_fuselage_weights_pspec.png")

        ##

        # Wing weights

        wing_center_weight_ref = ac_ref.wing.center.weight
        wing_inboard_weight_ref = ac_ref.wing.inboard.weight
        wing_outboard_weight_ref = ac_ref.wing.outboard.weight
        x_values = range(sigma_fc_max, stop=sigma_fc_min, length=N_points) / 1e3

        # Define colors for cases
        colors = [:black, :blue, :purple]

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            ac_list = result_dict[case_idx]["ac"]
            plot!(x_values, [ac.wing.center.weight for ac in ac_list] / wing_center_weight_ref, label="", color=color, linestyle=:solid)
            plot!(x_values, [ac.wing.inboard.weight for ac in ac_list] / wing_inboard_weight_ref, label="", color=color, linestyle=:dash)
            plot!(x_values, [ac.wing.outboard.weight for ac in ac_list] / wing_outboard_weight_ref, label="", color=color, linestyle=:dot)
        end

        # Custom legend

        x_dummy = [NaN]

        # Cases (colors)
        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]
            plot!(x_dummy, x_dummy, color=color, linestyle=:solid, label="Case $case_idx")
        end

        # Components (linestyles)
        plot!(x_dummy, x_dummy, color=:black, linestyle=:solid, label="Centre box")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dash, label="Inboard section")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dot, label="Outboard section")

        # Final plot settings
        plot!(
            xlabel="FC system specific power (kW/kg)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            # title="Normalised Wing component weights"
        )
        display(p)
        # savefig("normalized_wing_weights_pspec.png")

    end

    ###

    """
        plot_wmto_vs_sigma_fcs(
            ac_ref::aircraft, result_dict::Dict, N_points::Int,
            sigma_fc_min::Number, sigma_fc_max::Number,
        )
    
    Plot relative MTOW for varying FC system specific power.
    """
    function plot_wmto_vs_sigma_fcs(
        ac_ref::aircraft, result_dict::Dict, N_points::Int,
        sigma_fc_min::Number, sigma_fc_max::Number,
    )

        WMTO_ref = ac_ref.parg[igWMTO]
        x_values = range(sigma_fc_max, stop=sigma_fc_min, length=N_points) / 1e3

        p = plot(x_values, [ac.parg[igWMTO] for ac in result_dict[1]["ac"]] / WMTO_ref, color="black", linestyle=:solid, label="Case 1")
        plot!(x_values, [ac.parg[igWMTO] for ac in result_dict[21]["ac"]] / WMTO_ref, color="blue", linestyle=:solid, label="Case 22")
        plot!(x_values, [ac.parg[igWMTO] for ac in result_dict[22]["ac"]] / WMTO_ref, color="green", linestyle=:solid, label="Case 23")

        plot!(
            xlabel="FC system specific power (kW/kg)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            title="Normalised take-off weights"
        )
        display(p)
        # savefig("normalised_mto_weights_pspec.png")

    end

    ###

    """
        plot_weight_drag_vs_weight_split(
            ac_ref::aircraft, result_dict::Dict, N_points::Int,
            Weng_single_ref::Number, sigma_gt::Number, sigma_fc_wing_min::Number, sigma_fc_wing_max::Number,
        )
    
    Plot relative MTOW and component weights and drags for varying FC weight split between wing and nacelle.
    """
    function plot_weight_drag_vs_weight_split(
        ac_ref::aircraft, result_dict::Dict, N_points::Int,
        Weng_single_ref::Number, sigma_gt::Number, sigma_fc_wing_min::Number, sigma_fc_wing_max::Number,
    )

        # MTOWs
        
        WMTO_ref = ac_ref.parg[igWMTO]
        W_fc_wing_min = Weng_single_ref * sigma_gt / sigma_fc_wing_max
        W_fc_wing_max = Weng_single_ref * sigma_gt / sigma_fc_wing_min
        x_values = range(W_fc_wing_min, stop=W_fc_wing_max, length=N_points) / W_fc_wing_max

        p = plot(x_values, [ac.parg[igWMTO] for ac in result_dict[31]["ac"]] / WMTO_ref, color="black", linestyle=:solid, label="Case 31")
        plot!(x_values, [ac.parg[igWMTO] for ac in result_dict[32]["ac"]] / WMTO_ref, color="blue", linestyle=:solid, label="Case 32")

        plot!(
            xlabel="Fraction of weight in nacelle (-)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            title="Normalised take-off weights"
        )
        display(p)
        # savefig("normalised_mto_weights_split.png")

        ##

        # Aircraft weights

        fuselage_weight_ref = ac_ref.fuselage.weight
        wing_weight_ref = ac_ref.wing.weight
        htail_weight_ref = ac_ref.htail.weight
        vtail_weight_ref = ac_ref.vtail.weight

        # Define colors and linestyles for cases and components
        colors = [:black, :blue, :red, :green, :purple]
        case_indices = [31, 32]

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            ac_list = result_dict[case_idx]["ac"]
            Wpointload_fuselage_list = result_dict[case_idx]["Wpointload_fuselage"]

            plot!(x_values, [ac.fuselage.weight - Wpointload_fuselage_list[j] for (j, ac) in enumerate(ac_list)] / fuselage_weight_ref, color=color, linestyle=:solid, label="", yscale=:log10)
            plot!(x_values, [ac.wing.weight for ac in ac_list] / wing_weight_ref, color=color, linestyle=:dash, label="", yscale=:log10)
            plot!(x_values, [ac.htail.weight + ac.vtail.weight for ac in ac_list] / (htail_weight_ref + vtail_weight_ref), color=color, linestyle=:dashdot, label="", yscale=:log10)
        end

        yticks = ([1, 2, 3, 5, 10], string.([1, 2, 3, 5, 10]))
        plot!(yticks=yticks)

        # Custom legend

        x_dummy = [NaN]

        # Cases (colors)
        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]
            plot!(x_dummy, x_dummy, color=color, linestyle=:solid, label="Case $case_idx")
        end

        # Components (linestyles)
        plot!(x_dummy, x_dummy, color=:black, linestyle=:solid, label="Fuselage")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dash, label="Wing")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Tail")

        # Final plot settings
        plot!(
            xlabel="Fraction of weight in nacelle (-)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            # title="Normalised component weights"
            dpi=400
        )
        display(p)
        # savefig("normalised_aircraft_weights_split.png")

        ##

        # Aircraft drags

        # Define colors and linestyles for cases and components
        colors = [:black, :blue, :red, :green, :purple]

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            CLSref = ac_ref.wing.layout.S * ac_ref.para[iaCL, ipcruise1]
            CDSref = ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1]
            CDfuseSref = ac_ref.wing.layout.S * ac_ref.para[iaCDfuse, ipcruise1]
            CDiSref = ac_ref.wing.layout.S * ac_ref.para[iaCDi, ipcruise1]
            CDwingSref = ac_ref.wing.layout.S * ac_ref.para[iaCDwing, ipcruise1]
            CDhtailSref = ac_ref.wing.layout.S * ac_ref.para[iaCDhtail, ipcruise1]
            CDvtailSref = ac_ref.wing.layout.S * ac_ref.para[iaCDvtail, ipcruise1]

            ac_list = result_dict[case_idx]["ac"]
            
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDfuse, ipcruise1] for ac in ac_list] / CDfuseSref, color=color, linestyle=:solid, label="", yscale=:log10)
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDwing, ipcruise1] for ac in ac_list] / CDwingSref, color=color, linestyle=:dash, label="", yscale=:log10)
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDhtail, ipcruise1] + ac.wing.layout.S * ac.para[iaCDvtail, ipcruise1] for ac in ac_list] / CDhtailSref, color=color, linestyle=:dashdot, label="", yscale=:log10)

        end

        yticks = ([1, 2, 3, 5, 10], string.([1, 2, 3, 5, 10]))  # choose ticks on case-by-case basis
        # yticks = ([1, 1.1, 1.2, 1.3, 1.5], string.([1, 1.1, 1.2, 1.3, 1.5]))  # choose ticks on case-by-case basis
        plot!(yticks=yticks)

        # Custom legend

        x_dummy = [NaN]

        # Cases (colors)
        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]
            plot!(x_dummy, x_dummy, color=color, linestyle=:solid, label="Case $case_idx")
        end

        # Components (linestyles)
        plot!(x_dummy, x_dummy, color=:black, linestyle=:solid, label="Fuselage")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dash, label="Wing")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Tail")

        # Final plot settings
        plot!(
            xlabel="Fraction of weight in nacelle (-)",
            ylabel="Normalised drag (-)",
            legend=:outerright,
            # title="Normalised component drags"
            dpi=400
        )
        display(p)
        # savefig("normalised_aircraft_drags_split.png")

    end

    ###

    """
        plot_wmto_vs_weight_position(
            ac_ref::aircraft, result_dict::Dict, N_points_4::Int, N_points_5::Int,
        )
    
    Plot relative MTOW and component weights and drags for varying FC weight position along fuselage and number of engines.
    """
    function plot_wmto_vs_weight_position(
        ac_ref::aircraft, result_dict::Dict, N_points_4::Int, N_points_5::Int,
    )
        # MTOWs

        WMTO_ref = ac_ref.parg[igWMTO]
        x_values_4 = range(0, stop=1, length=N_points_4)
        x_values_5 = range(0, stop=1, length=N_points_5)

        # p = plot(x_values_4, [ac.parg[igWMTO] for ac in result_dict[4]["ac"]] / WMTO_ref, color="black", linestyle=:solid, label="Case 4, f(% pressure shell)")
        p = plot!(x_values_5, [ac.parg[igWMTO] for ac in result_dict[5]["ac"]] / WMTO_ref, color="blue", linestyle=:solid, label="Case 5, f(# engines)")

        plot!(
            # xlabel="Position of weight along cylinder (-)",
            xlabel="Number of engines (-)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            title="Normalised take-off weights"
        )
        display(p)
        # savefig("normalised_mto_weights.png")

        ##

        # Aircraft weights

        fuselage_weight_ref = ac_ref.fuselage.weight
        wing_weight_ref = ac_ref.wing.weight
        htail_weight_ref = ac_ref.htail.weight
        vtail_weight_ref = ac_ref.vtail.weight

        # Define colors and linestyles for cases and components
        colors = [:black, :blue, :red, :green, :purple]
        case_indices = [5]
        x_values = x_values_5

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            ac_list = result_dict[case_idx]["ac"]
            Wpointload_fuselage_list = result_dict[case_idx]["Wpointload_fuselage"]

            plot!(x_values, [ac.fuselage.weight - Wpointload_fuselage_list[j] for (j, ac) in enumerate(ac_list)] / fuselage_weight_ref, color=color, linestyle=:solid, label="", yscale=:log10)
            plot!(x_values, [ac.wing.weight for ac in ac_list] / wing_weight_ref, color=color, linestyle=:dash, label="", yscale=:log10)
            plot!(x_values, [ac.htail.weight + ac.vtail.weight for ac in ac_list] / (htail_weight_ref + vtail_weight_ref), color=color, linestyle=:dashdot, label="", yscale=:log10)
        end

        yticks = ([1, 1.1, 1.2, 1.3, 1.5], string.([1, 1.1, 1.2, 1.3, 1.5]))
        plot!(yticks=yticks)

        # Custom legend

        x_dummy = [NaN]

        # Cases (colors)
        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]
            plot!(x_dummy, x_dummy, color=color, linestyle=:solid, label="Case $case_idx")
        end

        # Components (linestyles)
        plot!(x_dummy, x_dummy, color=:black, linestyle=:solid, label="Fuselage")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dash, label="Wing")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Tail")

        # Final plot settings
        plot!(
            # xlabel="Position of weight along cylinder (-)",
            xlabel="Number of engines (-)",
            ylabel="Normalised weight (-)",
            legend=:outerright,
            # title="Normalised component weights"
            dpi=400
        )
        display(p)
        savefig("normalised_aircraft_weights.png")

        ##

        # Aircraft drags

        # Define colors and linestyles for cases and components
        colors = [:black, :blue, :red, :green, :purple]

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            CLSref = ac_ref.wing.layout.S * ac_ref.para[iaCL, ipcruise1]
            CDSref = ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1]
            CDfuseSref = ac_ref.wing.layout.S * ac_ref.para[iaCDfuse, ipcruise1]
            CDiSref = ac_ref.wing.layout.S * ac_ref.para[iaCDi, ipcruise1]
            CDwingSref = ac_ref.wing.layout.S * ac_ref.para[iaCDwing, ipcruise1]
            CDhtailSref = ac_ref.wing.layout.S * ac_ref.para[iaCDhtail, ipcruise1]
            CDvtailSref = ac_ref.wing.layout.S * ac_ref.para[iaCDvtail, ipcruise1]

            ac_list = result_dict[case_idx]["ac"]
            
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDfuse, ipcruise1] for ac in ac_list] / CDfuseSref, color=color, linestyle=:solid, label="", yscale=:log10)
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDwing, ipcruise1] for ac in ac_list] / CDwingSref, color=color, linestyle=:dash, label="", yscale=:log10)
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDhtail, ipcruise1] + ac.wing.layout.S * ac.para[iaCDvtail, ipcruise1] for ac in ac_list] / CDhtailSref, color=color, linestyle=:dashdot, label="", yscale=:log10)

        end

        # yticks = ([1, 2, 3, 5, 10], string.([1, 2, 3, 5, 10]))  # choose ticks on case-by-case basis
        yticks = ([1, 1.1, 1.2, 1.3, 1.5], string.([1, 1.1, 1.2, 1.3, 1.5]))  # choose ticks on case-by-case basis
        plot!(yticks=yticks)

        # Custom legend

        x_dummy = [NaN]

        # Cases (colors)
        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]
            plot!(x_dummy, x_dummy, color=color, linestyle=:solid, label="Case $case_idx")
        end

        # Components (linestyles)
        plot!(x_dummy, x_dummy, color=:black, linestyle=:solid, label="Fuselage")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dash, label="Wing")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Tail")

        # Final plot settings
        plot!(
            # xlabel="Position of weight along cylinder (-)",
            xlabel="Number of engines (-)",
            ylabel="Normalised drag (-)",
            legend=:outerright,
            # title="Normalised component drags"
            dpi=400
        )
        display(p)
        savefig("normalised_aircraft_drags.png")

    end

    ###

    """
        plot_drag_vs_sigma_fcs(
            ac_refs::Vector{Any}, result_dict::Dict, case_indices::Vector{Int}, N_points::Int,
            sigma_fc_min::Number, sigma_fc_max::Number,
        )
    
    Plot relative drags of fuselage, wing, and tail for regional and narrowbody aircraft.
    """
    function plot_drag_vs_sigma_fcs(
        ac_refs::Vector{Any}, result_dict::Dict, case_indices::Vector{Int}, N_points::Int,
        sigma_fc_min::Number, sigma_fc_max::Number,
    )

        # Aircraft drags

        # Define colors and linestyles for cases and components
        colors = [:black, :blue, :red, :green, :purple]

        p = plot()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            CLSref = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCL, ipcruise1]
            CDSref = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCD, ipcruise1]
            CDfuseSref = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCDfuse, ipcruise1]
            CDiSref = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCDi, ipcruise1]
            CDwingSref = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCDwing, ipcruise1]
            CDhtailSref = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCDhtail, ipcruise1]
            CDvtailSref = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCDvtail, ipcruise1]

            x_values = range(sigma_fc_max, stop=sigma_fc_min, length=N_points) / 1e3

            ac_list = result_dict[case_idx]["ac"]
            
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDfuse, ipcruise1] for ac in ac_list] / CDfuseSref, color=color, linestyle=:solid, label="", yscale=:log10)
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDwing, ipcruise1] for ac in ac_list] / CDwingSref, color=color, linestyle=:dash, label="", yscale=:log10)
            plot!(x_values, [ac.wing.layout.S * ac.para[iaCDhtail, ipcruise1] + ac.wing.layout.S * ac.para[iaCDvtail, ipcruise1] for ac in ac_list] / CDhtailSref, color=color, linestyle=:dashdot, label="", yscale=:log10)

        end

        # yticks = ([1, 2, 3, 5, 10], string.([1, 2, 3, 5, 10]))  # choose ticks on case-by-case basis
        yticks = ([1, 1.1, 1.2, 1.3, 1.5, 2.0], string.([1, 1.1, 1.2, 1.3, 1.5, 2.0]))  # choose ticks on case-by-case basis
        plot!(yticks=yticks)

        # Custom legend

        x_dummy = [NaN]

        # Cases (colors)
        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]
            plot!(x_dummy, x_dummy, color=color, linestyle=:solid, label="Case $case_idx")
        end

        # Components (linestyles)
        plot!(x_dummy, x_dummy, color=:black, linestyle=:solid, label="Fuselage")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dash, label="Wing")
        plot!(x_dummy, x_dummy, color=:black, linestyle=:dashdot, label="Tail")

        # Final plot settings
        plot!(
            xlabel="FC system specific power (kW/kg)",
            ylabel="Normalised drag (-)",
            legend=:outerright,
            # title="Normalised component drags"
            dpi=400
        )
        display(p)
        # savefig("normalised_aircraft_drags_pspec.png")

    end

    ###

    """
        scatter_weight_drag_segment(
            ac_refs::Vector{Any}, result_dict::Dict, case_indices::Vector{Int},
        )
    
    Plot overall weights and drags of regional and narrowbody aircraft,
    both normalised to their respective references, and relative to each other.
    """
    function scatter_weight_drag_segment(
        ac_refs::Vector{Any}, result_dict::Dict, case_indices::Vector{Int},
    )

        # Regional and narrowbody vs respective reference

        # Define colors and linestyles for cases and components
        colors = vcat(fill(:black, length(case_indices) ÷ 2), fill(:red, length(case_indices) ÷ 2))

        p = scatter()

        for (i, case_idx) in enumerate(case_indices)
            color = colors[i % length(colors) + 1]

            WMTO_ref = ac_refs[i].parg[igWMTO]
            CDSref = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCD, ipcruise1]

            ac_list = result_dict[case_idx]["ac"]

            scatter!(
                [ac.parg[igWMTO] for ac in ac_list] / WMTO_ref,
                [ac.wing.layout.S * ac.para[iaCD, ipcruise1] for ac in ac_list] / CDSref,
                color=color, marker=:circle, label="", markerstrokecolor=color, markerstrokewidth=0.5
            )
        end

        # Custom legend

        x_dummy = [NaN]   # invisible point

        # Components (linestyles)
        scatter!(x_dummy, x_dummy, color=:black, markerstrokecolor=:black, markerstrokewidth=0.5, label="Narrowbody")
        scatter!(x_dummy, x_dummy, color=:red, markerstrokecolor=:red, markerstrokewidth=0.5, label="Regional")

        # Set equal limits for axes
        x_min, x_max = xlims(p)
        y_min, y_max = ylims(p)
        xlims!(p, min(x_min, y_min), max(x_max, y_max))
        ylims!(p, xlims(p))

        # Set equal aspect ratio for axes
        scatter!(aspect_ratio=:equal)

        # Final plot settings
        scatter!(
            xlabel="Normalised weight (-)",
            ylabel="Normalised drag (-)",
            legend=:outerright,
            title="Regional and narrowbody vs respective reference",
            dpi=400
        )
        display(p)
        # savefig("normalised_scatter_regional_and_narrow.png")

        ##

        # Regional relative to narrowbody

        # Define colors and linestyles for cases and components
        colors = fill(:black, length(case_indices) ÷ 2)

        p = scatter()

        for (i, case_idx_i) in enumerate(case_indices[1: length(case_indices) ÷ 2])
            color = colors[i % length(colors) + 1]

            j = i + length(case_indices) ÷ 2
            case_idx_j = case_indices[j]

            WMTO_ref_narrowbody = ac_refs[i].parg[igWMTO]
            WMTO_ref_regional = ac_refs[j].parg[igWMTO]
            CDSref_narrowbody = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCD, ipcruise1]
            CDSref_regional = ac_refs[j].wing.layout.S * ac_refs[j].para[iaCD, ipcruise1]

            ac_list_narrowbody = result_dict[case_idx_i]["ac"]
            ac_list_regional = result_dict[case_idx_j]["ac"]

            scatter!(
                [ac.parg[igWMTO] for ac in ac_list_regional] / [ac.parg[igWMTO] for ac in ac_list_narrowbody],
                [ac.wing.layout.S * ac.para[iaCD, ipcruise1] for ac in ac_list_regional] / [ac.wing.layout.S * ac.para[iaCD, ipcruise1] for ac in ac_list_narrowbody],
                color=color, marker=:circle, label="", markerstrokecolor=color, markerstrokewidth=0.5
            )
        end

        # Set equal limits for axes
        x_min, x_max = xlims(p)
        y_min, y_max = ylims(p)
        xlims!(p, min(x_min, y_min), max(x_max, y_max))
        ylims!(p, xlims(p))

        # Set equal aspect ratio for axes
        scatter!(aspect_ratio=:equal)

        # Final plot settings
        scatter!(
            xlabel="Normalised weight (-)",
            ylabel="Normalised drag (-)",
            title="Regional relative to narrowbody",
            dpi=400
        )
        display(p)
        # savefig("scatter_regional_rel_narrow.png")

        ##

        # Relative regional relative to relative narrowbody

        # Define colors and linestyles for cases and components
        colors = fill(:black, length(case_indices) ÷ 2)

        p = scatter()

        for (i, case_idx_i) in enumerate(case_indices[1: length(case_indices) ÷ 2])
            color = colors[i % length(colors) + 1]

            j = i + length(case_indices) ÷ 2
            case_idx_j = case_indices[j]

            WMTO_ref_narrowbody = ac_refs[i].parg[igWMTO]
            WMTO_ref_regional = ac_refs[j].parg[igWMTO]
            CDSref_narrowbody = ac_refs[i].wing.layout.S * ac_refs[i].para[iaCD, ipcruise1]
            CDSref_regional = ac_refs[j].wing.layout.S * ac_refs[j].para[iaCD, ipcruise1]

            ac_list_narrowbody = result_dict[case_idx_i]["ac"]
            ac_list_regional = result_dict[case_idx_j]["ac"]

            scatter!(
                [ac.parg[igWMTO] / WMTO_ref_regional for ac in ac_list_regional] / [ac.parg[igWMTO] / WMTO_ref_narrowbody for ac in ac_list_narrowbody],
                [ac.wing.layout.S * ac.para[iaCD, ipcruise1] / CDSref_regional for ac in ac_list_regional] / [ac.wing.layout.S * ac.para[iaCD, ipcruise1] / CDSref_narrowbody for ac in ac_list_narrowbody],
                color=color, marker=:circle, label="", markerstrokecolor=color, markerstrokewidth=0.5
            )
        end

        # Set equal limits for axes
        x_min, x_max = xlims(p)
        y_min, y_max = ylims(p)
        xlims!(p, min(x_min, y_min), max(x_max, y_max))
        ylims!(p, xlims(p))

        # Set equal aspect ratio for axes
        scatter!(aspect_ratio=:equal)

        # Final plot settings
        scatter!(
            xlabel="Normalised weight (-)",
            ylabel="Normalised drag (-)",
            title="Relative regional relative to relative narrowbody",
            dpi=400
        )
        display(p)
        # savefig("scatter_rel_regional_rel_rel_narrow.png")

    end

end # module PointLoadPlots

###

# Stacked area plots (not used currently)

# stacked = [fuselage_weight_list wing_weight_list htail_weight_list vtail_weight_list]
# cumulative = cumsum(stacked, dims=2)
# p = plot(x_values, cumulative[:,1],
#     xlabel=xlabel, ylabel="Structural weight (N)",
#     label="Fuselage", lw=0, fill_between=(0, cumulative[:,1]))
# plot!(x_values, cumulative[:,2], label="Wing", lw=0, fill_between=(cumulative[:,1], cumulative[:,2]))
# plot!(x_values, cumulative[:,3], label="Horizontal tail", lw=0, fill_between=(cumulative[:,2], cumulative[:,3]))
# plot!(x_values, cumulative[:,4], label="Vertical tail", lw=0, fill_between=(cumulative[:,3], cumulative[:,4]))

# stacked = [fuselage_shell_weight_list fuselage_bendingmaterial_h_weight_list fuselage_bendingmaterial_v_weight_list]
# cumulative = cumsum(stacked, dims=2)
# p = plot(x_values, cumulative[:,1],
#     xlabel=xlabel, ylabel="Structural weight (N)",
#     label="Shell", lw=0, fill_between=(0, cumulative[:,1]))
# plot!(x_values, cumulative[:,2], label="Horizontal bending material", lw=0, fill_between=(cumulative[:,1], cumulative[:,2]))
# plot!(x_values, cumulative[:,3], label="Vertical bending material", lw=0, fill_between=(cumulative[:,2], cumulative[:,3]))

# stacked = [wing_center_weight_list wing_inboard_weight_list wing_outboard_weight_list]
# cumulative = cumsum(stacked, dims=2)
# p = plot(x_values, cumulative[:,1],
#     xlabel=xlabel, ylabel="Structural weight(N)",
#     label="Centre box", lw=0, fill_between=(0, cumulative[:,1]))
# plot!(x_values, cumulative[:,2], label="Inboard section", lw=0, fill_between=(cumulative[:,1], cumulative[:,2]))
# plot!(x_values, cumulative[:,3], label="Outboard section", lw=0, fill_between=(cumulative[:,2], cumulative[:,3]))

# p = TASOPT.MomentShear(ac)
# p = TASOPT.plot_details(ac)
