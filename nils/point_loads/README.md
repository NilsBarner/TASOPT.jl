# Nils' notes on implementation of fuselage- and wing-housed fuel cell
# propulsion system structural analysis

As of 03.09.2025, I implemented the ability to add a vector point load
based on PointLoad mutable structure from /src/structures/loads.jl of
arbitrary component magnitudes and location anywhere inside the fuselage
(for now I have placed it at the location of the APU, so rear end). This
should simulate a (partially) fuselage-housed fuel cell propulsion system.

To model a wing-mounted fuel cell system, I have decided to go with a
slighty different approach, namely to add a `custom_weight_delta` term
inside `tfweightwrap!()` to the value returned by the TASOPT.jl engine
weight estimation method. The reason for that is that a separate point load
requires lots of additional code, e.g. in `balance_aircraft!` of balance.jl,
to correctly account for moments etc., so I deemed it easier to just
modify the weight of the engine itself.

It would obviously be desirable, for consistency, to also treat the fuselage-
housed fuel cell propulsion system in this way, however, when I tried setting
`engine_location = "fuselage"` in e.g. default_input.toml, this did not
appear to be fully implemented (e.g. the stickfig plot still had the engine
under the wing), so I decided to add a separate point load for now.