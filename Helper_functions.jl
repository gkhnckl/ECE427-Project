############## FUNCTIONS TO ADD NETWORK ELEMENTS #################

# Function to add a generator to the network
function add_gen!(
    network,
    gen_bus::Int,
    pg::Float64,
    qmax::Float64, qmin::Float64,
    pmax::Float64=0.0, pmin::Float64=0.0;
    qg::Float64=0.0,
    gen_status::Int=1,
    mu_pmin::Float64=0.0, mu_pmax::Float64=0.0, mu_qmin::Float64=0.0, mu_qmax::Float64=0.0,
    pc1::Float64=0.0, pc2::Float64=0.0,
    qc1min::Float64=0.0, qc1max::Float64=0.0, qc2min::Float64=0.0, qc2max::Float64=0.0,
    ramp_agc::Float64=0.0, ramp_10::Float64=0.0, ramp_30::Float64=0.0, ramp_q::Float64=0.0,
    apf::Float64=0.0,
    model::Int=1,
    startup::Float64=0.0, shutdown::Float64=0.0,
    ncost::Int=0, cost::Array{Float64}=Float64[],
    xcoord::Float64=NaN, ycoord::Float64=NaN
)
# Find the highest index among the generators already in the network so that the new one is added in order
gen_indices = parse.(Int64, keys(network["gen"]))
new_gen_index = maximum(gen_indices) + 1

# Set the voltage setpoint to match the bus voltage
vg = network["bus"]["$gen_bus"]["vm"]

if pmax < pg
    pmax = pg
end

# Define a dictionary with all the values PowerModels expects for a generator
new_gen_dict = Dict{String, Any}(
    "index"      => new_gen_index,
    "gen_bus"    => gen_bus,
    "pg"         => pg,
    "qmax"       => qmax,
    "qmin"       => qmin,
    "qg"         => qg,
    "pmax"       => pmax,
    "pmin"       => pmin,
    "vg"         => vg,
    "mbase"      => network["baseMVA"],
    "source_id"  => Any["gen", new_gen_index],
    "gen_status" => gen_status,
    "mu_pmin"    => mu_pmin,
    "mu_pmax"    => mu_pmax,
    "mu_qmin"    => mu_qmin,
    "mu_qmax"    => mu_qmax,
    "pc1"        => pc1,
    "pc2"        => pc2,
    "qc1min"     => qc1min,
    "qc2min"     => qc2min,
    "qc1max"     => qc1max,
    "qc2max"     => qc2max,
    "ramp_agc"   => ramp_agc,
    "ramp_q"     => ramp_q,
    "ramp_10"    => ramp_10,
    "ramp_30"    => ramp_30,
    "apf"        => apf,
    "ncost"      => ncost,
    "model"      => model,
    "shutdown"   => shutdown,
    "startup"    => startup,
    "cost"       => cost,
    "xcoord_1"   => xcoord,
    "ycoord_1"   => ycoord
)

# Add this dictionary to the network data under the "gen" key
network["gen"]["$new_gen_index"] = new_gen_dict

println("Added a generator with index $new_gen_index at bus $gen_bus.")
end


# Function to add a branch to the network
function add_branch!(
    network,
    f_bus::Int, t_bus::Int,
    br_r::Float64, br_x::Float64,
    shunt_b::Float64,
    rate_a::Float64,
    angmin::Float64, angmax::Float64,
    transformer::Bool;
    br_status::Int=1,
    tap::Float64=1.0,
    shift::Float64=0.0,
    rate_b::Float64=-1.0, rate_c::Float64=-1.0,
    pf::Float64=0.0, pt::Float64=0.0,
    qf::Float64=0.0, qt::Float64=0.0,
    mu_angmin::Float64=0.0, mu_angmax::Float64=0.0,
    mu_sf::Float64=0.0, mu_st::Float64=0.0
)
# Find the highest index among the branches already in the network so that the new one is added in order
branch_indices = parse.(Int64, keys(network["branch"]))
new_branch_index = maximum(branch_indices) + 1

# Set the emergency thermal limits as the same as the normal limits if they were not set
if rate_b < 0
    rate_b = rate_a
end
if rate_c < 0
    rate_c = rate_a
end

b_fr = shunt_b/2
b_to = shunt_b/2

g_fr = 0.0
g_to = 0.0

# Define a dictionary with all the values PowerModels expects for a branch
new_branch_dict = Dict{String, Any}(
    "index"=>new_branch_index,
    "f_bus"=>f_bus,
    "t_bus"=>t_bus,
    "source_id"=>Any["branch", new_branch_index],
    "br_status"=>br_status,
    "br_r"=>br_r,
    "br_x"=>br_x,
    "b_fr"=>b_fr,
    "g_fr"=>g_fr,
    "b_to"=>b_to,
    "g_to"=>g_to,
    "rate_a"=>rate_a,
    "rate_b"=>rate_b,
    "rate_c"=>rate_c,
    "angmin"=>angmin,
    "angmax"=>angmax,
    "transformer"=>transformer,
    "tap"=>tap,
    "shift"=>shift,
    "pf"=>pf,
    "pt"=>pt,
    "qf"=>qf,
    "qt"=>qt,
    "mu_angmin"=>mu_angmin,
    "mu_angmax"=>mu_angmax,
    "mu_sf"=>mu_sf,
    "mu_st"=>mu_st
)

# Add this dictionary to the network data under the "branch" key
network["branch"]["$new_branch_index"] = new_branch_dict

println("Added a branch with index $new_branch_index from bus $f_bus to $t_bus.")
end

# Function to add a bus to the network
function add_bus!(
    network,
    bus_type::Int,
    vbase_kv::Float64,
    vmin::Float64, vmax::Float64,
    xcoord::Float64, ycoord::Float64;
    bus_name::String="",
    area::Int=1, zone::Int=1,
    va::Float64=0.0, vm::Float64=1.0,
    mu_vmax::Float64=0.0, mu_vmin::Float64=0.0,
    lam_p::Float64=0.0, lam_q::Float64=0.0,
    index::Int64=-1
)
# If no index is given, find the highest index among the buses already in the network
# so that the new one is added in order
if index < 0
    bus_indices = parse.(Int64, keys(network["bus"]))
    new_bus_index = maximum(bus_indices) + 1
else
    new_bus_index = index
end

# Define a dictionary with all the values PowerModels expects for a bus
new_bus_dict = Dict{String, Any}(
    "index"=>new_bus_index,
    "bus_i"=>new_bus_index,
    "bus_type"=>bus_type,
    "base_kv"=>vbase_kv,
    "vmin"=>vmin,
    "vmax"=>vmax,
    "area"=>area,
    "zone"=>zone,
    "name"=>bus_name,
    "source_id"=>Any["bus", new_bus_index],
    "xcoord_1"=>xcoord,
    "ycoord_1"=>ycoord,
    "va"=>va,
    "vm"=>vm,
    "mu_vmax"=>mu_vmax,
    "mu_vmin"=>mu_vmin,
    "lam_p"=>lam_p,
    "lam_q"=>lam_q
)

# Add this dictionary to the network data under the "bus" key
network["bus"]["$new_bus_index"] = new_bus_dict

println("Added a bus with index $new_bus_index.")
return new_bus_index
end


# Function to add a load to the network
function add_load!(
    network,
    bus::Int,
    pd::Float64, qd::Float64;
    status::Int=1
)
# Find the highest index among the loads already in the network so that the new one is added in order
load_indices = parse.(Int64, keys(network["load"]))
new_load_index = maximum(load_indices) + 1

# Define a dictionary with all the values PowerModels expects for a load
new_load_dict = Dict{String, Any}(
    "index"=>new_load_index,
    "load_bus"=>bus,
    "source_id"=>Any["bus", bus],
    "pd"=>pd,
    "qd"=>qd,
    "status"=>status
)

# Add this dictionary to the network data under the "load" key
network["load"]["$new_load_index"] = new_load_dict

println("Added a load with index $new_load_index at bus $bus.")
end

# Function to scale loads in the network
function scale_loads!(
    network,
    scale_factor::Float64;
    load=[],
    bus=[]
)

loads_scaled = []
bus_loads_scaled = []
# Scale the loads with the given indices, or if bus indices are given, at those buses.
# If no load of bus arguments are given, scale all loads in the system.
if length(load) > 0
    for (idx, l) in network["load"]
        if parse(Int,idx) in load
            l["pd"] *= scale_factor
            l["qd"] *= scale_factor
            push!(loads_scaled, idx)
            push!(bus_loads_scaled, l["load_bus"])
        end
    end
elseif length(bus) > 0
    for (idx, l) in network["load"]
        if l["load_bus"] in bus
            l["pd"] *= scale_factor
            l["qd"] *= scale_factor
            push!(loads_scaled, idx)
            push!(bus_loads_scaled, l["load_bus"])
        end
    end
else # if neither load nor bus is specified, scale all loads
    for (idx, l) in network["load"]
        l["pd"] *= scale_factor
        l["qd"] *= scale_factor
        push!(loads_scaled, idx)
        push!(bus_loads_scaled, l["load_bus"])
    end
end

println("Scaled load number $loads_scaled, connected at bus(es) $bus_loads_scaled, by the scaling factor $scale_factor.")
end

############ PRINTING FUNCTIONS ################

# Function to print the active and reactive power output of generators
function print_gen_dispatch(network_data)

for i in keys(network_data["gen"])
    # active and reactive power generation
    p_gen = network_data["gen"][i]["pg"]*network_data["gen"][i]["mbase"]
    q_gen = network_data["gen"][i]["qg"]*network_data["gen"][i]["mbase"]
    # maximum allowable active and reactive power generation
    p_max = network_data["gen"][i]["pmax"]*network_data["gen"][i]["mbase"]
    q_max = network_data["gen"][i]["qmax"]*network_data["gen"][i]["mbase"]
    # bus location of the generator
    bus = network_data["bus"][string(network_data["gen"][i]["gen_bus"])]["name"]

    # printing the results
    if p_gen>p_max
        println("OVERLOAD!!!!! Gen $i at bus $bus: OVERLOAD!!!!!") 
    elseif p_gen>=0.9*p_max
        println("Gen $i at bus $bus: HIGH LOAD")
    else
        println("Gen $i at bus $bus:")
    end
    println("Active power: $p_gen (max is $p_max)")
    println("Reactive power: $q_gen (max is $q_max)")
    println()
end

end


# Printing the branch data

function print_power_flows(network_data)

for i in keys(network_data["branch"])
    # Power flow at from-bus
    pf = network_data["branch"][i]["pf"]*network_data["baseMVA"]
    qf = network_data["branch"][i]["qf"]*network_data["baseMVA"]
    sf = sqrt(pf^2 + qf^2)

    # Power flow at to-bus
    pt = network_data["branch"][i]["pt"]*network_data["baseMVA"]
    qt = network_data["branch"][i]["qt"]*network_data["baseMVA"]
    st = sqrt(pt^2 + qt^2)

    # Maximum power flow
    smax = network_data["branch"][i]["rate_a"]*network_data["baseMVA"]
    sflow = max(sf, st)

    # From-bus and to-bus
    f_bus = network_data["bus"][string(network_data["branch"][i]["f_bus"])]["name"]
    t_bus = network_data["bus"][string(network_data["branch"][i]["t_bus"])]["name"]

    if sflow>smax
        println("OVERLOAD!!!!! Line $i from $f_bus to $t_bus: OVERLOAD!!!!!") 
    elseif sflow>=0.9*smax
        println("Line $i from $f_bus to $t_bus: HIGH LOAD")
    else
        println("Line $i from $f_bus to $t_bus:")
    end
    println("Power flow: $sflow (max is $smax)")
    println()
end

end

# Printing the branch data

function print_voltage_magnitudes(network_data)

for i in keys(network_data["bus"])
    
    # Read in voltage data
    vm = network_data["bus"][i]["vm"]
    vmax = network_data["bus"][i]["vmax"]
    vmin = network_data["bus"][i]["vmin"]
    
    if vm>vmax
        println("OVERVOLTAGE!!!!! Bus $i: Voltage magnitude is $vm OVERVOLTAGE!!!!!") 
    elseif vm<vmin
        println("UNDERVOLTAGE!!!!! Bus $i: Voltage magnitude is $vm UNDERVOLTAGE!!!!!") 
    else
        println("Bus $i: Voltage magnitude is $vm")
    end
end

end

# Function to plot Wisconsin map
function plot_WI(;size=(800,600))
# download from  https://github.com/deldersveld/topojson/blob/master/countries/us-states/WI-55-wisconsin-counties.json
file_path = "./WI-55-wisconsin-counties.json"
WI = VegaDatasets.VegaJSONDataset(JSON.parsefile(file_path), Path(file_path))

(width,height) = size
p = @vlplot(
     mark={
         :geoshape,
         fill=:lightgray,
         stroke=:white
     },
     data={
         values= WI,
         format={
             type=:topojson,
             feature=:cb_2015_wisconsin_county_20m
         }
     },
     projection={type=:albersUsa},
 )
 return p
end

# Function to plot power flow in network
function plot_powerflow(case)
plot1 = powerplot!(plot_WI(), case,
            bus_data=:vm,
            bus_data_type=:quantitative,
            gen_data=:pg,
            gen_data_type=:quantitative,
            branch_data=:pt,
            branch_data_type=:quantitative,
            branch_color=["black","black","red"],
            gen_color=["black","black","red"],
bus_color=colorscheme2array(ColorSchemes.winter),
fixed=true,
)

plot1.layer[2]["transform"] = Dict{String, Any}[
Dict("calculate"=>"abs(datum.pt)/datum.rate_a*100", "as"=>"branch_Percent_Loading"),
Dict("calculate"=>"abs(datum.pt)", "as"=>"BranchPower")
]
plot1.layer[2]["encoding"]["color"]["field"]="branch_Percent_Loading"
plot1.layer[2]["encoding"]["color"]["title"]="Branch Utilization %"
plot1.layer[2]["encoding"]["color"]["scale"]["domain"]=[0,100]

plot1.layer[4]["encoding"]["color"]["title"]="Voltage magnitude"
plot1.layer[4]["encoding"]["size"]=Dict("value"=>100)

plot1.layer[5]["transform"] = Dict{String, Any}[
Dict("calculate"=>"datum.pg/datum.pmax*100", "as"=>"gen_Percent_Loading"),
Dict("calculate"=>"datum.pg", "as"=>"GenPower")
]
plot1.layer[5]["encoding"]["color"]["field"]="gen_Percent_Loading"
plot1.layer[5]["encoding"]["color"]["scale"]["domain"]=[0,100]
plot1.layer[5]["encoding"]["color"]["title"]="Gen Utilization %"

plot1.layer[4]["encoding"]["color"]["legend"]=Dict("orient"=>"bottom", "offset"=>-20)
plot1.layer[2]["encoding"]["color"]["legend"]=Dict("orient"=>"right", "offset"=>50)
plot1.layer[5]["encoding"]["color"]["legend"]=Dict("orient"=>"right")

@set! plot1.resolve.scale.size=:independent
@set! plot1.resolve.scale.color=:shared

return electrondisplay(Experimental.cartesian2geo!(plot1))
end


function solve_power_flow!(network; kwargs...)

result = PowerModels.compute_ac_pf(network; kwargs...)


# Check if the power flow found a solution
#if "bus" in keys(result)
if result["termination_status"]==true
    println("The solver found an AC power flow solution.")
    success = true
    # Apply the solution to the network data
    PowerModels.update_data!(network, result["solution"])
    
    # Also calculate and apply the branch flows
    flows = calc_branch_flow_ac(network)
    update_data!(network, flows)
else # if the power flow did not find a solution
    println("WARNING! The solver did not find an AC power flow solution! The values of the power flow variables in the network are not a power flow solution.")
    success = false
    PowerModels.logger_config!("warn")
    result = PowerModels.compute_ac_pf(network; kwargs...)
    PowerModels.logger_config!("error")
end

# Return a Boolean indicating whether the solver successfully found a solution
return success
end