#=

mathematical billiards bounce sequence generator

run from the command line, format is:
$ julia billiard_bounce_sequences.jl "[[x1,y1],[x2,y2],[x3,y3]]" low_theta:theta_step:high_theta low_x:x_step:high_x sequence_length

the [xi,yi] correspond to vertices of the triangular billiard table;
by convention, (1) [x1,y1] is at the origin (=[0.0,0.0]); (2) the triangle's longest side (call it the base side) is the segment between [x1,y1] and [x2,y2], and it must be on the x-axis; (3) the other two sides are facing "up" (in +y halfplane)

the x range and theta range describe a discretized window in the (x,theta) phase space, corresponding to the total set of initial conditions (trajectories) we're interested in the bounce sequences of;
by convention, x values correspond to trajectory foot points on the triangle's base side, in interval (0.0,x2) (see note below on x foot point limits), and angles point upward (in (0,pi))

eg, for 1,2,7 triangle:
$ julia billiard_bounce_sequences.jl "[[0.0,0.0],[2.0,0.0],[0.6180339887498950,0.4490279765795854]]" 2.17975:0.000061966:2.241716 0.25568:0.00029656:0.55224 250

to save the results, in header area set out_file to some appropriate name (eg out_file = "/home/bounces.txt");
the out_file will have on each line a bounce sequence; bounce sequences will be ordered as follows:
<angle 1><position 1>, <angle 1><position 2>,...,<angle 1><position end>, <angle 2><position 1>, ...
where angle i and position j go over their respective ranges

note: including trajectories with foot points (x values) too close to the limits 0.0 and x2 (ie too close to the x-axis base side's endpoints) can create errors (eg throw error "ERROR: LoadError: BoundsError:..."); the foot points need to be sufficiently "inside" (0.0,x2) so that after the small displacement to get the footpoint just "above" the x-axis (see eps_dsp parameter), the trajectory starting point is still inside the triangle

=#

using Printf
#using PlotlyJS	# if wanted, for visualizing or debugging
using DynamicalBilliards
using DelimitedFiles

out_file = ""	# set this to a file that will contain the bounce sequences

eps_dsp = 1e-9	# small displacement to get initial particle point off base wall of triangle


# (optional) this allows plotting the billiard boundaries, including for Sinai billiard
function parse_gen_billiard_type(data_in,f)
	if typeof(data_in)==FiniteWall{Float64} || typeof(data_in)==InfiniteWall{Float64}
		p_1 = data_in.sp
		p_2 = data_in.ep
		x_vals = Float64[]
		y_vals = Float64[]
		push!(x_vals,p_1[1])
		push!(x_vals,p_2[1])
		push!(y_vals,p_1[2])
		push!(y_vals,p_2[2])
		return x_vals, y_vals	
	elseif typeof(data_in)== Disk{Float64}
		# f will be an integer 1 thru 4 direction it bulges toward--1 is "east", 2 is "north, 3 is "west", 4 is "south"; if it's -1, then whole circle
		if f==-1
			theta = LinRange(0,2*pi,500)
		else
			theta = LinRange((f-1)*pi/2,(f+1)*pi/2,500)
		end
		out_vec = data_in.c[1] .+ (data_in.r)*sin.(theta), data_in.c[2] .+ (data_in.r)*cos.(theta)
		return out_vec
	end
end

# this checks vertex set--of form [[x1,y1],[x2,y2],[x3,y3]]
function tri_len_proc(tri_vrt)
    side_len_A = sqrt((tri_vrt[2][1]-tri_vrt[1][1])^2 +
                  (tri_vrt[2][2]-tri_vrt[1][2])^2)
    side_len_B = sqrt((tri_vrt[3][1]-tri_vrt[2][1])^2 +
                  (tri_vrt[3][2]-tri_vrt[2][2])^2)
    side_len_C = sqrt((tri_vrt[1][1]-tri_vrt[3][1])^2 +
                  (tri_vrt[1][2]-tri_vrt[3][2])^2)
    if side_len_A < side_len_B || side_len_A < side_len_C
        println("triangle side A is not the longest")
        exit()
	end
end

# wall collision index propogator, which returns which obstacles (wall indices) were hit;
# note this always seems to start with index 0 (zero)--this seems like a null index, since the obstacles seem always indexed starting at 1; so can discard the 1st wall index (if it's 0)
function bill_prop_vobs(billiard,angle,position,time_tot)

	# create ball / particle; displace slightly in direction of angle to get off base wall;
	# Particle expects as first 2 arguments the x and y value for particle position; 
	# Particle's 3rd argument is trajectory angle;
	# Particle automatically assigns norm of velocity (speed) to 1
	pt = Particle(position+eps_dsp*cos(angle),eps_dsp*sin(angle),angle)

	# iterate
	time_col, obst_ind = visited_obstacles!(pt,billiard,time_tot)	# if fed an integer for time_tot, it will run through that number of boundary collisions / bounces

	return obst_ind

end

# (optional) time series propogator, which returns the total path of the billiard, under some sampling frequency, which can get interspersed by collision events
function bill_prop_timeseries(billiard, angle, time_tot, dt_val)

	# create ball / particle; displace slightly in direction of angle to get off base wall;
	# Particle expects as first 2 arguments the x and y value for particle position; 
	# Particle's 3rd argument is trajectory angle;
	# Particle automatically assigns norm of velocity (speed) to 1
	pt = Particle(position+eps_dsp*cos(angle),eps_dsp*sin(angle),angle)

	# iterate; timeseries returns the whole time series (path) of the particle
	xt,yt,vxt,vyt,t = timeseries!(pt,billiard,time_tot; dt=dt_val)

	return xt, yt, t	# for timeseries!

end

# (optional) "evolve" propogator, which goes from collision to collision, registering the [x,y] of each collision event
function bill_prop_evolve(billiard, angle, time_tot)

	# create ball / particle; displace slightly in direction of angle to get off base wall;
	# Particle expects as first 2 arguments the x and y value for particle position; 
	# Particle's 3rd argument is trajectory angle;
	# Particle automatically assigns norm of velocity (speed) to 1
	pt = Particle(position+eps_dsp*cos(angle),eps_dsp*sin(angle),angle)

	# iterate; evolve returns the time, position, and velocity at the collision points only
	ct,poss,vels = evolve!(pt,billiard,time_tot)	# if fed an integer for time_tot, it will run through that number of boundary collisions / bounces

	return poss	# for evolve!

end

# parse the series of integer wall collision indices into a symbol string
function os_parse(wall_hits)
	symbs = ["1", "2", "3"] # assume index #1 is always the longest wall of the triangle
	str_out = "";
	for ii in wall_hits[2:end]
		str_out *= symbs[ii]
	end
	return str_out
end



###
# main routine
###

# command line invocation: $ julia billiard_bounce_sequences <tri_vertices> <ang_vec> <pos_vec> <t_max>
# tri_vertices is 2-dim array of vertices; eg "[[0.0,0.0],[3.0,0.0],[1.0,1.0]]" (must be in float format);
# ang_vec will be a range of angles to try; eg 0.1:pi/10:pi/2;
# pos_vec will be a range of positions of the ball along the base side (on x-axis);
# t_max is the number of wall collisions (symbols) to generate for each initial condition (should be an integer)


## parse input arguments

if length(ARGS)!=4
	println("problem with input arguments in billiard_driver")
	exit()
end

tmp_arg = ARGS[1]
tmp_arg = replace(tmp_arg,"],["=>";")
tmp_arg = replace(tmp_arg,"["=>"")
tmp_arg = replace(tmp_arg,"]"=>"")	# now have x1,y1;x2,y2;x3,y3
tmp_cp = split(tmp_arg,";")
tri_vertices = []
for pair in tmp_cp
	tmp_xy = split(pair,",")
	push!(tri_vertices,[parse(Float64,tmp_xy[1]),parse(Float64,tmp_xy[2])])
end

ang_arr = Float64[]
ang_ran = split(ARGS[2],":")
for ii in eachindex(ang_ran)
	push!(ang_arr,parse(Float64,ang_ran[ii]))
end

pos_arr = Float64[]
pos_ran = split(ARGS[3],":");
for ii in eachindex(pos_ran)
	push!(pos_arr,parse(Float64,pos_ran[ii]))
end

t_max = parse(Int64,ARGS[4])

if length(ang_arr)==3
	ang_vec = Array(ang_arr[1]:ang_arr[2]:ang_arr[3])
elseif length(ang_arr)==1
	ang_vec = ang_arr[1]
else
	println("problem with inputs")
	exit()
end

if length(pos_arr)==3
	pos_vec = Array(pos_arr[1]:pos_arr[2]:pos_arr[3])
elseif length(pos_arr)==1
	pos_vec = pos_arr[1]
else
	println("problem with inputs")
	exit()
end

@printf("%d angle elements\n", length(ang_vec))
@printf("%d position elements\n", length(pos_vec))


## initialize for generating bounces sequences over (x,theta) discretized ranges

tri_len_proc(tri_vertices)

# (optional) for full timeseries and/or "evolve":
dt_val = 0.05 # note, this can be also tuned based on ts_parse picking too many wall intersections (can check plots)--if it's always picking a time point on the boundary, can lower dt_val to make the mesh a bit finer--this collision-point tendency was from the combo of how ts_parse looks for max index of time point at or below desired, and the timeseries function filling the time vector by overriding the regular dt schedule w/ time points it hits a wall (then resuming dt spacing till next wall hit)
xt = Array{Float64}
yt = Array{Float64}
poss = Array{Float64,2}

obs_vec = Array{Int64}	# vector of obstacle collisions
out_vec_str = String[] # a vector of strings

# create billiard table; note this creates a complete billiard table object (no need to call Billiard() on the result)
billiard = billiard_vertices(tri_vertices)


## main bounce sequence generation loop

for ang in ang_vec
	for pos in pos_vec

		# run propogator
		#xt, yt, t_in = bill_prop_timeseries(billiard,ang_vec[1],t_max,dt_val)	# if wanting full ball track, timeseries
		obs_vec = bill_prop_vobs(billiard,ang,pos,t_max)	# just returns integer indicies of wall collisions; this always(?) starts with a null wall index (of 0) (all actual walls are indexed starting from 1)
		str_out = os_parse(obs_vec)
		push!(out_vec_str,str_out)
	end
end


if out_file == ""
	@printf("no output file set; bounce sequence output will not be saved")
else
	writedlm(out_file,out_vec_str)
end

















