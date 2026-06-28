abstract type AbstractStar end
abstract type AbstractStarWithSpot <: AbstractStar end

"""
	struct Star <: AbstractStar

Simple star with radius `R` solar radii, mass `M` solar masses, and temperature `T` kelvins.

# Fields

- `name :: String`     Star name
- `R :: Float64`       Star radius in Solar radii
- `M :: Float64`       Star mass in Solar masses
- `T :: Float64`       Star Temperature in kelvins
- `v_eq :: Float64`    Star escape velocity in km/s
- `v_esc :: Float64`   Star equatorial rotational velocity in km/s

# Constructors

    Star(name :: AbstractString, R :: Real, M :: Real, T :: Real, v_eq :: Real)

Creates star from arguments. `v_esc` is computed automatically from radius and mass.

    Star(name :: AbstractString)
	
Loads star with name if it has been saved before. To save star use `savestar(star :: AbstractStar)`

# Supertypes Hierarchy

    Star <: AbstractStar
"""
struct Star <: AbstractStar
	name :: String   # Star name
	R :: Float64     # Star radius in Solar radii
	M :: Float64     # Star mass in Solar masses
	T :: Float64     # Star Temperature in kelvins
	v_eq :: Float64  # Star escape velocity in km/s
	v_esc :: Float64 # Star equatorial rotational velocity in km/s
	function Star(name :: AbstractString, R, M, T, v_eq)
		v_esc = √(2*G*M☉/R☉*M/R)*1e-5
		new(name, Float64(R), Float64(M), Float64(T), Float64(v_eq), Float64(v_esc))
	end
	function Star(name)
		loadstar(name)
	end
end

"""
    starradiiincm(star :: AbstractStar)

Returns `star` radius in centimeters
"""
starradiiincm(star :: AbstractStar) = star.R*R☉

"""
    starmassing(star :: AbstractStar)

Returns `star` mass in grams
"""
starmassing(star :: AbstractStar) = star.M*M☉


"""
    starvescincms(star :: AbstractStar)

Returns `star` escape velocity in centimeters per second
"""
starvescincms(star :: AbstractStar) = star.v_esc*1e5

"""
    starvescincms(star :: AbstractStar)

Returns `star` equatorial rotational velocity in centimeters per second
"""
starveqincms(star :: AbstractStar) = star.v_eq*1e5

"""
    struct MagnetosphereSpotStar <: AbstractStarWithSpot

Star with hotspots from dipole magnetosphere: two strips between two polar angles 
on the star surface (one strip for each pole).

# Fields

Standart fields of `Star` + hot spot parameteres
- `T_spot :: Float64` Spot temperature [K]
- `θ_1 :: Float64` polar angle of spot border that is closer to the pole
- `θ_2 :: Float64` polar angle of spot border that is further from the pole

# Constructors

    MagnetosphereSpotStar(name :: AbstractString, R, M, T, v_eq, T_spot, θ_1, θ_2)

Creates star from arguments. `v_esc` is computed automatically from radius and mass.

    MagnetosphereSpotStar(star :: AbstractStar, T_spot, r_mi, r_mo)

Returns new star with the same parameteres as `star` but with dipole magnetosphere hotspots.
`θ_1` and `θ_2` are computed from `r_mi` and `r_mo` assuming they are given in star radii:
- `θ_1 = asin(√(1/r_mo))`
- `θ_2 = asin(√(1/r_mi))`

`MagnetosphereSpotStar(name :: AbstractString, T_spot, r_mi, r_mo)`

Loads star (see `Star` docs) and adds spots on it

    MagnetosphereSpotStarFromMdot(star :: AbstractStar, Ṁ, r_mi, r_mo)

Spot termperature is computed from accretion temp `Ṁ` (given in [M_sun/yr]). Assumes that all kinetic energy
is radiated away as blackbody radiation.


# Supertypes Hierarchy

    MagnetosphereSpotStar <: AbstractStarWithSpot <: AbstractStar
"""
struct MagnetosphereSpotStar <: AbstractStarWithSpot
	name :: String
	R :: Float64
	M :: Float64
	T :: Float64
	v_eq :: Float64
	v_esc :: Float64

	T_spot :: Float64
	θ_1 :: Float64
	θ_2 :: Float64

	function MagnetosphereSpotStar(name :: AbstractString, R, M, T, v_eq, T_spot, θ_1, θ_2)
		v_esc = √(2*G*M☉/R☉*M/R)*1e-5
		new(name, R, M, T, v_eq, v_esc, T_spot, min(θ_1, θ_2), max(θ_1, θ_2))
	end

	function MagnetosphereSpotStar(name :: AbstractString, T_spot, θ_1, θ_2)
		star = loadstar(name)
		MagnetosphereSpotStar(star, T_spot, θ_1, θ_2)
	end

	function MagnetosphereSpotStar(star :: AbstractStar, T_spot, θ_1, θ_2)
		name = star.name
		R = star.R; M = star.M; T = star.T
		v_eq = star.v_eq
		MagnetosphereSpotStar(name, Float64.([R, M, T, v_eq, T_spot, θ_1, θ_2])...)
	end
end

"""

    MagnetosphereSpotStarFromMdot(star :: AbstractStar, Ṁ, r_mi, r_mo)

Returns `MagnetosphereSpotStar`. Spot termperature is computed from accretion temp `Ṁ` 
(given in [M_sun/yr]). Assumes that all kinetic energy is radiated away as blackbody radiation.

"""
function MagnetosphereSpotStarFromMdot(star :: AbstractStar, Ṁ, r_mi, r_mo)
	T_spot = calcmagspottemperature(star, Ṁ, r_mi, r_mo)
	MagnetosphereSpotStar(star, T_spot, asin(√(1/r_mi)), asin(√(1/r_mo)))
end

planckfunction(ν, T) =  2*h*ν^3/c^2/(exp(h*ν/(kB*T)) - 1)
starcontinuum(star :: AbstractStar, ν) = planckfunction(ν, star.T)
spotcontinuum(star :: AbstractStarWithSpot, ν) = planckfunction(ν, star.T_spot)

function pictureplanecontinuum(star :: AbstractStar, ν, x, y, star_axis)
	starcontinuum(star, ν)
end

function pointraycontinuum(star :: AbstractStar, ν, r, α, β, star_axis)
	starcontinuum(star, ν)
end

function isthispointonspot(star :: AbstractStar, x, y, star_axis)
	return false
end

function isthispointonspot(star :: MagnetosphereSpotStar, x, y, star_axis)
	if star.T_spot < star.T
		return false # there is no spot 
	end
	z = real(√(1 - x^2 - y^2 + 0im))
	star_pos = SA[x, y, z]
	θ = acos(abs(sum(star_pos .* star_axis)))
	midspot = (star.θ_1 + star.θ_2)/2
	spotwidth = abs(star.θ_1 - star.θ_2)
	if abs(θ - midspot) < spotwidth/2
		return true
	else
		return false
	end
end

function isthisdirectiononspot(star :: AbstractStarWithSpot, z :: Real, α :: Real, β :: Real, 
												star_axis :: AbstractArray{T, 1}) where {T <: Real}
	if z ≈ 1.0
		z = 1.0 + 10*eps(1.0)
	end
	if z*sin(α) < 1e-6
		d = z - z^2*sin(α)^2
	else
		d = -z*cos(α) - real(√(1 - z^2*sin(α)^2 + 0im))
	end
	R = d*sin(α)
	x, y = R*sin(β), -R*cos(β)
	return isthispointonspot(star, x, y, star_axis)
end

function isthisdirectiononspot(star :: AbstractStarWithSpot, r::Real, θ::Real, α::Real, β::Real)
	star_axis = SA[0, sin(θ), cos(θ)]
	return isthisdirectiononspot(star, r, α, β, star_axis)
end

function pictureplanecontinuum(star :: AbstractStarWithSpot, ν, x, y, star_axis)
	if isthispointonspot(star, x, y, star_axis)
		return spotcontinuum(star, ν)
	else
		return starcontinuum(star, ν)
	end
end

function pointraycontinuum(star :: AbstractStarWithSpot, ν :: Real, z :: Real, α :: Real, β :: Real,
											 star_axis :: Array{T, 1}) where {T <: Real}
	if isthisdirectiononspot(star, z, α, β, star_axis)
		return spotcontinuum(star, ν)
	else
		return starcontinuum(star, ν)
	end
end

function pointraycontinuum(star :: AbstractStarWithSpot, ν :: Real, r :: Real, θ :: Real, α :: Real, β :: Real)
	if isthisdirectiononspot(star, r, θ, α, β)
		return spotcontinuum(star, ν)
	else
		return starcontinuum(star, ν)
	end
end

function calcmagspottemperature(star :: AbstractStar, M_dot, r_mi, r_mo)
    θ_1 = asin(√(1/r_mo)); θ_2 = asin(√(1/r_mi))
    T_spot = (1e10*M☉/year_seconds/R☉^2*(star.v_esc)^2/2*M_dot*(1 - 1/2*(r_mi + r_mo)/(r_mi*r_mo))
						/(star.R^2*4*π*abs(cos(θ_2)-cos(θ_1))*σ) + star.T^4)^0.25
    T_spot = T_spot
    return T_spot
end

corotationradius(star :: AbstractStar) = ((star.v_esc^2/star.v_eq^2/2)^(1/3))

Base.length(star :: AbstractStar) = 1
Base.iterate(star :: AbstractStar) = (star, nothing)
Base.iterate(star :: AbstractStar, n :: Nothing) = nothing

genstardata(star :: AbstractStar) = (
			@sprintf("# R[R☉] M[M☉] T[K] v_eq[km/s]\n %5.2f %5.2f %5.5g %6.2f\n", 
										        star.R, star.M, star.T, star.v_eq))

genspotdata(star :: MagnetosphereSpotStar) = (
				@sprintf("# MagnetosphereSpot\n# θ_1[deg] θ_2[deg] T_spot[K]\n#%8.1f %8.1f %9.0f\n",
													star.θ_1/π*180, star.θ_2/π*180, star.T_spot))

function loadstar(name :: AbstractString; dir = "stars")
	star_file = dir*'/'*name*'/'*name*".dat"
	try
		global star_file_lines
		star_file_lines = readlines(star_file)
	catch err
		if err isa SystemError
			print("Can't find file "*name*" in directoty "*star_dir)
		end
	end
	data_line = star_file_lines[2]
	data_strings = split(data_line)
	R, M, T, v_eq = parse.(Float64, data_strings)
	star = Star(name, R, M, T, v_eq)
	return star
end

function savestar(star :: AbstractStar; dir = "stars")
	star_file = dir*'/'*star.name*'/'*star.name*".dat"
	mkpath(dir*'/'*star.name)
	star_out = open(star_file, "w")
	data = genstardata(star)
	print(star_out, data)
	close(star_out)
end
