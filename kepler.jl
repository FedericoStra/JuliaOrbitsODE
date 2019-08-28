using StaticArrays
using ReferenceFrameRotations: angle_to_dcm

# units: position [km], velocity [km/2]
const μ = 398600.4418

@inline typetol(::Type{Float16}) = 1e-2
@inline typetol(::Type{Float32}) = 1e-4
@inline typetol(::Type{Float64}) = 1e-10
@inline typetol(::Type{BigFloat}) = 1e-50

function M2f(e, M, tol=nothing)
    T = promote_type(typeof(e), typeof(M), Float64)
    if isnothing(tol)
        tol = typetol(T)
    else
        tol::Float64
    end
    M = mod2pi(M)
    E = (M > π) ? M - e : M + e
    sin_E, cos_E = sincos(E)
    while ( abs(E - e*sin_E - M) > tol )
        E = E - (E - e*sin_E - M) / (1 - e*cos_E)
        sin_E, cos_E = sincos(E)
    end
    E = mod2pi(E)
    sin_Eo2, cos_Eo2 = sincos(E/2)
    mod2pi( 2atan(sqrt(1+e)*sin_Eo2, sqrt(1-e)*cos_Eo2) )
end

function f2M(e, f)
    sin_fo2, cos_fo2 = sincos(f/2)
    E = 2atan(sqrt(1-e)*sin_fo2, sqrt(1+e)*cos_fo2)
    mod2pi(E - e*sin(E))
end

function rv_to_coe(r0::AbstractVector{<:Real}, v0::AbstractVector{<:Real})
    r_i = SVector{3}(map(big, r0))
    v_i = SVector{3}(map(big, v0))

    @inbounds begin

        # Position and velocity vector norms and auxiliary dot products.
        r2 = dot(r_i,r_i)
        v2 = dot(v_i,v_i)

        r  = sqrt(r2)
        v  = sqrt(v2)

        rv = dot(r_i,v_i)

        # The type of `r` will be the type of the orbit elements.
        T   = typeof(r)
        Tm0 = T(μ)

        # Angular momentum vector.
        h_i = cross( r_i, v_i )
        h   = norm(h_i)

        # Vector that points to the right ascension of the ascending node (RAAN).
        n_i = SVector{3}(0,0,1) × h_i
        n   = norm(n_i)

        # Eccentricity vector.
        e_i = ( (v2 - Tm0/r)*r_i - rv*v_i )/Tm0

        # Orbit energy.
        ξ = v2/2 - Tm0/r

        # Eccentricity
        # ============

        ecc = norm(e_i)

        # Semi-major axis
        # ===============

        if ecc < 1.0
            a = -Tm0/(2ξ)
        else
            error("Could not convert the provided Cartesian values to Kepler elements.\n" *
                  "The computed eccentricity was not in the interval [0,1.)");
        end

        # Inclination
        # ===========

        cos_i = h_i[3]/h
        cos_i = abs(cos_i) > 1 ? sign(cos_i) : cos_i
        i     = acos(cos_i)

        # Check the type of the orbit to account for special cases
        # ======================================================================

        # Equatorial
        # ----------------------------------------------------------------------

        if n == 0

            # Right Ascension of the Ascending Node.
            # ======================================

            Ω = T(0)

            # Equatorial and elliptical
            # ------------------------------------------------------------------

            if ecc > 0

                # Argument of Perigee
                # ===================

                cos_ω = e_i[1]/ecc
                cos_ω = abs(cos_ω) > 1 ? sign(cos_ω) : cos_ω
                ω     = acos(cos_ω)

                (e_i[3] < 0) && (ω = T(2π) - ω)

                # True anomaly
                # ============

                cos_v = dot(e_i,r_i)/(ecc*r)
                cos_v = abs(cos_v) > 1 ? sign(cos_v) : cos_v
                v     = acos(cos_v)

                (rv < 0) && (v = T(2π) - v)

            # Equatorial and circular
            # ------------------------------------------------------------------

            else
                # Argument of Perigee
                # ===================

                ω = T(0)

                # True anomaly
                # ============

                cos_v = r_i[1]/r
                cos_v = abs(cos_v) > 1 ? sign(cos_v) : cos_v
                v     = acos(cos_v)

                (r_i[2] < 0) && (v = T(2π) - v)
            end

        # Inclined
        # ----------------------------------------------------------------------
        else

            # Right Ascension of the Ascending Node.
            # ======================================

            cos_Ω = n_i[1]/n
            cos_Ω = abs(cos_Ω) > 1 ? sign(cos_Ω) : cos_Ω
            Ω     = acos(cos_Ω)

            (n_i[2] < 0) && (Ω = T(2π) - Ω)

            # Circular and inclined
            # ------------------------------------------------------------------

            if ecc == 0.

                # Argument of Perigee
                # ===================

                ω = T(0)

                # True anomaly
                # ============

                cos_v = dot(n_i,r_i)/(n*r)
                cos_v = abs(cos_v) > 1 ? sign(cos_v) : cos_v
                v     = acos(cos_v)

                (r_i[3] < 0) && (v = T(2π) - v)
            else

                # Argument of Perigee
                # ===================

                cos_ω = dot(n_i,e_i)/(n*ecc)
                cos_ω = abs(cos_ω) > 1 ? sign(cos_ω) : cos_ω
                ω     = acos(cos_ω)

                (e_i[3] < 0) && (ω = T(2π) - ω)

                # True anomaly
                # ============

                cos_v = dot(e_i,r_i)/(ecc*r)
                cos_v = abs(cos_v) > 1 ? sign(cos_v) : cos_v
                v     = acos(cos_v)

                (rv < 0) && (v = T(2π) - v)
            end
        end
    end

    # Return the Keplerian elements.
    # ==============================

    return (a, ecc, i, Ω, ω, v)
end

function coe_to_rv(a::Real, e::Real, i::Real, Ω::Real, ω::Real, f::Real)
    # Check eccentricity.
    if !(0 <= e < 1)
        throw(ArgumentError("Eccentricity must be in the interval [0,1)."))
    end

    # Auxiliary variables.
    sin_f, cos_f = sincos(f)

    # Compute the geocentric distance.
    r = a*(1-e^2)/(1+e*cos_f)

    # Compute the position vector in the orbit plane, defined as:
    #   - The X axis points towards the perigee;
    #   - The Z axis is perpendicular to the orbital plane (right-hand);
    #   - The Y axis completes a right-hand coordinate system.
    r_o = SVector{3}(r*cos_f, r*sin_f, 0)

    # Compute the velocity vector in the orbit plane without perturbations.
    n = sqrt(μ/a^3)
    v_o = n*a/sqrt(1-e^2)*SVector{3}(-sin_f, e+cos_f, 0)

    # Compute the matrix that rotates the orbit reference frame into the
    # inertial reference frame.
    Dio = angle_to_dcm(-ω, -i, -Ω, :ZXZ)

    # Compute the position and velocity represented in the inertial frame.
    r_i = Dio*r_o
    v_i = Dio*v_o

    vcat(r_i, v_i)
end

function coe_to_rv(a::Real, e::Real, i::Real, Ω::Real, ω::Real, f0::Real, t::Real)
    M0 = f2M(e, f0)
    M = M0 + t * sqrt(μ/a^3)
    f = M2f(e, M)
    coe_to_rv(a, e, i, Ω, ω, f)
end

function kepler_solver(r0, v0)
    coe = rv_to_coe(r0, v0)
    function (t)
        coe_to_rv(coe..., t)
    end
end
