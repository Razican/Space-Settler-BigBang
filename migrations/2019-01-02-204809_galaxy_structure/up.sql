-- This creates the tables for the galaxy information.

CREATE TYPE star_class AS ENUM ('O', 'B', 'A', 'F', 'G', 'K', 'M', 'black_hole', 'neutron_star', 'quark_star', 'white_dwarf');

-- TODO: Star systems with multiple stars

-- Table of stars
CREATE TABLE stars (
    id SERIAL NOT NULL PRIMARY KEY,
    orb_radius INTEGER NOT NULL CHECK (orb_radius > 0), -- Orbit radius, in 1/100 light years
    class star_class NOT NULL, -- Star class
    mass DOUBLE PRECISION NOT NULL CHECK (radius > 0), -- Mass of the star, in kg
    radius DOUBLE PRECISION NOT NULL CHECK (radius > 0), -- Radius of the star, in meters
    temperature INTEGER CHECK (temperature IS NULL OR temperature > 0), -- Effective surface temperature of the star, in 1/1,000 Kelvin
    CHECK ((class = 'black_hole' AND temperature IS NULL) OR (class != 'black_hole' AND temperature IS NOT NULL))
);

-- Type of planets
CREATE TYPE planet_type AS ENUM ('rocky', 'gaseous');

-- Table of planets
CREATE TABLE planets (
    id BIGSERIAL NOT NULL PRIMARY KEY,
    star_id INTEGER NOT NULL REFERENCES stars(id) ON DELETE CASCADE, -- Star the planet is orbiting
    position SMALLINT NOT NULL CHECK (position > 0), -- Position of the planet in the solar system
    orb_ecc DOUBLE PRECISION NOT NULL CHECK (orb_ecc >= 0), -- Eccentricity of the planet's orbit (0 - 1 for closed orbits)
    orb_sma DOUBLE PRECISION NOT NULL CHECK (orb_sma > 0), -- Semimajor axis of the planet's orbit, in meters
    orb_incl DOUBLE PRECISION NOT NULL CHECK (orb_incl >= 0 AND orb_incl <= 2*PI()), -- Inclination of the planet's orbit, in radians
    orb_lan DOUBLE PRECISION NOT NULL CHECK (orb_lan >= 0 AND orb_lan <= 2*PI()), -- Longitude of the ascending node of the planet's orbit, in radians
    orb_arg_p DOUBLE PRECISION NOT NULL CHECK (orb_arg_p >= 0 AND orb_arg_p <= 2*PI()), -- Argument of the periapsis of the planet's orbit, in radians
    orb_m0 DOUBLE PRECISION NOT NULL CHECK (orb_m0 >= 0 AND orb_m0 <= 2*PI()), -- Mean anomaly of the planet's orbit, in radians
    orb_period DOUBLE PRECISION NOT NULL CHECK (orb_period > 0), -- Planet's orbital period, in seconds
    atm_pressure DOUBLE PRECISION CHECK (atm_pressure >= 0), -- Pressure of the planet's atmosphere
    atm_h2o DOUBLE PRECISION CHECK (atm_h2o >= 0 AND atm_h2o <= 1), -- Water percentage of the planet's atmosphere
    atm_co2 DOUBLE PRECISION CHECK (atm_co2 >= 0 AND atm_co2 <= 1), -- Carbon dioxide percentage of the planet's atmosphere
    atm_co DOUBLE PRECISION CHECK (atm_co >= 0 AND atm_co <= 1), -- Carbon monoxide percentage of the planet's atmosphere
    atm_n2 DOUBLE PRECISION CHECK (atm_n2 >= 0 AND atm_n2 <= 1), -- Nitrogen percentage of the planet's atmosphere
    atm_o2 DOUBLE PRECISION CHECK (atm_o2 >= 0 AND atm_o2 <= 1), -- Oxygen percentage of the planet's atmosphere
    atm_ar DOUBLE PRECISION CHECK (atm_ar >= 0 AND atm_ar <= 1), -- Argon percentage of the planet's atmosphere
    atm_so2 DOUBLE PRECISION CHECK (atm_so2 >= 0 AND atm_so2 <= 1), -- Sulfure dioxide percentage of the planet's atmosphere
    atm_ne DOUBLE PRECISION CHECK (atm_ne >= 0 AND atm_ne <= 1), -- Neon percentage of the planet's atmosphere
    atm_ch4 DOUBLE PRECISION CHECK (atm_ch4 >= 0 AND atm_ch4 <= 1), -- Methane percentage of the planet's atmosphere
    atm_he DOUBLE PRECISION CHECK (atm_he >= 0 AND atm_he <= 1), -- Helium percentage of the planet's atmosphere
    ax_tilt DOUBLE PRECISION NOT NULL CHECK (ax_tilt >= 0 AND ax_tilt <= 2*PI()), -- Axial tilt of the planet, in radians
    rot_period DOUBLE PRECISION NOT NULL CHECK (rot_period > 0), -- Rotation period of the planet, in seconds
    planet_type planet_type NOT NULL, -- Planet type (rocky or gaseous)
    surf_fresh_water DOUBLE PRECISION CHECK (surf_fresh_water >= 0 AND surf_fresh_water <= 1), -- Percentage of fresh water in the surface of the planet
    surf_ocean_water DOUBLE PRECISION CHECK (surf_ocean_water >= 0 AND surf_ocean_water <= 1), -- Percentage of ocean water in the surface of the planet
    surf_snow DOUBLE PRECISION CHECK (surf_snow >= 0 AND surf_snow <= 1), -- Percentage of snow in the surface of the planet
    surf_land DOUBLE PRECISION CHECK (surf_land >= 0 AND surf_land <= 1), -- Percentage of land in the surface of the planet
    -- TODO: surface vegetation
    -- TODO: crust element distribution
    -- TODO: life distribution
    bond_albedo DOUBLE PRECISION NOT NULL CHECK (bond_albedo >= 0 AND bond_albedo <= 1), -- The bond albedo of the planet, from 0 to 1
    geometric_albedo DOUBLE PRECISION NOT NULL CHECK (geometric_albedo >= 0 AND geometric_albedo <= 1), -- The geometric albedo of the planet, from 0 to 1
    mass DOUBLE PRECISION NOT NULL CHECK (mass > 0), -- The mass of the planet, in kg
    radius DOUBLE PRECISION NOT NULL CHECK (radius > 0), -- The radius of the planet, in meters
    -- TODO: temperature as a fixed point number
    eff_temp DOUBLE PRECISION NOT NULL CHECK (eff_temp > 0), -- The effective temperature of the planet, in Kelvin
    min_temp DOUBLE PRECISION NOT NULL CHECK (min_temp > 0), -- The minimum average temperature of the planet, in Kelvin
    max_temp DOUBLE PRECISION NOT NULL CHECK (max_temp > 0), -- The maximum average temperature of the planet, in Kelvin
    avg_temp DOUBLE PRECISION NOT NULL CHECK (avg_temp > 0), -- The average temperature of the planet, in Kelvin
    habitable BOOLEAN NOT NULL, -- Wether the planet is habitable
    CHECK (atm_h2o + atm_co2 + atm_co + atm_n2 + atm_o2 + atm_ar + atm_so2 + atm_ne + atm_ch4 + atm_he > 0.99999 AND
           atm_h2o + atm_co2 + atm_co + atm_n2 + atm_o2 + atm_ar + atm_so2 + atm_ne + atm_ch4 + atm_he < 1.00001 AND
           surf_fresh_water + surf_ocean_water + surf_snow + surf_land > 0.99999 AND
           surf_fresh_water + surf_ocean_water + surf_snow + surf_land < 1.00001 AND
           max_temp >= avg_temp AND
           min_temp <= avg_temp
    )
);