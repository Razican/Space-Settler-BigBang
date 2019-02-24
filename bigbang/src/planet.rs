//! Planet module
//!
//! This is the `planet` module. This module contains all the structures, enumerations and
//! implementations needed to define a planet.

use std::{
    f64::consts::{FRAC_PI_6, FRAC_PI_8, PI},
    fmt,
};

use rand::{thread_rng, Rng};
use sps_db::planets as db;

use crate::{consts::*, star::Star, utils::*};

const ALBEDO_TEMP_EPSILON: f64 = 0.005;

/// Planet structure
///
/// This structure defines a planet and contains all the structures needed for a correct definition
/// and representation of itself.
#[derive(Debug, Clone, Copy)]
pub struct Planet {
    orbit: Orbit,
    atmosphere: Option<Atmosphere>,
    planet_type: Type,
    surface: Option<Surface>,
    // crust: Crust,
    // life: Life,
    bond_albedo: f64,
    geometric_albedo: f64,
    mass: f64,
    radius: f64,
    eff_temp: f64,
    min_temp: f64,
    max_temp: f64,
    avg_temp: f64,
    // habitable: bool
}

/// Planet type enumeration
///
/// This enumeration contains the planet types available when creating a new planet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Type {
    /// This represents a rocky planet like Earth. Y has a surface where objects and even life can
    /// be found.
    Rocky,
    /// This represents a gaseous planet, like Jupiter. It has no known surface, and in the case of
    /// having it it would be so deep that it supposes nothing for the atmospheric configuration or
    /// habitability.
    Gaseous,
}

impl fmt::Display for Type {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Type::Rocky => write!(f, "rocky"),
            Type::Gaseous => write!(f, "gaseous"),
        }
    }
}

impl Into<db::PlanetType> for Type {
    fn into(self) -> db::PlanetType {
        match self {
            Type::Rocky => db::PlanetType::Rocky,
            Type::Gaseous => db::PlanetType::Gaseous,
        }
    }
}

/// Orbit structure
///
/// This structure contains all the needed information to define an orbit, and as an extra, it
/// contains the axial tilt of the rotation of the body and the rotation period of the body. It also
/// includes a reference to the star being orbited. It also contains the orbital period, even if it
/// can be calculated with the rest of the parameters.
#[derive(Debug, Clone, Copy)]
pub struct Orbit {
    position: u8,
    ecc: f64,
    sma: f64,
    incl: f64,
    lan: f64,
    arg_p: f64,
    m0: f64,
    period: f64,
    ax_tilt: f64,
    rot_period: f64,
}

/// Atmosphere structure
///
/// This structure defines an atmosphere with all the needed parameters and chemical elements. This
/// will be used not only as a mere curiosity but to calculate greenhouse effect and habitability of
/// the body.
#[derive(Debug, Clone, Copy)]
pub struct Atmosphere {
    pressure: f64,
    h2o: f64,
    co2: f64,
    co: f64,
    n2: f64,
    o2: f64,
    ar: f64,
    so2: f64,
    ne: f64,
    ch4: f64,
    he: f64,
}

/// Surface structure
///
/// This structure represents the composition of the surface of the planet. Gaseous planets will not
/// have it, while rocky for ones it will contain the surface composition in 3 variables (water,
/// snow and land) in percentages.
#[derive(Debug, Clone, Copy)]
pub struct Surface {
    fresh_water: f64,
    ocean_water: f64,
    snow: f64,
    land: f64,
}

/// Crust structure
///
/// This structure represents the crust composition of the planet, time in kg per element. It will
/// only contain the composition of the crust of the planet, that will allow for the players to dig
/// and take resources.
// pub struct Crust {
// }
/// Life structure
///
/// This structure represents the life in the planet.
// pub struct Life {
// }

impl Planet {
    /// Constructs a new `Planet`.
    ///
    /// It creates a random planet taking into account real planet statistics. It requires the
    /// reference to parent star and the two values of the Titius–Bode law, along with the order
    /// of the planet in the solar system and the last body's semi-major axis in the solar system.
    ///
    /// # Examples
    ///
    /// ```
    /// use sps_creation_core::{star::Star, planet::Planet};
    ///
    /// let star = Star::new();
    ///
    /// let num_bodies = star.generate_num_bodies();
    /// let (tb_m, tb_n) = star.generate_titius_bode(num_bodies);
    ///
    /// if num_bodies > 0 {
    ///     let planet = Planet::new(&star, tb_m, tb_n, 1, 0_f64);
    /// }
    /// ```
    pub fn new(star: &Star, m: f64, n: f64, position: u8, last_sm_a: f64) -> Self {
        let orbit = Self::generate_orbit(star, m, n, position, last_sm_a);

        let planet_type = Self::generate_type(
            orbit.sma(),
            star.calculate_luminosity()
                .expect("black holes should not have planets"),
        );
        let mut surface = if planet_type == Type::Rocky {
            Some(Self::generate_surface(None))
        } else {
            None
        };
        let (mass, radius) = Self::generate_properties(planet_type);

        let atmosphere = if planet_type == Type::Rocky {
            Some(Self::generate_atmosphere(Self::calculate_surface_gravity(
                mass, radius,
            )))
        } else {
            None
        };
        let (mut bond_albedo, mut geometric_albedo) =
            Self::calculate_albedo(planet_type, atmosphere.as_ref(), surface.as_ref());

        let mut eff_temp = Self::calculate_t_eff(star, orbit.sma(), bond_albedo);
        let greenhouse = if planet_type == Type::Rocky {
            Self::calculate_greenhouse(atmosphere.as_ref())
        } else {
            1_f64
        };
        let mut avg_temp = eff_temp * greenhouse;
        let (mut min_temp, mut max_temp) = if planet_type == Type::Rocky {
            Self::calculate_surface_temp(avg_temp, atmosphere.as_ref().unwrap().pressure(), &orbit)
        } else {
            (avg_temp, avg_temp)
        };

        if planet_type == Type::Rocky {
            loop {
                let before_temp = eff_temp;

                surface = Some(Self::generate_surface(Some((
                    min_temp,
                    max_temp,
                    atmosphere.as_ref().unwrap(),
                ))));
                let (new_bond_albedo, new_geometric_albedo) =
                    Self::calculate_albedo(planet_type, atmosphere.as_ref(), surface.as_ref());

                eff_temp = Self::calculate_t_eff(star, orbit.sma(), new_bond_albedo);
                avg_temp = eff_temp * greenhouse;
                let (new_min_temp, new_max_temp) = Self::calculate_surface_temp(
                    avg_temp,
                    atmosphere.as_ref().unwrap().pressure(),
                    &orbit,
                );

                bond_albedo = new_bond_albedo;
                geometric_albedo = new_geometric_albedo;
                min_temp = new_min_temp;
                max_temp = new_max_temp;

                if (1_f64 - eff_temp / before_temp).abs() < ALBEDO_TEMP_EPSILON {
                    break;
                }
            }
        }

        Self {
            orbit,
            atmosphere,
            surface,
            planet_type,
            bond_albedo,
            geometric_albedo,
            mass,
            radius,
            eff_temp,
            min_temp,
            max_temp,
            avg_temp,
        }
    }

    /// Gets the database insertable object for the star.
    pub fn insertable(&self) -> db::NewPlanet {
        // use std::i32;

        // debug_assert!(
        //     self.orbit_radius < i32::max_value() as u32,
        //     "star orbit radius too big"
        // );
        // debug_assert!(self.mass > 0_f64, "star mass must be positive");
        // debug_assert!(self.radius > 0_f64, "star radius must be positive");
        // if let Some(temp) = self.temperature {
        //     debug_assert!(
        //         temp <= i32::max_value() as u32,
        //         "star temperature is too big"
        //     );
        // }

        let (surf_fresh_water, surf_ocean_water, surf_snow, surf_land) =
            if let Some(surf) = self.surface {
                (
                    Some(surf.fresh_water),
                    Some(surf.ocean_water),
                    Some(surf.snow),
                    Some(surf.land),
                )
            } else {
                (None, None, None, None)
            };

        db::NewPlanet {
            orbit,      // TODO
            atmosphere, // TODO
            surf_fresh_water,
            surf_ocean_water,
            surf_snow,
            surf_land,
            planet_type: self.planet_type.into(),
            bond_albedo: self.bond_albedo,
            geometric_albedo: self.geometric_albedo,
            mass: self.mass,
            radius: self.radius,
            eff_temp: self.eff_temp,
            min_temp: self.min_temp,
            max_temp: self.max_temp,
            avg_temp: self.avg_temp,
        }
    }

    /// Get `Orbit`
    ///
    /// Gets the orbit information of the planet.
    pub fn orbit(&self) -> &Orbit {
        &self.orbit
    }

    /// Get `Atmosphere`
    ///
    /// Gets the atmosphere information of the planet.
    pub fn atmosphere(&self) -> Option<&Atmosphere> {
        self.atmosphere.as_ref()
    }

    /// Get `Surface`
    ///
    /// Gets the surface information of the planet.
    pub fn surface(&self) -> Option<&Surface> {
        self.surface.as_ref()
    }

    /// Get planet type
    ///
    /// Gets the type of the planet.
    pub fn planet_type(&self) -> &Type {
        &self.planet_type
    }

    /// Get Bond albedo
    ///
    /// Gets the Bond albedo of the planet. The Bond albedo is the reflectivity of the planet, or in
    /// other words, the percentage of light reflected by the planet, from 0 to 1.
    pub fn bond_albedo(&self) -> f64 {
        self.bond_albedo
    }

    /// Get geometric albedo
    ///
    /// Gets the geometric albedo of the planet. The geometric albedo is the reflectivity of the
    /// planet, but only taking into account the light that reflects directly back to the emission,
    /// or in other words, how bright an observer would see the planet if looking at it directly
    /// from the light source, from 0 to 1.
    pub fn geometric_albedo(&self) -> f64 {
        self.geometric_albedo
    }

    /// Get mass
    ///
    /// Gets the mass of the planet in *kg*.
    pub fn mass(&self) -> f64 {
        self.mass
    }

    /// Get radius
    ///
    /// Gets the radius of the planet, in meters (*m*).
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Calculates density
    ///
    /// Calculates the density of the planet, in *kg/m³*.
    pub fn calculate_density(&self) -> f64 {
        self.mass / self.calculate_volume()
    }

    /// Calculate volume
    ///
    /// Calculates the volume of the planet, in *m³*.
    pub fn calculate_volume(&self) -> f64 {
        4_f64 / 3_f64 * PI * self.radius.powi(3)
    }

    /// Calculate greenhouse effect
    ///
    /// Calculates the greenhouse effect of the planet, as a multiplier.
    pub fn greenhouse(&self) -> f64 {
        self.avg_temp / self.eff_temp
    }

    /// Get effective temperature
    ///
    /// Gets the effective temperature of the planet, in Kelvin (*K*).
    pub fn eff_temp(&self) -> f64 {
        self.eff_temp
    }

    /// Get minimum temperature
    ///
    /// Gets the minimum temperature of the planet, in Kelvin (*K*).
    pub fn min_temp(&self) -> f64 {
        self.min_temp
    }

    /// Get average temperature
    ///
    /// Gets the average temperature of the planet, in Kelvin (*K*).
    pub fn avg_temp(&self) -> f64 {
        self.avg_temp
    }

    /// Get maximum temperature
    ///
    /// Gets the maximum temperature of the planet, in Kelvin (*K*).
    pub fn max_temp(&self) -> f64 {
        self.max_temp
    }

    /// Calculate surface gravity
    ///
    /// Calculates the surface gravity of the planet, in *m/s²*.
    pub fn surface_gravity(&self) -> f64 {
        Self::calculate_surface_gravity(self.mass, self.radius)
    }

    /// Check if is earth twin
    ///
    /// Checks if the properties of the planet are similar to the ones in Earth
    pub fn is_earth_twin(&self) -> bool {
        self.planet_type == Type::Rocky
            && self.mass < 2.5_f64 * EARTH_MASS
            && self.mass > 0.5_f64 * EARTH_MASS
            && self.radius < 1.5_f64 * EARTH_RADIUS
            && self.radius > 0.6_f64 * EARTH_RADIUS
            && self.min_temp > 200_f64
            && self.min_temp < 280_f64
            && self.avg_temp > 270_f64
            && self.avg_temp < 310_f64
            && self.max_temp > 275_f64
            && self.max_temp < 330_f64
            && self.atmosphere.unwrap().pressure() > 40_000_f64
            && self.atmosphere.unwrap().pressure() < 150_000_f64
            && self.atmosphere.unwrap().o2() < 0.35_f64
            && self.atmosphere.unwrap().o2() > 0.15_f64
            && self.atmosphere.unwrap().co2() < 0.1_f64
            && self.atmosphere.unwrap().co() < 0.01_f64
    }

    /// Returns wether the planet is habitable.
    pub fn is_habitable(&self) -> bool {
        self.planet_type == Type::Rocky
            && self.surface.unwrap().ocean_water() > 0.1
            && self.surface.unwrap().fresh_water() > 0.01
            && self.surface.unwrap().snow() > 0.01
            && self.surface.unwrap().land() > 0.1
            && self.atmosphere.unwrap().pressure() < 200_000_f64
            && self.surface_gravity() > 3_f64
            && self.surface_gravity() < 15_f64
    }

    /// Check Roche limit
    ///
    /// Checks if the Roche limit for the planet is correct.
    pub fn is_roche_ok(&self, star: &Star) -> bool {
        let rigid = self.radius * (3_f64 * star.mass() / self.mass).powf(1_f64 / 3_f64);
        let fluid = 2.455 * self.radius * (star.mass() / self.mass).powf(1_f64 / 3_f64);

        let roche_limit = match self.planet_type {
            Type::Rocky => thread_rng().gen_range(rigid, rigid * 0.3 + fluid * 0.7),
            Type::Gaseous => thread_rng().gen_range(rigid * 0.3 + fluid * 0.7, fluid),
        };

        self.orbit().periapsis() > roche_limit
    }

    // ----------  generators  ----------

    /// Generate orbit
    ///
    /// Generates the orbit of the planet taking into account the Titius-Bode law, the last planet's
    /// semimajor axis and the position in the system.
    fn generate_orbit(st: &Star, m: f64, n: f64, position: u8, last_sma: f64) -> Orbit {
        let mut sma = (m * f64::from(position) - n).exp() * AU;
        sma = thread_rng().gen_range(sma * 0.9, sma * 1.15);

        if sma < last_sma * 1.15 {
            sma = last_sma * thread_rng().gen_range(1.15_f64, 1.25_f64);
        }

        let ecc = if sma / AU < st.mass() / (SUN_MASS * 2_f64) {
            thread_rng().gen_range(0.05_f64, 0.3_f64)
        } else if sma / AU < st.mass() / SUN_MASS * 30_f64 {
            thread_rng().gen_range(0_f64, 0.15_f64)
        } else {
            thread_rng().gen_range(0_f64, 0.1_f64)
        };

        let period = 2_f64 * PI * (sma.powi(3) / (G * st.mass())).sqrt();

        let incl = thread_rng().gen_range(0_f64, PI / 18_f64);
        let lan = thread_rng().gen_range(0_f64, 2_f64 * PI);
        let arg_p = thread_rng().gen_range(0_f64, 2_f64 * PI);
        let m0 = thread_rng().gen_range(0_f64, 2_f64 * PI);

        let (ax_tilt, rot_period) = Self::generate_rotation(st, sma, period);

        Orbit {
            position,
            ecc,
            sma,
            incl,
            lan,
            arg_p,
            m0,
            period,
            ax_tilt,
            rot_period,
        }
    }

    /// Generate rotation
    ///
    /// Generates the rotation of the planet taking into account tidal lock and orbital resonance.
    fn generate_rotation(st: &Star, sma: f64, orb_period: f64) -> (f64, f64) {
        let tidal_lock = (st.mass() / SUN_MASS).sqrt() / 2_f64;

        if sma / AU > tidal_lock {
            let ax_tilt = if thread_rng().gen() {
                thread_rng().gen_range(0_f64, PI)
            } else {
                thread_rng().gen_range(FRAC_PI_8, FRAC_PI_6) // 15° - 30°
            };

            let rot_period = if ax_tilt > PI / 2_f64 {
                if orb_period < 50_000_f64 {
                    -thread_rng().gen_range(orb_period * 0.8, orb_period - 1_f64)
                } else {
                    -thread_rng().gen_range(
                        50_000_f64,
                        if orb_period < 25_000_000_f64 {
                            orb_period - 1_f64
                        } else {
                            25_000_000_f64
                        },
                    )
                }
            } else if orb_period < 18_000_f64 {
                thread_rng().gen_range(orb_period * 0.8, orb_period - 1_f64)
            } else {
                thread_rng().gen_range(
                    18_000_f64,
                    if orb_period < 180_000_f64 {
                        orb_period - 1_f64
                    } else {
                        180_000_f64
                    },
                )
            };

            (ax_tilt, rot_period)
        } else if sma > tidal_lock.sqrt() / 3_f64 {
            let ax_tilt = thread_rng().gen_range(0_f64, 0.017454_f64); // 0° - 1°
                                                                       // Resonance
            let rot_period = orb_period * 2_f64 / f64::from(thread_rng().gen_range(3, 6));

            (ax_tilt, rot_period)
        } else {
            // Tidal lock
            (0_f64, orb_period)
        }
    }

    /// Generate planet type
    ///
    /// Generates the Type of the planet depending on star and the SMa of the orbit.
    fn generate_type(sma: f64, luminosity: f64) -> Type {
        if sma / (AU * 2_f64) < (luminosity / SUN_LUMINOSITY).sqrt() {
            if thread_rng().gen_range(0, 2) == 0 {
                Type::Gaseous
            } else {
                Type::Rocky
            }
        } else if luminosity > 1.923e+27_f64 && // 5*SUN_LUMINOSITY
            sma/AU < (luminosity/SUN_LUMINOSITY).sqrt()*50_f64
        {
            if thread_rng().gen_range(0, 3) == 0 {
                Type::Rocky
            } else {
                Type::Gaseous
            }
        } else if sma / AU < luminosity / SUN_LUMINOSITY * 200_f64 {
            if thread_rng().gen_range(0, 5) == 0 {
                Type::Gaseous
            } else {
                Type::Rocky
            }
        } else {
            Type::Rocky
        }
    }

    /// Generate atmosphere
    ///
    /// Generates a random atmosphere that can be mostly nitrogen, CO₂ or oxygen.
    fn generate_atmosphere(gravity: f64) -> Atmosphere {
        let pressure = if gravity < 5_f64 && thread_rng().gen_range(0, 5) == 0 {
            thread_rng().gen_range(0_f64, 0.01_f64) * gravity
        } else if gravity < 5_f64 && thread_rng().gen_range(0, 2) == 0 {
            thread_rng().gen_range(0_f64, 10_000_f64) * gravity / 10_f64
        } else if gravity < 15_f64 {
            thread_rng().gen_range(0_f64, 200_000_000_f64) * gravity / 10_f64
        } else if thread_rng().gen_range(0, 5) == 0 {
            thread_rng().gen_range(0_f64, 15_000_000_f64)
        } else {
            thread_rng().gen_range(0_f64, 3_000_000_f64)
        };

        let mut left = 1_f64;
        let mut co2 = 0_f64;
        let mut n2 = 0_f64;
        let mut o2 = 0_f64;

        if pressure > 0_f64 && thread_rng().gen_range(0, 2) == 0 {
            co2 = thread_rng().gen_range(0.75_f64, 0.99_f64);
            left -= co2;
            n2 = thread_rng().gen_range(0_f64, left);
            left -= n2;
            o2 = thread_rng().gen_range(0_f64, left);
            left -= o2;
        } else if pressure > 0_f64 {
            n2 = thread_rng().gen_range(0.5_f64, 0.95_f64);
            left -= n2;
            if thread_rng().gen_range(0, 2) == 0 {
                co2 = thread_rng().gen_range(0.004_f64, left);
                left -= co2;
                o2 = thread_rng().gen_range(0_f64, left);
                left -= o2;
            } else {
                o2 = thread_rng().gen_range(0.004_f64, left);
                left -= o2;
                co2 = thread_rng().gen_range(0_f64, left);
                left -= co2;
            }
        }

        let ar = thread_rng().gen_range(0_f64, left);
        left -= ar;
        let ne = thread_rng().gen_range(0_f64, left);
        left -= ne;
        let co = thread_rng().gen_range(0_f64, left);
        left -= co;
        let so2 = thread_rng().gen_range(0_f64, left);
        left -= so2;
        let ch4 = thread_rng().gen_range(0_f64, left);
        left -= ch4;
        let he = thread_rng().gen_range(0_f64, left);
        left -= he;
        let h2o = left;

        Atmosphere {
            pressure,
            h2o,
            co2,
            co,
            n2,
            o2,
            ar,
            so2,
            ne,
            ch4,
            he,
        }
    }

    /// Generate properties
    ///
    /// This function generates the basic properties ob the planet. The mass and the radius.
    fn generate_properties(planet_type: Type) -> (f64, f64) {
        match planet_type {
            Type::Rocky => {
                let radius = thread_rng().gen_range(2e+6_f64, 15e+6_f64); // m

                let density = if radius < 75e+5_f64 {
                    thread_rng().gen_range(1_500_f64, 6_000_f64) // kg/m³
                } else {
                    thread_rng().gen_range(5_000_f64, 13_000_f64) // kg/m³
                };

                let mut mass = 4_f64 * PI * radius.powi(3) * density / 3_f64; // kg
                if mass > 2e+25_f64 {
                    mass = thread_rng().gen_range(9e+24_f64, 2e+25_f64) // kg
                }

                (mass, radius)
            }
            Type::Gaseous => {
                let radius = thread_rng().gen_range(2e+7_f64, 1.5e+8_f64); // m

                let mut mass = (radius / 1e+3_f64).powf(1.3) * 1.445e+21 - 5e+26; // kg
                mass = thread_rng().gen_range(mass / 5_f64, mass * 5_f64);

                if mass > 1e+28 && thread_rng().gen_range(0, 3001) != 0 {
                    mass /= 10_f64;
                }

                (mass, radius)
            }
        }
    }

    /// Generate surface
    ///
    /// This function generates properties of the surface of the planet. The first time will be
    /// random, while the next ones will receive the minimum temperature, the average temperature
    /// and the maximum temperature, along with the atmosphere, to check the pressure in the water
    /// phase diagram.
    fn generate_surface(base: Option<(f64, f64, &Atmosphere)>) -> Surface {
        if base.is_some() {
            let (min_temp, max_temp, atm) = base.unwrap();
            let mut left = 1_f64;

            let ocean = if can_water_be_liquid(min_temp, max_temp, atm.pressure()) {
                if thread_rng().gen_range(0, 3) != 0 {
                    0_f64
                } else if thread_rng().gen_range(0, 5) == 0 {
                    1_f64
                } else {
                    thread_rng().gen_range(0_f64, 1_f64)
                }
            } else {
                0_f64
            };
            left -= ocean;

            let snow = if can_water_be_ice(min_temp, atm.pressure()) && left > 0_f64 {
                if ocean > 0_f64 {
                    thread_rng().gen_range(0_f64, left)
                } else if thread_rng().gen_range(0, 2) == 0 {
                    thread_rng().gen_range(0_f64, 1_f64)
                } else {
                    thread_rng().gen_range(0.9_f64, 1_f64)
                }
            } else {
                0_f64
            };
            left -= snow;

            let land = if left > 0_f64 {
                if ocean > 0_f64 {
                    thread_rng().gen_range(left - 0.05_f64, left)
                } else {
                    left
                }
            } else {
                0_f64
            };

            let fresh_water = if ocean > 0_f64 && left > 0_f64 {
                thread_rng().gen_range(0_f64, if left > 0.1_f64 { 0.1_f64 } else { left })
            } else {
                0_f64
            };

            Surface::new(fresh_water, ocean, snow, land)
        } else {
            Surface::new(0.0177, 0.6903, 0.0584, 0.2336)
        }
    }

    // ----------  calculators  ---------

    fn calculate_surface_gravity(mass: f64, radius: f64) -> f64 {
        G * mass / radius.powi(2)
    }

    /// Calculate atmosphere albedo
    ///
    /// This function calculates the albedo produced by the atmosphere. It will not be the final
    /// albedo, since it needs the surface albedo to do the final calculation depending on the
    /// atmospheric pressure.
    fn calculate_atmosphere_albedo(atm: &Atmosphere) -> (f64, f64) {
        let bond = (atm.co2() + atm.co()) * thread_rng().gen_range(0.5_f64, 0.7_f64)
            + atm.n2() * thread_rng().gen_range(0.2_f64, 0.3_f64)
            + atm.ch4() * thread_rng().gen_range(0.15_f64, 0.25_f64)
            + atm.o2() * thread_rng().gen_range(0.25_f64, 0.5_f64)
            + atm.h2o() * thread_rng().gen_range(0.2_f64, 0.5_f64)
            + (atm.so2() + atm.ne() + atm.he() + atm.ar()) * thread_rng().gen_range(0_f64, 0.9_f64);

        (bond, thread_rng().gen_range(bond * 0.85, bond * 1.15))
    }

    /// Calculate surface albedo
    ///
    /// This function calculates the albedo produced by the surface. It will not be the final albedo
    /// since it needs the atmospheric albedo and the atmospheric pressure to be able to calculate
    /// the final one.
    fn calculate_surface_albedo(surface: &Surface) -> (f64, f64) {
        let geom = (surface.ocean_water() + surface.fresh_water())
            * (thread_rng().gen_range(0.05_f64, 0.15_f64))
            * surface.snow()
            * (thread_rng().gen_range(0.5_f64, 1.5_f64))
            * surface.land()
            * (thread_rng().gen_range(0.1_f64, 0.6_f64));

        let bond = (surface.ocean_water() + surface.fresh_water())
            * (thread_rng().gen_range(0.1_f64, 0.2_f64))
            * surface.snow()
            * (thread_rng().gen_range(0.6_f64, 0.999_f64))
            * surface.land()
            * (thread_rng().gen_range(0.05_f64, 0.4_f64));
        (bond, geom)
    }

    /// Calculate final albedo
    ///
    /// This function calculates the final albedo for the body. This albedo will be calculated
    /// depending on the contribution of the atmosphere to the final albedo.
    fn calculate_albedo(
        planet_type: Type,
        atm: Option<&Atmosphere>,
        surface: Option<&Surface>,
    ) -> (f64, f64) {
        match planet_type {
            Type::Rocky => {
                let (surface_bond, surface_geom) = Self::calculate_surface_albedo(surface.unwrap());
                let (atm_bond, atm_geom) = Self::calculate_atmosphere_albedo(atm.unwrap());

                let pressure = atm.unwrap().pressure();
                if pressure <= 250_f64 {
                    (surface_bond, surface_geom)
                } else if pressure <= 750_f64 {
                    (
                        surface_bond * 0.75 + atm_bond * 0.25,
                        surface_geom * 0.75 + atm_geom * 0.25,
                    )
                } else if pressure <= 1_250_f64 {
                    (
                        surface_bond * 0.50 + atm_bond * 0.50,
                        surface_geom * 0.50 + atm_geom * 0.50,
                    )
                } else if pressure <= 2_000_f64 {
                    (
                        surface_bond * 0.25 + atm_bond * 0.75,
                        surface_geom * 0.25 + atm_geom * 0.75,
                    )
                } else {
                    (atm_bond, atm_geom)
                }
            }
            Type::Gaseous => {
                let bond = thread_rng().gen_range(0.25_f64, 0.4_f64);
                let geometric = thread_rng().gen_range(0.35_f64, 0.55_f64);

                (bond, geometric)
            }
        }
    }

    /// Calculate effective temperature
    ///
    /// This function calculates the effective temperature of the planet using real calculations
    /// given the star, the semimajor axis of the orbit and the albedo of the planet.
    fn calculate_t_eff(st: &Star, sma: f64, albedo: f64) -> f64 {
        let solar_constant = st
            .calculate_luminosity()
            .expect("black holes should not have planets")
            / (4_f64 * PI * sma.powi(2)); // W/m²

        (solar_constant * (1_f64 - albedo) / (4_f64 * BOLTZ)).powf(1_f64 / 4_f64)
    }

    /// Calculate greenhouse effect
    ///
    /// This function calculates the greenhouse effect of the planet. It uses an approximation
    /// taking into account the atmosphere composition and pressure.
    fn calculate_greenhouse(atm: Option<&Atmosphere>) -> f64 {
        let atmosphere = atm.unwrap();

        1_f64
            + (atmosphere.co2().powi(6) / 835_f64
                + atmosphere.h2o().sqrt() / 250_f64
                + atmosphere.ch4().powf(0.25_f64) / 1_000_f64)
                * atmosphere.pressure().sqrt()
    }

    /// Calculate surface temperature
    ///
    /// This function calculates the maximum and minimum surface temperatures of the planet. It uses
    /// an approximation taking into account the orbit and the atmospheric pressure.
    fn calculate_surface_temp(avg_temp: f64, atm_pressure: f64, orbit: &Orbit) -> (f64, f64) {
        let f1 = orbit.calculate_day().powf(1_f64 / 8_f64) / 8.25_f64;
        let min_temp = avg_temp
            * (1_f64
                - (if f1 > 1_f64 { 0.99999_f64 } else { f1 })
                    * (1_f64 - atm_pressure.powf(1_f64 / 3_f64) / 300_f64)
                    * (1_f64 - orbit.ecc())
                    * ((orbit.ax_tilt() - PI) / PI).abs());

        let max_temp = avg_temp
            * (1_f64
                + orbit.calculate_day().powf(1_f64 / 2.1_f64)
                    * (1_f64 - (atm_pressure * 1000_f64).powf(1_f64 / 2.3_f64) / 28_000_f64)
                    * (1_f64 - orbit.ecc().powi(4))
                    * ((orbit.ax_tilt() - PI) / PI).abs()
                    / 1_150_f64);

        (min_temp, max_temp)
    }
}

impl Orbit {
    /// Get position
    ///
    /// Gets the position of the orbit in the solar system. For example, for Earth would be 3.
    pub fn position(&self) -> u8 {
        self.position
    }

    /// Get eccentricity
    ///
    /// Gets the eccentricity of the orbit. Since all planets will be in closed orbits, the
    /// eccentricity will be between 0 and 1.
    pub fn ecc(&self) -> f64 {
        self.ecc
    }

    /// Get semimajor axis
    ///
    /// Gets the semimajor axis of the orbit, in meters (*m*).
    pub fn sma(&self) -> f64 {
        self.sma
    }

    /// Get inclination
    ///
    /// Gets the inclination of the orbit, in radians (*rad*).
    pub fn incl(&self) -> f64 {
        self.incl
    }

    /// Get longitude of ascending node
    ///
    /// Gets the longitude of the ascending node of the orbit, in radians (*rad*).
    pub fn lan(&self) -> f64 {
        self.lan
    }

    /// Get argument of periapsis
    ///
    /// Gets the argument of the periapsis of the orbit, in radians (*rad*).
    pub fn arg_p(&self) -> f64 {
        self.arg_p
    }

    /// Get mean anomaly
    ///
    /// Gets the mean anomaly of the orbit at the beginning of the universe, in radians (*rad*).
    pub fn anomaly(&self) -> f64 {
        self.m0
    }

    /// Get orbital period
    ///
    /// Gets the period of the orbit, in seconds (*s*).
    pub fn orb_period(&self) -> f64 {
        self.period
    }

    /// Get apoapsis
    ///
    /// Gets the apoapsis of the orbit, in meters (*m*).
    pub fn apoapsis(&self) -> f64 {
        self.sma * (1_f64 + self.ecc)
    }

    /// Get periapsis
    ///
    /// Gets the periapsis of the orbit, in meters (*m*).
    pub fn periapsis(&self) -> f64 {
        self.sma * (1_f64 - self.ecc)
    }

    /// Get axial tilt
    ///
    /// Gets the axial tilt of the rotation of the body, in radians (*rad*).
    pub fn ax_tilt(&self) -> f64 {
        self.ax_tilt
    }

    /// Get rotation period
    ///
    /// Gets the sidereal rotation period of the body, in seconds (*s*).
    pub fn rot_period(&self) -> f64 {
        self.rot_period
    }

    /// Calculate day duration
    ///
    /// Calculates the day duration of the body, in seconds (*s*).
    pub fn calculate_day(&self) -> f64 {
        if (self.ax_tilt > 1.308997_f64 && self.ax_tilt < 1.832596_f64)
            || (self.ax_tilt > 4.450589_f64 && self.ax_tilt < 4.974188_f64)
        {
            self.orb_period()
        } else if self.rot_period > 0_f64 {
            if (self.orb_period() - self.rot_period).abs() < self.orb_period() * 0.001 {
                0_f64
            } else {
                self.rot_period / (1_f64 - self.rot_period / self.orb_period())
            }
        } else if self.rot_period < 0_f64 {
            self.rot_period.abs() / (1_f64 + self.rot_period.abs() / self.orb_period())
        } else {
            unreachable!()
        }
    }
}

impl Atmosphere {
    /// Get pressure
    ///
    /// Gets the pressure of the atmosphere in Pascals (*Pa*).
    pub fn pressure(&self) -> f64 {
        self.pressure
    }

    /// Get H₂O
    ///
    /// Gets the percentage of H₂O (water vapour) in the atmosphere, from 0 to 1.
    pub fn h2o(&self) -> f64 {
        self.h2o
    }

    /// Get CO₂
    ///
    /// Gets the percentage of CO₂ (carbon dioxide) in the atmosphere, from 0 to 1.
    pub fn co2(&self) -> f64 {
        self.co2
    }

    /// Get CO
    ///
    /// Gets the percentage of CO (carbon monoxide) in the atmosphere, from 0 to 1.
    pub fn co(&self) -> f64 {
        self.co
    }

    /// Get N₂
    ///
    /// Gets the percentage of N₂ (nitrogen) in the atmosphere, from 0 to 1.
    pub fn n2(&self) -> f64 {
        self.n2
    }

    /// Get O₂
    ///
    /// Gets the percentage of O₂ (oxygen) in the atmosphere, from 0 to 1.
    pub fn o2(&self) -> f64 {
        self.o2
    }

    /// Get Ar
    ///
    /// Gets the percentage of Ar (argon) in the atmosphere, from 0 to 1.
    pub fn ar(&self) -> f64 {
        self.ar
    }

    /// Get SO₂
    ///
    /// Gets the percentage of SO₂ (sulfur dioxide) in the atmosphere, from 0 to 1.
    pub fn so2(&self) -> f64 {
        self.so2
    }

    /// Get Ne
    ///
    /// Gets the percentage of Ne (neon) in the atmosphere, from 0 to 1.
    pub fn ne(&self) -> f64 {
        self.ne
    }

    /// Get CH₄
    ///
    /// Gets the percentage of CH₄ (methane) in the atmosphere, from 0 to 1.
    pub fn ch4(&self) -> f64 {
        self.ch4
    }

    /// Get He
    ///
    /// Gets the percentage of He (helium) in the atmosphere, from 0 to 1.
    pub fn he(&self) -> f64 {
        self.he
    }
}

impl Surface {
    /// Constructs a new `Surface` structure.
    ///
    /// It creates a new surface structure with the needed information for representing the surface
    /// of the planet.
    fn new(fresh_water: f64, ocean_water: f64, snow: f64, land: f64) -> Self {
        Self {
            fresh_water,
            ocean_water,
            snow,
            land,
        }
    }

    /// Get fresh water
    ///
    /// Gets the percentage of fresh water on the surface of the planet. From 0 to 1.
    pub fn fresh_water(&self) -> f64 {
        self.fresh_water
    }

    /// Get ocean water
    ///
    /// Gets the percentage of ocean water on the surface of the planet. From 0 to 1.
    pub fn ocean_water(&self) -> f64 {
        self.ocean_water
    }

    /// Get snow
    ///
    /// Gets the percentage of the surface of the planet covered in snow. From 0 to 1.
    pub fn snow(&self) -> f64 {
        self.snow
    }

    /// Get land
    ///
    /// Gets the percentage of land on the surface of the planet. From 0 to 1.
    pub fn land(&self) -> f64 {
        self.land
    }
}

#[cfg(test)]
mod tests {
    use super::{super::star::Star, Atmosphere, Orbit, Planet, Surface, Type};
    use crate::consts::EARTH_ATM_PRESSURE;
    use std::f64::EPSILON;

    #[test]
    fn it_orbit_getters() {
        let orb = Orbit {
            position: 3,
            ecc: 0.5,
            sma: 150e+9,
            incl: 1.5,
            lan: 1.2,
            arg_p: 1.3,
            m0: 1.4,
            period: 31_558_118.4,
            ax_tilt: 1.1,
            rot_period: 80_600_f64,
        };

        assert!(orb.ecc() >= 0.5 - EPSILON && orb.ecc() <= 0.5 + EPSILON);
        assert!(orb.sma() >= 150e+9 - EPSILON && orb.sma() <= 150e+9 + EPSILON);
        assert!(orb.incl() >= 1.5 - EPSILON && orb.incl() <= 1.5 + EPSILON);
        assert!(orb.lan() >= 1.2 - EPSILON && orb.lan() <= 1.2 + EPSILON);
        assert!(orb.arg_p() >= 1.3 - EPSILON && orb.arg_p() <= 1.3 + EPSILON);
        assert!(orb.anomaly() >= 1.4 - EPSILON && orb.anomaly() <= 1.4 + EPSILON);
        assert!(
            orb.orb_period() >= 31_558_118.4 - EPSILON
                && orb.orb_period() <= 31_558_118.4 + EPSILON
        );
        assert!(orb.ax_tilt() >= 1.1 - EPSILON && orb.ax_tilt() <= 1.1 + EPSILON);
        assert!(
            orb.rot_period() >= 80_600_f64 - EPSILON && orb.rot_period() <= 80_600_f64 + EPSILON
        );
    }

    #[test]
    fn it_atm_getters() {
        let atm = Atmosphere {
            pressure: EARTH_ATM_PRESSURE,
            h2o: 0.01,
            co2: 0.0397,
            co: 0_f64,
            n2: 78.084,
            o2: 20.946,
            ar: 0.9340,
            so2: 0.1,
            ne: 0.00181,
            ch4: 0.00017,
            he: 0.00052,
        };

        assert!(
            atm.pressure() >= EARTH_ATM_PRESSURE - EPSILON
                && atm.pressure() <= EARTH_ATM_PRESSURE + EPSILON
        );
        assert!(atm.h2o() >= 0.01 - EPSILON && atm.h2o() <= 0.01 + EPSILON);
        assert!(atm.co2() >= 0.0397 - EPSILON && atm.co2() <= 0.0397 + EPSILON);
        assert!(atm.co() >= 0_f64 - EPSILON && atm.co() <= 0_f64 + EPSILON);
        assert!(atm.n2() >= 78.084 - EPSILON && atm.n2() <= 78.084 + EPSILON);
        assert!(atm.o2() >= 20.946 - EPSILON && atm.o2() <= 20.946 + EPSILON);
        assert!(atm.ar() >= 0.9340 - EPSILON && atm.ar() <= 0.9340 + EPSILON);
        assert!(atm.so2() >= 0.1 - EPSILON && atm.so2() <= 0.1 + EPSILON);
        assert!(atm.ne() >= 0.00181 - EPSILON && atm.ne() <= 0.00181 + EPSILON);
        assert!(atm.ch4() >= 0.00017 - EPSILON && atm.ch4() <= 0.00017 + EPSILON);
        assert!(atm.he() >= 0.00052 - EPSILON && atm.he() <= 0.00052 + EPSILON);
    }

    #[test]
    fn it_surface_getters() {
        let surface = Surface::new(0.0177, 0.6903, 0.0584, 0.2336);

        assert!(
            surface.fresh_water() >= 0.0177 - EPSILON && surface.fresh_water() <= 0.0177 + EPSILON
        );
        assert!(
            surface.ocean_water() >= 0.6903 - EPSILON && surface.ocean_water() <= 0.6903 + EPSILON
        );
        assert!(surface.snow() >= 0.0584 - EPSILON && surface.snow() <= 0.0584 + EPSILON);
        assert!(surface.land() >= 0.2336 - EPSILON && surface.land() <= 0.2336 + EPSILON);
    }

    #[test]
    fn it_planet_getters() {
        let star = Star::new();
        let orb = Orbit {
            position: 3,
            ecc: 0.5,
            sma: 150e+9,
            incl: 1.5,
            lan: 1.2,
            arg_p: 1.3,
            m0: 1.4,
            period: 31_558_118.4,
            ax_tilt: 1.1,
            rot_period: 80_600_f64,
        };
        let atm = Atmosphere {
            pressure: EARTH_ATM_PRESSURE,
            h2o: 0.01,
            co2: 0.0397,
            co: 0_f64,
            n2: 78.084,
            o2: 20.946,
            ar: 0.9340,
            so2: 0.1,
            ne: 0.00181,
            ch4: 0.00017,
            he: 0.00052,
        };
        let surface = super::Surface::new(0.0177, 0.6903, 0.0584, 0.2336);

        let planet = Planet {
            orbit: orb,
            atmosphere: Some(atm),
            surface: Some(surface),
            planet_type: Type::Rocky,
            bond_albedo: 0.306,
            geometric_albedo: 0.367,
            mass: 5.9726e+24,
            radius: 6.371e+6,
            eff_temp: 254.3367460856,
            min_temp: 183.95,
            max_temp: 329.85,
            avg_temp: 289.15,
        };

        assert!(
            planet.atmosphere().unwrap().pressure() >= EARTH_ATM_PRESSURE - EPSILON
                && planet.atmosphere().unwrap().pressure() <= EARTH_ATM_PRESSURE + EPSILON
        );
        assert!(
            planet.surface().unwrap().ocean_water() >= 0.6903 - EPSILON
                && planet.surface().unwrap().ocean_water() <= 0.6903 + EPSILON
        );
        assert_eq!(&Type::Rocky, planet.planet_type());
        assert!(planet.bond_albedo() >= 0.306 - EPSILON && planet.bond_albedo() <= 0.306 + EPSILON);
        assert!(
            planet.geometric_albedo() >= 0.367 - EPSILON
                && planet.geometric_albedo() <= 0.367 + EPSILON
        );
        assert!(planet.mass() >= 5.9726e+24 - EPSILON && planet.mass() <= 5.9726e+24 + EPSILON);
        assert!(planet.radius() >= 6.371e+6 - EPSILON && planet.radius() <= 6.371e+6 + EPSILON);
        assert!(planet.is_roche_ok(&star));
        assert!(planet.min_temp() >= 183.95 - EPSILON && planet.min_temp() <= 183.95 + EPSILON);
        assert!(planet.max_temp() >= 329.85 - EPSILON && planet.max_temp() <= 329.85 + EPSILON);
        assert!(planet.avg_temp() >= 289.15 - EPSILON && planet.avg_temp() <= 289.15 + EPSILON);
    }

    #[test]
    fn it_volume() {
        let orb = Orbit {
            position: 3,
            ecc: 0.5,
            sma: 150e+9,
            incl: 1.5,
            lan: 1.2,
            arg_p: 1.3,
            m0: 1.4,
            period: 31_558_118.4,
            ax_tilt: 1.1,
            rot_period: 80_600_f64,
        };
        let atm = Atmosphere {
            pressure: EARTH_ATM_PRESSURE,
            h2o: 0.01,
            co2: 0.0397,
            co: 0_f64,
            n2: 78.084,
            o2: 20.946,
            ar: 0.9340,
            so2: 0.1,
            ne: 0.00181,
            ch4: 0.00017,
            he: 0.00052,
        };
        let surface = super::Surface::new(0.0177, 0.6903, 0.0584, 0.2336);

        let planet = Planet {
            orbit: orb,
            atmosphere: Some(atm),
            surface: Some(surface),
            planet_type: Type::Rocky,
            bond_albedo: 0.306,
            geometric_albedo: 0.367,
            mass: 5.9726e+24,
            radius: 6.371e+6,
            eff_temp: 254.3367460856,
            min_temp: 183.95,
            max_temp: 329.85,
            avg_temp: 289.15,
        };

        assert!(
            10.8321e+20 * 0.999 < planet.calculate_volume()
                && 10.8321e+20 * 1.001 > planet.calculate_volume()
        );
    }

    #[test]
    fn it_density() {
        let orb = Orbit {
            position: 3,
            ecc: 0.5,
            sma: 150e+9,
            incl: 1.5,
            lan: 1.2,
            arg_p: 1.3,
            m0: 1.4,
            period: 31_558_118.4,
            ax_tilt: 1.1,
            rot_period: 80_600_f64,
        };
        let atm = Atmosphere {
            pressure: EARTH_ATM_PRESSURE,
            h2o: 0.01,
            co2: 0.0397,
            co: 0_f64,
            n2: 78.084,
            o2: 20.946,
            ar: 0.9340,
            so2: 0.1,
            ne: 0.00181,
            ch4: 0.00017,
            he: 0.00052,
        };
        let surface = Surface::new(0.0177, 0.6903, 0.0584, 0.2336);

        let planet = Planet {
            orbit: orb,
            atmosphere: Some(atm),
            surface: Some(surface),
            planet_type: Type::Rocky,
            bond_albedo: 0.306,
            geometric_albedo: 0.367,
            mass: 5.9726e+24,
            radius: 6.371e+6,
            eff_temp: 254.3367460856,
            min_temp: 183.95,
            max_temp: 329.85,
            avg_temp: 289.15,
        };

        assert!(
            5_514_f64 * 0.999 < planet.calculate_density()
                && 5_514_f64 * 1.001 > planet.calculate_density()
        );
    }
}
