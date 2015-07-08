//! Planet module
//!
//! This is the `planet` module. This module contains all the structures, enumerations and
//! implementations needed to define a planet.

extern crate rand;

use std::f64::consts::PI;
use std::fmt;

use self::rand::Rng;

use consts::*;
use star::Star;

/// Planet structure
///
/// This structure defines a planet and contains all the structures needed for a correct definition
/// and representation of itself.
pub struct Planet<'p> {
    orbit: Orbit<'p>,
    atmosphere: Option<Atmosphere>,
    planet_type: PlanetType,
    // soil: Soil,
    // life: Life,
    albedo: f64,
    mass: f64,
    radius: f64,
    min_temp: f64,
    max_temp: f64,
    avg_temp: f64
    // habitable: bool,
    // double_planet: bool,
}

/// PlanetType enumeration
///
/// This enumeration contains the planet types available when creating a new planet.
#[derive(PartialEq)]
pub enum PlanetType {
    /// This represents a rocky planet like Earth. Y has a surface where objects and even life can
    /// be found.
    Rocky,
    /// This represents a gaseous planet, like Jupiter. It has no known surface, and in the case of
    /// having it it would be so deep that it supposes nothing for the atmospheric configuration or
    /// habitability.
    Gaseous,
}

impl fmt::Debug for PlanetType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            PlanetType::Rocky => write!(f, "rocky"),
            PlanetType::Gaseous => write!(f, "gaseous"),
        }
    }
}

/// Orbit structure
///
/// This structure contains all the needed information to define an orbit, and as an extra, it
/// contains the axial tilt of the rotation of the body and the rotation period of the body. It also
/// includes a reference to the star being orbited. It also contains the orbital period, even if it
/// can be calculated with the rest of the parameters.
pub struct Orbit<'o> {
    star: &'o Star,
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
pub struct Atmosphere {
    pressure: f64,
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

// struct Soil {
// }

// struct Life {
// }

impl<'p>  Planet<'p> {
    /// Constructs a new `Planet`.
    ///
    /// It creates a random planet taking into account real planet statistics. It requires the
    /// reference to parent star and the two values of the Titius–Bode law, along with the order
    /// of the planet in the solar system and the last body's semi-major axis in the solar system.
    ///
    /// # Examples
    ///
    /// ```
    /// use star::Star;
    /// use planet::Planet;
    ///
    /// let st = Star::new(0, 1);
    ///
    /// let num_bodies = star.calculate_num_bodies();
    /// let (tb_m, tb_n) = star.calculate_titius_bode(num_bodies);
    ///
    /// if num_bodies > 0 {
    ///     let planet = Planet::new(&star, tb_m, tb_n, 1);
    /// }
    /// ```
    pub fn new(st: &'p Star, m: f64, n: f64, position: u8, last_sm_a: f64) -> Planet {
        let orb = Planet::generate_orbit(st, m, n, position, last_sm_a);

        let planet_type = Planet::generate_type(orb.get_sma(), st.get_luminosity());
        let atm = if planet_type == PlanetType::Rocky {Some(Planet::generate_atmosphere())} else {None};
        let (mass, radius) = Planet::generate_properties(&planet_type);
        let alb = 0.306;
        let min_temp = 0_f64;
        let max_temp = 0_f64;
        let avg_temp = 0_f64;

        Planet {orbit: orb, atmosphere: atm, planet_type: planet_type, albedo: alb, mass: mass,
            radius: radius, min_temp: min_temp, max_temp: max_temp, avg_temp: avg_temp}
    }

    /// Get `Orbit`
    ///
    /// Gets the orbit information of the planet.
    pub fn get_orbit(&self) -> &Orbit {
        &self.orbit
    }

    /// Get `Atmosphere`
    ///
    /// Gets the atmosphere information of the planet.
    pub fn get_atmosphere(&self) -> Option<&Atmosphere> {
        self.atmosphere.as_ref()
    }

    /// Get planet type
    ///
    /// Gets the type of the planet.
    pub fn get_type(&self) -> &PlanetType {
        &self.planet_type
    }

    /// Get albedo
    ///
    /// Gets the albedo of the planet. The albedo is the reflectivity of the planet, or in other
    /// words, the percentage of light reflected by the planet, from 0 to 1.
    pub fn get_albedo(&self) -> f64 {
        self.albedo
    }

    /// Get mass
    ///
    /// Gets the mass of the planet in *kg*.
    pub fn get_mass(&self) -> f64 {
        self.mass
    }

    /// Get radius
    ///
    /// Gets the radius of the planet, in meters (*m*).
    pub fn get_radius(&self) -> f64 {
        self.radius
    }

    /// Get density
    ///
    /// Gets the density of the planet, in *kg/m³*.
    pub fn get_density(&self) -> f64 {
        self.mass/self.get_volume()
    }

    /// Get volume
    ///
    /// Gets the volume of the planet, in *m³*.
    pub fn get_volume(&self) -> f64 {
        4_f64/3_f64*PI*self.radius.powi(3)
    }

    /// Get minimum temperature
    ///
    /// Gets the minimum temperature of the planet, in Kelvin (*K*).
    pub fn get_min_temp(&self) -> f64 {
        self.min_temp
    }

    /// Get maximum temperature
    ///
    /// Gets the maximum temperature of the planet, in Kelvin (*K*).
    pub fn get_max_temp(&self) -> f64 {
        self.max_temp
    }

    /// Get average temperature
    ///
    /// Gets the average temperature of the planet, in Kelvin (*K*).
    pub fn get_avg_temp(&self) -> f64 {
        self.avg_temp
    }

    /// Generate orbit
    ///
    /// Generates the orbit of the planet taking into account the Titius-Bode law, the last planet's
    /// semimajor axis and the position in the system.
    fn generate_orbit(st: &Star, m: f64, n: f64, position: u8, last_sma: f64) -> Orbit {
        let mut sma = (m*(position as f64) - n).exp()*AU;
        sma = rand::thread_rng().gen_range(sma*0.9, sma*1.1);

        if sma < last_sma*1.15 {
            sma = last_sma*rand::thread_rng().gen_range(1.15_f64, 1.25_f64);
        }

        let ecc = if sma/AU < st.get_mass()/(SUN_MASS*2_f64) {
            rand::thread_rng().gen_range(0.05_f64, 0.25_f64)
        } else {
            rand::thread_rng().gen_range(0_f64, 0.1_f64)
        };

        let period = 2_f64*PI*(sma.powi(3)/(G*st.get_mass())).sqrt();

        let incl = rand::thread_rng().gen_range(0_f64, PI/18_f64);
        let lan = rand::thread_rng().gen_range(0_f64, 2_f64*PI);
        let arg_p = rand::thread_rng().gen_range(0_f64, 2_f64*PI);
        let m0 = rand::thread_rng().gen_range(0_f64, 2_f64*PI);

        let (ax_tilt, rot_period) = Planet::generate_rotation(st, sma, period);

        Orbit::new(st, position, ecc, sma, incl, lan, arg_p, m0, period, ax_tilt, rot_period)
    }

    /// Generate rotation
    ///
    /// Generates the rotation of the planet taking into account tidal lock and orbital resonance.
    fn generate_rotation(st: &Star, sma: f64, orb_period: f64) -> (f64, f64) {
        let tidal_lock =  (st.get_mass()/SUN_MASS).sqrt()/2_f64;

        if sma/AU > tidal_lock {
            let ax_tilt = if rand::thread_rng().gen_range(0, 1) != 0 {
                rand::thread_rng().gen_range(0.349_f64, 0.5236_f64) // 20° - 30°
            } else {
                rand::thread_rng().gen_range(0_f64, PI)
            };

            let rot_period = if ax_tilt > PI/2_f64 {
                if orb_period < 50_000_f64 {
                    -rand::thread_rng().gen_range(orb_period*0.8, orb_period-1_f64)
                } else {
                    -rand::thread_rng().gen_range(50_000_f64, if orb_period < 25_000_000_f64
                        {orb_period-1_f64} else {25_000_000_f64})
                }
            } else {
                if orb_period < 18_000_f64 {
                    rand::thread_rng().gen_range(orb_period*0.8, orb_period-1_f64)
                } else {
                    rand::thread_rng().gen_range(18_000_f64, if orb_period < 180_000_f64
                        {orb_period-1_f64} else {180_000_f64})
                }
            };

            (ax_tilt, rot_period)
        } else if sma > tidal_lock.sqrt()/3_f64 {
            let ax_tilt = rand::thread_rng().gen_range(0_f64, 0.017454_f64); // 0° - 1°
            let rot_period = orb_period*2_f64/(rand::thread_rng().gen_range(3, 6) as f64); // Resonance

            (ax_tilt, rot_period)
        } else { // Tidal lock
            (0_f64, orb_period)
        }
    }

    /// Generate planet type
    ///
    /// Generates the PlanetType of the planet depending on star and the SMa of the orbit.
    fn generate_type(sma: f64, luminosity: f64) -> PlanetType {
        if sma/(AU*2_f64) < (luminosity/SUN_LUMINOSITY).sqrt() {
            if rand::thread_rng().gen_range(0, 2) == 0 {PlanetType::Gaseous} else {PlanetType::Rocky}
        } else if luminosity > 1.923e+27_f64 && // 5*SUN_LUMINOSITY
            sma/AU < luminosity/SUN_LUMINOSITY.sqrt()*50_f64 {
            if rand::thread_rng().gen_range(0, 4) == 0 {PlanetType::Rocky} else {PlanetType::Gaseous}
        } else if sma/AU < luminosity/SUN_LUMINOSITY*200_f64 {
            if rand::thread_rng().gen_range(0, 2) == 0 {PlanetType::Gaseous} else {PlanetType::Rocky}
        } else {
            PlanetType::Rocky
        }
    }

    /// Generate atmosphere
    ///
    /// Generates a random atmosphere that can be mostly nitrogen, CO₂ or oxygen.
    fn generate_atmosphere() -> Atmosphere {
        let pressure = if rand::thread_rng().gen_range(0, 5) == 0 {
                rand::thread_rng().gen_range(0_f64, 1e-5_f64)
            } else if rand::thread_rng().gen_range(0, 2) == 0 {
                rand::thread_rng().gen_range(0_f64, 100_f64)
            } else if rand::thread_rng().gen_range(0, 5) != 0 {
                rand::thread_rng().gen_range(0_f64, 3_000_f64)
            } else {
                rand::thread_rng().gen_range(0_f64, 100_000_f64)
            };

        let mut left = 1_f64;
        let mut co2 = 0_f64;
        let mut n2 = 0_f64;
        let mut o2 = 0_f64;

        if pressure > 0_f64 && rand::thread_rng().gen_range(0, 2) == 0 {
            co2 = rand::thread_rng().gen_range(0.75_f64, 0.99_f64);
            left -= co2;
            n2 = rand::thread_rng().gen_range(0_f64, left);
            left -= n2;
            o2 = rand::thread_rng().gen_range(0_f64, left);
            left -= o2;
        } else if pressure > 0_f64 {
            n2 = rand::thread_rng().gen_range(0.5_f64, 0.95_f64);
            left -= n2;
            if rand::thread_rng().gen_range(0, 2) == 0 {
                co2 = rand::thread_rng().gen_range(0.004_f64, left);
                left -= co2;
                o2 = rand::thread_rng().gen_range(0_f64, left);
                left -= o2;
            } else {
                o2 = rand::thread_rng().gen_range(0.004_f64, left);
                left -= o2;
                co2 = rand::thread_rng().gen_range(0_f64, left);
                left -= co2;
            }

        }

        let ar = rand::thread_rng().gen_range(0_f64, left);
        left -= ar;
        let ne = rand::thread_rng().gen_range(0_f64, left);
        left -= ne;
        let co = rand::thread_rng().gen_range(0_f64, left);
        left -= co;
        let so2 = rand::thread_rng().gen_range(0_f64, left);
        left -= so2;
        let ch4 = rand::thread_rng().gen_range(0_f64, left);
        left -= ch4;
        let he = left;

        Atmosphere::new(pressure, co2, co, n2, o2, ar, so2, ne, ch4, he)

        // TODO: greenhouse effect
    }

    /// Generate properties
    ///
    /// This function generates the basic properties ob the planet. The mass and the radius.
    fn generate_properties(planet_type: &PlanetType) -> (f64, f64) {
        match *planet_type {
            PlanetType::Rocky => {
                let radius = rand::thread_rng().gen_range(2e+6_f64, 15e+6_f64); // m

                let density = if radius < 75e+5_f64 {
                    rand::thread_rng().gen_range(2_700_f64, 6_500_f64) // kg/m³
                } else {
                    rand::thread_rng().gen_range(5_500_f64, 15_000_f64) // kg/m³
                };

                let mut mass = 4_f64*PI*radius.powi(3)*density/3_f64; // kg
                if mass > 9e+25_f64 {
                    mass = rand::thread_rng().gen_range(5e+25_f64, 9e+25_f64) // kg
                }

                (mass, radius)
            },
            PlanetType::Gaseous => {
                let radius = rand::thread_rng().gen_range(2e+7_f64, 1.5e+8_f64); // m

                let mut mass = (radius/1e+3_f64).powf(1.3)*1.445e+21-5e+26; // kg
                mass = rand::thread_rng().gen_range(mass/5_f64, mass*5_f64);

                if mass > 1e+28 && rand::thread_rng().gen_range(0, 3001) != 0 {
                    mass /= 10_f64;
                }

                (mass, radius)
            }
        }
    }
}

impl<'o> Orbit<'o> {
    /// Constructs a new `Orbit`.
    ///
    /// It creates a new orbit structure with all the needed parameters for complete representation.
    fn new(star: &'o Star, position: u8, ecc: f64, sma: f64, incl: f64, lan: f64, arg_p: f64,
        m0: f64, period: f64, ax_tilt: f64, rot_period: f64) -> Orbit {
        Orbit {star: star, position: position, ecc: ecc, sma: sma, incl: incl, lan: lan,
            arg_p: arg_p, m0: m0, period: period, ax_tilt: ax_tilt, rot_period: rot_period}
    }

    /// Get `Star`
    ///
    /// Gets the Star being orbited.
    pub fn get_star(&self) -> &Star {
        self.star
    }

    /// Get position
    ///
    /// Gets the position of the orbit in the solar system. For example, for Earth would be 3.
    pub fn get_position(&self) -> u8 {
        self.position
    }

    /// Get eccentricity
    ///
    /// Gets the eccentricity of the orbit. Since all planets will be in closed orbits, the
    /// eccentricity will be between 0 and 1.
    pub fn get_ecc(&self) -> f64 {
        self.ecc
    }

    /// Get semimajor axis
    ///
    /// Gets the semimajor axis of the orbit, in meters (*m*).
    pub fn get_sma(&self) -> f64 {
        self.sma
    }

    /// Get inclination
    ///
    /// Gets the inclination of the orbit, in radians (*rad*).
    pub fn get_incl(&self) -> f64 {
        self.incl
    }

    /// Get longitude of ascending node
    ///
    /// Gets the longitude of the ascending node of the orbit, in radians (*rad*).
    pub fn get_lan(&self) -> f64 {
        self.lan
    }

    /// Get argument of periapsis
    ///
    /// Gets the argument of the periapsis of the orbit, in radians (*rad*).
    pub fn get_arg_p(&self) -> f64 {
        self.arg_p
    }

    /// Get mean anomaly
    ///
    /// Gets the mean anomaly of the orbit at the beginning of the universe, in radians (*rad*).
    pub fn get_anomaly(&self) -> f64 {
        self.m0
    }

    /// Get orbital period
    ///
    /// Gets the period of the orbit, in seconds (*s*).
    pub fn get_orb_period(&self) -> f64 {
        self.period
    }

    /// Get apoapsis
    ///
    /// Gets the apoapsis of the orbit, in meters (*m*).
    pub fn get_apoapsis(&self) -> f64 {
        self.sma*(1_f64+self.ecc)
    }

    /// Get periapsis
    ///
    /// Gets the periapsis of the orbit, in meters (*m*).
    pub fn get_periapsis(&self) -> f64 {
        self.sma*(1_f64-self.ecc)
    }

    /// Get axial tilt
    ///
    /// Gets the axial tilt of the rotation of the body, in radians (*rad*).
    pub fn get_ax_tilt(&self) -> f64 {
        self.ax_tilt
    }

    /// Get rotation period
    ///
    /// Gets the sidereal rotation period of the body, in seconds (*s*).
    pub fn get_rot_period(&self) -> f64 {
        self.rot_period
    }

    /// Get day
    ///
    /// Gets the day length of the body, in seconds (*s*).
    pub fn get_day(&self) -> f64 {
        if self.ax_tilt > 1.3962634 && self.ax_tilt < 1.3962634 {
            self.get_orb_period()
        } else {
            if self.rot_period > 0_f64 {
                if (self.get_orb_period() - self.rot_period).abs() < self.get_orb_period() * 0.001 {
                    0_f64
                } else {
                    self.rot_period/(1_f64-self.rot_period/self.get_orb_period())
                }
            } else if self.rot_period < 0_f64 {
                self.rot_period.abs()/(1_f64+self.rot_period.abs()/self.get_orb_period())
            } else {unreachable!()}
        }
    }
}

impl Atmosphere {
    /// Constructs a new `Atmosphere` structure.
    ///
    /// It creates a new atmosphere structure with all the percentages of the composition and its
    /// pressure.
    fn new(pressure: f64, co2: f64, co: f64, n2: f64, o2: f64, ar: f64, so2: f64,
        ne: f64, ch4: f64, he: f64) -> Atmosphere {
        Atmosphere {pressure: pressure, co2: co2, co: co, n2: n2, o2: o2, ar: ar,
            so2: so2, ne: ne, ch4: ch4, he: he}
    }

    /// Get pressure
    ///
    /// Gets the pressure of the atmosphere in Pascals (*Pa*).
    pub fn get_pressure(&self) -> f64 {
        self.pressure
    }

    /// Get CO₂
    ///
    /// Gets the percentage of CO₂ (carbon dioxide) in the atmosphere, from 0 to 1.
    pub fn get_co2(&self) -> f64 {
        self.co2
    }

    /// Get CO
    ///
    /// Gets the percentage of CO (carbon monoxide) in the atmosphere, from 0 to 1.
    pub fn get_co(&self) -> f64 {
        self.co
    }

    /// Get N₂
    ///
    /// Gets the percentage of N₂ (nitrogen) in the atmosphere, from 0 to 1.
    pub fn get_n2(&self) -> f64 {
        self.n2
    }

    /// Get O₂
    ///
    /// Gets the percentage of O₂ (oxygen) in the atmosphere, from 0 to 1.
    pub fn get_o2(&self) -> f64 {
        self.o2
    }

    /// Get Ar
    ///
    /// Gets the percentage of Ar (argon) in the atmosphere, from 0 to 1.
    pub fn get_ar(&self) -> f64 {
        self.ar
    }

    /// Get SO₂
    ///
    /// Gets the percentage of SO₂ (sulfur dioxide) in the atmosphere, from 0 to 1.
    pub fn get_so2(&self) -> f64 {
        self.so2
    }

    /// Get Ne
    ///
    /// Gets the percentage of Ne (neon) in the atmosphere, from 0 to 1.
    pub fn get_ne(&self) -> f64 {
        self.ne
    }

    /// Get CH₄
    ///
    /// Gets the percentage of CH₄ (methane) in the atmosphere, from 0 to 1.
    pub fn get_ch4(&self) -> f64 {
        self.ch4
    }

    /// Get He
    ///
    /// Gets the percentage of He (helium) in the atmosphere, from 0 to 1.
    pub fn get_he(&self) -> f64 {
        self.he
    }
}

#[cfg(test)]
mod tests {
    use super::Planet;
    use super::PlanetType;
    use super::super::star::Star;

    #[test]
    fn it_orbit_getters() {
        let st = Star::new(2, 0);

        let orb = super::Orbit::new(&st, 3, 0.5_f64, 150e+9_f64, 1.5_f64, 1.2_f64, 1.3_f64, 1.4_f64,
            31_558_118.4_f64, 1.1_f64, 80_600_f64);

        assert_eq!(3, orb.get_star().get_id());
        assert_eq!(0.5_f64, orb.get_ecc());
        assert_eq!(150e+9_f64, orb.get_sma());
        assert_eq!(1.5_f64, orb.get_incl());
        assert_eq!(1.2_f64, orb.get_lan());
        assert_eq!(1.3_f64, orb.get_arg_p());
        assert_eq!(1.4_f64, orb.get_anomaly());
        assert_eq!(31_558_118.4_f64, orb.get_orb_period());
        assert_eq!(1.1_f64, orb.get_ax_tilt());
        assert_eq!(80_600_f64, orb.get_rot_period());
    }

    #[test]
    fn it_atm_getters() {
        let atm = super::Atmosphere::new(101325_f64, 0.0397_f64, 0_f64, 78.084_f64, 20.946_f64,
            0.9340_f64, 0.1_f64, 0.00181_f64, 0.00017_f64, 0.00052_f64);

        assert_eq!(101325_f64,  atm.get_pressure());
        assert_eq!(0.0397_f64,  atm.get_co2());
        assert_eq!(0_f64,       atm.get_co());
        assert_eq!(78.084_f64,  atm.get_n2());
        assert_eq!(20.946_f64,  atm.get_o2());
        assert_eq!(0.9340_f64,  atm.get_ar());
        assert_eq!(0.1_f64,     atm.get_so2());
        assert_eq!(0.00181_f64, atm.get_ne());
        assert_eq!(0.00017_f64, atm.get_ch4());
        assert_eq!(0.00052_f64, atm.get_he());
    }

    #[test]
    fn it_parameter_test() {
        let st = Star::new(2, 0);
        let pl = Planet::new(&st, 0.0183, 1.0643, 3, 0_f64);

        assert_eq!(3, pl.get_orbit().get_star().get_id());
    }

    #[test]
    fn it_planet_getters() {
        let st = Star::new(4, 6);
        let orb = super::Orbit::new(&st, 3, 0.5_f64, 150e+9_f64, 1.5_f64, 1.2_f64, 1.3_f64, 1.4_f64,
            31_558_118.4_f64, 1.1_f64, 80_600_f64);
        let atm = super::Atmosphere::new(101325_f64, 0.0397_f64, 0_f64, 78.084_f64, 20.946_f64,
            0.9340_f64, 0.1_f64, 0.00181_f64, 0.00017_f64, 0.00052_f64);

        let planet = Planet {orbit: orb, atmosphere: Some(atm), planet_type: PlanetType::Rocky,
            albedo: 0.306_f64, mass: 5.9726e+24_f64, radius: 6.371e+6_f64, min_temp: 183.95_f64,
            max_temp: 329.85_f64, avg_temp: 289.15_f64};

        assert_eq!(5, planet.get_orbit().get_star().get_id());
        assert_eq!(101325_f64, planet.get_atmosphere().unwrap().get_pressure());
        assert_eq!(&PlanetType::Rocky, planet.get_type());
        assert_eq!(0.306_f64, planet.get_albedo());
        assert_eq!(5.9726e+24_f64, planet.get_mass());
        assert_eq!(6.371e+6_f64, planet.get_radius());
        assert_eq!(183.95_f64, planet.get_min_temp());
        assert_eq!(329.85_f64, planet.get_max_temp());
        assert_eq!(289.15_f64, planet.get_avg_temp());
    }

    #[test]
    fn it_volume() {
        let st = Star::new(4, 6);
        let orb = super::Orbit::new(&st, 3, 0.5_f64, 150e+9_f64, 1.5_f64, 1.2_f64, 1.3_f64, 1.4_f64,
            31_558_118.4_f64, 1.1_f64, 80_600_f64);
        let atm = super::Atmosphere::new(101325_f64, 0.0397_f64, 0_f64, 78.084_f64, 20.946_f64,
            0.9340_f64, 0.1_f64, 0.00181_f64, 0.00017_f64, 0.00052_f64);

        let planet = Planet {orbit: orb, atmosphere: Some(atm), planet_type: PlanetType::Rocky,
            albedo: 0.306_f64, mass: 5.9726e+24_f64, radius: 6.371e+6_f64, min_temp: 183.95_f64,
            max_temp: 329.85_f64, avg_temp: 289.15_f64};

        assert!(10.8321e+20_f64*0.999 < planet.get_volume() && 10.8321e+20_f64*1.001 > planet.get_volume());
    }

    #[test]
    fn it_density() {
        let st = Star::new(4, 6);
        let orb = super::Orbit::new(&st, 3, 0.5_f64, 150e+9_f64, 1.5_f64, 1.2_f64, 1.3_f64, 1.4_f64,
            31_558_118.4_f64, 1.1_f64, 80_600_f64);
        let atm = super::Atmosphere::new(101325_f64, 0.0397_f64, 0_f64, 78.084_f64, 20.946_f64,
            0.9340_f64, 0.1_f64, 0.00181_f64, 0.00017_f64, 0.00052_f64);

        let planet = Planet {orbit: orb, atmosphere: Some(atm), planet_type: PlanetType::Rocky,
            albedo: 0.306_f64, mass: 5.9726e+24_f64, radius: 6.371e+6_f64, min_temp: 183.95_f64,
            max_temp: 329.85_f64, avg_temp: 289.15_f64};

        assert!(5_514_f64*0.999 < planet.get_density() && 5_514_f64*1.001 > planet.get_density());
    }
}
