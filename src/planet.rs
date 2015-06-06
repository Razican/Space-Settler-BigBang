//! Planet implementation

use star::Star;
use std::f64::consts::PI;
use super::consts::*;

/// Orbit structure
struct Orbit<'o> {
    star: &'o Star,     // Orbiting star
    ecc: f64,           // Eccentricity
    sMa: f64,           // Semimajor axis (in meters)
    incl: f64,          // Inclination (in radians)
    lan: f64,           // Longitude of ascending node (in radians)
    argP: f64,          // Argument of periapsis (in radians)
    m0: f64,            // Mean anomaly at the creation of the universe (in radians)
    ax_tilt: f64,       // Axial tilt
    rot_period: f64,    // Rotation period in seconds
}

/// Atmosphere structure
struct Atmosphere {
    pressure: f64,  // Pressure in Pascals
    co2: f64,       // Carbon dioxide (percentage)
    co: f64,        // Carbon monoxide (percentage)
    n2: f64,        // Nitrogen (percentage)
    o2: f64,        // Oxygen (percentage)
    ar: f64,        // Argon (percentage)
    so2: f64,       // Sulfur dioxide (percentage)
    ne: f64,        // Neon (percentage)
    ch4: f64,       // Methane (percentage)
    he: f64,        // Helium (percentage)
}

// struct Soil {
// }

// struct Life {
// }

/// Planet structure
pub struct Planet<'p> {
    orbit: Orbit<'p>,
    atmosphere: Atmosphere,
    // soil: Soil,
    // life: Life,
    albedo: f32,
    mass: f64,
    radius: f64,
    min_temp: f64,
    max_temp: f64,
    avg_temp: f64
    // habitable: bool,
    // double_planet: bool,
}

impl<'p>  Planet<'p> {
    /// Constructs a new `Planet`.
    ///
    /// It creates a random planet taking into account real planet statistics. It requires the
    /// reference to parent star and the two values of the titius bode law, along with the order
    /// of the planet in the solar system.
    ///
    /// # Examples
    ///
    /// ```
    /// use star::Star;
    /// use planet::Planet;
    ///
    /// let st = Star::new(0, 1);
    ///
    /// let pl = Planet::new(&st, 0.0183, 1.0643, 3);
    /// ```
    pub fn new(st: &'p Star, n: f32, m: f32, pos: u32) -> Planet { // maybe position could  u8, must test performance

        let orb = Orbit::new(st, 0.01671022, 149.60e+9_f64, 0_f64, -0.196535244, 1.79676742,
            1.75343369, 0.409105177_f64, 86_164.2_f64);
        let atm = Atmosphere::new(101325_f64, 0.0397, 0_f64, 78.084, 20.946, 0.9340, 0_f64, 0.00181,
            0.000179, 0.000524);
        let alb = 0.306;
        let mass = 0_f64;
        let rad = 0_f64;
        let min_temp = 0_f64;
        let max_temp = 0_f64;
        let avg_temp = 0_f64;

        Planet {orbit: orb, atmosphere: atm, albedo: alb, mass: mass, radius: rad,
            min_temp: max_temp, max_temp: max_temp, avg_temp: avg_temp}
    }

    /// Gets the orbit information of the planet.
    pub fn get_orbit(&self) -> &Orbit {
        &self.orbit
    }

    /// Gets the atmosphere information of the planet.
    pub fn get_atmosphere(&self) -> &Atmosphere {
        &self.atmosphere
    }

    /// Gets the albedo of the planet.
    pub fn get_albedo(&self) -> f32 {
        self.albedo
    }

    /// Gets the mass of the planet in kg.
    pub fn get_mass(&self) -> f64 {
        self.mass
    }

    /// Gets the radius of the planet in meters.
    pub fn get_radius(&self) -> f64 {
        self.radius
    }

    /// Gets the density of the planet in kg/m³.
    pub fn get_density(&self) -> f64 {
        self.mass/self.get_volume()
    }

    /// Gets the minimum temperature of the planet in Kelvin.
    pub fn get_min_temp(&self) -> f64 {
        self.min_temp
    }

    /// Gets the maximum temperature of the planet in Kelvin.
    pub fn get_max_temp(&self) -> f64 {
        self.max_temp
    }

    /// Gets the average temperature of the planet in Kelvin.
    pub fn get_avg_temp(&self) -> f64 {
        self.avg_temp
    }

    /// Gets the volume of the planet in m³
    pub fn get_volume(&self) -> f64 {
        4_f64/3_f64*PI*self.radius.powi(3)
    }
}

impl<'o> Orbit<'o> {
    /// Constructs a new `Orbit`.
    ///
    /// It creates a new orbit structure with all the needed parameters for complete representation.
    pub fn new(star: &'o Star, ecc: f64, sMa: f64, incl: f64, lan: f64, argP: f64, m0: f64,
        ax_tilt: f64, rot_period: f64) -> Orbit {
        Orbit {star: star, ecc: ecc, sMa: sMa, incl: incl, lan: lan, argP: argP, m0: m0,
            ax_tilt: ax_tilt, rot_period: rot_period}
    }

    pub fn get_star(&self) -> &Star {
        self.star
    }

    /// Gets the eccentricity of the orbit.
    pub fn get_ecc(&self) -> f64 {
        self.ecc
    }

    /// Gets the semimajor axis of the orbit in meters.
    pub fn get_sMa(&self) -> f64 {
        self.sMa
    }

    /// Gets the inclination of the orbit in radians.
    pub fn get_incl(&self) -> f64 {
        self.incl
    }

    /// Gets the longitude of the ascending node of the orbit in radians.
    pub fn get_lan(&self) -> f64 {
        self.lan
    }

    /// Gets the argument of the periapsis of the orbit in radians.
    pub fn get_argP(&self) -> f64 {
        self.argP
    }

    /// Gets the mean anomaly of the orbit at the beginning of the universe in radians.
    pub fn get_anomaly(&self) -> f64 {
        self.m0
    }

    /// Gets the period of the orbit in seconds.
    pub fn get_orb_period(&self) -> f64 {
        2_f64*PI*(self.sMa.powi(3)/(G*self.star.get_mass())).sqrt()
    }

    /// Gets the apoapsis ob the orbit in meters.
    pub fn get_apoapsis(&self) -> f64 {
        self.sMa*(1_f64+self.ecc)
    }

    /// Gets the periapsis ob the orbit in meters.
    pub fn get_periapsis(&self) -> f64 {
        self.sMa*(1_f64-self.ecc)
    }

    /// Gets the axial tilt of the rotation of the body in radians.
    pub fn get_ax_tilt(&self) -> f64 {
        self.ax_tilt
    }

    /// Gets the sidereal rotation period of the body in seconds.
    pub fn get_rot_period(&self) -> f64 {
        self.rot_period
    }

    /// Gets the day length of the body in seconds.
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
    pub fn new(pressure: f64, co2: f64, co: f64, n2: f64, o2: f64, ar: f64, so2: f64,
        ne: f64, ch4: f64, he: f64) -> Atmosphere {
        Atmosphere {pressure: pressure, co2: co2, co: co, n2: n2, o2: o2, ar: ar,
            so2: so2, ne: ne, ch4: ch4, he: he}
    }

    /// Gets the pressure of the atmosphere in Pascals.
    pub fn get_pressure(&self) -> f64 {
        self.pressure
    }

    /// Gets the percentage of CO₂ (carbon dioxide) in the atmosphere.
    pub fn get_co2(&self) -> f64 {
        self.co2
    }

    /// Gets the percentage of CO in the atmosphere.
    pub fn get_co(&self) -> f64 {
        self.co
    }

    /// Gets the percentage of N₂ (nitrogen) in the atmosphere.
    pub fn get_n2(&self) -> f64 {
        self.n2
    }

    /// Gets the percentage of O₂ (oxygen) in the atmosphere.
    pub fn get_o2(&self) -> f64 {
        self.o2
    }

    /// Gets the percentage of Ar (argon) in the atmosphere.
    pub fn get_ar(&self) -> f64 {
        self.ar
    }

    /// Gets the percentage of SO₂ (sulfur dioxide) in the atmosphere.
    pub fn get_so2(&self) -> f64 {
        self.so2
    }

    /// Gets the percentage of Ne (neon) in the atmosphere.
    pub fn get_ne(&self) -> f64 {
        self.ne
    }

    /// Gets the percentage of CH₄ (methane) in the atmosphere.
    pub fn get_ch4(&self) -> f64 {
        self.ch4
    }

    /// Gets the percentage of He (helium) in the atmosphere.
    pub fn get_he(&self) -> f64 {
        self.he
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::star::Star;

    #[test]
    fn it_orbit_getters() {
        let st = Star::new(2, 0);

        let orb = super::Orbit {star: &st, ecc: 0.5_f64, sMa: 150e+9_f64, incl: 1.5_f64, lan: 1.2_f64,
            argP: 1.3_f64, m0: 1.4_f64, ax_tilt: 1.1_f64, rot_period: 80_600_f64};

        assert_eq!(3, orb.get_star().get_id());
        assert_eq!(0.5_f64, orb.get_ecc());
        assert_eq!(150e+9_f64, orb.get_sMa());
        assert_eq!(1.5_f64, orb.get_incl());
        assert_eq!(1.2_f64, orb.get_lan());
        assert_eq!(1.3_f64, orb.get_argP());
        assert_eq!(1.4_f64, orb.get_anomaly());
        assert_eq!(1.1_f64, orb.get_ax_tilt());
        assert_eq!(80_600_f64, orb.get_rot_period());
    }

    // TODO test day calculation and orbit period

    #[test]
    fn it_atm_getters() {
        let atm = super::Atmosphere {pressure: 101325_f64, co2: 0.0397_f64, co: 0_f64,
            n2: 78.084_f64, o2: 20.946_f64, ar: 0.9340_f64, so2: 0.1_f64, ne: 0.00181_f64,
            ch4: 0.00017_f64, he: 0.00052_f64};

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

    // TODO Planet getters tests

    #[test]
    fn it_star_id() {
        let st = Star::new(2, 0);
        let pl = Planet::new(&st, 0.0183, 1.0643, 3);

        assert_eq!(3, pl.get_orbit().get_star().get_id());
    }

    // TODO Planet creation tests
}
