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
    // habitable: bool,
    // double_planet: bool,
}

impl<'p>  Planet<'p> {
    /// Constructs a new `Planet`.
    ///
    /// It creates a random planet taking into account real planet statistics. It requires the
    /// reference to parent star.
    ///
    /// # Examples
    ///
    /// ```
    /// use star::Star;
    /// use planet::Planet;
    ///
    /// let st = Star::new(0, 1);
    ///
    /// let pl = Planet::new(&st, 0.0183, 1.0643);
    /// ```
    pub fn new(st: &'p Star, n: f32, m: f32) -> Planet {

        let orb = Orbit::new(st, 0.01671022, 149.60e+9_f64, 0_f64, -0.196535244, 1.79676742,
            1.75343369, 0.409105177_f64, 86_164.2_f64);
        let atm = Atmosphere::new(101325_f64, 0.0397, 0_f64, 78.084, 20.946, 0.9340, 0_f64, 0.00181,
            0.000179, 0.000524);
        let alb = 0.306;

        Planet {orbit: orb, atmosphere: atm, albedo: alb}
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

    /// Gets the eccentricity of the orbit
    pub fn get_ecc(&self) -> f64 {
        self.ecc
    }

    /// Gets the semimajor axis of the orbit in meters
    pub fn get_sMa(&self) -> f64 {
        self.sMa
    }

    /// Gets the inclination of the orbit in radians
    pub fn get_incl(&self) -> f64 {
        self.incl
    }

    /// Gets the longitude of the ascending node of the orbit in radians
    pub fn get_lan(&self) -> f64 {
        self.lan
    }

    /// Gets the argument of the periapsis of the orbit in radians
    pub fn get_argP(&self) -> f64 {
        self.argP
    }

    /// Gets the mean anomaly of the orbit at the beginning of the universe in radians
    pub fn get_anomaly(&self) -> f64 {
        self.m0
    }

    /// Gets the period of the orbit in seconds
    pub fn get_orb_period(&self) -> f64 {
        2_f64*PI*(self.sMa.powi(3)/(G*self.star.get_mass())).sqrt()
    }

    /// Gets the apoapsis ob the orbit in meters
    pub fn get_apoapsis(&self) -> f64 {
        self.sMa*(1_f64+self.ecc)
    }

    /// Gets the periapsis ob the orbit in meters
    pub fn get_periapsis(&self) -> f64 {
        self.sMa*(1_f64-self.ecc)
    }

    /// Gets the axial tilt of the rotation of the body in radians
    pub fn get_ax_tilt(&self) -> f64 {
        self.ax_tilt
    }

    /// Gets the sidereal rotation period of the body in seconds
    pub fn get_rot_period(&self) -> f64 {
        self.rot_period
    }

    /// Gets the day length of the body in seconds
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

    /// Gets the pressure of the atmosphere in Pascals
    pub fn get_pressure(&self) -> f64 {
        self.pressure
    }

    /// Gets the percentage of CO₂ (carbon dioxide) in the atmosphere
    pub fn get_co2(&self) -> f64 {
        self.co2
    }

    /// Gets the percentage of CO in the atmosphere
    pub fn get_co(&self) -> f64 {
        self.co
    }

    /// Gets the percentage of N₂ (nitrogen) in the atmosphere
    pub fn get_n2(&self) -> f64 {
        self.n2
    }

    /// Gets the percentage of O₂ (oxygen) in the atmosphere
    pub fn get_o2(&self) -> f64 {
        self.o2
    }

    /// Gets the percentage of Ar (argon) in the atmosphere
    pub fn get_ar(&self) -> f64 {
        self.ar
    }

    /// Gets the percentage of SO₂ (sulfur dioxide) in the atmosphere
    pub fn get_so2(&self) -> f64 {
        self.so2
    }

    /// Gets the percentage of Ne (neon) in the atmosphere
    pub fn get_ne(&self) -> f64 {
        self.ne
    }

    /// Gets the percentage of CH₄ (methane) in the atmosphere
    pub fn get_ch4(&self) -> f64 {
        self.ch4
    }

    /// Gets the percentage of He (helium) in the atmosphere
    pub fn get_he(&self) -> f64 {
        self.he
    }
}
