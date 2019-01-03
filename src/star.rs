//! Star module
//!
//! This is the `star` module. This module contains all the structures, enumerations and
//! implementations needed to define a star.

extern crate rand;

use std::{f64::consts::PI, fmt};

use self::rand::Rng;

use crate::{consts::*, planet::Planet};

/// Star class enumeration
///
/// This enumeration contains all the star types available in this universe creation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Class {
    /// This represents a black hole, an structure so dense that not even the light can escape. Due
    /// to the high radiation, no planets can survive near one of these structures. Their mass is
    /// between 2 and 20 solar masses. Approximately the 0.03333% of the stars in the galaxy will be
    /// of this type. Its effective temperature will be represented as 0K.
    BlackHole,
    /// This represents a neutron star. Neutron stars are really compact stars, that due to their
    /// extreme gravitational forces not even atomic nuclei can survive, and the protons are
    /// degraded to neutrons. Their mass is between 1.1 and 2.2 solar masses. Due to the high
    /// radiation, no planet can survive near a neutron star. Some neutron stars can be seen as
    /// pulsars. Approximately the 0.03333% of the stars in the galaxy will be of this type.
    NeutronStar,
    /// This represents a quark star. Quark stars are hypothetical stars that are even more compact
    /// than neutron stars but they have not collapsed into black holes. They have between 2 and 3
    /// solar masses. They have such big pressures that even the quarks that are part of the
    /// neutrons get free. Approximately the 0.00001% of the stars in the galaxy will be of this
    /// type, so most of the galaxies will not even have one.
    QuarkStar,
    /// This represents a white dwarf. A white dwarf is a compact star that is usually the remnant
    /// of a Sun-type star. They are very dense and hot, and their size is inversely proportional to
    /// the cube of their mass, which can be between 0.2 and 1.39 solar masses. This last mass is
    /// the Chandrasekhar limit. There are not many planets around quark stars, and never more than
    /// 2. Approximately the 0.03332% of the stars in the galaxy will be of this type.
    WhiteDwarf,
    /// This represents an O type star. This stars are hot and blue colored. Their surface
    /// temperature can be as high as 52,000K. Approximately the 0.00003% of the stars will be of
    /// this type, so many galaxies will not even have one. No planet can survive to the high
    /// radiation coming from the star.
    O,
    /// This represents a B type star. This stars are hot and blue-white colored. Their surface
    /// temperature is usually between 10,000K and 30,000K and their mass between 2.1 and 16 solar
    /// masses. Approximately the 0.13% of the stars will be of this type. No more than 2 planets
    /// will survive to the radiations, and only in small B type stars.
    B,
    /// This represents an A type of star. This stars are white stars. Their surface temperature is
    /// usually between 7,500K and 10,000K and their mass between 1.4 and 2.1 solar masses.
    /// Approximately the 0.6% of the stars will be of this type. They will have up to 3 planets
    /// orbiting them.
    A,
    /// This represents an F type of star. This stars are yellow-white stars. Their surface
    /// temperature is usually between 6,000K and 7,500K and their mass between 1.04 and 1.4 solar
    /// masses. Approximately the 3% of the stars will be of this type. They will have up to 8
    /// planets orbiting them.
    F,
    /// This represents an G type of star, the same type as our Sun. This stars are yellow stars.
    /// Their surface temperature is usually between 5,200K and 6,000K and their mass between 0.8
    /// and 1.04 solar masses. Approximately the 7.6% of the stars will be of this type. They will
    /// have up to 10 planets orbiting them.
    G,
    /// This represents an K type of star. This stars are orange stars. Their surface temperature is
    /// usually between 3,700K and 5,200K and their mass between 0.45 and 0.8 solar masses.
    /// Approximately the 12.1% of the stars will be of this type. They will have up to 10 planets
    /// orbiting them.
    K,
    /// This represents an M type of star. This stars are red stars. Their surface temperature is
    /// usually between 2,400K and 3,700K and their mass between 0.08 and 0.45 solar masses.
    /// Approximately the 76.45% of the stars will be of this type. They will have up to 9 planets
    /// orbiting them.
    M,
}

impl fmt::Display for Class {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Class::BlackHole => write!(f, "black hole"),
            Class::NeutronStar => write!(f, "neutron star"),
            Class::QuarkStar => write!(f, "quark star"),
            Class::WhiteDwarf => write!(f, "white dwarf"),
            Class::O => write!(f, "O"),
            Class::B => write!(f, "B"),
            Class::A => write!(f, "A"),
            Class::F => write!(f, "F"),
            Class::G => write!(f, "G"),
            Class::K => write!(f, "K"),
            Class::M => write!(f, "M"),
        }
    }
}

/// Star structure
///
/// This structure represents all the parameters of the given star.
#[derive(Debug, Clone)]
pub struct Star {
    id: u64,
    orb_radius: i32,
    class: Class,
    mass: f64,
    radius: f64,
    temperature: i32,
    planets: Vec<Planet>,
}

impl Star {
    /// Constructs a new `Star`.
    ///
    /// It creates the star generating its parameters using statistical information from the known
    /// universe. It has two parameters: the first is the ID of the last star created and the
    /// second one is the ID of the current galaxy.
    ///
    /// # Examples
    ///
    /// ```
    /// use star::Star;
    ///
    /// let st = Star::new(0, 1);
    /// ```
    pub fn new(last_id: u64) -> Star {
        let orbit = rand::thread_rng().gen_range(20_000, 30_001);
        let star_class = Star::generate_class();
        let (star_mass, star_radius, star_temperature) = Star::generate_properties(star_class);

        Star {
            id: last_id + 1,
            orb_radius: orbit,
            class: star_class,
            mass: star_mass,
            radius: star_radius,
            temperature: star_temperature,
            planets: Vec::new(),
        }
    }

    /// Get ID
    ///
    /// Gets the ID of the star.
    pub fn get_id(&self) -> u64 {
        self.id
    }

    /// Get orbit
    ///
    /// Gets the orbit radius around the center of the galaxy, in 1/100 light years (*ly*).
    pub fn get_orbit(&self) -> i32 {
        self.orb_radius
    }

    /// Get class
    ///
    /// Gets the class of the star.
    pub fn get_class(&self) -> &Class {
        &self.class
    }

    /// Get mass
    ///
    /// Gets the mass of the star, in *kg*.
    pub fn get_mass(&self) -> f64 {
        self.mass
    }

    /// Get radius
    ///
    /// Gets the radius of the star, in meters (*m*).
    pub fn get_radius(&self) -> f64 {
        self.radius
    }

    /// Get density
    ///
    /// Gets the density of the star, in *kg/m³*.
    pub fn get_density(&self) -> f64 {
        self.mass / self.get_volume()
    }

    /// Get volume
    ///
    /// Gets the volume of the star, in *m³*.
    pub fn get_volume(&self) -> f64 {
        4_f64 / 3_f64 * PI * self.radius.powi(3)
    }

    /// Get temperature
    ///
    /// Gets the effective surface temperature of the star, in 1/1,000 Kelvin (*K*).
    pub fn get_temperature(&self) -> i32 {
        self.temperature
    }

    /// Get luminosity
    ///
    /// Gets the luminosity of the star, in watts (*W*).
    pub fn get_luminosity(&self) -> f64 {
        4_f64 * PI * self.radius.powi(2) * BOLTZ * f64::from(self.temperature).powi(4)
    }

    // ----------  generators  ----------

    /// Generate class
    ///
    /// Generates the class of the star taking into account real star proportions in the known
    /// universe.
    fn generate_class() -> Class {
        let prob = rand::thread_rng().gen_range(1, 10_000_001);

        // Generate star class based on probabilities.
        match prob {
            0...3_333 => Class::BlackHole,
            3_334...6_666 => Class::NeutronStar,
            6_667 => Class::QuarkStar,
            6_668...10_000 => Class::WhiteDwarf,
            10_001...10_003 => Class::O,
            10_004...23_000 => Class::B,
            23_001...83_000 => Class::A,
            83_001...383_000 => Class::F,
            383_001...1_143_000 => Class::G,
            1_143_001...2_353_000 => Class::K,
            2_353_001...10_000_000 => Class::M,
            _ => unreachable!(),
        }
    }

    /// Generate properties
    ///
    /// Generates the basic properties (mass, radius and effective surface temperature) of the star
    /// taking into account the star class.
    fn generate_properties(class: Class) -> (f64, f64, i32) {
        use crate::consts::{SUN_MASS, SUN_RADIUS};

        match class {
            Class::BlackHole => {
                let mass = if rand::thread_rng().gen_range(0, 2) == 0 {
                    // 5 - 10 solar masses
                    rand::thread_rng().gen_range(5_f64 * SUN_MASS, 10_f64 * SUN_MASS)
                } else {
                    // 2 - 20 solar masses
                    rand::thread_rng().gen_range(2_f64 * SUN_MASS, 20_f64 * SUN_MASS)
                };

                let radius = 2_f64 * G * mass / C.powi(2); // m
                let temperature = 0; // N/A

                (mass, radius, temperature)
            }
            Class::NeutronStar => {
                // 1.1 - 2.2 solar masses
                let mass = rand::thread_rng().gen_range(1.1 * SUN_MASS, 2.2 * SUN_MASS);

                let radius = rand::thread_rng().gen_range(11_000_f64, 15_000_f64); // m
                let temperature = rand::thread_rng().gen_range(100_000, 10_000_001); // Kelvin

                (mass, radius, temperature)
            }
            Class::QuarkStar => {
                // 2 - 3 solar masses
                let mass = rand::thread_rng().gen_range(2_f64 * SUN_MASS, 3_f64 * SUN_MASS);

                let radius = rand::thread_rng().gen_range(11_000_f64, 15_000_f64); // m
                let temperature = rand::thread_rng().gen_range(8_000_000, 100_000_001); // Kelvin

                (mass, radius, temperature)
            }
            Class::WhiteDwarf => {
                // 0.2 - 1.3 solar masses
                let mass = rand::thread_rng().gen_range(0.2 * SUN_MASS, 1.3 * SUN_MASS);

                let radius = 0.78_e+7_f64
                    * ((CH_LIMIT / mass).powf(2_f64 / 3_f64)
                        - (mass / CH_LIMIT).powf(2_f64 / 3_f64))
                    .sqrt();

                let temperature = if rand::thread_rng().gen_range(0, 4) == 0 {
                    rand::thread_rng().gen_range(4_000, 150_001) // Kelvin
                } else {
                    rand::thread_rng().gen_range(6_000, 30_000) // Kelvin
                };

                (mass, radius, temperature)
            }
            Class::O => {
                // 15 - 90 solar masses
                let mass = rand::thread_rng().gen_range(15_f64 * SUN_MASS, 90_f64 * SUN_MASS);

                // 6.6 - 40 solar radius
                let radius = rand::thread_rng().gen_range(6.6 * SUN_RADIUS, 40_f64 * SUN_RADIUS);
                let temperature = rand::thread_rng().gen_range(30_000, 52_001); // Kelvin

                (mass, radius, temperature)
            }
            Class::B => {
                // 2.1 - 16 solar masses
                let mass = rand::thread_rng().gen_range(2.1 * SUN_MASS, 16_f64 * SUN_MASS);

                // 1.8 - 6.6 solar radius
                let radius = rand::thread_rng().gen_range(1.8 * SUN_RADIUS, 6.6 * SUN_RADIUS);
                let temperature = rand::thread_rng().gen_range(10_000, 30_001); // Kelvin

                (mass, radius, temperature)
            }
            Class::A => {
                // 1.4 - 2.1 solar masses
                let mass = rand::thread_rng().gen_range(1.4 * SUN_MASS, 2.1 * SUN_MASS);

                // 1.4 - 1.8 solar radius
                let radius = rand::thread_rng().gen_range(1.4 * SUN_RADIUS, 1.8 * SUN_RADIUS);
                let temperature = rand::thread_rng().gen_range(7_500, 10_001); // Kelvin

                (mass, radius, temperature)
            }
            Class::F => {
                // 1.04 - 1.4 solar masses
                let mass = rand::thread_rng().gen_range(1.04 * SUN_MASS, 1.4 * SUN_MASS);

                // 1.15 - 1.4 solar radius
                let radius = rand::thread_rng().gen_range(1.15 * SUN_RADIUS, 1.4 * SUN_RADIUS);
                let temperature = rand::thread_rng().gen_range(6_000, 7_501); // Kelvin

                (mass, radius, temperature)
            }
            Class::G => {
                // 0.8 - 1.04 solar masses
                let mass = rand::thread_rng().gen_range(0.8 * SUN_MASS, 1.04 * SUN_MASS);

                // 0.96 - 1.15 solar radius
                let radius = rand::thread_rng().gen_range(0.96 * SUN_RADIUS, 1.15 * SUN_RADIUS);
                let temperature = rand::thread_rng().gen_range(5_200, 6_001); // Kelvin

                (mass, radius, temperature)
            }
            Class::K => {
                // 0.45 - 0.8 solar masses
                let mass = rand::thread_rng().gen_range(0.45 * SUN_MASS, 0.8 * SUN_MASS);

                // 0.7 - 0.96 solar radius
                let radius = rand::thread_rng().gen_range(0.7 * SUN_RADIUS, 0.96 * SUN_RADIUS);
                let temperature = rand::thread_rng().gen_range(3_700, 5_201); // Kelvin

                (mass, radius, temperature)
            }
            Class::M => {
                let mass = if rand::thread_rng().gen_range(0, 4) == 0 {
                    // 0.1 - 0.45 solar masses
                    rand::thread_rng().gen_range(0.1 * SUN_MASS, 0.45 * SUN_MASS)
                } else {
                    // 0.1 - 0.25 solar masses
                    rand::thread_rng().gen_range(0.1 * SUN_MASS, 0.25 * SUN_MASS)
                };

                let radius = if rand::thread_rng().gen_range(0, 4) == 0 {
                    // solar_masses-0.05 - 0.7 solar radius
                    rand::thread_rng()
                        .gen_range((mass / SUN_MASS - 0.05) * SUN_RADIUS, 0.7 * SUN_RADIUS)
                } else {
                    // solar_masses-0.05 - 0.5 solar radius
                    rand::thread_rng()
                        .gen_range((mass / SUN_MASS - 0.05) * SUN_RADIUS, 0.5 * SUN_RADIUS)
                };

                let temperature = rand::thread_rng().gen_range(2_500, 3_701); // Kelvin

                (mass, radius, temperature)
            }
        }
    }

    /// Generate number of bodies
    ///
    /// Generates a random number based on the star class and luminosity, that will be the number of
    /// bodies in the solar system.
    pub fn generate_num_bodies(&self) -> u8 {
        match self.class {
            Class::O | Class::BlackHole | Class::NeutronStar | Class::QuarkStar => 0,
            Class::WhiteDwarf | Class::B => {
                if self.get_luminosity() < 300_f64 * SUN_LUMINOSITY {
                    rand::thread_rng().gen_range(0, 3)
                } else {
                    0
                }
            }
            Class::A => rand::thread_rng().gen_range(0, 4),
            Class::F => rand::thread_rng().gen_range(0, 9),
            Class::G => {
                if rand::thread_rng().gen_range(0, 11) == 0 {
                    rand::thread_rng().gen_range(0, 4)
                } else {
                    rand::thread_rng().gen_range(4, 11)
                }
            }
            Class::K => {
                if rand::thread_rng().gen_range(0, 16) == 0 {
                    rand::thread_rng().gen_range(0, 5)
                } else {
                    rand::thread_rng().gen_range(5, 13)
                }
            }
            Class::M => {
                if rand::thread_rng().gen_range(0, 21) == 0 {
                    rand::thread_rng().gen_range(0, 7)
                } else {
                    rand::thread_rng().gen_range(7, 16)
                }
            }
        }
    }

    /// Generate Titius-Bode
    ///
    /// Generates Titius-Bode law *m* and *n* parameters. It supposes that bodies > 0, or it will
    /// panic.
    pub fn generate_titius_bode(&self, bodies: u8) -> (f64, f64) {
        if bodies > 0 {
            let luminosity = self.get_luminosity();
            let m = if luminosity < 0.01 * SUN_LUMINOSITY {
                rand::thread_rng().gen_range(0.004_f64, 0.01_f64)
            } else if luminosity < 0.1 * SUN_LUMINOSITY {
                rand::thread_rng().gen_range(0.0075_f64, 0.025_f64)
            } else if luminosity < 0.5 * SUN_LUMINOSITY {
                rand::thread_rng().gen_range(0.015_f64, 0.1_f64)
            } else if luminosity < 6_f64 * SUN_LUMINOSITY {
                if rand::thread_rng().gen_range(0, (10_f64 - luminosity / SUN_LUMINOSITY) as u8) > 0
                {
                    rand::thread_rng().gen_range(0.01_f64, 0.1_f64)
                } else {
                    rand::thread_rng().gen_range(0.025_f64, 0.6_f64)
                }
            } else if luminosity < 10_f64 * SUN_LUMINOSITY {
                rand::thread_rng().gen_range(
                    0.1 * luminosity / SUN_LUMINOSITY,
                    0.3 * luminosity / SUN_LUMINOSITY,
                )
            } else if luminosity < 30_f64 * SUN_LUMINOSITY {
                if luminosity * 0.15 > 3_f64 {
                    rand::thread_rng().gen_range(0.1 * luminosity / SUN_LUMINOSITY, 3_f64)
                } else {
                    rand::thread_rng().gen_range(
                        0.1 * luminosity / SUN_LUMINOSITY,
                        0.15 * luminosity / SUN_LUMINOSITY,
                    )
                }
            } else {
                rand::thread_rng().gen_range(2.5_f64, 5_f64)
            };

            let n = if m > 0.035 {
                m.sqrt() * 1.7 + 0.165
            } else {
                0.017 / m
            };

            (m, rand::thread_rng().gen_range(0.8 * n, 1.2 * n))
        } else {
            panic!("Tried to generate Titius-Bode parameters for a 0 body system!")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        super::consts::{C, CH_LIMIT, G, SUN_MASS, SUN_RADIUS},
        Class, Star,
    };
    use std::f64::EPSILON;

    #[test]
    fn it_star_getters() {
        let st = Star {
            id: 2,
            orb_radius: 26_000,
            class: Class::G,
            mass: 1.9885e+30,
            radius: 6.96e+8,
            temperature: 5_778,
            planets: Vec::new(),
        };

        assert_eq!(2, st.get_id());
        assert_eq!(26_000, st.get_orbit());
        assert_eq!(&Class::G, st.get_class());
        assert!(st.get_mass() >= 1.9885e+30 - EPSILON && st.get_mass() <= 1.9885e+30 + EPSILON);
        assert!(st.get_radius() >= 6.96e+8 - EPSILON && st.get_radius() <= 6.96e+8 + EPSILON);
        assert_eq!(5_778, st.get_temperature());
    }

    #[test]
    fn it_star_parameters() {
        let st = Star::new(1);

        assert_eq!(2, st.get_id());
    }

    #[test]
    fn it_luminosity() {
        let st = Star {
            id: 2,
            orb_radius: 26_000,
            class: Class::G,
            mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64,
            temperature: 5_778,
            planets: Vec::new(),
        };

        assert!(
            st.get_luminosity() > 384.6_e+24_f64 * 0.999
                && st.get_luminosity() < 384.6_e+24_f64 * 1.001
        );
    }

    #[test]
    fn it_volume() {
        let st = Star {
            id: 2,
            orb_radius: 26_000,
            class: Class::G,
            mass: 1.9885_e+30_f64,
            radius: 6.96_e+8_f64,
            temperature: 5_778,
            planets: Vec::new(),
        };

        assert!(
            st.get_volume() > 1.412_e+27_f64 * 0.999 && st.get_volume() < 1.412_e+27_f64 * 1.001
        );
    }

    #[test]
    fn it_density() {
        let st = Star {
            id: 2,
            orb_radius: 26_000,
            class: Class::G,
            mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64,
            temperature: 5_778,
            planets: Vec::new(),
        };

        assert!(st.get_density() > 1_408_f64 * 0.999 && st.get_density() < 1_408_f64 * 1.001);
    }

    #[test]
    fn it_bh_properties() {
        // Black hole test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::BlackHole);
            assert!(mass >= 2_f64 * SUN_MASS && mass <= 20_f64 * SUN_MASS);
            let cmp_radius = 2_f64 * G * mass / C.powi(2);
            assert!(cmp_radius >= radius - EPSILON && cmp_radius <= radius + EPSILON);
            assert_eq!(0, temperature);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert_eq!(0, st.generate_num_bodies());
        }
    }

    #[test]
    fn it_ns_properties() {
        // Neutron star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::NeutronStar);
            assert!(mass >= 1.1_f64 * SUN_MASS && mass <= 2.2_f64 * SUN_MASS);
            assert!(radius >= 11_000_f64 && radius <= 15_000_f64);
            assert!(radius > 2_f64 * G * mass / C.powi(2));
            assert!(temperature >= 10_000 && temperature <= 10_000_000);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert_eq!(0, st.generate_num_bodies());
        }
    }

    #[test]
    fn it_qs_properties() {
        // Quark star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::QuarkStar);
            assert!(mass >= 2_f64 * SUN_MASS && mass <= 3_f64 * SUN_MASS);
            assert!(radius >= 11_000_f64 && radius <= 15_000_f64);
            assert!(radius > 2_f64 * G * mass / C.powi(2));
            assert!(temperature >= 8_000_000 && temperature <= 100_000_000);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert_eq!(0, st.generate_num_bodies());
        }
    }

    #[test]
    fn it_wd_properties() {
        // White dwarf test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::WhiteDwarf);
            assert!(mass >= 0.2_f64 * SUN_MASS && mass <= 1.3_f64 * SUN_MASS);
            let cmp_radius = 0.78e+7_f64
                * ((CH_LIMIT / mass).powf(2_f64 / 3_f64) - (mass / CH_LIMIT).powf(2_f64 / 3_f64))
                    .sqrt();
            assert!(cmp_radius <= radius + EPSILON && cmp_radius >= radius - EPSILON);
            assert!(temperature >= 4_000 && temperature <= 150_000);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert!(st.generate_num_bodies() <= 2);
        }
    }

    #[test]
    fn it_o_properties() {
        // O star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::O);
            assert!(mass >= 15_f64 * SUN_MASS && mass <= 90_f64 * SUN_MASS);
            assert!(radius >= 6.6_f64 * SUN_RADIUS && radius <= 40_f64 * SUN_RADIUS);
            assert!(temperature >= 30_000 && temperature <= 52_000);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert_eq!(0, st.generate_num_bodies());
        }
    }

    #[test]
    fn it_b_properties() {
        // B star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::B);
            assert!(mass >= 2.1_f64 * SUN_MASS && mass <= 16_f64 * SUN_MASS);
            assert!(radius >= 1.8_f64 * SUN_RADIUS && radius <= 6.6_f64 * SUN_RADIUS);
            assert!(temperature >= 10_000 && temperature <= 30_000);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert!(st.generate_num_bodies() <= 2);
        }
    }

    #[test]
    fn it_a_properties() {
        // A star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::A);
            assert!(mass >= 1.4 * SUN_MASS && mass <= 2.1 * SUN_MASS);
            assert!(radius >= 1.4 * SUN_RADIUS && radius <= 1.8 * SUN_RADIUS);
            assert!(temperature >= 7_500 && temperature <= 10_000);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert!(st.generate_num_bodies() <= 3);
        }
    }

    #[test]
    fn it_f_properties() {
        // F star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::F);
            assert!(mass >= 1.04_f64 * SUN_MASS && mass <= 1.4_f64 * SUN_MASS);
            assert!(radius >= 1.15_f64 * SUN_RADIUS && radius <= 1.4_f64 * SUN_RADIUS);
            assert!(temperature >= 6_000 && temperature <= 7_500);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert!(st.generate_num_bodies() <= 8);
        }
    }

    #[test]
    fn it_g_properties() {
        // G star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::G);
            assert!(mass >= 0.8_f64 * SUN_MASS && mass <= 1.04_f64 * SUN_MASS);
            assert!(radius >= 0.96_f64 * SUN_RADIUS && radius <= 1.15_f64 * SUN_RADIUS);
            assert!(temperature >= 5_200 && temperature <= 6_000);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert!(st.generate_num_bodies() <= 10);
        }
    }

    #[test]
    fn it_k_properties() {
        // K star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::K);
            assert!(mass >= 0.45_f64 * SUN_MASS && mass <= 0.8_f64 * SUN_MASS);
            assert!(radius >= 0.7_f64 * SUN_RADIUS && radius <= 0.96_f64 * SUN_RADIUS);
            assert!(temperature >= 3_700 && temperature <= 5_200);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert!(st.generate_num_bodies() <= 10);
        }
    }

    #[test]
    fn it_m_properties() {
        // M star test
        for _i in 1..10_000 {
            let (mass, radius, temperature) = Star::generate_properties(Class::M);
            assert!(mass >= 0.1_f64 * SUN_MASS && mass <= 0.45_f64 * SUN_MASS);
            assert!(radius <= 0.7_f64 * SUN_RADIUS);
            assert!(temperature >= 2_400 && temperature <= 3_700);

            let st = Star {
                id: 2,
                orb_radius: 26_000,
                class: Class::NeutronStar,
                mass,
                radius,
                temperature,
                planets: Vec::new(),
            };

            assert!(st.generate_num_bodies() <= 9);
        }
    }

    #[test]
    #[should_panic]
    fn it_generate_titius_bode_fail() {
        let st = Star {
            id: 2,
            orb_radius: 26_000,
            class: Class::G,
            mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64,
            temperature: 5_778,
            planets: Vec::new(),
        };

        st.generate_titius_bode(0);
    }
}
