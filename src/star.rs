//! Star implementation

extern crate rand;

use std::f64::consts::PI;

use self::rand::Rng;

use super::consts::*;

#[derive(Debug, PartialEq)]
enum StarClass {
    BlackHole,
    NeutronStar,
    QuarkStar,
    WhiteDwarf,
    O,
    B,
    A,
    F,
    G,
    K,
    M,
}

/// Star structure
pub struct Star {
    id: u64,
    galaxy_id: u32,
    orb_radius: u32,
    class: StarClass,
    mass: f64,
    radius: f64,
    temp: u32,
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
    pub fn new(last_id:u64, galaxy: u32) -> Star {
        let orbit = rand::thread_rng().gen_range(20000, 30001);
        let star_class = Star::generate_class();
        let (star_mass, star_radius, star_temp) = Star::generate_properties(&star_class);

        Star {id: last_id+1, galaxy_id: galaxy, orb_radius: orbit, class: star_class,
            mass: star_mass, radius: star_radius, temp: star_temp}
    }

    /// Gets the ID of the star.
    pub fn get_id(&self) -> u64 {
        self.id
    }

    /// Gets the ID of the galaxy of the star.
    pub fn get_galaxy_id(&self) -> u32 {
        self.galaxy_id
    }

    /// Gets the distance to the center of the galaxy in light years.
    pub fn get_orbit(&self) -> u32 {
        self.orb_radius
    }

    /// Gets the star class.
    pub fn get_class(&self) -> &StarClass {
        &self.class
    }

    /// Gets the mass of the star in kg.
    pub fn get_mass(&self) -> f64 {
        self.mass
    }

    /// Gets the radius of the star in meters.
    pub fn get_radius(&self) -> f64 {
        self.radius
    }

    /// Gets the density of the star in kg/m³.
    pub fn get_density(&self) -> f64 {
        self.mass/self.get_volume()
    }

    /// Gets the volume of the star in m³
    pub fn get_volume(&self) -> f64 {
        4_f64/3_f64*PI*self.radius.powi(3)
    }

    /// Gets the effective surface temperature of the star in Kelvin.
    pub fn get_temp(&self) -> u32 {
        self.temp
    }

    /// Gets the luminosity of the star in Watts.
    pub fn get_luminosity(&self) -> f64 {
        4_f64*PI*self.radius.powi(2)*BOLTZ*(self.temp as f64).powi(4) // W
    }

    // ----------  generators  ----------

    fn generate_class() -> StarClass {
        let prob = rand::thread_rng().gen_range(1, 10000001);

        // Generate star class based on probabilities.
        match prob {
            0...3_333               => StarClass::BlackHole,
            3_334...6_666           => StarClass::NeutronStar,
            6_667                   => StarClass::QuarkStar,
            6_668...10_000          => StarClass::WhiteDwarf,
            10_001...10_003         => StarClass::O,
            10_004...23_000         => StarClass::B,
            23_001...83_000         => StarClass::A,
            83_001...383_000        => StarClass::F,
            383_001...1_143_000     => StarClass::G,
            1_143_001...2_353_000   => StarClass::K,
            2_353_001...10_000_000  => StarClass::M,
            _ => unreachable!()
        }
    }

    fn generate_properties(class: &StarClass) -> (f64, f64, u32) {
        match *class {
            StarClass::BlackHole => {
                    let mass = if rand::thread_rng().gen_range(0, 2) == 0 { // TODO gaussian distribution
                        rand::thread_rng().gen_range(9.942e+30_f64, 1.9884e+31_f64) // 5 - 10 solar masses
                    } else {
                        rand::thread_rng().gen_range(3.9768e+30_f64, 3.9768e+31_f64) // 2 - 20 solar masses
                    };

                    let radius = 2_f64*G*mass/C.powi(2); // m
                    let temp = 0; // N/A

                    (mass, radius, temp)
                },
            StarClass::NeutronStar => {
                    let mass = rand::thread_rng().gen_range(2.18724e+30_f64, 4.37448e+30_f64); // 1.1 - 2.2 solar masses

                    let radius = rand::thread_rng().gen_range(11_000_f64, 15_000_f64); // m
                    let temp = rand::thread_rng().gen_range(100_000, 10_000_001); // Kelvin

                    (mass, radius, temp)
                },
            StarClass::QuarkStar => {
                    let mass = rand::thread_rng().gen_range(3.9768e+30_f64, 5.9652e+30_f64); // 2 - 3 solar masses

                    let radius = rand::thread_rng().gen_range(11_000_f64, 15_000_f64); // m
                    let temp = rand::thread_rng().gen_range(8_000_000, 100_000_001); // Kelvin (N/A)

                    (mass, radius, temp)
                },
            StarClass::WhiteDwarf => {
                    let mass = rand::thread_rng().gen_range(3.9768e+29_f64, 2.58492e+30_f64); // 0.2 - 1.3 solar masses

                    let radius = 0.78e+7_f64*((CH_LIMIT/mass).powf(2_f64/3_f64) - (mass/CH_LIMIT).powf(2_f64/3_f64)).sqrt();

                    let temp = if rand::thread_rng().gen_range(0,4) == 0 {
                        rand::thread_rng().gen_range(4_000, 150_001) // Kelvin
                    } else {
                        rand::thread_rng().gen_range(6_000, 30_000) // Kelvin
                    };

                    (mass, radius, temp)
                },
            StarClass::O => {
                    let mass = rand::thread_rng().gen_range(2.9826e+31_f64, 1.78956e+32_f64); // 15 - 90 solar masses

                    let radius = rand::thread_rng().gen_range(4.5936e+9_f64, 2.784e+10_f64); // 6.6 - 40 solar radius
                    let temp = rand::thread_rng().gen_range(30_000, 52_001); // Kelvin

                    (mass, radius, temp)
                },
            StarClass::B => {
                    let mass = rand::thread_rng().gen_range(4.17564e+30_f64, 3.18144e+31_f64); // 2.1 - 16 solar masses

                    let radius = rand::thread_rng().gen_range(1.2528e+9_f64, 4.5936e+9_f64); // 1.8 - 6.6 solar radius
                    let temp = rand::thread_rng().gen_range(10_000, 30_001); // Kelvin

                    (mass, radius, temp)
                },
            StarClass::A => {
                    let mass = rand::thread_rng().gen_range(2.78376e+30_f64, 4.17564e+30_f64); // 1.4 - 2.1 solar masses

                    let radius = rand::thread_rng().gen_range(9.744e+8_f64, 1.2528e+9_f64); // 1.4 - 1.8 solar radius
                    let temp = rand::thread_rng().gen_range(7_500, 10_001); // Kelvin

                    (mass, radius, temp)
                },
            StarClass::F => {
                    let mass = rand::thread_rng().gen_range(2.067936e+30_f64, 2.78376e+30_f64); // 1.04 - 1.4 solar masses

                    let radius = rand::thread_rng().gen_range(8.004e+8_f64, 9.744e+8_f64); // 1.15 - 1.4 solar radius
                    let temp =rand::thread_rng().gen_range(6_000, 7_501); // Kelvin

                    (mass, radius, temp)
                },
            StarClass::G => {
                    let mass = rand::thread_rng().gen_range(1.59072e+30_f64, 2.067936e+30_f64); // 0.8 - 1.04 solar masses

                    let radius = rand::thread_rng().gen_range(6.6816e+8_f64, 8.004e+8_f64); // 0.96 - 1.15 solar radius
                    let temp = rand::thread_rng().gen_range(5_200, 6_001); // Kelvin

                    (mass, radius, temp)
                },
            StarClass::K => {
                    let mass = rand::thread_rng().gen_range(8.9478e+29_f64, 1.59072e+30_f64); // 0.45 - 0.8 solar masses

                    let radius = rand::thread_rng().gen_range(4.872e+8_f64, 6.6816e+8_f64); // 0.7 - 0.96 solar radius
                    let temp = rand::thread_rng().gen_range(3_700, 5_201); // Kelvin

                    (mass, radius, temp)
                },
            StarClass::M => {
                    let mass = if rand::thread_rng().gen_range(0, 4) == 0 {
                        rand::thread_rng().gen_range(1.9884e+29_f64, 8.9478e+29_f64) // 0.1 - 0.45 solar masses
                    } else {
                        rand::thread_rng().gen_range(1.9884e+29_f64, 4.971e+29_f64) // 0.1 - 0.25 solar masses
                    };

                    let radius = if rand::thread_rng().gen_range(0, 4) == 0 {
                        rand::thread_rng().gen_range((mass/SUN_MASS-0.05)*SUN_RADIUS, 4.872e+8_f64) // solar_masses-0.05 - 0.7 solar radius
                    } else {
                        rand::thread_rng().gen_range((mass/SUN_MASS-0.05)*SUN_RADIUS, 3.48e+8_f64) // solar_masses-0.05 - 0.5 solar radius
                    };

                    let temp = rand::thread_rng().gen_range(2_500, 3_701); // Kelvin

                    (mass, radius, temp)
                },
        }
    }

    /// Generates a random number based on the star class and luminosity, that will be the number of
    /// bodies in the solar system.
    pub fn generate_num_bodies(&self)  -> u8 {
        match self.class {
            StarClass::O | StarClass::BlackHole |
            StarClass::NeutronStar | StarClass::QuarkStar => 0,
            StarClass::WhiteDwarf | StarClass::B => {
                if self.get_luminosity() < 300_f64*SUN_LUMINOSITY {
                    rand::thread_rng().gen_range(0, 3)
                } else {0}
            },
            StarClass::A => rand::thread_rng().gen_range(0, 4),
            StarClass::F => rand::thread_rng().gen_range(0, 9),
            StarClass::G => {
                if rand::thread_rng().gen_range(0, 11) != 0 {
                    rand::thread_rng().gen_range(4, 11)
                } else {
                    rand::thread_rng().gen_range(0, 4)
                }
            },
            StarClass::K => {
                if rand::thread_rng().gen_range(0, 16) != 0 {
                    rand::thread_rng().gen_range(4, 11)
                } else {
                    rand::thread_rng().gen_range(0, 4)
                }
            },
            StarClass::M => {
                if rand::thread_rng().gen_range(0, 21) != 0 {
                    rand::thread_rng().gen_range(4, 10)
                } else {
                    rand::thread_rng().gen_range(0, 5)
                }
            }
        }
    }

    /// Generates Titius-Bode law *m* and *n* parameters. It supposes that bodies > 0, or it will
    /// panic.
    pub fn generate_titius_bode(&self, bodies: u8) -> (f64, f64) {
        if bodies > 0 {
            let luminosity = self.get_luminosity();
            let m = if luminosity < 0.01*SUN_LUMINOSITY {
                rand::thread_rng().gen_range(0.004_f64, 0.01_f64)
            } else if luminosity < 0.1*SUN_LUMINOSITY {
                rand::thread_rng().gen_range(0.0075_f64, 0.025_f64)
            } else if luminosity < 0.5*SUN_LUMINOSITY {
                rand::thread_rng().gen_range(0.015_f64, 0.1_f64)
            } else if luminosity < 6_f64*SUN_LUMINOSITY {
                if rand::thread_rng().gen_range(0, (10_f64-luminosity/SUN_LUMINOSITY) as u8) > 0 {
                    rand::thread_rng().gen_range(0.01_f64, 0.1_f64)
                } else {
                    rand::thread_rng().gen_range(0.025_f64, 0.6_f64)
                }
            } else if luminosity < 10_f64*SUN_LUMINOSITY {
                rand::thread_rng().gen_range(0.1*luminosity/SUN_LUMINOSITY, 0.3*luminosity/SUN_LUMINOSITY)
            } else if luminosity < 30_f64*SUN_LUMINOSITY {
                if luminosity*0.15 > 3_f64 {
                    rand::thread_rng().gen_range(0.1*luminosity/SUN_LUMINOSITY, 3_f64)
                } else {
                    rand::thread_rng().gen_range(0.1*luminosity/SUN_LUMINOSITY, 0.15*luminosity/SUN_LUMINOSITY)
                }
            } else {
                rand::thread_rng().gen_range(2.5_f64, 5_f64)
            };

            let n = if m > 0.035 {
                m.sqrt()*1.7+0.165
            } else {
                0.017/m
            };

            (m, rand::thread_rng().gen_range(0.8*n, 1.2*n))
        } else {panic!("Tried to generate Titius-Bode parameters for 0 body system!")}
    }
}

#[cfg(test)]
mod tests {
    use super::Star;
    use super::StarClass;
    use super::super::consts::*;

    #[test]
    fn it_star_getters() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::G, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};

        assert_eq!(2, st.get_id());
        assert_eq!(5, st.get_galaxy_id());
        assert_eq!(26_000, st.get_orbit());
        assert_eq!(&StarClass::G, st.get_class());
        assert_eq!(1.9885e+30_f64, st.get_mass());
        assert_eq!(6.96e+8_f64, st.get_radius());
        assert_eq!(5_778, st.get_temp());
    }

    #[test]
    fn it_star_parameters() {
        let st = Star::new(1, 3);

        assert_eq!(2, st.get_id());
        assert_eq!(3, st.get_galaxy_id());
    }

    #[test]
    fn it_luminosity() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::G,
            mass: 1.9885e+30_f64, radius: 6.96e+8_f64, temp: 5_778};

        assert!(st.get_luminosity() > 384.6e+24_f64*0.999 && st.get_luminosity() < 384.6e+24_f64*1.001);
    }

    #[test]
    fn it_volume() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::G,
            mass: 1.9885e+30_f64, radius: 6.96e+8_f64, temp: 5_778};

        assert!(st.get_volume() > 1.412e+27_f64*0.999 && st.get_volume() < 1.412e+27_f64*1.001);
    }

    #[test]
    fn it_density() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::G,
            mass: 1.9885e+30_f64, radius: 6.96e+8_f64, temp: 5_778};

        assert!(st.get_density() > 1_408_f64*0.999 && st.get_density() < 1_408_f64*1.001);
    }

    #[test]
    fn it_bh_properties() { // Black hole test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::BlackHole);
            assert!(mass >= 2_f64*SUN_MASS && mass <= 20_f64*SUN_MASS);
            assert_eq!(2_f64*G*mass/C.powi(2), radius);
            assert_eq!(0, temp);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert_eq!(0, st.generate_num_bodies());
        }
    }

    #[test]
    fn it_ns_properties() { // Neutron star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::NeutronStar);
            assert!(mass >= 1.1_f64*SUN_MASS && mass <= 2.2_f64*SUN_MASS);
            assert!(radius >= 11_000_f64 && radius <= 15_000_f64);
            assert!(radius > 2_f64*G*mass/C.powi(2));
            assert!(temp >= 10_000 && temp <= 10_000_000);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert_eq!(0, st.generate_num_bodies());
        }
    }

    #[test]
    fn it_qs_properties() { // Quark star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::QuarkStar);
            assert!(mass >= 2_f64*SUN_MASS && mass <= 3_f64*SUN_MASS);
            assert!(radius >= 11_000_f64 && radius <= 15_000_f64);
            assert!(radius > 2_f64*G*mass/C.powi(2));
            assert!(temp >= 8_000_000 && temp <= 100_000_000);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert_eq!(0, st.generate_num_bodies());
        }
    }

    #[test]
    fn it_wd_properties() { // White dwarf test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::WhiteDwarf);
            assert!(mass >= 0.2_f64*SUN_MASS && mass <= 1.3_f64*SUN_MASS);
            assert_eq!(0.78e+7_f64*((CH_LIMIT/mass).powf(2_f64/3_f64) - (mass/CH_LIMIT).powf(2_f64/3_f64)).sqrt(), radius);
            assert!(temp >= 4_000 && temp <= 150_000);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert!(st.generate_num_bodies() >= 0 && st.generate_num_bodies() <= 2);
        }
    }

    #[test]
    fn it_o_properties() { // O star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::O);
            assert!(mass >= 15_f64*SUN_MASS && mass <= 90_f64*SUN_MASS);
            assert!(radius >= 6.6_f64*SUN_RADIUS && radius <= 40_f64*SUN_RADIUS);
            assert!(temp >= 30_000 && temp <= 52_000);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert_eq!(0, st.generate_num_bodies());
        }
    }

    #[test]
    fn it_b_properties() { // B star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::B);
            assert!(mass >= 2.1_f64*SUN_MASS && mass <= 16_f64*SUN_MASS);
            assert!(radius >= 1.8_f64*SUN_RADIUS && radius <= 6.6_f64*SUN_RADIUS);
            assert!(temp >= 10_000 && temp <= 30_000);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert!(st.generate_num_bodies() >= 0 && st.generate_num_bodies() <= 2);
        }
    }

    #[test]
    fn it_a_properties() { // A star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::A);
            assert!(mass >= 1.4_f64*SUN_MASS && mass <= 2.1_f64*SUN_MASS);
            assert!(radius >= 1.4_f64*SUN_RADIUS && radius <= 1.8_f64*SUN_RADIUS);
            assert!(temp >= 7_500 && temp <= 10_000);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert!(st.generate_num_bodies() >= 0 && st.generate_num_bodies() <= 3);
        }
    }

    #[test]
    fn it_f_properties() { // F star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::F);
            assert!(mass >= 1.04_f64*SUN_MASS && mass <= 1.4_f64*SUN_MASS);
            assert!(radius >= 1.15_f64*SUN_RADIUS && radius <= 1.4_f64*SUN_RADIUS);
            assert!(temp >= 6_000 && temp <= 7_500);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert!(st.generate_num_bodies() >= 0 && st.generate_num_bodies() <= 8);
        }
    }

    #[test]
    fn it_g_properties() { // G star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::G);
            assert!(mass >= 0.8_f64*SUN_MASS && mass <= 1.04_f64*SUN_MASS);
            assert!(radius >= 0.96_f64*SUN_RADIUS && radius <= 1.15_f64*SUN_RADIUS);
            assert!(temp >= 5_200 && temp <= 6_000);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert!(st.generate_num_bodies() >= 0 && st.generate_num_bodies() <= 10);
        }
    }

    #[test]
    fn it_k_properties() { // K star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::K);
            assert!(mass >= 0.45_f64*SUN_MASS && mass <= 0.8_f64*SUN_MASS);
            assert!(radius >= 0.7_f64*SUN_RADIUS && radius <= 0.96_f64*SUN_RADIUS);
            assert!(temp >= 3_700 && temp <= 5_200);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert!(st.generate_num_bodies() >= 0 && st.generate_num_bodies() <= 10);
        }
    }

    #[test]
    fn it_m_properties() { // M star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(&StarClass::M);
            assert!(mass >= 0.1_f64*SUN_MASS && mass <= 0.45_f64*SUN_MASS);
            assert!(radius <= 0.7_f64*SUN_RADIUS);
            assert!(temp >= 2_400 && temp <= 3_700);

            let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::NeutronStar,
                mass: mass, radius: radius, temp: temp};

            assert!(st.generate_num_bodies() >= 0 && st.generate_num_bodies() <= 9);
        }
    }

    #[test]
    #[should_panic]
    fn it_generate_titius_bode_fail() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: StarClass::G,
            mass: 1.9885e+30_f64, radius: 6.96e+8_f64, temp: 5_778};

        st.generate_titius_bode(0);
    }
}
