//! Star implementation

extern crate rand;

use self::rand::Rng;
use std::f64::consts::PI;
use super::consts::*;

/// Star structure
pub struct Star {
    id: u64,
    galaxy_id: u32,
    orb_radius: u32,
    class: u32,
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
        let (star_mass, star_radius, star_temp) = Star::generate_properties(star_class);

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
    pub fn get_class(&self) -> u32 {
        self.class
    }

    /// Gets the human readable star class.
    pub fn get_class_str(&self) -> &'static str { // Use enum
        match self.class {
            0   => "Black hole",
            1   => "Neutron star",
            2   => "Quark star",
            3   => "White dwarf",
            4   => "O",
            5   => "B",
            6   => "A",
            7   => "F",
            8   => "G",
            9   => "K",
            10  => "M",
            _ => unreachable!()
        }
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

    fn generate_class() -> u32 { // Maybe here could be u8, must test performance.
        let prob = rand::thread_rng().gen_range(1, 10000001);

        // Numbers have to change meaning
        match prob {
            0...3_333               => 0, // Black hole
            3_334...6_666           => 1, // Neutron star
            6_667                   => 2, // Quark star
            6_668...10_000          => 3, // White dwarf
            10_001...10_003         => 4, // Class O
            10_004...23_000         => 5, // Class B
            23_001...83_000         => 6, // Class A
            83_001...383_000        => 7, // Class F
            383_001...1_143_000     => 8, // Class G
            1_143_001...2_353_000   => 9, // Class K
            2_353_001...10_000_000  => 10, // Class M
            _ => unreachable!()
        }
    }

    fn generate_properties(class: u32) -> (f64, f64, u32) {
        match class {
            0   => { // Black hole
                    let mass = if rand::thread_rng().gen_range(0, 2) == 0 { // TODO gaussian distribution
                        rand::thread_rng().gen_range(9.942e+30_f64, 1.9884e+31_f64) // 5 - 10 solar masses
                    } else {
                        rand::thread_rng().gen_range(3.9768e+30_f64, 3.9768e+31_f64) // 2 - 20 solar masses
                    };

                    let radius = 2_f64*G*mass/C.powi(2); // m
                    let temp = 0; // N/A

                    (mass, radius, temp)
                },
            1   => { // Neutron star
                    let mass = rand::thread_rng().gen_range(2.18724e+30_f64, 4.37448e+30_f64); // 1.1 - 2.2 solar masses

                    let radius = rand::thread_rng().gen_range(11_000_f64, 15_000_f64); // m
                    let temp = rand::thread_rng().gen_range(100_000, 10_000_000); // Kelvin

                    (mass, radius, temp)
                },
            2   => { // Quark star
                    let mass = rand::thread_rng().gen_range(3.9768e+30_f64, 5.9652e+30_f64); // 2 - 3 solar masses

                    let radius = rand::thread_rng().gen_range(11_000_f64, 15_000_f64); // m
                    let temp = rand::thread_rng().gen_range(8_000_000, 100_000_000); // Kelvin (N/A)

                    (mass, radius, temp)
                },
            3   => { // White dwarf
                    let mass = rand::thread_rng().gen_range(9.942e+29_f64, 1.59072e+31_f64); // 0.5 - 8 solar masses

                    let radius = (3_f64*mass/(4.0e+9_f64*PI)).powf(1_f64/3_f64);

                    let temp = if rand::thread_rng().gen_range(0,4) == 0 {
                        rand::thread_rng().gen_range(4_000, 150_000) // Kelvin
                    } else {
                        rand::thread_rng().gen_range(6_000, 30_000) // Kelvin
                    };

                    (mass, radius, temp)
                },
            4   => { // Class O
                    let mass = rand::thread_rng().gen_range(2.9826e+31_f64, 1.78956e+32_f64); // 15 - 90 solar masses

                    let radius = rand::thread_rng().gen_range(4.5936e+9_f64, 2.784e+10_f64); // 6.6 - 40 solar radius
                    let temp = rand::thread_rng().gen_range(30_000, 56_000); // Kelvin

                    (mass, radius, temp)
                },
            5   => { // Class B
                    let mass = rand::thread_rng().gen_range(4.17564e+30_f64, 3.18144e+31_f64); // 2.1 - 16 solar masses

                    let radius = rand::thread_rng().gen_range(1.2528e+9_f64, 4.5936e+9_f64); // 1.8 - 6.6 solar radius
                    let temp = rand::thread_rng().gen_range(10_000, 30_000); // Kelvin

                    (mass, radius, temp)
                },
            6   => { // Class A
                    let mass = rand::thread_rng().gen_range(2.78376e+30_f64, 4.17564e+30_f64); // 1.4 - 2.1 solar masses

                    let radius = rand::thread_rng().gen_range(9.744e+8_f64, 1.2528e+9_f64); // 1.4 - 1.8 solar radius
                    let temp = rand::thread_rng().gen_range(7_500, 10_000); // Kelvin

                    (mass, radius, temp)
                },
            7   => { // Class F
                    let mass = rand::thread_rng().gen_range(2.067936e+30_f64, 2.78376e+30_f64); // 1.04 - 1.4 solar masses

                    let radius = rand::thread_rng().gen_range(8.004e+8_f64, 9.744e+8_f64); // 1.15 - 1.4 solar radius
                    let temp =rand::thread_rng().gen_range(6_000, 7_500); // Kelvin

                    (mass, radius, temp)
                },
            8   => { // Class G
                    let mass = rand::thread_rng().gen_range(1.59072e+30_f64, 2.067936e+30_f64); // 0.8 - 1.04 solar masses

                    let radius = rand::thread_rng().gen_range(6.6816e+8_f64, 8.004e+8_f64); // 0.96 - 1.15 solar radius
                    let temp = rand::thread_rng().gen_range(5_200, 6_000); // Kelvin

                    (mass, radius, temp)
                },
            9   => { // Class K
                    let mass = rand::thread_rng().gen_range(8.9478e+29_f64, 1.59072e+30_f64); // 0.45 - 0.8 solar masses

                    let radius = rand::thread_rng().gen_range(4.872e+8_f64, 6.6816e+8_f64); // 0.7 - 0.96 solar radius
                    let temp = rand::thread_rng().gen_range(3_700, 5_200); // Kelvin

                    (mass, radius, temp)
                },
            10  => { // Class M
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

                    let temp = rand::thread_rng().gen_range(2_500, 3_700); // Kelvin

                    (mass, radius, temp)
                },
            _ => unreachable!()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Star;
    use super::super::consts::*;

    #[test]
    fn it_star_getters() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 8, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};

        assert_eq!(2, st.get_id());
        assert_eq!(5, st.get_galaxy_id());
        assert_eq!(26_000, st.get_orbit());
        assert_eq!(8, st.get_class());
        assert_eq!("G", st.get_class_str());
        assert_eq!(1.9885e+30_f64, st.get_mass());
        assert_eq!(6.96e+8_f64, st.get_radius());
        assert_eq!(5_778, st.get_temp());
    }

    #[test]
    fn it_class_str() {
        let st_bh = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 0, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("Black hole", st_bh.get_class_str());

        let st_ns = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 1, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("Neutron star", st_ns.get_class_str());

        let st_qs = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 2, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("Quark star", st_qs.get_class_str());

        let st_wd = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 3, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("White dwarf", st_wd.get_class_str());

        let st_o = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 4, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("O", st_o.get_class_str());

        let st_b = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 5, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("B", st_b.get_class_str());

        let st_a = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 6, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("A", st_a.get_class_str());

        let st_f = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 7, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("F", st_f.get_class_str());

        let st_g = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 8, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("G", st_g.get_class_str());

        let st_k = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 9, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("K", st_k.get_class_str());

        let st_m = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 10, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};
        assert_eq!("M", st_m.get_class_str());
    }

    #[test]
    fn it_star_parameters() {
        let st = Star::new(1, 3);

        assert_eq!(2, st.get_id());
        assert_eq!(3, st.get_galaxy_id());
    }

    #[test]
    fn it_luminosity() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 8, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};

        assert!(384.6e+24_f64*0.999 < st.get_luminosity() && 384.6e+24_f64*1.001 > st.get_luminosity());
    }

    #[test]
    fn it_volume() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 8, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};

        assert!(1.412e+27_f64*0.999 < st.get_volume() && 1.412e+27_f64*1.001 > st.get_volume());
    }

    #[test]
    fn it_density() {
        let st = Star {id: 2, galaxy_id: 5, orb_radius: 26_000, class: 8, mass: 1.9885e+30_f64,
            radius: 6.96e+8_f64, temp: 5_778};

        assert!(1_408_f64*0.999 < st.get_density() && 1_408_f64*1.001 > st.get_density());
    }

    #[test]
    fn it_bh_properties() { // Black hole test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(0);
            assert!(mass < 20_f64*SUN_MASS && mass > 2_f64*SUN_MASS);
            assert_eq!(2_f64*G*mass/C.powi(2), radius);
            assert_eq!(0, temp);
        }
    }

    #[test]
    fn it_ns_properties() { // Neutron star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(1);
            assert!(mass < 2.2_f64*SUN_MASS && mass > 1.1_f64*SUN_MASS);
            assert!(radius > 11_000_f64 && radius < 15_000_f64);
            assert!(radius > 2_f64*G*mass/C.powi(2));
            assert!(temp > 10_000 && temp < 10_000_000);
        }
    }

    #[test]
    fn it_qs_properties() { // Quark star test
        for _i in 1..10_000 {
            let (mass, radius, temp) = Star::generate_properties(2);
            assert!(mass < 3_f64*SUN_MASS && mass > 2_f64*SUN_MASS);
            assert!(radius > 11_000_f64 && radius < 15_000_f64);
            assert!(radius > 2_f64*G*mass/C.powi(2));
            assert!(temp > 8_000_000 && temp < 100_000_000);
        }
    }
}
