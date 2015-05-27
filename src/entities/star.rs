extern crate rand;

use self::rand::Rng;
use std::f64::consts::PI;
use super::super::consts::*;

const G: f64 = 6.67428e-11_f64; // m³·kg⁻¹·s⁻²
const C: f64 = 299_792_458_f64; // m·s⁻¹
const BOLTZ: f64 = 5.670373E-8_f64; // W·m⁻²·K⁻⁴

pub struct Star {
    id: i32,
    galaxy_id: i32,
    orb_radius: i32,
    class: usize,
    mass: f64,
    radius: f64,
    temp: u32,
}

impl Star {
    pub fn new(last_id:i32, galaxy: i32) -> Star {
        let orbit = rand::thread_rng().gen_range(20000, 30001);
        let star_class = Star::generate_class();
        let (star_mass, star_radius, star_temp) = Star::generate_properties(star_class);

        Star {id: last_id+1, galaxy_id: galaxy, orb_radius: orbit, class: star_class,
            mass: star_mass, radius: star_radius, temp: star_temp}
    }

    pub fn get_id(&self) -> i32 {
        self.id
    }

    pub fn get_galaxy_id(&self) -> i32 {
        self.galaxy_id
    }

    pub fn get_orbit(&self) -> i32 {
        self.orb_radius
    }

    pub fn get_class(&self) -> usize {
        self.class
    }

    pub fn get_class_str(&self) -> &'static str {
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

    pub fn get_mass(&self) -> f64 {
        self.mass
    }

    pub fn get_radius(&self) -> f64 {
        self.radius
    }

    pub fn get_density(&self) -> f64 {
        self.mass/self.get_volume()
    }

    pub fn get_temp(&self) -> u32 {
        self.temp
    }

    pub fn get_luminosity(&self) -> f64 {
        4_f64*PI*self.radius.powi(2)*BOLTZ*(self.temp as f64).powi(4) // W
    }

    // ----------  generators  ----------

    fn generate_class() -> usize { // Maybe here could be u8, must test performance
        let prob = rand::thread_rng().gen_range(1, 10000001);

        // Numbers have to change meaning
        match prob {
            0...3333            => 0, // Black hole
            3334...6666         => 1, // Neutron star
            6667                => 2, // Quark star
            6668...10000        => 3, // White dwarf
            10001...10003       => 4, // Class O
            10004...23000       => 5, // Class B
            23001...83000       => 6, // Class A
            83001...383000      => 7, // Class F
            383001...1143000    => 8, // Class G
            1143001...2353000   => 9, // Class K
            2353001...10000000  => 10, // Class M
            _ => unreachable!()
        }
    }

    fn generate_properties(class: usize) -> (f64, f64, u32) {
        match class {
            0   => { // Black hole
                    let mass = if rand::thread_rng().gen_range(0, 2) == 0 { // TODO gaussian distribution
                        rand::thread_rng().gen_range(9.942e+30_f64, 1.9884e+31_f64) // 5 - 10 solar masses
                    } else {
                        rand::thread_rng().gen_range(1.9884e+31_f64, 3.9768e+31_f64) // 3 - 20 solar masses
                    };

                    let radius = 2_f64*G*mass/C.powi(2); // m
                    let temp = 0; // N/A

                    (mass, radius, temp)
                },
            1   => { // Neutron star
                    let mass = rand::thread_rng().gen_range(2.743992e+30_f64, 3.9768e+30_f64); // 1.38 - 2 solar masses

                    let radius = rand::thread_rng().gen_range(11_000_f64, 13_000_f64); // m
                    let temp = rand::thread_rng().gen_range(100_000, 10_000_000); // Kelvin

                    (mass, radius, temp)
                },
            2   => { // Quark star
                    let mass = rand::thread_rng().gen_range(2.743992e+30_f64, 1.9884e+31_f64); // 2 - 3 solar masses

                    let radius = rand::thread_rng().gen_range(11_000_f64, 13_000_f64); // m
                    let temp = 0; // N/A

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

    // ----------  private functions  ----------

    fn get_volume(&self) -> f64 {
        4_f64/3_f64*PI*self.radius.powi(3)
    }
    // TODO
}
