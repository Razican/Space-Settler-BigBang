//! Galaxy creator for Space Settler
//!
//! This program creates a new galaxy for Space Settler. It has a statistically correct sample of
//! stars, and in those stars it creates complete solar systems (still in development). It creates
//! stars regarding to proportion shown [here](http://adsabs.harvard.edu/abs/2001JRASC..95...32L).
//! Then, the planets are created thinking on habitability, but following as far as possible all the
//! physical laws.

extern crate rand;

pub mod planet;
pub mod star;
pub mod consts;
pub mod utils;

use std::env;

use self::rand::Rng;

use star::{Star, StarClass};
use planet::{Planet, PlanetType};
use consts::*;
use utils::*;

fn main() {
    let mut args = env::args();
    let arg = args.nth(1);

    if arg.is_some() {
        match arg.unwrap().as_ref() {
            "-e" | "--earth" => {
                loop {
                    let new_star: Star = Star::new(0, 1);
                    let num_bodies = new_star.generate_num_bodies();

                    if num_bodies > 1 {
                        let (new_tb_m, new_tb_n) = new_star.generate_titius_bode(num_bodies);
                        let new_planet: Planet = Planet::new(&new_star, new_tb_m, new_tb_n,
                            rand::thread_rng().gen_range(1, num_bodies), 0_f64);

                        if new_planet.is_roche_ok() && new_planet.is_earth_twin() {
                            print_star(&new_star);
                            print_planet(&new_planet);
                            break;
                        }
                    }
                }
            },
            "-g" | "--create-galaxy" => {
                let galaxy_id_arg = args.next();
                if galaxy_id_arg.is_some() {
                    let galaxy_id_parsed = galaxy_id_arg.unwrap().parse::<u32>();
                    if galaxy_id_parsed.is_err() {
                        println!("Please, provide a valid galaxy ID.")
                    } else {
                        let mut stats = Stats::new();
                        let galaxy_id = galaxy_id_parsed.unwrap();
                        let star_count = rand::thread_rng().gen_range(950_000, 1_005_000);
                        println!("Creating a galaxy of {} stars...", star_count);
                        println!("");

                        for i in 0..star_count {
                            print!("Creating star {} and its planets...\r", i);

                            let star = Star::new(i, galaxy_id);
                            stats.add_star(&star);
                            let num_bodies = star.generate_num_bodies();

                            if num_bodies > 0 {
                                let (tb_m, tb_n) = star.generate_titius_bode(num_bodies);
                                let mut last_distance = 0_f64;

                                for g in 1..num_bodies {
                                    let planet = Planet::new(&star, tb_m, tb_n, g, last_distance);
                                    stats.add_planet(&planet);
                                    last_distance = planet.get_orbit().get_sma();
                                }
                            }
                        }

                        println!("Finished the creation of the galaxy.");
                        print_stats(&stats);
                    }
                } else {
                    println!("You must provide the galaxy ID.");
                }
            },
            _ => show_help()
        }
    } else {
        show_help();
    }
}

fn show_help() {
    println!("You must include the action you want. Currently you have these options:");
    println!("-e, --earth\t\tGenerates lots of random stars and planets until it finds an earth twin.");
    println!("-g, --create-galaxy\tGenerates a complete new galaxy and shows statistics (not yet functional).");
}

fn print_star(star: &Star) {
    println!("Created a new star:");

    println!("\tID: {}", star.get_id());
    println!("\tGalaxy: {}", star.get_galaxy_id());
    println!("\tOrbit: {} light years", star.get_orbit());
    println!("\tClass: {:?}", star.get_class());
    println!("\tMass: {:.2} M☉", star.get_mass()/consts::SUN_MASS);
    println!("\tRadius: {:.2} R☉", star.get_radius()/consts::SUN_RADIUS);
    println!("\tDensity: {:.2} kg/m³", star.get_density());
    println!("\tTemperature: {} K", star.get_temperature());
    println!("\tLuminosity: {:.4} suns", star.get_luminosity()/consts::SUN_LUMINOSITY);
    println!("");
}

fn print_planet(planet: &Planet) {
    println!("\nCreated a new Planet:");

    println!("\tStar ID: {}", planet.get_orbit().get_star().get_id());
    println!("\tPosition in solar system: {}", planet.get_orbit().get_position());
    println!("\tPlanet type: {:?}", planet.get_type());
    println!("\tBond albedo: {:.2}%", planet.get_bond_albedo()*100_f64);
    println!("\tGeometric albedo: {:.2}%", planet.get_geometric_albedo()*100_f64);
    if planet.get_mass() > 50_f64*consts::EARTH_MASS {
        println!("\tMass: {:.3} Mj", planet.get_mass()/consts::JUPITER_MASS);
        println!("\tRadius: {:.3} Rj", planet.get_radius()/consts::JUPITER_RADIUS);
    } else {
        println!("\tMass: {:.3} M⊕", planet.get_mass()/consts::EARTH_MASS);
        println!("\tRadius: {:.3} R⊕", planet.get_radius()/consts::EARTH_RADIUS);
    }
    println!("\tDensity: {} kg/m³", planet.get_density());
    if planet.get_type() == &PlanetType::Rocky {
        println!("\tMinimum temperature: {:.2}°C", kelvin_to_celsius(planet.get_min_temp()));
        println!("\tAverage temperature: {:.2}°C", kelvin_to_celsius(planet.get_avg_temp()));
        println!("\tMaximum temperature: {:.2}°C", kelvin_to_celsius(planet.get_max_temp()));
    }
    println!("\tOrbit:");
    println!("\t\tSemimajor axis: {:.3} AU", planet.get_orbit().get_sma()/AU);
    println!("\t\tEccentricity: {:.4}", planet.get_orbit().get_ecc());
    println!("\t\tApoapsis: {:.3} AU", planet.get_orbit().get_apoapsis()/AU);
    println!("\t\tPeriapsis: {:.3} AU", planet.get_orbit().get_periapsis()/AU);
    println!("\t\tOrbital period: {:.3} days", planet.get_orbit().get_orb_period()/planet.get_orbit().get_day());
    println!("\t\tDay length: {:.3} hours", planet.get_orbit().get_day()/3_600_f64);
    println!("\t\tAxial tilt: {:.3}°", rad_to_deg(planet.get_orbit().get_ax_tilt()));
    if planet.get_atmosphere().is_some() {
        println!("\tAtmosphere:");
        println!("\t\tPressure: {:.2} Pa", planet.get_atmosphere().unwrap().get_pressure());
        println!("\t\tCarbon dioxide (CO₂): {:.2}%", planet.get_atmosphere().unwrap().get_co2()*100_f64);
        println!("\t\tCarbon monoxide (CO): {:.2}%", planet.get_atmosphere().unwrap().get_co()*100_f64);
        println!("\t\tNitrogen (N₂): {:.2}%", planet.get_atmosphere().unwrap().get_n2()*100_f64);
        println!("\t\tOxygen (O₂): {:.2}%", planet.get_atmosphere().unwrap().get_o2()*100_f64);
        println!("\t\tArgon (Ar): {:.2}%", planet.get_atmosphere().unwrap().get_ar()*100_f64);
        println!("\t\tSulfur dioxide (SO₂): {:.2}%", planet.get_atmosphere().unwrap().get_so2()*100_f64);
        println!("\t\tNeon (Ne): {:.2}%", planet.get_atmosphere().unwrap().get_ne()*100_f64);
        println!("\t\tMethane (CH₄): {:.2}%", planet.get_atmosphere().unwrap().get_ch4()*100_f64);
        println!("\t\tHelium (He): {:.2}%", planet.get_atmosphere().unwrap().get_he()*100_f64);
    }
    println!("");
}

fn print_stats(st: &Stats) {
    println!("Statistics:");
    println!("Total stars: {}", st.get_total_stars());
    println!("\tBlack holes: {} ({:.5}%)", st.get_black_holes(), st.get_black_hole_percent()*100_f64);
    println!("\tNeutron stars: {} ({:.5}%)", st.get_neutron_stars(), st.get_neutron_star_percent()*100_f64);
    println!("\tQuark stars: {} ({:.5}%)", st.get_quark_stars(), st.get_quark_star_percent()*100_f64);
    println!("\tWhite dwarfs: {} ({:.5}%)", st.get_white_dwarfs(), st.get_white_dwarf_percent()*100_f64);
    println!("");
    println!("\tO type stars: {} ({:.5}%)", st.get_o_stars(), st.get_o_star_percent()*100_f64);
    println!("\tB type stars: {} ({:.5}%)", st.get_b_stars(), st.get_b_star_percent()*100_f64);
    println!("\tA type stars: {} ({:.5}%)", st.get_a_stars(), st.get_a_star_percent()*100_f64);
    println!("\tF type stars: {} ({:.5}%)", st.get_f_stars(), st.get_f_star_percent()*100_f64);
    println!("\tG type stars: {} ({:.5}%)", st.get_g_stars(), st.get_g_star_percent()*100_f64);
    println!("\tK type stars: {} ({:.5}%)", st.get_k_stars(), st.get_k_star_percent()*100_f64);
    println!("\tM type stars: {} ({:.5}%)", st.get_m_stars(), st.get_m_star_percent()*100_f64);

    println!("Total planets: {}", st.get_total_planets());
    println!("\tRocky planets: {} ({:.2}%)", st.get_rocky_planets(), st.get_rocky_planet_percent()*100_f64);
    println!("\t\tSuper earths: {}", st.get_super_earths());
    println!("\t\tSmall planets: {}", st.get_small_planets());
    println!("\t\tEarth twins (not necessarily the only habitable ones): {}", st.get_earth_twins());
    println!("\tGaseous planets: {} ({:.2}%)", st.get_gaseous_planets(), st.get_gaseous_planet_percent()*100_f64);
    println!("\t\tHot jupiters: {}", st.get_hot_jupiters());
    println!("\tMini-neptunes (some gaseous, some rocky, some in between): {}", st.get_mini_neptunes());

    println!("Records:");
    println!("\tMinimum temperature: {:.2}K ({:.2}°C)", st.get_min_temp(), kelvin_to_celsius(st.get_min_temp()));
    println!("\tMaximum temperature: {:.2}K ({:.2}°C)", st.get_max_temp(), kelvin_to_celsius(st.get_max_temp()));
    println!("");

    println!("\tMinimum sMa: {:.4} AU", st.get_min_sma()/AU);
    println!("\tMaximum sMa: {:.2} AU", st.get_max_sma()/AU);

    println!("\tShortest year: {:.2} hours", st.get_min_orb_period()/3600_f64);
    println!("\tLongest year: {:.3} years", st.get_max_orb_period()/(3600_f64*24_f64*365.256363004));

    println!("\tShortest day: {:.2} hours", st.get_min_day()/3600_f64);
    println!("\tLongest day: {:.2} years", st.get_max_day()/(3600_f64*24_f64*365.256363004));

    println!("");
}

struct Stats {
    // Stars
    black_holes: u32,
    neutron_stars: u32,
    quark_stars: u32,
    white_dwarfs: u32,
    o_stars: u32,
    b_stars: u32,
    a_stars: u32,
    f_stars: u32,
    g_stars: u32,
    k_stars: u32,
    m_stars: u32,
    // Planets
    gaseous: u32,
    rocky: u32,
    small_planets: u32,
    super_earths: u32,
    earth_twins: u32,
    hot_jupiters: u32,
    mini_neptunes: u32,
    //Records
    min_temp: f64,
    max_temp: f64,
    min_sma: f64,
    max_sma: f64,
    min_orb_period: f64,
    max_orb_period: f64,
    min_day: f64,
    max_day: f64,
}

impl Stats {
    pub fn new() -> Stats {
        Stats {black_holes: 0, neutron_stars: 0, quark_stars: 0, white_dwarfs: 0, o_stars: 0,
            b_stars: 0, a_stars: 0, f_stars: 0, g_stars: 0, k_stars: 0, m_stars: 0, gaseous: 0,
            rocky: 0, small_planets: 0, super_earths: 0, earth_twins: 0, hot_jupiters: 0,
            mini_neptunes: 0, min_temp: 0_f64, max_temp: 0_f64, min_sma: 0_f64, max_sma: 0_f64,
            min_orb_period: 0_f64, max_orb_period: 0_f64, min_day: 0_f64, max_day: 0_f64}
    }

    pub fn add_star(&mut self, st: &Star) {
        match *st.get_class() {
            StarClass::BlackHole => self.black_holes += 1,
            StarClass::NeutronStar => self.neutron_stars += 1,
            StarClass::QuarkStar => self.quark_stars += 1,
            StarClass::WhiteDwarf => self.white_dwarfs += 1,
            StarClass::O => self.o_stars += 1,
            StarClass::B => self.b_stars += 1,
            StarClass::A => self.a_stars += 1,
            StarClass::F => self.f_stars += 1,
            StarClass::G => self.g_stars += 1,
            StarClass::K => self.k_stars += 1,
            StarClass::M => self.m_stars += 1,
        }
    }

    pub fn add_planet(&mut self, planet: &Planet) {
        if self.min_sma > planet.get_orbit().get_sma() || self.min_sma == 0_f64 {
            self.min_sma = planet.get_orbit().get_sma();
        }
        if self.max_sma < planet.get_orbit().get_sma() {
            self.max_sma = planet.get_orbit().get_sma();
        }
        if self.min_orb_period > planet.get_orbit().get_orb_period() || self.min_orb_period == 0_f64 {
            self.min_orb_period = planet.get_orbit().get_orb_period();
        }
        if self.max_orb_period < planet.get_orbit().get_orb_period() {
            self.max_orb_period = planet.get_orbit().get_orb_period();
        }
        if self.min_day > planet.get_orbit().get_day() || self.min_day == 0_f64 {
            self.min_day = planet.get_orbit().get_day();
        }
        if self.max_day < planet.get_orbit().get_day() {
            self.max_day = planet.get_orbit().get_day();
        }
        if planet.get_radius() > 2_f64*EARTH_RADIUS && planet.get_radius() < 4_f64*EARTH_RADIUS {
            self.mini_neptunes += 1;
        }
        if planet.is_earth_twin() {
            self.earth_twins +=1;
        }
        match *planet.get_type() {
            PlanetType::Rocky => {
                self.rocky += 1;
                if planet.get_radius() > 0.8*EARTH_RADIUS && planet.get_radius() < 1.25*EARTH_RADIUS {
                    self.super_earths += 1;
                } else if planet.get_radius() < 0.8*EARTH_RADIUS {
                    self.small_planets += 1;
                }

                if planet.get_min_temp() < 0_f64 {
                    print_planet(planet);
                }

                if self.min_temp > planet.get_min_temp() || self.min_temp == 0_f64 {
                    self.min_temp = planet.get_min_temp();
                }
                if self.max_temp < planet.get_max_temp() {
                    self.max_temp = planet.get_max_temp();
                }
            },
            PlanetType::Gaseous => {
                self.gaseous += 1;
                if planet.get_eff_temp() > 700_f64 {
                    self.hot_jupiters += 1;
                }
            },
        }
    }

    // Star getters

    pub fn get_total_stars(&self) -> u32 {
        self.black_holes + self.neutron_stars + self.quark_stars + self.white_dwarfs + self.o_stars +
        self.b_stars + self.a_stars + self.f_stars + self.g_stars + self.k_stars + self.m_stars
    }

    pub fn get_black_hole_percent(&self) -> f64 {
        (self.black_holes as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_black_holes(&self) -> u32 {
        self.black_holes
    }

    pub fn get_neutron_star_percent(&self) -> f64 {
        (self.neutron_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_neutron_stars(&self) -> u32 {
        self.neutron_stars
    }

    pub fn get_quark_star_percent(&self) -> f64 {
        (self.quark_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_quark_stars(&self) -> u32 {
        self.quark_stars
    }

    pub fn get_white_dwarf_percent(&self) -> f64 {
        (self.white_dwarfs as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_white_dwarfs(&self) -> u32 {
        self.white_dwarfs
    }

    pub fn get_o_star_percent(&self) -> f64 {
        (self.o_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_o_stars(&self) -> u32 {
        self.o_stars
    }

    pub fn get_b_star_percent(&self) -> f64 {
        (self.b_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_b_stars(&self) -> u32 {
        self.b_stars
    }

    pub fn get_a_star_percent(&self) -> f64 {
        (self.a_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_a_stars(&self) -> u32 {
        self.a_stars
    }

    pub fn get_f_star_percent(&self) -> f64 {
        (self.f_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_f_stars(&self) -> u32 {
        self.f_stars
    }

    pub fn get_g_star_percent(&self) -> f64 {
        (self.g_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_g_stars(&self) -> u32 {
        self.g_stars
    }

    pub fn get_k_star_percent(&self) -> f64 {
        (self.k_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_k_stars(&self) -> u32 {
        self.k_stars
    }

    pub fn get_m_star_percent(&self) -> f64 {
        (self.m_stars as f64) / (self.get_total_stars() as f64)
    }

    pub fn get_m_stars(&self) -> u32 {
        self.m_stars
    }

    // Planet getters

    pub fn get_total_planets(&self) -> u32 {
        self.gaseous+self.rocky
    }

    pub fn get_gaseous_planet_percent(&self) -> f64 {
        (self.gaseous as f64) / (self.get_total_planets() as f64)
    }

    pub fn get_gaseous_planets(&self) -> u32 {
        self.gaseous
    }

    pub fn get_hot_jupiters(&self) -> u32 {
        self.hot_jupiters
    }

    pub fn get_mini_neptunes(&self) -> u32 {
        self.mini_neptunes
    }

    pub fn get_rocky_planet_percent(&self) -> f64 {
        (self.rocky as f64) / (self.get_total_planets() as f64)
    }

    pub fn get_rocky_planets(&self) -> u32 {
        self.rocky
    }

    pub fn get_super_earths(&self) -> u32 {
        self.super_earths
    }

    pub fn get_earth_twins(&self) -> u32 {
        self.earth_twins
    }

    pub fn get_small_planets(&self) -> u32 {
        self.small_planets
    }

    // Record getters

    pub fn get_min_temp(&self) -> f64 {
        self.min_temp
    }

    pub fn get_max_temp(&self) -> f64 {
        self.max_temp
    }

    pub fn get_min_sma(&self) -> f64 {
        self.min_sma
    }

    pub fn get_max_sma(&self) -> f64 {
        self.max_sma
    }

    pub fn get_min_orb_period(&self) -> f64 {
        self.min_orb_period
    }

    pub fn get_max_orb_period(&self) -> f64 {
        self.max_orb_period
    }

    pub fn get_min_day(&self) -> f64 {
        self.min_day
    }

    pub fn get_max_day(&self) -> f64 {
        self.max_day
    }
}
