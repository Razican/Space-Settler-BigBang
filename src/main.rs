//! Galaxy creator for Space Settler
//!
//! This program creates a new galaxy for Space Settler. It has a statistically correct sample of
//! stars, and in those stars it creates complete solar systems (still in development). It creates
//! stars regarding to proportion shown [here](http://adsabs.harvard.edu/abs/2001JRASC..95...32L).
//! Then, the planets are created thinking on habitability, but following as far as possible all the
//! physical laws.

mod planet;
mod star;
mod consts {
	pub const SUN_MASS: f64 = 1.9884E+30; // kg
	pub const SUN_RADIUS: f64 = 6.96E+8; // m
	pub const SUN_LUMINOSITY: f64 = 3.846E+26; // W

    pub const G: f64 = 6.67428e-11_f64; // m³·kg⁻¹·s⁻²
    pub const C: f64 = 299_792_458_f64; // m·s⁻¹
    pub const BOLTZ: f64 = 5.670373E-8_f64; // W·m⁻²·K⁻⁴
}

use std::f64::consts::PI;
use star::Star;
use planet::Planet;

pub fn rad_to_deg(rads: f64) -> f64 {
    rads*180_f64/PI
}

pub fn deg_to_rad(degs: f64) -> f64 {
    degs*PI/180_f64
}

fn main() {
	let sun = Star::new(0, 1);
	println!("Created a new star:");

    println!("\tID: {}", sun.get_id());
    println!("\tGalaxy: {}", sun.get_galaxy_id());
    println!("\tOrbit: {} light years", sun.get_orbit());
    println!("\tClass: {}", sun.get_class_str());
    println!("\tMass: {} suns", sun.get_mass()/consts::SUN_MASS);
    println!("\tRadius: {} suns", sun.get_radius()/consts::SUN_RADIUS);
    println!("\tDensity: {} kg/m³", sun.get_density());
    println!("\tTemperature: {} K", sun.get_temp());
    println!("\tLuminosity: {} suns", sun.get_luminosity()/consts::SUN_LUMINOSITY);

    let earth = Planet::new(&sun, 0.0183, 1.0643, 3);
    println!("\nCreated a new Planet:");

    println!("\tStar ID: {}", earth.get_orbit().get_star().get_id());
    println!("\tAlbedo: {}", earth.get_albedo());
    println!("\tOrbit:");
    println!("\t\tSemimajor axis: {:e} meters", earth.get_orbit().get_sMa());
    println!("\t\tEccentricity: {}", earth.get_orbit().get_ecc());
    println!("\t\tApoapsis: {:e} meters", earth.get_orbit().get_apoapsis());
    println!("\t\tPeriapsis: {:e} meters", earth.get_orbit().get_periapsis());
    println!("\t\tOrbital period: {} days", earth.get_orbit().get_orb_period()/earth.get_orbit().get_day());
    println!("\t\tDay length: {} hours", earth.get_orbit().get_day()/3_600_f64);
    println!("\t\tAxial tilt: {}°", rad_to_deg(earth.get_orbit().get_ax_tilt()));
    println!("\tAtmosphere:");
    println!("\t\tPressure: {} Pa", earth.get_atmosphere().get_pressure());
    println!("\t\tCarbon dioxide (CO₂): {}%", earth.get_atmosphere().get_co2());
    println!("\t\tCarbon monoxide (CO): {}%", earth.get_atmosphere().get_co());
    println!("\t\tNitrogen (N₂): {}%", earth.get_atmosphere().get_n2());
    println!("\t\tOxygen (O₂): {}%", earth.get_atmosphere().get_o2());
    println!("\t\tArgon (Ar): {}%", earth.get_atmosphere().get_ar());
    println!("\t\tSulfur dioxide (SO₂): {}%", earth.get_atmosphere().get_so2());
    println!("\t\tNeon (Ne): {}%", earth.get_atmosphere().get_ne());
    println!("\t\tMethane (CH₄): {}%", earth.get_atmosphere().get_ch4());
    println!("\t\tHelium (He): {}%", earth.get_atmosphere().get_he());
}
