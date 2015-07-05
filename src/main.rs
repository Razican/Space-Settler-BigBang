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
	pub const SUN_MASS: f64 = 1.9884e+30; // kg
	pub const SUN_RADIUS: f64 = 6.96e+8; // m
	pub const SUN_LUMINOSITY: f64 = 3.846e+26; // W

    pub const G: f64 = 6.67428e-11; // m³·kg⁻¹·s⁻²
    pub const C: f64 = 299_792_458_f64; // m·s⁻¹
    pub const BOLTZ: f64 = 5.670373e-8; // W·m⁻²·K⁻⁴

    pub const CH_LIMIT: f64 = 2.765e+30; // kg
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
	let star = Star::new(0, 1);
	println!("Created a new star:");

    println!("\tID: {}", star.get_id());
    println!("\tGalaxy: {}", star.get_galaxy_id());
    println!("\tOrbit: {} light years", star.get_orbit());
    println!("\tClass: {:?}", star.get_class());
    println!("\tMass: {} suns", star.get_mass()/consts::SUN_MASS);
    println!("\tRadius: {} suns", star.get_radius()/consts::SUN_RADIUS);
    println!("\tDensity: {} kg/m³", star.get_density());
    println!("\tTemperature: {} K", star.get_temp());
    println!("\tLuminosity: {} suns", star.get_luminosity()/consts::SUN_LUMINOSITY);

    let num_bodies = star.calculate_num_bodies();
    let (tb_m, tb_n) = star.calculate_titius_bode(num_bodies);

    if num_bodies > 0 {
        let planet = Planet::new(&star, tb_m, tb_n, 1);
        println!("\nCreated a new Planet:");

        println!("\tStar ID: {}", planet.get_orbit().get_star().get_id());
        println!("\tAlbedo: {}", planet.get_albedo());
        println!("\tOrbit:");
        println!("\t\tSemimajor axis: {:e} meters", planet.get_orbit().get_sm_a());
        println!("\t\tEccentricity: {}", planet.get_orbit().get_ecc());
        println!("\t\tApoapsis: {:e} meters", planet.get_orbit().get_apoapsis());
        println!("\t\tPeriapsis: {:e} meters", planet.get_orbit().get_periapsis());
        println!("\t\tOrbital period: {} days", planet.get_orbit().get_orb_period()/planet.get_orbit().get_day());
        println!("\t\tDay length: {} hours", planet.get_orbit().get_day()/3_600_f64);
        println!("\t\tAxial tilt: {}°", rad_to_deg(planet.get_orbit().get_ax_tilt()));
        println!("\tAtmosphere:");
        println!("\t\tPressure: {} Pa", planet.get_atmosphere().get_pressure());
        println!("\t\tCarbon dioxide (CO₂): {}%", planet.get_atmosphere().get_co2());
        println!("\t\tCarbon monoxide (CO): {}%", planet.get_atmosphere().get_co());
        println!("\t\tNitrogen (N₂): {}%", planet.get_atmosphere().get_n2());
        println!("\t\tOxygen (O₂): {}%", planet.get_atmosphere().get_o2());
        println!("\t\tArgon (Ar): {}%", planet.get_atmosphere().get_ar());
        println!("\t\tSulfur dioxide (SO₂): {}%", planet.get_atmosphere().get_so2());
        println!("\t\tNeon (Ne): {}%", planet.get_atmosphere().get_ne());
        println!("\t\tMethane (CH₄): {}%", planet.get_atmosphere().get_ch4());
        println!("\t\tHelium (He): {}%", planet.get_atmosphere().get_he());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use super::consts::*;

    #[test]
    fn it_main() {
        super::main();
    }

    #[test]
    fn it_rad_to_deg() {
        assert_eq!(360_f64, rad_to_deg(2_f64*PI));
    }

    #[test]
    fn it_deg_to_rad() {
        assert_eq!(2_f64*PI, deg_to_rad(360_f64));
    }

    #[test]
    fn it_constants() {
        assert_eq!(1.9884E+30_f64, SUN_MASS);
        assert_eq!(6.96E+8_f64, SUN_RADIUS);
        assert_eq!(3.846E+26_f64, SUN_LUMINOSITY);

        assert_eq!(6.67428e-11_f64, G);
        assert_eq!(299_792_458_f64, C);
        assert_eq!(5.670373E-8_f64, BOLTZ)
    }
}
