//! Galaxy creator for Space Settler
//!
//! This program creates a new galaxy for Space Settler. It has a statistically correct sample of
//! stars, and in those stars it creates complete solar systems (still in development). It creates
//! stars regarding to proportion shown [here](http://adsabs.harvard.edu/abs/2001JRASC..95...32L).
//! Then, the planets are created thinking on habitability, but following as far as possible all the
//! physical laws.

pub mod planet;
pub mod star;
pub mod consts;
pub mod utils;

use star::Star;
use planet::Planet;
use consts::*;
use utils::*;

fn main() {
	let star = Star::new(0, 1);
	println!("Created a new star:");

    println!("\tID: {}", star.get_id());
    println!("\tGalaxy: {}", star.get_galaxy_id());
    println!("\tOrbit: {} light years", star.get_orbit());
    println!("\tClass: {:?}", star.get_class());
    println!("\tMass: {} M☉", star.get_mass()/consts::SUN_MASS);
    println!("\tRadius: {} R☉", star.get_radius()/consts::SUN_RADIUS);
    println!("\tDensity: {} kg/m³", star.get_density());
    println!("\tTemperature: {} K", star.get_temperature());
    println!("\tLuminosity: {} suns", star.get_luminosity()/consts::SUN_LUMINOSITY);

    let num_bodies = star.generate_num_bodies();

    if num_bodies > 2 {
        let (tb_m, tb_n) = star.generate_titius_bode(num_bodies);
        let planet = Planet::new(&star, tb_m, tb_n, 3, 0_f64);
        println!("\nCreated a new Planet:");

        println!("\tStar ID: {}", planet.get_orbit().get_star().get_id());
        println!("\tPosition in solar system: {}", planet.get_orbit().get_position());
        println!("\tPlanet type: {:?}", planet.get_type());
        println!("\tAlbedo: {}", planet.get_albedo());
        if planet.get_mass() > 50_f64*consts::EARTH_MASS {
            println!("\tMass: {} Mj", planet.get_mass()/consts::JUPITER_MASS);
            println!("\tRadius: {} Rj", planet.get_radius()/consts::JUPITER_RADIUS);
        } else {
            println!("\tMass: {} M⊕", planet.get_mass()/consts::EARTH_MASS);
            println!("\tRadius: {} R⊕", planet.get_radius()/consts::EARTH_RADIUS);
        }
        println!("\tDensity: {} kg/m³", planet.get_density());
        println!("\tOrbit:");
        println!("\t\tSemimajor axis: {} AU", planet.get_orbit().get_sma()/AU);
        println!("\t\tEccentricity: {}", planet.get_orbit().get_ecc());
        println!("\t\tApoapsis: {} AU", planet.get_orbit().get_apoapsis()/AU);
        println!("\t\tPeriapsis: {} AU", planet.get_orbit().get_periapsis()/AU);
        println!("\t\tOrbital period: {} days", planet.get_orbit().get_orb_period()/planet.get_orbit().get_day());
        println!("\t\tDay length: {} hours", planet.get_orbit().get_day()/3_600_f64);
        println!("\t\tAxial tilt: {}°", rad_to_deg(planet.get_orbit().get_ax_tilt()));
        if planet.get_atmosphere().is_some() {
            println!("\tAtmosphere:");
            println!("\t\tPressure: {} Pa", planet.get_atmosphere().unwrap().get_pressure());
            println!("\t\tCarbon dioxide (CO₂): {}%", planet.get_atmosphere().unwrap().get_co2()*100_f64);
            println!("\t\tCarbon monoxide (CO): {}%", planet.get_atmosphere().unwrap().get_co()*100_f64);
            println!("\t\tNitrogen (N₂): {}%", planet.get_atmosphere().unwrap().get_n2()*100_f64);
            println!("\t\tOxygen (O₂): {}%", planet.get_atmosphere().unwrap().get_o2()*100_f64);
            println!("\t\tArgon (Ar): {}%", planet.get_atmosphere().unwrap().get_ar()*100_f64);
            println!("\t\tSulfur dioxide (SO₂): {}%", planet.get_atmosphere().unwrap().get_so2()*100_f64);
            println!("\t\tNeon (Ne): {}%", planet.get_atmosphere().unwrap().get_ne()*100_f64);
            println!("\t\tMethane (CH₄): {}%", planet.get_atmosphere().unwrap().get_ch4()*100_f64);
            println!("\t\tHelium (He): {}%", planet.get_atmosphere().unwrap().get_he()*100_f64);
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_main() {
        super::main();
    }
}
