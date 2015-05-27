mod entities;
mod consts {
	pub const SUN_MASS: f64 = 1.9884E+30; // kg
	pub const SUN_RADIUS: f64 = 6.96E+8; // m
	pub const SUN_LUMINOSITY: f64 = 3.846E+26; // W
}

use entities::star::*;

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
}
