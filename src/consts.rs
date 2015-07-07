//! Constants
//!
//! This module contains all the needed constants in the galaxy creation.

/// The mass of the Sun, in *kg*
pub const SUN_MASS: f64 = 1.9884e+30;
/// The radius of the Sun, in *m*
pub const SUN_RADIUS: f64 = 6.96e+8;
/// The luminosity of the Sun, in *W*
pub const SUN_LUMINOSITY: f64 = 3.846e+26;

/// *G*, the gravitational constant, in *m³·kg⁻¹·s⁻²*
pub const G: f64 = 6.67428e-11;
/// *c*, the speed of light in vacuum, in *m·s⁻¹*
pub const C: f64 = 299_792_458_f64;
/// *σ*, the Stefan–Boltzmann constant, in *W·m⁻²·K⁻⁴*
pub const BOLTZ: f64 = 5.670373e-8;

/// The Chandrasekhar limit, in *kg*
pub const CH_LIMIT: f64 = 2.765e+30;

/// The astronomical unit (*au* or *AU*), the average distance from earth to the sun, in *m*
pub const AU: f64 = 149_597_870_700_f64;

#[cfg(test)]
mod tests {
    use super::*;
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
