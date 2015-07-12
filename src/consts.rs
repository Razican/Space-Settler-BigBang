//! Constants
//!
//! This module contains all the needed constants in the galaxy creation.

/// The mass of the Sun (M<sub>☉</sub>), in *kg*.
pub const SUN_MASS: f64 = 1.9884e+30;
/// The radius of the Sun (R<sub>☉</sub>), in meters (*m*).
pub const SUN_RADIUS: f64 = 6.96e+8;
/// The luminosity of the Sun (L<sub>☉</sub>), in watts (*W*).
pub const SUN_LUMINOSITY: f64 = 3.846e+26;

/// The mass of the Earth (M<sub>⊕</sub>), in *kg*.
pub const EARTH_MASS: f64 = 5.9726e+24;
/// The radius of the Earth (R<sub>⊕</sub>), in meters (*m*).
pub const EARTH_RADIUS: f64 = 6_371_000_f64;
/// The atmospheric pressure at sea level in Earth, in pascals (*Pa*).
pub const EARTH_ATM_PRESSURE: f64 = 101_325_f64;

/// The mass of Jupiter (M<sub>j</sub>), in *kg*.
pub const JUPITER_MASS: f64 = 1.898e+27;
/// The radius of Jupiter (R<sub>j</sub>), in meters (*m*).
pub const JUPITER_RADIUS: f64 = 69_911_000_f64;

/// *G*, the gravitational constant, in *m³·kg⁻¹·s⁻²*.
pub const G: f64 = 6.67428e-11;
/// *c*, the speed of light in vacuum, in *m·s⁻¹*.
pub const C: f64 = 299_792_458_f64;
/// *σ*, the Stefan–Boltzmann constant, in *W·m⁻²·K⁻⁴*.
pub const BOLTZ: f64 = 5.670373e-8;

/// The Chandrasekhar limit, in *kg*.
pub const CH_LIMIT: f64 = 2.765e+30;

/// The astronomical unit (*au* or *AU*), the average distance from earth to the sun, in meters
/// (*m*).
pub const AU: f64 = 149_597_870_700_f64;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_constants() {
        assert_eq!(1.9884E+30_f64, SUN_MASS);
        assert_eq!(6.96E+8_f64, SUN_RADIUS);
        assert_eq!(3.846E+26_f64, SUN_LUMINOSITY);

        assert_eq!(5.9726e+24_f64, EARTH_MASS);
        assert_eq!(6_371_000_f64, EARTH_RADIUS);

        assert_eq!(1.898e+27_f64, JUPITER_MASS);
        assert_eq!(69_911_000_f64, JUPITER_RADIUS);

        assert_eq!(6.67428e-11_f64, G);
        assert_eq!(299_792_458_f64, C);
        assert_eq!(5.670373E-8_f64, BOLTZ);

        assert_eq!(2.765e+30_f64, CH_LIMIT);

        assert_eq!(149_597_870_700_f64, AU);
    }
}
