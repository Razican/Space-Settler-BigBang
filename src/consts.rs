//! Constants
//!
//! This module contains all the needed constants in the galaxy creation.
//! Information retrieved from <https://ssd.jpl.nasa.gov/?planet_phys_par>.

/// The mass of the Sun (M<sub>☉</sub>), in *kg*.
pub const SUN_MASS: f64 = 1.9885e+30;
/// The radius of the Sun (R<sub>☉</sub>), in meters (*m*).
pub const SUN_RADIUS: f64 = 6.957e+8;
/// The luminosity of the Sun (L<sub>☉</sub>), in watts (*W*).
pub const SUN_LUMINOSITY: f64 = 3.828e+26;

/// The mass of the Earth (M<sub>⊕</sub>), in *kg*.
pub const EARTH_MASS: f64 = 5.97237e+24;
/// The radius of the Earth (R<sub>⊕</sub>), in meters (*m*).
pub const EARTH_RADIUS: f64 = 6_371_008.4;
/// The atmospheric pressure at sea level in Earth, in pascals (*Pa*).
pub const EARTH_ATM_PRESSURE: f64 = 101_325_f64;
/// The surface gravity in Earth, in *m/s²*.
pub const EARTH_GRAVITY: f64 = 9.80665;

/// The mass of Jupiter (M<sub>j</sub>), in *kg*.
pub const JUPITER_MASS: f64 = 1.898_187e+27;
/// The radius of Jupiter (R<sub>j</sub>), in meters (*m*).
pub const JUPITER_RADIUS: f64 = 69_911_000_f64;

/// *G*, the gravitational constant, in *m³·kg⁻¹·s⁻²*.
pub const G: f64 = 6.6740831e-11;
/// *c*, the speed of light in vacuum, in *m·s⁻¹*.
pub const C: f64 = 299_792_458_f64;
/// *σ*, the Stefan–Boltzmann constant, in *W·m⁻²·K⁻⁴*.
pub const BOLTZ: f64 = 5.67036713e-8;

/// The Chandrasekhar limit, in *kg*.
pub const CH_LIMIT: f64 = 2.765e+30;

/// The astronomical unit (*au* or *AU*), the average distance from earth to the sun, in meters
/// (*m*).
pub const AU: f64 = 149_597_870_700_f64;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::EPSILON;

    #[test]
    fn it_constants() {
        assert!(SUN_MASS <= 1.9885e+30 + EPSILON && SUN_MASS >= 1.9885e+30 - EPSILON);
        assert!(SUN_RADIUS <= 6.957e+8 + EPSILON && SUN_RADIUS >= 6.957e+8 - EPSILON);
        assert!(SUN_LUMINOSITY <= 3.828e+26 + EPSILON && SUN_LUMINOSITY >= 3.828e+26 - EPSILON);

        assert!(EARTH_MASS <= 5.97237e+24 + EPSILON && EARTH_MASS >= 5.97237e+24 - EPSILON);
        assert!(
            EARTH_RADIUS <= 6_371_008.4_f64 + EPSILON && EARTH_RADIUS >= 6_371_008.4_f64 - EPSILON
        );

        assert!(JUPITER_MASS <= 1.898_187e+27 + EPSILON && JUPITER_MASS >= 1.898_187e+27 - EPSILON);
        assert!(
            JUPITER_RADIUS <= 69_911_000_f64 + EPSILON
                && JUPITER_RADIUS >= 69_911_000_f64 - EPSILON
        );

        assert!(G <= 6.6740831e-11 + EPSILON && G >= 6.6740831e-11 - EPSILON);
        assert!(C <= 299_792_458_f64 + EPSILON && C >= 299_792_458_f64 - EPSILON);
        assert!(BOLTZ <= 5.67036713e-8 + EPSILON && BOLTZ >= 5.67036713e-8 - EPSILON);

        assert!(CH_LIMIT <= 2.765e+30 + EPSILON && CH_LIMIT >= 2.765e+30 - EPSILON);

        assert!(AU <= 149_597_870_700_f64 + EPSILON && AU >= 149_597_870_700_f64 - EPSILON);
    }
}
