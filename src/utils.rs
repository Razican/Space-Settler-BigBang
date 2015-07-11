//! Utilities for Space Settler Galaxy Creator
//!
//! This module contains all the common utilities needed during the creation process of the galaxy.

use std::f64::consts::PI;

/// Converts radians to degrees
///
/// # Examples
///
/// ```
/// use std::f64::consts::PI;
///
/// let rad = 2_f64*PI;
/// let deg = rad_to_deg(rad);
///
/// // Now deg should be 360
/// # assert_eq!(360_f64, deg);
/// ```
pub fn rad_to_deg(rads: f64) -> f64 {
    rads*180_f64/PI
}

/// Converts degrees to radians
///
/// # Examples
///
/// ```
/// let deg = 360_f64;
/// let rad = deg_to_rad(deg);
///
/// // Now rad should be 2*PI
/// # use std::f64::consts::PI;
/// # assert_eq!(2_f64*PI, rad);
/// ```
pub fn deg_to_rad(degs: f64) -> f64 {
    degs*PI/180_f64
}

/// Converts Celsius degrees (*°C*)to Kelvin (*K*)
///
/// # Examples
///
/// ```
/// let celsius = 200_f64;
/// let kelvin = celsius_to_kelvin(celsius);
///
/// // Now kelvin should be 473.15 K
/// # assert_eq!(473.15_f64, kelvin);
/// ```
pub fn celsius_to_kelvin(celsius: f64) -> f64 {
    celsius+273.15_f64
}

/// Converts Kelvins (*K*) to Celsius (*°C*)
///
/// # Examples
///
/// ```
/// let kelvin = 273.15_f64;
/// let celsius = kelvin_to_celsius(kelvin);
///
/// // Now kelvin should be 0°C
/// # assert_eq!(0_f64, celsius);
/// ```
pub fn kelvin_to_celsius(kelvin: f64) -> f64 {
    kelvin-273.15_f64
}

#[cfg(test)]
mod tests {
	use std::f64::consts::PI;
	use super::*;

    #[test]
    fn it_rad_to_deg() {
        assert_eq!(360_f64, rad_to_deg(2_f64*PI));
    }

    #[test]
    fn it_deg_to_rad() {
        assert_eq!(2_f64*PI, deg_to_rad(360_f64));
    }

    #[test]
    fn it_celsius_to_kelvin() {
        assert_eq!(473.15_f64, celsius_to_kelvin(200_f64));
    }

    #[test]
    fn it_kelvin_to_celsius() {
        assert_eq!(0_f64, kelvin_to_celsius(273.15_f64));
    }
}
