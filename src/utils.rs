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

/// Checks if water can be liquid
///
/// # Examples
///
/// ```
/// let min_temp = 260_f64;
/// let max_temp = 290_f64;
/// let pressure = 101_325_f64;
///
/// let can_be_liquid = can_water_be_liquid(min_temp, max_temp, pressure);
///
/// // Now can_be_liquid should be true
/// # assert_eq!(true, can_be_liquid);
/// ```
pub fn can_water_be_liquid(min_temp: f64, max_temp: f64, pressure: f64) -> bool {
    pressure > 611.657  &&
    pressure > get_water_vaporize_pressure(min_temp) && // not to vaporize in low temperatures
    pressure < get_water_melt_pressure(max_temp) // Not to ice in high temperatures
}

/// Checks if water can be ice
///
/// # Examples
///
/// ```
/// let min_temp = 260_f64;
/// let max_temp = 290_f64;
/// let pressure = 101_325_f64;
///
/// let can_be_ice = can_water_be_ice(min_temp, pressure);
///
/// // Now can_be_ice should be true
/// # assert_eq!(true, can_be_ice);
/// ```
pub fn can_water_be_ice(min_temp: f64, pressure: f64) -> bool {
    if pressure < 611.657 {
        pressure > get_water_sublimation_pressure(min_temp)
    } else {
        min_temp < get_water_melt_pressure(min_temp)
    }
}

fn get_water_vaporize_pressure(temp: f64) -> f64 {
    -2836.5744*temp.powi(-2) -6028.076559/temp +19.54263612 -0.02737830188*temp
    +1.6261698e-5*temp.powi(2)+7.022905e-10*temp.powi(3)
    -1.8680009e-13*temp.powi(4)+2.7150305*temp.ln()
}

fn get_water_sublimation_pressure(temp: f64) -> f64 {
    (-5723.265/temp+9.550426-0.00728332*temp+3.53068*temp.ln()).exp()
}

fn get_water_melt_pressure(temp: f64) -> f64 {
    -395.2*((temp/273.16).powi(9)-1_f64)*1e+6
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
