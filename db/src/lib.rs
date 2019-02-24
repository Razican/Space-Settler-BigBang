#![recursion_limit = "1024"]
#![allow(proc_macro_derive_resolution_fallback)]

#[macro_use]
extern crate diesel;

pub mod planets;
pub mod schema;
pub mod stars;
pub mod users;

pub use crate::{planets::PlanetDao, stars::StarDao, users::UserDao};
