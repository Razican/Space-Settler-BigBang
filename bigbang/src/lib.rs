//! Space Settler big bang creation core.

#![forbid(anonymous_parameters)]
#![warn(clippy::pedantic)]
#![deny(
    clippy::all,
    variant_size_differences,
    unused_results,
    unused_qualifications,
    unused_import_braces,
    unsafe_code,
    trivial_numeric_casts,
    trivial_casts,
    missing_docs,
    missing_debug_implementations,
    missing_copy_implementations,
    unused_extern_crates
)]
#![allow(clippy::unreadable_literal)]

pub mod consts;
pub mod planet;
pub mod star;
pub mod utils;
