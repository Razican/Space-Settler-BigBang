//! Program entry point.
#![feature(proc_macro_hygiene)]
// Lints
#![forbid(anonymous_parameters)]
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
    unused_extern_crates,
    missing_debug_implementations,
    missing_copy_implementations
)]
#![warn(clippy::pedantic)]
#![allow(clippy::non_ascii_literal)]

use dotenv::dotenv;
use rocket::routes;
use sps_web_core::*;

/// Program entry point.
fn main() {
    let _ = dotenv().ok();

    let _ = rocket::ignite()
        .attach(SpsDb::fairing())
        .mount(
            "/",
            routes![
                index,
                favicon,
                android_config,
                windows_config,
                images,
                css,
                js
            ],
        )
        .launch();
}
