//! Implementation of the web page.
#![feature(proc_macro_hygiene, decl_macro)]
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

#[macro_use]
extern crate rocket;
#[macro_use]
extern crate rocket_contrib;

mod compress;
mod templates;

use std::path::{Path, PathBuf};

use rocket::{http::ContentType, response::NamedFile};
use rocket_contrib::databases::diesel;

use crate::{
    compress::{CompressedFile, CompressedTemplate},
    templates::{footer, header, home, HeaderContext},
};

/// Space Settler database
#[database("sps_db")]
#[derive(Debug)]
pub struct SpsDb(diesel::PgConnection);

/// Index page.
#[get("/")]
pub fn index() -> CompressedTemplate {
    let header_context = HeaderContext::new("Space Settler - Home");

    let header_markup = header(&header_context);
    let footer_markup = footer();

    CompressedTemplate::new(home(&header_markup, &footer_markup))
}

/// Main favicon.
#[get("/favicon.ico")]
pub fn favicon() -> Option<NamedFile> {
    NamedFile::open("dist/fav/favicon.ico").ok()
}

/// Android app configuration.
#[get("/fav/browserconfig.xml")]
pub fn android_config() -> Option<CompressedFile> {
    CompressedFile::new("dist/fav/browserconfig.xml", ContentType::XML).ok()
}

/// Windows app configuration.
#[get("/fav/manifest.json")]
pub fn windows_config() -> Option<CompressedFile> {
    CompressedFile::new("dist/fav/manifest.json", ContentType::JSON).ok()
}

/// Favicons files.
#[get("/fav/<file..>", rank = 2)]
pub fn favicons(file: PathBuf) -> Option<NamedFile> {
    NamedFile::open(Path::new("dist/fav/").join(file)).ok()
}

/// Image files.
#[get("/img/<file..>", rank = 2)]
pub fn images(file: PathBuf) -> Option<NamedFile> {
    NamedFile::open(Path::new("dist/img/").join(file)).ok()
}

// TODO: do not build in release mode.
/// CSS file.
#[get("/css/<file..>")]
pub fn css(file: PathBuf) -> Option<CompressedFile> {
    CompressedFile::new(Path::new("dist/css").join(file), ContentType::CSS).ok()
}

/// Gets a javascript file.
#[get("/js/<file..>")]
pub fn js(file: PathBuf) -> Option<CompressedFile> {
    CompressedFile::new(Path::new("dist/js").join(file), ContentType::JavaScript).ok()
}
