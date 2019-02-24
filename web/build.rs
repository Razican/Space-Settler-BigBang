//! GitWorld build script.

use std::{env, fs, path::Path, process::Command};

use failure::{Error, ResultExt};

fn main() {
    // Clean the `dist` folder.
    clean_dist().expect("could not clean the `dist` folder.");

    // Install NPM dependencies:
    npm_install().expect("there was an error installing NPM dependencies. Is npm installed?");

    // Generating Webpack resources:
    webpack_generate().expect("there was an error running Webpack");
}

/// Cleans the `dist` folder.
fn clean_dist() -> Result<(), Error> {
    let path_css = Path::new("dist/css");
    if path_css.exists() {
        fs::remove_dir_all(path_css)?;
    }
    let path_js = Path::new("dist/js");
    if path_js.exists() {
        fs::remove_dir_all(path_js)?;
    }
    Ok(())
}

/// Installs NPM packages.
fn npm_install() -> Result<(), Error> {
    let output = Command::new("npm").arg("install").output()?;
    if !output.status.success() {
        let stdout = String::from_utf8(output.stdout).context("error reading `stdout`")?;
        panic!("Error installing NPM dependencies: {}", stdout);
    }
    Ok(())
}

/// Generates Webpack resources.
fn webpack_generate() -> Result<(), Error> {
    let mut command = Command::new("npx");
    command.arg("webpack");
    match env::var("PROFILE") {
        Ok(ref var) if var == "release" => {
            command.env("NODE_ENV", "production");
        }
        _ => {
            command.env("NODE_ENV", "development");
        }
    }
    let output = command.output()?;
    if !output.status.success() {
        let stdout = String::from_utf8(output.stdout).context("error reading `stdout`")?;
        panic!("Error generating Webpack resources: {}", stdout);
    }
    Ok(())
}
