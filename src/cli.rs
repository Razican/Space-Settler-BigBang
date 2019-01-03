use clap::{App, Arg};

/// Generates the command line interface.
pub fn generate() -> App<'static, 'static> {
    App::new("Space Settler Big Bang")
        .version(crate_version!())
        .author("Razican <razican@protonmail.ch>")
        .about("Creates a new galaxy for Space Settler")
        .arg(
            Arg::with_name("earth")
            .long("earth").short("e")
                .help("Just search for an earth-like planet.")
                .required_unless("galaxy")
                .conflicts_with("galaxy")
                .conflicts_with("stars")
                .takes_value(false)
        ).arg(
            Arg::with_name("galaxy")
            .long("galaxy").short("g")
                .help("Create a complete galaxy.")
                .required_unless("earth")
                .conflicts_with("earth")
                .takes_value(false)
        ).arg(
            Arg::with_name("stars")
            .long("stars").short("s")
                .help("Approximate number of stars in the galaxy (it will generate between 90% and 110% of the given number).")
                .conflicts_with("earth")
                .takes_value(true).default_value("100000")
        ).arg(
            Arg::with_name("threads")
            .long("threads").short("t")
                .help("Number of threads to use to generate the galaxy, by default, the number of virtual CPUs.")
                .conflicts_with("earth")
                .takes_value(true)
        ).arg(
            Arg::with_name("verbose")
            .long("verbose").short("v")
                .help("Adds extra information to the output.")
                .takes_value(false)
        )
}
