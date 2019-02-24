//! Module containing structures and methods to manage planets.

use crate::schema::planets::{self, all_columns, columns, dsl, table};
use diesel::{
    deserialize::{self, FromSql},
    dsl::Select,
    insert_into,
    pg::Pg,
    prelude::*,
    query_builder::InsertStatement,
    serialize::{self, IsNull, Output, ToSql},
};
use std::io::Write;

/// Planet information structure
#[derive(Identifiable, Queryable, Associations, PartialEq, Debug)]
#[belongs_to(crate::stars::Star)]
#[table_name = "planets"]
pub struct Planet {
    id: i32,
    position: i16,
    orb_ecc: f64,
    orb_sma: f64,
    orb_incl: f64,
    orb_lan: f64,
    orb_arg_p: f64,
    orb_m0: f64,
    orb_period: f64,
    atm_pressure: Option<f64>,
    atm_h2o: Option<f64>,
    atm_co2: Option<f64>,
    atm_co: Option<f64>,
    atm_n2: Option<f64>,
    atm_o2: Option<f64>,
    atm_ar: Option<f64>,
    atm_so2: Option<f64>,
    atm_ne: Option<f64>,
    atm_ch4: Option<f64>,
    atm_he: Option<f64>,
    ax_tilt: f64,
    rot_period: f64,
    planet_type: PlanetType,
    surf_fresh_water: Option<f64>,
    surf_ocean_water: Option<f64>,
    surf_snow: Option<f64>,
    surf_land: Option<f64>,
    // crust: Crust,
    // life: Life,
    bond_albedo: f64,
    geometric_albedo: f64,
    mass: f64,
    radius: f64,
    eff_temp: f64,
    min_temp: f64,
    max_temp: f64,
    avg_temp: f64,
    // habitable: bool
}

/// New planet creation.
#[derive(Debug, Clone, Copy, Insertable)]
#[table_name = "planets"]
pub struct NewPlanet {
    pub position: i16,
    pub orb_ecc: f64,
    pub orb_sma: f64,
    pub orb_incl: f64,
    pub orb_lan: f64,
    pub orb_arg_p: f64,
    pub orb_m0: f64,
    pub orb_period: f64,
    pub atm_pressure: Option<f64>,
    pub atm_h2o: Option<f64>,
    pub atm_co2: Option<f64>,
    pub atm_co: Option<f64>,
    pub atm_n2: Option<f64>,
    pub atm_o2: Option<f64>,
    pub atm_ar: Option<f64>,
    pub atm_so2: Option<f64>,
    pub atm_ne: Option<f64>,
    pub atm_ch4: Option<f64>,
    pub atm_he: Option<f64>,
    pub ax_tilt: f64,
    pub rot_period: f64,
    pub planet_type: PlanetType,
    pub surf_fresh_water: Option<f64>,
    pub surf_ocean_water: Option<f64>,
    pub surf_snow: Option<f64>,
    pub surf_land: Option<f64>,
    // crust: Crust,
    // life: Life,
    pub bond_albedo: f64,
    pub geometric_albedo: f64,
    pub mass: f64,
    pub radius: f64,
    pub eff_temp: f64,
    pub min_temp: f64,
    pub max_temp: f64,
    pub avg_temp: f64,
    // habitable: bool
}

/// Data access object for planets.
pub struct PlanetDao {}

type AllColumns = (
    columns::id,
    columns::star_id,
    columns::position,
    columns::orb_ecc,
    columns::orb_sma,
    columns::orb_incl,
    columns::orb_lan,
    columns::orb_arg_p,
    columns::orb_m0,
    columns::orb_period,
    columns::atm_pressure,
    columns::atm_h2o,
    columns::atm_co2,
    columns::atm_co,
    columns::atm_n2,
    columns::atm_o2,
    columns::atm_ar,
    columns::atm_so2,
    columns::atm_ne,
    columns::atm_ch4,
    columns::atm_he,
    columns::ax_tilt,
    columns::rot_period,
    columns::planet_type,
    columns::surf_fresh_water,
    columns::surf_ocean_water,
    columns::surf_snow,
    columns::surf_land,
    columns::bond_albedo,
    columns::geometric_albedo,
    columns::mass,
    columns::radius,
    columns::eff_temp,
    columns::min_temp,
    columns::max_temp,
    columns::avg_temp,
    columns::habitable,
);

impl PlanetDao {
    /// Selects all columns from the planets.
    pub fn all() -> Select<dsl::planets, AllColumns> {
        use crate::schema::planets::dsl::*;
        planets.select(all_columns)
    }

    /// Inserts a planet or multiple planets in the database.
    pub fn insert<P>(pl: P) -> InsertStatement<table, P::Values>
    where
        P: Insertable<dsl::planets>,
    {
        use crate::schema::planets::dsl::*;
        insert_into(planets).values(pl)
    }
}

/// Mapping of the SQL `planet_type` type to the `PlanetType` enumeration.
#[derive(SqlType)]
#[postgres(type_name = "planet_type")]
pub struct PlanetTypeMapping;

/// Enumeration representing the type of planet.
#[derive(Debug, Clone, Copy, PartialEq, FromSqlRow, AsExpression)]
#[sql_type = "PlanetTypeMapping"]
pub enum PlanetType {
    Rocky,
    Gaseous,
}

impl FromSql<PlanetTypeMapping, Pg> for PlanetType {
    fn from_sql(bytes: Option<&[u8]>) -> deserialize::Result<Self> {
        match not_none!(bytes) {
            b"rocky" => Ok(PlanetType::Rocky),
            b"gaseous" => Ok(PlanetType::Gaseous),
            _ => Err("unrecognized enum variant".into()),
        }
    }
}

impl ToSql<PlanetTypeMapping, Pg> for PlanetType {
    fn to_sql<W>(&self, out: &mut Output<W, Pg>) -> serialize::Result
    where
        W: Write,
    {
        match *self {
            PlanetType::Rocky => out.write_all(b"rocky")?,
            PlanetType::Gaseous => out.write_all(b"gaseous")?,
        }
        Ok(IsNull::No)
    }
}
