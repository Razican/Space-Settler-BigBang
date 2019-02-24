//! Module containing structures and methods to manage stars.

use crate::schema::stars::{self, all_columns, columns, dsl, table};
use diesel::{
    associations::HasTable,
    deserialize::{self, FromSql},
    dsl::Select,
    insert_into,
    pg::Pg,
    prelude::*,
    query_builder::{DeleteStatement, InsertStatement, IntoUpdateTarget},
    serialize::{self, IsNull, Output, ToSql},
};
use std::io::Write;

// TODO: unsigned
/// Star information structure.
#[derive(Identifiable, Queryable, PartialEq, Debug, Copy, Clone)]
#[table_name = "stars"]
pub struct Star {
    id: i32,
    orb_radius: i32,
    class: StarClass,
    mass: f64,
    radius: f64,
    temperature: Option<i32>,
}

impl Star {
    /// Gets the ID of the star in the database.
    pub fn id(&self) -> u32 {
        self.id as u32
    }

    /// Gets the radius of the orbit of the star in the galaxy, in 1/100 light years.
    pub fn orb_radius(&self) -> u32 {
        self.orb_radius as u32
    }

    /// Gets the class of the star.
    pub fn class(&self) -> StarClass {
        self.class
    }

    /// Gets the mass of the star, in *kg*.
    pub fn mass(&self) -> f64 {
        self.mass
    }

    /// Gets the radius of the star, in *m*.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Gets the temperature of the star, in 1/1,000 *K*.
    pub fn temperature(&self) -> Option<u32> {
        match self.temperature {
            Some(t) => Some(t as u32),
            None => None,
        }
    }
}

/// New star creation.
#[derive(Debug, Clone, Copy, Insertable)]
#[table_name = "stars"]
pub struct NewStar {
    pub orb_radius: i32,
    pub class: StarClass,
    pub mass: f64,
    pub radius: f64,
    pub temperature: Option<i32>,
}

/// Data access object for stars.
pub struct StarDao {}

type AllColumns = (
    columns::id,
    columns::orb_radius,
    columns::class,
    columns::mass,
    columns::radius,
    columns::temperature,
);

impl StarDao {
    /// Selects all columns from the stars.
    pub fn all() -> Select<dsl::stars, AllColumns> {
        use crate::schema::stars::dsl::*;
        stars.select(all_columns)
    }

    /// Inserts a star or multiple stars in the database.
    pub fn insert<S>(st: S) -> InsertStatement<table, S::Values>
    where
        S: Insertable<dsl::stars>,
    {
        use crate::schema::stars::dsl::*;
        insert_into(stars).values(st)
    }

    /// Deletes all stars from the database.
    pub fn delete_all() -> DeleteStatement<
        <crate::schema::stars::dsl::stars as HasTable>::Table,
        <crate::schema::stars::dsl::stars as IntoUpdateTarget>::WhereClause,
    > {
        use crate::schema::stars::dsl::*;
        diesel::delete(stars)
    }
}

/// Mapping of the SQL `star_class` type to the `StarClass` enumeration.
#[derive(SqlType)]
#[postgres(type_name = "star_class")]
pub struct StarClassMapping;

/// Enumeration representing the class of a star.
#[derive(Debug, Clone, Copy, PartialEq, FromSqlRow, AsExpression)]
#[sql_type = "StarClassMapping"]
pub enum StarClass {
    O,
    B,
    A,
    F,
    G,
    K,
    M,
    BlackHole,
    NeutronStar,
    QuarkStar,
    WhiteDwarf,
}

impl FromSql<StarClassMapping, Pg> for StarClass {
    fn from_sql(bytes: Option<&[u8]>) -> deserialize::Result<Self> {
        match not_none!(bytes) {
            b"O" => Ok(StarClass::O),
            b"B" => Ok(StarClass::B),
            b"A" => Ok(StarClass::A),
            b"F" => Ok(StarClass::F),
            b"G" => Ok(StarClass::G),
            b"K" => Ok(StarClass::K),
            b"M" => Ok(StarClass::M),
            b"black_hole" => Ok(StarClass::BlackHole),
            b"neutron_star" => Ok(StarClass::NeutronStar),
            b"quark_star" => Ok(StarClass::QuarkStar),
            b"white_dwarf" => Ok(StarClass::WhiteDwarf),
            _ => Err("unrecognized enum variant".into()),
        }
    }
}

impl ToSql<StarClassMapping, Pg> for StarClass {
    fn to_sql<W>(&self, out: &mut Output<W, Pg>) -> serialize::Result
    where
        W: Write,
    {
        match *self {
            StarClass::O => out.write_all(b"O")?,
            StarClass::B => out.write_all(b"B")?,
            StarClass::A => out.write_all(b"A")?,
            StarClass::F => out.write_all(b"F")?,
            StarClass::G => out.write_all(b"G")?,
            StarClass::K => out.write_all(b"K")?,
            StarClass::M => out.write_all(b"M")?,
            StarClass::BlackHole => out.write_all(b"black_hole")?,
            StarClass::NeutronStar => out.write_all(b"neutron_star")?,
            StarClass::QuarkStar => out.write_all(b"quark_star")?,
            StarClass::WhiteDwarf => out.write_all(b"white_dwarf")?,
        }
        Ok(IsNull::No)
    }
}
