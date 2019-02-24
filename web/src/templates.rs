//! Templates for the web.

mod footer;
mod header;
mod home;

pub use self::footer::template as footer;
pub use self::header::template as header;
pub use self::home::template as home;

/// Header context structure.
#[derive(Debug)]
pub struct HeaderContext {
    /// Page title.
    title: String,
}

impl HeaderContext {
    /// Creates a new header context.
    pub fn new<T: Into<String>>(title: T) -> Self {
        Self {
            title: title.into(),
        }
    }
}
