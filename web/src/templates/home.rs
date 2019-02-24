//! Panel template.

use maud::{html, Markup};

/// Panel template.
pub fn template(header: &Markup, footer: &Markup) -> Markup {
    html! {
        // Add the header markup to the page
        (header)
        body {
            header {
                h1 {
                    "Hello World"
                };
            };

            // Add the footer markup to the page
            (footer)
        }
    }
}
