//! Footer template.

use maud::{html, Markup};

/// Footer template.
pub fn template() -> Markup {
    html! {
        footer {
            p {
                "Copyright © GitWorld 2018"
            }
        }
        script src="/js/bundle.js" defer? {};
    }
}
