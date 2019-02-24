//! Header template.

use maud::{html, Markup, DOCTYPE};

use super::HeaderContext;

/// Header template.
pub fn template(context: &HeaderContext) -> Markup {
    html! {
        (DOCTYPE)
        head {
            meta charset="utf-8";
            title { (context.title) }
            link rel="stylesheet" type="text/css" href="/css/main.css";
        }
    }
}
