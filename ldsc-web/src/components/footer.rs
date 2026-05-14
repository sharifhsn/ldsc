//! Sticky footer — LDLink's blue `#2a71a5` agency footer, but with
//! project links instead of HHS / FOIA / Vulnerability Disclosure.

use leptos::prelude::*;

#[component]
pub fn Footer() -> impl IntoView {
    view! {
        <footer class="app-footer">
            <div class="footer-grid">
                <div>
                    <h6>"Project"</h6>
                    <div>"ldsc-rs · v0.5.0"</div>
                    <div style="opacity: 0.85;">
                        "Rust port of Bulik-Sullivan et al."
                    </div>
                </div>
                <div>
                    <h6>"Code"</h6>
                    <div>
                        <a href="https://github.com/sharifhsn/ldsc"
                           target="_blank"
                           rel="noopener noreferrer">
                            "github.com/sharifhsn/ldsc"
                        </a>
                    </div>
                    <div>
                        <a href="https://github.com/sharifhsn/ldsc/issues"
                           target="_blank"
                           rel="noopener noreferrer">"Report an issue"</a>
                    </div>
                </div>
                <div>
                    <h6>"Reference"</h6>
                    <div>
                        <a href="https://www.nature.com/articles/ng.3211"
                           target="_blank"
                           rel="noopener noreferrer">"LDSC paper (2015)"</a>
                    </div>
                    <div>
                        <a href="https://github.com/bulik/ldsc"
                           target="_blank"
                           rel="noopener noreferrer">"Original Python LDSC"</a>
                    </div>
                </div>
                <div>
                    <h6>"License"</h6>
                    <div>
                        <a href="https://www.gnu.org/licenses/gpl-3.0.html"
                           target="_blank"
                           rel="noopener noreferrer">"GPL-3.0"</a>
                    </div>
                </div>
            </div>
        </footer>
    }
}
