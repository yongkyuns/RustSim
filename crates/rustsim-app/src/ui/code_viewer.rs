//! Code viewer UI component with syntax highlighting.

use egui::{text::LayoutJob, Color32, FontId, ScrollArea, TextFormat, Ui};
use rustsim_codegen::EmbeddedGenerator;
use syntect::easy::HighlightLines;
use syntect::highlighting::{Style, ThemeSet};
use syntect::parsing::SyntaxSet;
use syntect::util::LinesWithEndings;

use crate::state::AppState;

/// Syntax highlighter state (cached for performance)
struct SyntaxHighlighter {
    syntax_set: SyntaxSet,
    theme_set: ThemeSet,
}

impl SyntaxHighlighter {
    fn new() -> Self {
        Self {
            syntax_set: SyntaxSet::load_defaults_newlines(),
            theme_set: ThemeSet::load_defaults(),
        }
    }

    fn highlight(&self, code: &str, dark_mode: bool) -> LayoutJob {
        let theme_name = if dark_mode {
            "base16-ocean.dark"
        } else {
            "base16-ocean.light"
        };

        let theme = &self.theme_set.themes[theme_name];
        let syntax = self
            .syntax_set
            .find_syntax_by_extension("rs")
            .unwrap_or_else(|| self.syntax_set.find_syntax_plain_text());

        let mut h = HighlightLines::new(syntax, theme);
        let mut job = LayoutJob::default();

        for line in LinesWithEndings::from(code) {
            match h.highlight_line(line, &self.syntax_set) {
                Ok(ranges) => {
                    for (style, text) in ranges {
                        let color = syntect_color_to_egui(style);
                        job.append(
                            text,
                            0.0,
                            TextFormat {
                                font_id: FontId::monospace(12.0),
                                color,
                                ..Default::default()
                            },
                        );
                    }
                }
                Err(_) => {
                    // Fallback: plain text
                    job.append(
                        line,
                        0.0,
                        TextFormat {
                            font_id: FontId::monospace(12.0),
                            color: Color32::GRAY,
                            ..Default::default()
                        },
                    );
                }
            }
        }

        job
    }
}

/// Convert syntect color to egui color
fn syntect_color_to_egui(style: Style) -> Color32 {
    Color32::from_rgb(style.foreground.r, style.foreground.g, style.foreground.b)
}

thread_local! {
    static HIGHLIGHTER: SyntaxHighlighter = SyntaxHighlighter::new();
}

/// Render the code viewer panel
pub fn render_code_viewer(ui: &mut Ui, state: &mut AppState) {
    ui.horizontal(|ui| {
        ui.heading("Generated Code");
        ui.separator();

        if ui.button("Copy to Clipboard").clicked() {
            let code = generate_code(state);
            ui.output_mut(|o| o.copied_text = code);
        }
    });

    ui.separator();

    if state.graph().nodes.is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label("Add blocks to see generated code");
        });
        return;
    }

    // Generate the code
    let code = generate_code(state);

    // Get dark mode from UI visuals
    let dark_mode = ui.visuals().dark_mode;

    // Create highlighted layout job
    let layout_job = HIGHLIGHTER.with(|h| h.highlight(&code, dark_mode));

    // Display with scroll area
    ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            // Use a label with the layout job for syntax highlighted text
            ui.add(egui::Label::new(layout_job).selectable(true));
        });
}

/// Generate code from the current graph
fn generate_code(state: &AppState) -> String {
    let generator = EmbeddedGenerator::new("Simulation");
    generator.generate(state.graph(), state.settings())
}
