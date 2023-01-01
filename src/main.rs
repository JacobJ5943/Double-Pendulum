use std::{
    collections::{HashMap, HashSet},
    f64::RADIX,
    sync::Arc,
};

use eframe::egui;
use egui::{Color32, Pos2, Stroke, Ui};

fn main() {
    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "My egui App",
        native_options,
        Box::new(|cc| Box::new(MyEguiApp::new(cc))),
    );
}

#[derive(Default, Debug)]
struct MyEguiApp {
    inner: Pos2,
    outer: Pos2,
    point: Pos2,
    is_draggin: bool,
    started: bool,
}

impl MyEguiApp {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            point: Pos2 { x: 1.0, y: 2.0 },
            is_draggin: false,
            inner: Pos2 { x: 0.0, y: 0.0 },
            outer: Pos2 { x: 0.0, y: 0.0 },
            started: false,
        }
    }
}

impl MyEguiApp {
    fn reset_button(&mut self, ui: &mut Ui, center: Pos2, total_radius: f32) {
        if ui.button("reset").clicked() {
            self.inner = Pos2 {
                x: center.x - 120.0,
                y: center.y,
            };
            self.outer = Pos2 {
                x: center.x - total_radius,
                y: center.y,
            };
        }
    }

    fn handle_primary_mouse_button_down(&mut self, ui: &mut Ui, center: Pos2, total_radius: f32) {
        if ui.input().pointer.primary_down() {
            if let Some(thing) = { ui.input().pointer.interact_pos() } {
                if self.is_draggin {
                    self.outer = thing;

                    // keep it in the greater radius
                    let distance =
                        ((center.x - thing.x).powi(2) + (center.y - thing.y).powi(2)).sqrt();
                    if distance >= total_radius {
                        let scaling = (total_radius / distance);
                        //dbg!(scaling);
                        let delta_x = center.x - thing.x;
                        let delta_y = center.y - thing.y;

                        self.outer = Pos2 {
                            x: center.x - delta_x * scaling,
                            y: center.y - delta_y * scaling,
                        };
                    }
                } else if intersect_circle(self.outer, 15.0, thing) {
                    self.outer = thing;
                    self.is_draggin = true;
                }
            }
        } else {
            self.is_draggin = false;
        }

        let distance =
            ((center.x - self.outer.x).powi(2) + (center.y - self.outer.y).powi(2)).sqrt();
        if let Some((p1, p2)) = dbg!(calculate_intersecting_points(
            center, self.outer, distance, 120.0, 120.0
        )) {
            let delta_p1 = ((p1.x - self.inner.x).powi(2) + (p1.y - self.inner.y).powi(2)).sqrt();
            let delta_p2 = ((p2.x - self.inner.x).powi(2) + (p2.y - self.inner.y).powi(2)).sqrt();
            let nan = Pos2 {
                x: f32::NAN,
                y: f32::NAN,
            };

            if p1.x.is_nan() && p2.x.is_nan() {
                // dbg!((p1,p2));
            } else {
                if delta_p1 > delta_p2 {
                    self.inner = p2;
                } else {
                    self.inner = p1;
                }
            }
        }
    }

    fn draw_pendulum(&mut self, ui: &mut Ui, center: Pos2) {
        ui.painter().circle(
            center,
            15.0,
            Color32::GREEN,
            Stroke::new(1.0, Color32::GOLD),
        );
        ui.painter()
            .circle(self.inner, 15.0, Color32::RED, Stroke::default());
        ui.painter()
            .circle(self.outer, 15.0, Color32::BLUE, Stroke::default());
        ui.painter()
            .line_segment([self.inner, center], Stroke::new(3.0, Color32::LIGHT_RED));
        ui.painter().line_segment(
            [self.outer, self.inner],
            Stroke::new(3.0, Color32::LIGHT_RED),
        );
    }
}
impl eframe::App for MyEguiApp {
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.label(format!("{},{}", self.point.x, self.point.y));
            let available_size = ui.available_size();
            let total_radius: f32 = 240.0;
            let center = Pos2 {
                x: available_size.x / 2.0,
                y: available_size.y / 2.0,
            };
            if !self.started {
                self.inner = Pos2 {
                    x: center.x - 120.0,
                    y: center.y,
                };
                self.outer = Pos2 {
                    x: center.x - total_radius,
                    y: center.y,
                };

                self.started = true;
            }
            // Find the center of the screen

            self.reset_button(ui, center, total_radius);

            self.handle_primary_mouse_button_down(ui, center, total_radius);

            self.draw_pendulum(ui, center);
        });
    }
}

/// this uses the sqrt(deltax ^2 + deltay^2) so if the delta is great enough an overlow may occur
fn intersect_circle(object_center: Pos2, object_radius: f32, clicked_pos: Pos2) -> bool {
    // This seems scuffed with both -
    object_radius
        >= ((object_center.x - clicked_pos.x).powi(2) + (object_center.y - clicked_pos.y).powi(2))
            .sqrt()
}

fn calculate_intersecting_points(
    B: Pos2,
    A: Pos2,
    d: f32,
    a_radius: f32,
    b_radius: f32,
) -> Option<(Pos2, Pos2)> {
    let ex = (B.x - A.x) / d;
    let ey = (B.y - A.y) / d;

    let x = (a_radius * a_radius - b_radius * b_radius + d * d) / (2.0 * d);
    let y = (a_radius * a_radius - x * x).sqrt();

    let p1 = Pos2 {
        x: A.x + x * ex - y * ey,
        y: A.y + x * ey + y * ex,
    };

    let p2 = Pos2 {
        x: A.x + x * ex + y * ey,
        y: A.y + x * ey - y * ex,
    };
    Some((p1, p2))
}
