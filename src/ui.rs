use std::f32::consts::PI;
use std::time::Duration;

use crate::{
    constants::{MASS_GRAVITY, MASS_MAX, MAX_FRAMERATE, MAX_GRANULARITY, MAX_STEPS_PER_TICK},
    math::{
        calculate_intersecting_points, calculate_theta_one_double_dot,
        calculate_theta_two_double_dot, cartesian_to_polar, intersect_circle, polar_to_cartesian,
        Polar2,
    },
};
use eframe::egui;
use egui::{Color32, Pos2, Stroke, Ui};

#[derive(Default, Debug)]
pub struct MyEguiApp {
    inner: Polar2,
    //inner_rod_length: f32,
    //outer_rod_length: f32,
    theta_one_dot: f32,
    theta_two_dot: f32,
    outer: Polar2,
    is_dragging: bool,
    is_dragging_p1: bool,
    is_dragging_p2: bool,
    inner_mass: f32,
    outer_mass: f32,
    counter: usize,
    steps_per_tick: usize,
    granularity: f32,
    gravity: f32,
    frame_rate: f32,
    playing: bool,
}

impl eframe::App for MyEguiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Inner rod length");
                ui.add(egui::Slider::new(&mut self.inner.radius, 0.1..=500.0));
            });
            ui.horizontal(|ui| {
                ui.label("Outer rod length");
                ui.add(egui::Slider::new(&mut self.outer.radius, 0.1..=500.0));
            });
            ui.horizontal(|ui| {
                ui.label("Inner point mass");
                ui.add(egui::Slider::new(&mut self.inner_mass, 0.1..=MASS_MAX));
            });
            ui.horizontal(|ui| {
                ui.label("Outer point mass");
                ui.add(egui::Slider::new(&mut self.outer_mass, 0.1..=MASS_MAX));
            });
            ui.horizontal(|ui| {
                ui.label("Gravity");
                ui.add(egui::Slider::new(&mut self.gravity, 0.1..=MASS_GRAVITY));
            });

            ui.horizontal(|ui| {
                ui.label("Steps per tick point mass");
                ui.add(egui::Slider::new(
                    &mut self.steps_per_tick,
                    1..=MAX_STEPS_PER_TICK,
                ));
            });
            ui.horizontal(|ui| {
                ui.label("Desired frame rate");
                ui.add(egui::Slider::new(&mut self.frame_rate, 1.0..=MAX_FRAMERATE));
            });

            ui.horizontal(|ui| {
                ui.label("Accuracy of simulation (lower more accurate)");
                ui.add(
                    egui::Slider::new(&mut self.granularity, 0.001..=MAX_GRANULARITY)
                        .fixed_decimals(3)
                        .step_by(0.0001),
                );
            });

            let available_size = ui.available_size();
            let center = Pos2 {
                x: available_size.x / 2.0,
                y: available_size.y / 2.0,
            };
            // Find the center of the screen

            self.reset_button(ui);

            if self.playing && ui.button("pause").clicked() {
                self.playing = false;
            }

            if !self.playing && ui.button("resume").clicked() {
                self.playing = true;
            }

            self.handle_primary_mouse_button_down(ui, center, self.inner.radius, self.outer.radius);

            if !self.is_dragging && self.playing {
                self.simulation_step();
                self.counter += 1;
            }

            ui.label(format!("tick_count:{}", self.counter));
            self.draw_pendulum(ui, center);
        });

        // roughly 1/144 seconds
        ctx.request_repaint_after(Duration::from_nanos(
            (1000000000.0 / self.frame_rate) as u64,
        ));
    }
}

impl MyEguiApp {
    pub fn new(
        _cc: &eframe::CreationContext<'_>,
        inner_rod_length: f32,
        outer_rod_length: f32,
    ) -> Self {
        Self {
            is_dragging: false,
            inner: Polar2 {
                radius: inner_rod_length,
                theta: PI,
            },
            outer: Polar2 {
                radius: outer_rod_length,
                theta: PI,
            },
            theta_one_dot: 0.0,
            theta_two_dot: 0.0,
            inner_mass: 15.0,
            outer_mass: 15.0,
            counter: 0,
            steps_per_tick: 10,
            granularity: 0.01,
            gravity: 9.8,
            frame_rate: 144.0,
            playing: true,
            is_dragging_p1: false,
            is_dragging_p2: false,
        }
    }
}

impl MyEguiApp {
    /// Resets the pendulums positions to straight up
    fn reset_button(&mut self, ui: &mut Ui) {
        if ui.button("reset").clicked() {
            self.inner = Polar2 {
                radius: self.inner.radius,
                theta: PI,
            };
            self.outer = Polar2 {
                radius: self.outer.radius,
                theta: PI,
            };
            self.theta_one_dot = 0.0;
            self.theta_two_dot = 0.0;
        }
    }

    /// The function responsible for allowing the user to move the points
    ///
    /// If a node is being grabbed it's position will be placed to as close to the cursor as it can.  
    /// The pendulums velocity is set to 0 each time this function is called and a has been and is dragging or is being clicked
    fn handle_primary_mouse_button_down(
        &mut self,
        ui: &mut Ui,
        center: Pos2,
        inner_rod_length: f32,
        outer_rod_length: f32,
    ) {
        if ui.input().pointer.primary_down() {
            if let Some(pointer_pos) = { ui.input().pointer.interact_pos() } {
                // This should probably be fixed
                if self.is_dragging {
                    if self.is_dragging_p2 {
                        let pointer_polar = cartesian_to_polar(pointer_pos, center); // Might be able to not have this be the center later
                        if pointer_polar.radius >= inner_rod_length + outer_rod_length {
                            self.outer.theta = pointer_polar.theta;
                            self.inner.theta = pointer_polar.theta;
                        } else {
                            self.outer = pointer_polar;
                            self.dragging_p2(center, outer_rod_length);
                        }
                    } else if self.is_dragging_p1 {
                        let radius = self.inner.radius;
                        self.inner = cartesian_to_polar(pointer_pos, center);
                        self.inner.radius = radius;
                    } else {
                        unreachable!();
                    }
                } else if intersect_circle(
                    polar_to_cartesian(&self.outer, polar_to_cartesian(&self.inner, center)),
                    15.0,
                    pointer_pos,
                ) {
                    let pointer_polar = cartesian_to_polar(pointer_pos, center); // Might be able to not have this be the center later
                    if pointer_polar.radius >= inner_rod_length + outer_rod_length {
                        self.outer.theta = pointer_polar.theta;
                        self.inner.theta = pointer_polar.theta;
                    } else {
                        self.outer = pointer_polar;
                        self.dragging_p2(center, outer_rod_length);
                    }
                    self.is_dragging = true;
                    self.is_dragging_p2 = true;
                } else if intersect_circle(
                    polar_to_cartesian(&self.inner, center),
                    15.0,
                    pointer_pos,
                ) {
                    self.is_dragging = true;
                    self.is_dragging_p1 = true;

                    let radius = self.inner.radius;
                    self.inner = cartesian_to_polar(pointer_pos, center);
                    self.inner.radius = radius;
                }
            }
        } else if self.is_dragging {
            self.is_dragging = false;
            self.is_dragging_p1 = false;
            self.is_dragging_p2 = false;
            self.theta_one_dot = 0.0;
            self.theta_two_dot = 0.0;
        }
    }

    /// At the moment this function relies on self.outer being a polar coordinate with the center being the center
    /// and not self.inner
    fn dragging_p2(&mut self, center: Pos2, outer_rod_length: f32) {
        // Update inner based on where outer is
        //ui.painter().circle(polar_to_cartesian(&self.outer, center), self.outer.radius, fill_color, stroke)
        let distance = self.outer.radius;
        if let Some((p1, p2)) = calculate_intersecting_points(
            center,
            polar_to_cartesian(&self.outer, center),
            distance,
            self.inner.radius,
            outer_rod_length, // might have outer and inner backwards
        ) {
            let inner = polar_to_cartesian(&self.inner, center);

            let delta_p1 = ((p1.x - inner.x).powi(2) + (p1.y - inner.y).powi(2)).sqrt();
            let delta_p2 = ((p2.x - inner.x).powi(2) + (p2.y - inner.y).powi(2)).sqrt();
            let _nan = Pos2 {
                x: f32::NAN,
                y: f32::NAN,
            };

            match (p1.x.is_nan(), p2.x.is_nan()) {
                (true, true) => {
                    self.outer = cartesian_to_polar(
                        polar_to_cartesian(&self.outer, center),
                        polar_to_cartesian(&self.inner, center),
                    );
                    self.outer.radius = outer_rod_length;
                }
                (true, false) | (false, true) => unreachable!(),
                (false, false) => {
                    if delta_p1 > delta_p2 {
                        self.inner = cartesian_to_polar(p2, center);
                    } else {
                        self.inner = cartesian_to_polar(p1, center);
                    }
                    self.outer = cartesian_to_polar(
                        polar_to_cartesian(&self.outer, center),
                        polar_to_cartesian(&self.inner, center),
                    );
                }
            }
        }
    }

    fn draw_pendulum(&mut self, ui: &mut Ui, center: Pos2) {
        let inner = polar_to_cartesian(&self.inner, center);
        let outer = polar_to_cartesian(&self.outer, inner);
        ui.painter().circle(
            center,
            15.0,
            Color32::GREEN,
            Stroke::new(1.0, Color32::GOLD),
        );
        ui.painter().circle(
            inner,
            60.0 * (self.inner_mass / MASS_MAX),
            Color32::RED,
            Stroke::default(),
        );
        ui.painter().circle(
            outer,
            60.0 * (self.outer_mass / MASS_MAX),
            Color32::BLUE,
            Stroke::default(),
        );
        ui.painter()
            .line_segment([inner, center], Stroke::new(3.0, Color32::LIGHT_RED));
        ui.painter()
            .line_segment([outer, inner], Stroke::new(3.0, Color32::LIGHT_RED));
    }

    /// Calculates the theta_double_dot values and use those to calculate and set new theta values.
    fn simulation_step(&mut self) {
        for _ in 0..self.steps_per_tick {
            let theta_two_double_dot = calculate_theta_two_double_dot(
                self.theta_one_dot,
                self.inner.theta,
                self.outer.theta,
                self.inner_mass,
                self.outer_mass,
                self.inner.radius,
                self.outer.radius,
                self.gravity,
                self.theta_two_dot,
            );
            let theta_one_double_dot = calculate_theta_one_double_dot(
                theta_two_double_dot,
                self.theta_two_dot,
                self.inner.theta,
                self.outer.theta,
                self.inner_mass,
                self.outer_mass,
                self.inner.radius,
                self.outer.radius,
                self.gravity,
            );

            // Now calculate the updated theta_dot_values
            self.theta_two_dot += theta_two_double_dot * self.granularity;
            self.theta_one_dot += theta_one_double_dot * self.granularity;

            let theta_two = self.outer.theta + self.theta_two_dot * self.granularity;
            let theta_one = self.inner.theta + self.theta_one_dot * self.granularity;

            self.inner.theta = theta_one;
            self.outer.theta = theta_two;
        }
    }
}
