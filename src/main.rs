use std::f32::consts::PI;

use eframe::egui;
use egui::{Color32, Pos2, Stroke, Ui};

fn main() {
    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "My egui App",
        native_options,
        Box::new(|cc| Box::new(MyEguiApp::new(cc, 120.0, 100.0))),
    );
}

#[derive(Default, Debug)]
struct MyEguiApp {
    inner: Polar2,
    inner_rod_length: f32,
    outer_rod_length: f32,
    theta_one_dot: f32,
    theta_two_dot: f32,
    outer: Polar2,
    point: Pos2,
    is_dragging: bool,
    inner_mass: f32,
    outer_mass: f32,
}

impl MyEguiApp {
    fn new(
        _cc: &eframe::CreationContext<'_>,
        inner_rod_length: f32,
        outer_rod_length: f32,
    ) -> Self {
        Self {
            point: Pos2 { x: 1.0, y: 2.0 },
            is_dragging: false,
            inner: Polar2 {
                radius: inner_rod_length,
                theta: PI,
            },
            outer: Polar2 {
                radius: outer_rod_length,
                theta: PI,
            },
            inner_rod_length,
            outer_rod_length,
            theta_one_dot: 0.0,
            theta_two_dot: 0.0,
            inner_mass: 1.0,
            outer_mass: 1.0,
        }
    }
}

///
/// This will return the polar coordinates as if y was going up.
/// That is theta is negative to what would be expected with y being positive going down the screen.
///
/// The reasoning for this is that the equations of motion I am taking use y positive going up
/// and it's been a while since I've done this level of math in school.  So I'm opting to just change
/// my coordinates
fn cartesian_to_polar(input: Pos2, center: Pos2) -> Polar2 {
    //https://brilliant.org/wiki/convert-cartesian-coordinates-to-polar/
    // I'm just going to keep operating on the normal cartesian plane with y being positive up
    // When the simulation breaks down I'll deal with it

    let delta_y = center.y - input.y;
    let delta_x = center.x - input.x;
    //dbg!((delta_y, delta_x));

    let mut theta = -f32::atan(delta_y / delta_x);
    if delta_x > 0.0 {
        theta += PI;
    }
    // if theta.is_nan() {
    //     if delta_y > 0.0 {
    //         theta = 3.0 * PI / 2.0;
    //     } else {
    //         theta = PI / 2.0;
    //     }
    // }
    Polar2 {
        radius: (delta_x.powi(2) + delta_y.powi(2)).sqrt(),
        theta,
    }
}

///
/// This will return the x,y coordinate with respect to the origin being the top left of the screen
/// with pos x going right
/// with pos y going down
///
/// This relies on theta of PI / 2.0 being y positive.
// the phrasing of this is weird
fn polar_to_cartesian(input: &Polar2, center: Pos2) -> Pos2 {
    Pos2 {
        x: center.x + (input.radius * f32::cos(input.theta)),
        y: center.y + -1.0 * (input.radius * f32::sin(input.theta)),
    }
}

impl MyEguiApp {
    fn reset_button(&mut self, ui: &mut Ui, inner_rod_length: f32, outer_rod_length: f32) {
        if ui.button("reset").clicked() {
            self.inner = Polar2 {
                radius: inner_rod_length,
                theta: PI,
            };
            self.outer = Polar2 {
                radius: outer_rod_length,
                theta: PI,
            };
            self.theta_one_dot = 0.0;
            self.theta_two_dot = 0.0;
        }
    }

    fn handle_primary_mouse_button_down(
        &mut self,
        ui: &mut Ui,
        center: Pos2,
        total_radius: f32,
        inner_rod_length: f32,
        outer_rod_length: f32,
    ) {
        let mut temp = Pos2 { x: 0.0, y: 0.0 };
        if ui.input().pointer.primary_down() {
            if let Some(thing) = { ui.input().pointer.interact_pos() } {
                //temp = polar_to_cartesian(&cartesian_to_polar(thing, center), center);
                //let thing = cartesian_to_polar(thing, center);
                if self.is_dragging {
                    self.outer = cartesian_to_polar(thing, center); // Might be able to not have this be the cneter later

                    if self.outer.radius > total_radius {
                        self.outer.radius = total_radius
                    }
                } else if _intersect_circle(polar_to_cartesian(&self.outer, polar_to_cartesian(&self.inner, center)), 15.0, thing) {
                    self.outer = cartesian_to_polar(thing, center); // Might be able to not have this be the cneter later
                    self.is_dragging = true;
                }

                if self.is_dragging {
                    // Update inner based on where outer is
                    let distance = self.outer.radius;
                    if let Some((p1, p2)) = calculate_intersecting_points(
                        center,
                        polar_to_cartesian(&self.outer, center),
                        distance,
                        inner_rod_length,
                        outer_rod_length // might have outer and inner backwards
                    ) {
                        let inner = polar_to_cartesian(&self.inner, center);
                        let delta_p1 = ((p1.x - inner.x).powi(2) + (p1.y - inner.y).powi(2)).sqrt();
                        let delta_p2 = ((p2.x - inner.x).powi(2) + (p2.y - inner.y).powi(2)).sqrt();
                        let _nan = Pos2 {
                            x: f32::NAN,
                            y: f32::NAN,
                        };

                        if p1.x.is_nan() && p2.x.is_nan() {
                            // dbg!((p1,p2));
                        } else if delta_p1 > delta_p2 {
                            self.inner = cartesian_to_polar(p2, center);
                        } else {
                            self.inner = cartesian_to_polar(p1, center);
                        }
                            self.outer = cartesian_to_polar(polar_to_cartesian(&self.outer, center), polar_to_cartesian(&self.inner, center))
                    }
                } // end update inner
            }
        } else {
            self.is_dragging = false;
            self.theta_one_dot = 0.0;
            self.theta_two_dot = 0.0;
        }
        ui.painter()
            .circle(temp, 10.0, Color32::WHITE, Stroke::default());
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
        ui.painter()
            .circle(inner, 15.0, Color32::RED, Stroke::default());
        ui.painter()
            .circle(outer, 15.0, Color32::BLUE, Stroke::default());
        ui.painter()
            .line_segment([inner, center], Stroke::new(3.0, Color32::LIGHT_RED));
        ui.painter()
            .line_segment([outer, inner], Stroke::new(3.0, Color32::LIGHT_RED));
    }
}
impl eframe::App for MyEguiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.label(format!("{},{}", self.point.x, self.point.y));
            let available_size = ui.available_size();
            let total_radius: f32 = self.inner_rod_length + self.outer_rod_length;
            let time_interval = 0.1; // This should have a slider or something
            let center = Pos2 {
                x: available_size.x / 2.0,
                y: available_size.y / 2.0,
            };
            // Find the center of the screen

            self.reset_button(ui, self.inner_rod_length, self.outer_rod_length);

            self.handle_primary_mouse_button_down(
                ui,
                center,
                total_radius,
                self.outer_rod_length,
                self.inner_rod_length,
            );

            self.draw_pendulum(ui, center);
        });
    }
}


/// this uses the sqrt(deltax ^2 + deltay^2) so if the delta is great enough an overflow may occur
fn _intersect_circle(object_center: Pos2, object_radius: f32, clicked_pos: Pos2) -> bool {
    // This seems scuffed with both -
    object_radius
        >= ((object_center.x - clicked_pos.x).powi(2) + (object_center.y - clicked_pos.y).powi(2))
            .sqrt()
}

#[derive(Default, Debug, PartialEq)]
struct Polar2 {
    radius: f32,
    theta: f32,
}

fn intersect_circle_polar(main_object: &Polar2, object_radius: f32, other_object: &Polar2) -> bool {
    // https://www.kristakingmath.com/blog/distance-between-polar-points
    //dbg!((main_object,other_object));
    object_radius
        >= (main_object.radius.powi(2) + other_object.radius.powi(2)
            - 2.0
                * main_object.radius
                * other_object.radius
                * f32::cos(main_object.theta - other_object.theta))
        .sqrt()
}

/// If there is too much overlap there nan will be returned.
// This is because of taking the square root of a negative when x^2 > a_radius^2
fn calculate_intersecting_points(
    b: Pos2,
    a: Pos2,
    d: f32,
    a_radius: f32,
    b_radius: f32,
) -> Option<(Pos2, Pos2)> {
    let ex = (b.x - a.x) / d;
    let ey = (b.y - a.y) / d;

    let x = (a_radius * a_radius - b_radius * b_radius + d * d) / (2.0 * d);
    // This line is the reason for the position not being updated if the outer circle is intersecting the center circle
    let y = (a_radius * a_radius - x * x).sqrt();

    let p1 = Pos2 {
        x: a.x + x * ex - y * ey,
        y: a.y + x * ey + y * ex,
    };

    let p2 = Pos2 {
        x: a.x + x * ex + y * ey,
        y: a.y + x * ey - y * ex,
    };
    Some((p1, p2))
}

#[cfg(test)]
mod tests {
    use std::f32::consts::PI;

    use egui::Pos2;

    use crate::{cartesian_to_polar, polar_to_cartesian, Polar2};

    #[test]
    fn test_polar_to_cartesian_offset_center() {
        // This center will always have a positive x and y
        let center = Pos2 { x: 7.0, y: 3.0 };
        let tolerance = 0.00005;
        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 1.0,
                theta: 0.0,
            },
            center,
            &Pos2 { x: 8.0, y: 3.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: 0.0,
            },
            center,
            &Pos2 { x: 9.0, y: 3.0 },
            tolerance,
        );
        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: PI / 2.0,
            },
            center,
            &Pos2 { x: 7.0, y: 1.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: PI,
            },
            center,
            &Pos2 { x: 5.0, y: 3.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: 3.0 * PI / 2.0,
            },
            center,
            &Pos2 { x: 7.0, y: 5.0 },
            tolerance,
        );
    }
    // Let's see how this works with floats
    // Might what to have a close enough so like
    #[test]
    fn test_polar_to_cartesian_center_origin() {
        let center = Pos2 { x: 0.0, y: 0.0 };
        let tolerance = 0.00005;
        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 1.0,
                theta: 0.0,
            },
            center,
            &Pos2 { x: 1.0, y: 0.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: 0.0,
            },
            center,
            &Pos2 { x: 2.0, y: 0.0 },
            tolerance,
        );
        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: PI / 2.0,
            },
            center,
            &Pos2 { x: 0.0, y: -2.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: PI,
            },
            center,
            &Pos2 { x: -2.0, y: 0.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: 3.0 * PI / 2.0,
            },
            center,
            &Pos2 { x: 0.0, y: 2.0 },
            tolerance,
        );
    }

    fn compare_pos2_with_tolerance(a: &Pos2, b: &Pos2, tolerance: f32) -> bool {
        compare_floats_with_variance(a.x, b.x, tolerance)
            && compare_floats_with_variance(a.y, b.y, tolerance)
    }
    fn compare_polar2_with_tolerance(a: &Polar2, b: &Polar2, tolerance: f32) -> bool {
        println!("---------");
        compare_floats_with_variance(a.radius, b.radius, tolerance)
            && compare_floats_with_variance(f32::sin(a.theta), f32::sin(b.theta), tolerance)
            && compare_floats_with_variance(f32::cos(a.theta), f32::cos(b.theta), tolerance)
    }

    fn test_cartesian_to_polar_helper(
        cartesian: Pos2,
        center: Pos2,
        expected_polar: Polar2,
        tolerance: f32,
    ) {
        let polar_coord = cartesian_to_polar(cartesian, center);
        assert!(
            compare_polar2_with_tolerance(&polar_coord, &expected_polar, tolerance),
            "{:?} -> {:?} did not match {:?}",
            cartesian,
            &polar_coord,
            expected_polar
        )
    }

    fn test_polar_to_cartesian_helper(
        polar: &Polar2,
        center: Pos2,
        expected_cartesian: &Pos2,
        tolerance: f32,
    ) {
        let cartesian_coord = polar_to_cartesian(polar, center);
        assert!(
            compare_pos2_with_tolerance(&cartesian_coord, expected_cartesian, tolerance),
            "{:?} -> {:?} did not match {:?}",
            polar,
            cartesian_coord,
            expected_cartesian
        )
    }

    fn compare_floats_with_variance(a: f32, b: f32, tolerance: f32) -> bool {
        if a > b {
            a - b <= tolerance
        } else {
            b - a <= tolerance
        }
    }

    #[test]
    fn test_polar_too_cartesian() {
        //todo!();
    }

    #[test]
    fn test_cartesian_to_polar() {
        // This center will always have a positive x and y
        let center = Pos2 { x: 7.0, y: 3.0 };
        let tolerance = 0.00005;
        test_cartesian_to_polar_helper(
            Pos2 { x: 8.0, y: 3.0 },
            center,
            Polar2 {
                radius: 1.0,
                theta: 0.0,
            },
            tolerance,
        );
        test_cartesian_to_polar_helper(
            Pos2 { x: 9.0, y: 3.0 },
            center,
            Polar2 {
                radius: 2.0,
                theta: 0.0,
            },
            tolerance,
        );
        test_cartesian_to_polar_helper(
            Pos2 { x: 7.0, y: 1.0 },
            center,
            Polar2 {
                radius: 2.0,
                theta: -PI / 2.0,
            },
            tolerance,
        );

        test_cartesian_to_polar_helper(
            Pos2 { x: 5.0, y: 3.0 },
            center,
            Polar2 {
                radius: 2.0,
                theta: PI,
            },
            tolerance,
        );

        test_cartesian_to_polar_helper(
            Pos2 { x: 7.0, y: 5.0 },
            center,
            Polar2 {
                radius: 2.0,
                theta: -3.0 * PI / 2.0,
            },
            tolerance,
        );
    }
}
