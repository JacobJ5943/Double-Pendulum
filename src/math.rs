use std::f32::consts::PI;

use egui::Pos2;

#[derive(Default, Debug, PartialEq)]
pub struct Polar2 {
    pub radius: f32,
    pub theta: f32,
}

///
/// This will return the polar coordinates as if y was going up.
/// That is theta is negative to what would be expected with y being positive going down the screen.
///
/// The reasoning for this is that the equations of motion I am taking use y positive going up
/// and it's been a while since I've done this level of math in school.  So I'm opting to just change
/// my coordinates
pub fn cartesian_to_polar(input: Pos2, center: Pos2) -> Polar2 {
    //https://brilliant.org/wiki/convert-cartesian-coordinates-to-polar/
    // I'm just going to keep operating on the normal cartesian plane with y being positive up
    // When the simulation breaks down I'll deal with it

    let delta_y = center.y - input.y;
    let delta_x = center.x - input.x;

    let mut theta = -f32::atan(delta_y / delta_x);
    if delta_x > 0.0 {
        theta += PI;
    }

    // If deltax is 0 then it's on the center and it's vertical
    // I should probably just not do the calculation if input.x == center.x
    if theta.is_nan() {
        if delta_y > 0.0 {
            // Theta is 0.0 because the polar coordinates are orientated such that theta = 0.0 is y positive or straight down on the screen
            theta = 0.0;
        } else {
            theta = PI;
        }
    }
    Polar2 {
        radius: (delta_x.powi(2) + delta_y.powi(2)).sqrt(),
        theta: theta + (PI / 2.0), // The addition is to rotate the pendulum 90 degrees to clockwise.
    }
}

///
/// This will return the x,y coordinate with respect to the origin being the top left of the screen
/// with pos x going right
/// with pos y going down
///
/// This relies on theta of PI / 2.0 being y positive.
// the phrasing of this is weird
pub fn polar_to_cartesian(input: &Polar2, center: Pos2) -> Pos2 {
    Pos2 {
        x: center.x + (input.radius * f32::cos(input.theta - (PI / 2.0))),
        y: center.y + -1.0 * (input.radius * f32::sin(input.theta - (PI / 2.0))),
    }
}

//  source https://planetcalc.com/8098/
/// If you point returned is NaN the other is guaranteed to be NaN as well
pub fn calculate_intersecting_points(
    pos_a: Pos2,
    pos_b: Pos2,
    distance: f32, // distance between pos_a and pos_b
    a_radius: f32,
    b_radius: f32,
) -> Option<(Pos2, Pos2)> {
    let a = (a_radius.powi(2) - b_radius.powi(2) + distance.powi(2)) / (2.0 * distance);
    let h = (a_radius.powi(2) - a.powi(2)).sqrt();
    let x_3 = pos_a.x + (a / distance) * (pos_b.x - pos_a.x);
    let y_3 = pos_a.y + (a / distance) * (pos_b.y - pos_a.y);

    let p1 = Pos2 {
        x: x_3 + (h / distance) * (pos_b.y - pos_a.y),
        y: y_3 - (h / distance) * (pos_b.x - pos_a.x),
    };

    let p2 = Pos2 {
        x: x_3 - (h / distance) * (pos_b.y - pos_a.y),
        y: y_3 + (h / distance) * (pos_b.x - pos_a.x),
    };

    Some((p1, p2))
}

/// Equation of motion
/// If this works
/// $\ddot{\theta_1}=-\dfrac{a}{k}\ddot{\theta_2}\cos(\theta_1-\theta_2)-\dfrac{a}{k}\dot{\theta_2^2}\sin(\theta_1-\theta_2)-\dfrac{f}{k}\sin(\theta_1)$
pub fn calculate_theta_one_double_dot(
    theta_two_double_dot: f32,
    theta_two_dot: f32,
    theta_one: f32,
    theta_two: f32,
    mass_one: f32,
    mass_two: f32,
    r_one: f32,
    r_two: f32,
    g: f32,
) -> f32 {
    let a = mass_two * r_two;
    let _b = mass_two * r_one;
    let _c = mass_two * g;
    let f = (mass_one + mass_two) * g;
    let k = (mass_one + mass_two) * r_one;

    -(a / k) * theta_two_double_dot * f32::cos(theta_one - theta_two)
        - ((a / k) * (theta_two_dot.powi(2)) * f32::sin(theta_one - theta_two))
        - (f / k) * f32::sin(theta_one)
}

/// Equation if this works
/// $\ddot{\theta_2}=(\dfrac{1}{a-\dfrac{ba}{k}\cos^2(\theta_1 - \theta_2)})(\dfrac{ab}{k}\dot{\theta_2}\sin(\theta_1-\theta_2)\cos(\theta_1-\theta_2)+\dfrac{fb}{k}\sin(\theta_1)\cos(\theta_1-\theta_2)+b\dot{\theta_1^2}\sin(\theta_1-\theta_2)-c\sin(\theta_2))$
pub fn calculate_theta_two_double_dot(
    theta_one_dot: f32,
    theta_one: f32,
    theta_two: f32,
    mass_one: f32,
    mass_two: f32,
    r_one: f32,
    r_two: f32,
    g: f32,
    theta_two_dot: f32,
) -> f32 {
    let a = mass_two * r_two;
    let b = mass_two * r_one;
    let c = mass_two * g;
    let f = (mass_one + mass_two) * g;
    let k = (mass_one + mass_two) * r_one;

    (1.0 / (a - ((b * a) / k) * (f32::cos(theta_one - theta_two).powi(2))))
        * ((((a * b) / k)
            * theta_two_dot
            * f32::sin(theta_one - theta_two)
            * f32::cos(theta_one - theta_two))
            + (((f * b) / k) * f32::sin(theta_one) * f32::cos(theta_one - theta_two)
                + b * (theta_one_dot.powi(2)) * f32::sin(theta_one - theta_two)
                - c * f32::sin(theta_two)))
}

/// this uses the sqrt(deltax ^2 + deltay^2) so if the delta is great enough an overflow may occur
pub fn intersect_circle(object_center: Pos2, object_radius: f32, clicked_pos: Pos2) -> bool {
    // This seems scuffed with both -
    object_radius
        >= ((object_center.x - clicked_pos.x).powi(2) + (object_center.y - clicked_pos.y).powi(2))
            .sqrt()
}

pub fn _intersect_circle_polar(
    main_object: &Polar2,
    object_radius: f32,
    other_object: &Polar2,
) -> bool {
    // https://www.kristakingmath.com/blog/distance-between-polar-points
    object_radius
        >= (main_object.radius.powi(2) + other_object.radius.powi(2)
            - 2.0
                * main_object.radius
                * other_object.radius
                * f32::cos(main_object.theta - other_object.theta))
        .sqrt()
}

#[cfg(test)]
mod tests {
    use std::f32::consts::PI;

    use egui::Pos2;

    use crate::math::{cartesian_to_polar, polar_to_cartesian, Polar2};

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
            &Pos2 { x: 7.0, y: 4.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: 0.0,
            },
            center,
            &Pos2 { x: 7.0, y: 5.0 },
            tolerance,
        );
        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: PI / 2.0,
            },
            center,
            &Pos2 { x: 9.0, y: 3.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: PI,
            },
            center,
            &Pos2 { x: 7.0, y: 1.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: 3.0 * PI / 2.0,
            },
            center,
            &Pos2 { x: 5.0, y: 3.0 },
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
            &Pos2 { x: 0.0, y: 1.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: 0.0,
            },
            center,
            &Pos2 { x: 0.0, y: 2.0 },
            tolerance,
        );
        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: PI / 2.0,
            },
            center,
            &Pos2 { x: 2.0, y: 0.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: PI,
            },
            center,
            &Pos2 { x: 0.0, y: -2.0 },
            tolerance,
        );

        test_polar_to_cartesian_helper(
            &Polar2 {
                radius: 2.0,
                theta: 3.0 * PI / 2.0,
            },
            center,
            &Pos2 { x: -2.0, y: 0.0 },
            tolerance,
        );
    }

    fn compare_pos2_with_tolerance(a: &Pos2, b: &Pos2, tolerance: f32) -> bool {
        compare_floats_with_variance(a.x, b.x, tolerance)
            && compare_floats_with_variance(a.y, b.y, tolerance)
    }
    fn compare_polar2_with_tolerance(a: &Polar2, b: &Polar2, tolerance: f32) -> bool {
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
    fn test_cartesian_to_polar() {
        // This center will always have a positive x and y
        let center = Pos2 { x: 7.0, y: 3.0 };
        let tolerance = 0.00005;
        test_cartesian_to_polar_helper(
            Pos2 { x: 8.0, y: 3.0 },
            center,
            Polar2 {
                radius: 1.0,
                theta: PI / 2.0,
            },
            tolerance,
        );
        test_cartesian_to_polar_helper(
            Pos2 { x: 9.0, y: 3.0 },
            center,
            Polar2 {
                radius: 2.0,
                theta: PI / 2.0,
            },
            tolerance,
        );
        test_cartesian_to_polar_helper(
            Pos2 { x: 7.0, y: 1.0 },
            center,
            Polar2 {
                radius: 2.0,
                theta: 0.0,
            },
            tolerance,
        );

        test_cartesian_to_polar_helper(
            Pos2 { x: 5.0, y: 3.0 },
            center,
            Polar2 {
                radius: 2.0,
                theta: 3.0 * PI / 2.0,
            },
            tolerance,
        );

        test_cartesian_to_polar_helper(
            Pos2 { x: 7.0, y: 5.0 },
            center,
            Polar2 {
                radius: 2.0,
                theta: PI,
            },
            tolerance,
        );
    }
}
