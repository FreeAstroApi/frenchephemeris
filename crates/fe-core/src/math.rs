// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Small fixed-size linear algebra and angle helpers shared by every module.

pub const PI: f64 = std::f64::consts::PI;
pub const TAU: f64 = std::f64::consts::TAU;
pub const DEG: f64 = PI / 180.0;
pub const ARCSEC: f64 = DEG / 3600.0;
pub const J2000: f64 = 2451545.0;
pub const DAYS_PER_CENTURY: f64 = 36525.0;
/// Speed of light, km/s (IAU).
pub const C_KM_S: f64 = 299_792.458;
pub const C_KM_DAY: f64 = C_KM_S * 86400.0;
/// Astronomical unit, km (IAU 2012).
pub const AU_KM: f64 = 149_597_870.700;

pub type Vec3 = [f64; 3];
/// Row-major 3x3 matrix: `m[row][col]`.
pub type Mat3 = [[f64; 3]; 3];

pub const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

#[inline]
pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline]
pub fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
pub fn mat_vec(m: &Mat3, v: Vec3) -> Vec3 {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Transpose times vector (inverse rotation).
#[inline]
pub fn mat_t_vec(m: &Mat3, v: Vec3) -> Vec3 {
    [
        m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
        m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
        m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2],
    ]
}

pub fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    r
}

pub fn transpose(m: &Mat3) -> Mat3 {
    [[m[0][0], m[1][0], m[2][0]], [m[0][1], m[1][1], m[2][1]], [m[0][2], m[1][2], m[2][2]]]
}

/// Rotation of the *frame* about x by angle `a` (SOFA `iauRx` convention).
pub fn rot_x(a: f64) -> Mat3 {
    let (s, c) = a.sin_cos();
    [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]]
}
/// Rotation of the frame about y (SOFA `iauRy`).
pub fn rot_y(a: f64) -> Mat3 {
    let (s, c) = a.sin_cos();
    [[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]]
}
/// Rotation of the frame about z (SOFA `iauRz`).
pub fn rot_z(a: f64) -> Mat3 {
    let (s, c) = a.sin_cos();
    [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
}

/// Normalize to [0, 2π).
#[inline]
pub fn norm_rad(a: f64) -> f64 {
    let r = a.rem_euclid(TAU);
    if r >= TAU {
        0.0
    } else {
        r
    }
}
/// Normalize to [0, 360).
#[inline]
pub fn norm_deg(a: f64) -> f64 {
    let r = a.rem_euclid(360.0);
    if r >= 360.0 {
        0.0
    } else {
        r
    }
}
/// Signed difference a-b in (-180, 180].
#[inline]
pub fn diff_deg(a: f64, b: f64) -> f64 {
    let d = norm_deg(a - b);
    if d > 180.0 {
        d - 360.0
    } else {
        d
    }
}

/// Cartesian -> (lon [0,2π), lat, r).
#[inline]
pub fn to_spherical(v: Vec3) -> (f64, f64, f64) {
    let rxy2 = v[0] * v[0] + v[1] * v[1];
    let r = (rxy2 + v[2] * v[2]).sqrt();
    let lon = if rxy2 == 0.0 { 0.0 } else { norm_rad(v[1].atan2(v[0])) };
    let lat = if r == 0.0 { 0.0 } else { v[2].atan2(rxy2.sqrt()) };
    (lon, lat, r)
}

/// Rates of (lon, lat, r) given position and velocity.
#[inline]
pub fn spherical_rates(p: Vec3, v: Vec3) -> (f64, f64, f64) {
    let rxy2 = p[0] * p[0] + p[1] * p[1];
    let r2 = rxy2 + p[2] * p[2];
    let r = r2.sqrt();
    let rxy = rxy2.sqrt();
    let dlon = if rxy2 == 0.0 { 0.0 } else { (p[0] * v[1] - p[1] * v[0]) / rxy2 };
    let dr = if r == 0.0 { 0.0 } else { dot(p, v) / r };
    let dlat = if rxy == 0.0 || r2 == 0.0 {
        0.0
    } else {
        (v[2] * rxy2 - p[2] * (p[0] * v[0] + p[1] * v[1])) / (r2 * rxy)
    };
    (dlon, dlat, dr)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rates_match_finite_difference() {
        let p = [1.0, 2.0, 0.5];
        let v = [-0.3, 0.1, 0.2];
        let h = 1e-6;
        let (l0, b0, r0) = to_spherical(sub(p, scale(v, h)));
        let (l1, b1, r1) = to_spherical(add(p, scale(v, h)));
        let (dl, db, dr) = spherical_rates(p, v);
        assert!(((l1 - l0) / (2.0 * h) - dl).abs() < 1e-8);
        assert!(((b1 - b0) / (2.0 * h) - db).abs() < 1e-8);
        assert!(((r1 - r0) / (2.0 * h) - dr).abs() < 1e-8);
    }
}
