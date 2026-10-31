// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! Charts for an observer away from Earth (cargo feature `space`): a site on the surface
//! of another planet or of the Moon, the centre of a body, or any point in space.
//!
//! - Orientation of each world: IAU WGCCRE 2015 rotation models (Archinal et al. 2018,
//!   Celest. Mech. Dyn. Astr. 130:22): pole right ascension and declination in the ICRF and
//!   prime meridian angle W, functions of TDB. Body-fixed frame `R3(W)·R1(90°−δ0)·R3(90°+α0)`.
//! - Sites: planetodetic latitude and east longitude on the IAU 2015 reference ellipsoid,
//!   height in km above it. East longitudes are positive for every world (right-handed,
//!   planetocentric convention), also for the IAU west-positive ones.
//! - Positions: the engine's pipeline (light time, solar deflection, aberration from the
//!   observer's barycentric velocity, rotation included) with the observer moved off Earth.
//! - Zodiac: [`Zodiac::Earth`] expresses everything in Earth's true ecliptic and equinox of
//!   date, as an ordinary chart. [`Zodiac::Native`] uses the world's own orbital plane
//!   (instantaneous heliocentric orbit normal; the Earth Moon barycentre's for the Moon) with
//!   0° at the world's vernal equinox, the ascending node of its equator (IAU north pole) on
//!   that plane.
//! - Houses: the house circle is the world's equator, oriented by its spin axis, so the
//!   Ascendant is always the rising point of the zodiac even on retrograde rotators (Venus,
//!   Uranus, Pluto). The house routines receive the sidereal angle, latitude and obliquity in
//!   that frame; cusps are then shifted to the chosen zodiac's zero point.

use crate::context::{Body, Engine};
use crate::ephem::{apparent, Corrections, Epoch, LightTime, ObsState, SSB, SUN};
use crate::houses::{self, Houses};
use crate::math::{
    add, cross, dot, mat_mul, mat_t_vec, mat_vec, norm, norm_deg, rot_x, rot_z, scale, to_spherical, Mat3,
    Vec3, AU_KM, DEG,
};
use crate::{Error, Result};

/// Days per Julian century.
const CENTURY: f64 = 36525.0;
const J2000: f64 = 2451545.0;
/// Earth rotation rate (IERS nominal 7.2921150e-5 rad/s), radians per day.
const EARTH_RATE_RAD_DAY: f64 = 7.292_115_0e-5 * 86400.0;
/// Step for numerical speeds, days.
const SPEED_STEP: f64 = 1e-3;

/// A body an observer can stand on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum World {
    Mercury,
    Venus,
    Earth,
    Moon,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
}

/// Orientation of a world at one instant, IAU convention, degrees.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orientation {
    /// Right ascension of the IAU north pole, ICRF.
    pub pole_ra: f64,
    /// Declination of the IAU north pole, ICRF.
    pub pole_dec: f64,
    /// Prime meridian angle W, measured along the equator from its ascending node on the
    /// ICRF equator.
    pub prime_meridian: f64,
    /// dW/dt, degrees per day (negative for retrograde rotators).
    pub rotation_rate: f64,
}

impl Orientation {
    /// Rotation from ICRF axes to the body-fixed frame.
    pub fn icrf_to_body(&self) -> Mat3 {
        mat_mul(
            &rot_z(self.prime_meridian * DEG),
            &mat_mul(&rot_x((90.0 - self.pole_dec) * DEG), &rot_z((90.0 + self.pole_ra) * DEG)),
        )
    }

    /// Unit vector of the IAU north pole, ICRF.
    pub fn pole(&self) -> Vec3 {
        let (sa, ca) = (self.pole_ra * DEG).sin_cos();
        let (sd, cd) = (self.pole_dec * DEG).sin_cos();
        [cd * ca, cd * sa, sd]
    }

    /// Unit vector of the spin axis (the pole about which the world turns counterclockwise).
    pub fn spin_axis(&self) -> Vec3 {
        let p = self.pole();
        if self.rotation_rate < 0.0 {
            scale(p, -1.0)
        } else {
            p
        }
    }
}

fn sind(x: f64) -> f64 {
    (x * DEG).sin()
}
fn cosd(x: f64) -> f64 {
    (x * DEG).cos()
}

impl World {
    pub const ALL: [World; 10] = [
        World::Mercury,
        World::Venus,
        World::Earth,
        World::Moon,
        World::Mars,
        World::Jupiter,
        World::Saturn,
        World::Uranus,
        World::Neptune,
        World::Pluto,
    ];

    pub fn name(self) -> &'static str {
        match self {
            World::Mercury => "mercury",
            World::Venus => "venus",
            World::Earth => "earth",
            World::Moon => "moon",
            World::Mars => "mars",
            World::Jupiter => "jupiter",
            World::Saturn => "saturn",
            World::Uranus => "uranus",
            World::Neptune => "neptune",
            World::Pluto => "pluto",
        }
    }

    /// Case-insensitive name.
    pub fn from_name(name: &str) -> Result<World> {
        let n = name.trim().to_ascii_lowercase();
        World::ALL
            .into_iter()
            .find(|w| w.name() == n)
            .ok_or_else(|| Error::InvalidArgument(format!("unknown world '{name}'")))
    }

    /// The world matching a body, if an observer can stand on it.
    pub fn from_body(body: Body) -> Option<World> {
        World::ALL.into_iter().find(|w| w.body() == body)
    }

    pub fn body(self) -> Body {
        match self {
            World::Mercury => Body::Mercury,
            World::Venus => Body::Venus,
            World::Earth => Body::Earth,
            World::Moon => Body::Moon,
            World::Mars => Body::Mars,
            World::Jupiter => Body::Jupiter,
            World::Saturn => Body::Saturn,
            World::Uranus => Body::Uranus,
            World::Neptune => Body::Neptune,
            World::Pluto => Body::Pluto,
        }
    }

    /// Pack body for the world's centre. The outer planets are their system barycentres,
    /// as for every other position of the engine.
    fn naif(self) -> i32 {
        self.body().naif().expect("worlds are physical bodies")
    }

    /// Pack body whose heliocentric orbit defines the native zodiac.
    fn orbit_naif(self) -> i32 {
        match self {
            World::Earth | World::Moon => crate::ephem::EMB,
            w => w.naif(),
        }
    }

    /// Reference ellipsoid: (equatorial radius, polar radius), km. IAU 2015, except WGS84
    /// for Earth, as the engine's terrestrial observers.
    pub fn radii_km(self) -> (f64, f64) {
        match self {
            World::Mercury => (2440.53, 2438.26),
            World::Venus => (6051.8, 6051.8),
            World::Earth => (6378.137, 6378.137 * (1.0 - 1.0 / 298.257_223_563)),
            World::Moon => (1737.4, 1737.4),
            World::Mars => (3396.19, 3376.20),
            World::Jupiter => (71492.0, 66854.0),
            World::Saturn => (60268.0, 54364.0),
            World::Uranus => (25559.0, 24973.0),
            World::Neptune => (24764.0, 24341.0),
            World::Pluto => (1188.3, 1188.3),
        }
    }

    /// IAU WGCCRE 2015 orientation at `jd_tdb`.
    pub fn orientation(self, jd_tdb: f64) -> Orientation {
        let d = jd_tdb - J2000;
        let t = d / CENTURY;
        let (ra, dec, w, rate) = match self {
            World::Mercury => {
                let m = |a: f64, b: f64| sind(a + b * d);
                let w = 329.5988 + 6.1385108 * d + 0.01067257 * m(174.7910857, 4.092335)
                    - 0.00112309 * m(349.5821714, 8.184670)
                    - 0.00011040 * m(164.3732571, 12.277005)
                    - 0.00002539 * m(339.1643429, 16.369340)
                    - 0.00000571 * m(153.9554286, 20.461675);
                (281.0103 - 0.0328 * t, 61.4155 - 0.0049 * t, w, 6.1385108)
            }
            World::Venus => (272.76, 67.16, 160.20 - 1.4813688 * d, -1.4813688),
            World::Earth => (-0.641 * t, 90.0 - 0.557 * t, 190.147 + 360.9856235 * d, 360.9856235),
            World::Moon => {
                let e = [
                    125.045 - 0.0529921 * d,
                    250.089 - 0.1059842 * d,
                    260.008 + 13.0120009 * d,
                    176.625 + 13.3407154 * d,
                    357.529 + 0.9856003 * d,
                    311.589 + 26.4057084 * d,
                    134.963 + 13.0649930 * d,
                    276.617 + 0.3287146 * d,
                    34.226 + 1.7484877 * d,
                    15.134 - 0.1589763 * d,
                    119.743 + 0.0036096 * d,
                    239.961 + 0.1643573 * d,
                    25.053 + 12.9590088 * d,
                ];
                let s = |i: usize| sind(e[i - 1]);
                let c = |i: usize| cosd(e[i - 1]);
                let ra = 269.9949 + 0.0031 * t - 3.8787 * s(1) - 0.1204 * s(2) + 0.0700 * s(3)
                    - 0.0172 * s(4)
                    + 0.0072 * s(6)
                    - 0.0052 * s(10)
                    + 0.0043 * s(13);
                let dec = 66.5392 + 0.0130 * t + 1.5419 * c(1) + 0.0239 * c(2) - 0.0278 * c(3)
                    + 0.0068 * c(4)
                    - 0.0029 * c(6)
                    + 0.0009 * c(7)
                    + 0.0008 * c(10)
                    - 0.0009 * c(13);
                let w = 38.3213 + 13.17635815 * d - 1.4e-12 * d * d + 3.5610 * s(1) + 0.1208 * s(2)
                    - 0.0642 * s(3)
                    + 0.0158 * s(4)
                    + 0.0252 * s(5)
                    - 0.0066 * s(6)
                    - 0.0047 * s(7)
                    - 0.0046 * s(8)
                    + 0.0028 * s(9)
                    + 0.0052 * s(10)
                    + 0.0040 * s(11)
                    + 0.0019 * s(12)
                    - 0.0044 * s(13);
                (ra, dec, w, 13.17635815)
            }
            World::Mars => {
                let ra = 317.269202 - 0.10927547 * t
                    + 0.000068 * sind(198.991226 + 19139.4819985 * t)
                    + 0.000238 * sind(226.292679 + 38280.8511281 * t)
                    + 0.000052 * sind(249.663391 + 57420.7251593 * t)
                    + 0.000009 * sind(266.183510 + 76560.6367950 * t)
                    + 0.419057 * sind(79.398797 + 0.5042615 * t);
                let dec = 54.432516 - 0.05827105 * t
                    + 0.000051 * cosd(122.433576 + 19139.9407476 * t)
                    + 0.000141 * cosd(43.058401 + 38280.8753272 * t)
                    + 0.000031 * cosd(57.663379 + 57420.7517205 * t)
                    + 0.000005 * cosd(79.476401 + 76560.6495004 * t)
                    + 1.591274 * cosd(166.325722 + 0.5042615 * t);
                let w = 176.049863
                    + 350.891982443297 * d
                    + 0.000145 * sind(129.071773 + 19140.0328244 * t)
                    + 0.000157 * sind(36.352167 + 38281.0473591 * t)
                    + 0.000040 * sind(56.668646 + 57420.9295360 * t)
                    + 0.000001 * sind(67.364003 + 76560.2552215 * t)
                    + 0.000001 * sind(104.792680 + 95700.4387578 * t)
                    + 0.584542 * sind(95.391654 + 0.5042615 * t);
                (ra, dec, w, 350.891982443297)
            }
            World::Jupiter => {
                let ja = 99.360714 + 4850.4046 * t;
                let jb = 175.895369 + 1191.9605 * t;
                let jc = 300.323162 + 262.5475 * t;
                let jd = 114.012305 + 6070.2476 * t;
                let je = 49.511251 + 64.3000 * t;
                let ra = 268.056595 - 0.006499 * t
                    + 0.000117 * sind(ja)
                    + 0.000938 * sind(jb)
                    + 0.001432 * sind(jc)
                    + 0.000030 * sind(jd)
                    + 0.002150 * sind(je);
                let dec = 64.495303
                    + 0.002413 * t
                    + 0.000050 * cosd(ja)
                    + 0.000404 * cosd(jb)
                    + 0.000617 * cosd(jc)
                    - 0.000013 * cosd(jd)
                    + 0.000926 * cosd(je);
                (ra, dec, 284.95 + 870.5360000 * d, 870.5360000)
            }
            World::Saturn => (40.589 - 0.036 * t, 83.537 - 0.004 * t, 38.90 + 810.7939024 * d, 810.7939024),
            World::Uranus => (257.311, -15.175, 203.81 - 501.1600928 * d, -501.1600928),
            World::Neptune => {
                let n = 357.85 + 52.316 * t;
                (
                    299.36 + 0.70 * sind(n),
                    43.46 - 0.51 * cosd(n),
                    249.978 + 541.1397757 * d - 0.48 * sind(n),
                    541.1397757,
                )
            }
            World::Pluto => (132.993, -6.163, 302.695 + 56.3625225 * d, 56.3625225),
        };
        Orientation { pole_ra: ra, pole_dec: dec, prime_meridian: norm_deg(w), rotation_rate: rate }
    }
}

/// A place on a world's surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Site {
    pub world: World,
    /// East longitude, degrees (positive east on every world).
    pub lon_deg: f64,
    /// Planetodetic latitude on the IAU reference ellipsoid, degrees.
    pub lat_deg: f64,
    /// Height above the reference ellipsoid, km.
    pub height_km: f64,
}

impl Site {
    pub fn new(world: World, lon_deg: f64, lat_deg: f64) -> Site {
        Site { world, lon_deg, lat_deg, height_km: 0.0 }
    }

    /// Body-fixed position (km) and outward ellipsoid normal (unit).
    fn body_fixed(&self) -> (Vec3, Vec3) {
        let (a, c) = self.world.radii_km();
        let (sl, cl) = (self.lat_deg * DEG).sin_cos();
        let (so, co) = (self.lon_deg * DEG).sin_cos();
        let n = a * a / (a * a * cl * cl + c * c * sl * sl).sqrt();
        let h = self.height_km;
        let r = [(n + h) * cl * co, (n + h) * cl * so, (n * c * c / (a * a) + h) * sl];
        (r, [cl * co, cl * so, sl])
    }
}

/// Where the chart is cast from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpaceObserver {
    /// On the surface of a world: positions, houses and angles.
    Surface(Site),
    /// At the centre of a world: positions only.
    Center(World),
    /// Any point: `offset_km` from the centre of pack body `center` (NAIF id: 0 solar system
    /// barycentre, 10 Sun, 399 Earth, ...), ICRF axes, moving with that body plus
    /// `velocity_km_day`. Positions only.
    Point { center: i32, offset_km: Vec3, velocity_km_day: Vec3 },
}

impl SpaceObserver {
    /// The world the observer stands on or at, when there is one.
    pub fn world(&self) -> Option<World> {
        match self {
            SpaceObserver::Surface(s) => Some(s.world),
            SpaceObserver::Center(w) => Some(*w),
            SpaceObserver::Point { center, offset_km, .. } if *offset_km == [0.0; 3] => {
                World::ALL.into_iter().find(|w| w.naif() == *center)
            }
            SpaceObserver::Point { .. } => None,
        }
    }
}

/// Reference plane and zero point for longitudes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Zodiac {
    /// Earth's true ecliptic and equinox of date: the usual tropical zodiac.
    #[default]
    Earth,
    /// The observer's world: its orbital plane, 0° at its own vernal equinox.
    Native,
}

/// A body seen from the observer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpacePosition {
    /// Longitude and latitude in the chosen zodiac, degrees.
    pub lon: f64,
    pub lat: f64,
    /// Distance from the observer, AU.
    pub dist_au: f64,
    /// Rates per day (degrees, degrees, AU).
    pub lon_speed: f64,
    pub lat_speed: f64,
    pub dist_speed: f64,
}

/// Houses cast on another world.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpaceHouses {
    /// Cusps and angles in the chosen zodiac.
    pub houses: Houses,
    /// Sidereal angle of the local meridian on the house circle (degrees), the analogue of
    /// the ARMC: measured along the world's equator from its node on the zodiac plane.
    pub sidereal_angle: f64,
    /// Latitude relative to the spin axis, degrees (the IAU latitude with its sign flipped on
    /// retrograde rotators).
    pub spin_latitude: f64,
    /// Angle between the spin axis and the zodiac pole, degrees (the world's obliquity
    /// against that zodiac; above 90° for retrograde rotators).
    pub obliquity: f64,
}

/// Orientation of a world as the engine uses it at one epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorldFrame {
    /// Rotation from GCRS axes to the body-fixed frame.
    pub to_body: Mat3,
    /// North pole (IAU convention), GCRS unit vector.
    pub pole: Vec3,
    /// Rotation rate about `pole`, radians per day (negative for retrograde rotators).
    pub rate: f64,
}

impl WorldFrame {
    /// Spin axis: the pole about which the world turns counterclockwise.
    pub fn spin_axis(&self) -> Vec3 {
        if self.rate < 0.0 {
            scale(self.pole, -1.0)
        } else {
            self.pole
        }
    }
}

/// Horizontal coordinates of a body for a surface observer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Horizontal {
    /// Azimuth from the world's (IAU) north towards east, degrees.
    pub azimuth: f64,
    /// Elevation above the ellipsoid tangent plane, degrees, no refraction.
    pub elevation: f64,
}

fn unit(v: Vec3) -> Vec3 {
    scale(v, 1.0 / norm(v))
}

impl Engine {
    /// Body-fixed frame of `world` at `ep`. Earth uses the engine's own Earth rotation
    /// (Greenwich apparent sidereal time on the true equator of date, polar motion
    /// neglected), as for ordinary charts; the other worlds their IAU 2015 model.
    pub fn world_frame(&self, ep: &Epoch, world: World) -> WorldFrame {
        if world == World::Earth {
            let to_body = mat_mul(&rot_z(ep.gast()), &ep.frames.bpn);
            WorldFrame { to_body, pole: to_body[2], rate: EARTH_RATE_RAD_DAY }
        } else {
            let o = world.orientation(ep.jd_tdb);
            WorldFrame { to_body: o.icrf_to_body(), pole: o.pole(), rate: o.rotation_rate * DEG }
        }
    }

    /// Epoch for a UT instant with this engine's Delta T and frame model.
    pub fn space_epoch(&self, jd_ut: f64) -> Epoch {
        self.prepare_ut(jd_ut, false).epoch
    }

    /// Barycentric state of the observer, as the position pipeline uses it.
    pub fn space_observer_state(&self, ep: &Epoch, obs: &SpaceObserver) -> Result<ObsState> {
        let pack = self.pack();
        let (pos, vel) = match obs {
            SpaceObserver::Surface(site) => {
                let (c, cv) = pack.state(site.world.naif(), SSB, ep.jd_tdb)?;
                let f = self.world_frame(ep, site.world);
                let (r_bf, _) = site.body_fixed();
                let r = mat_t_vec(&f.to_body, r_bf);
                let v = scale(cross(f.pole, r), f.rate);
                (add(c, r), add(cv, v))
            }
            SpaceObserver::Center(w) => pack.state(w.naif(), SSB, ep.jd_tdb)?,
            SpaceObserver::Point { center, offset_km, velocity_km_day } => {
                let (c, cv) =
                    if *center == SSB { ([0.0; 3], [0.0; 3]) } else { pack.state(*center, SSB, ep.jd_tdb)? };
                (add(c, *offset_km), add(cv, *velocity_km_day))
            }
        };
        let sun = pack.position(SUN, SSB, ep.jd_tdb)?;
        Ok(ObsState { pos, vel, sun, earthbound: true, barycentric: false })
    }

    /// Rotation from GCRS axes to the zodiac frame (x at 0°, z at the zodiac pole).
    pub fn space_zodiac_matrix(&self, ep: &Epoch, world: Option<World>, zodiac: Zodiac) -> Result<Mat3> {
        match zodiac {
            Zodiac::Earth => Ok(ep.frames.ecl_true),
            Zodiac::Native => {
                let w = world.ok_or_else(|| {
                    Error::InvalidArgument("the native zodiac needs an observer on or at a world".into())
                })?;
                let (r, v) = self.pack().state(w.orbit_naif(), SUN, ep.jd_tdb)?;
                let k = unit(cross(r, v));
                let x = unit(cross(self.world_frame(ep, w).pole, k));
                Ok([x, cross(k, x), k])
            }
        }
    }

    /// GCRS direction vector (scaled to the light-time distance, km) of `body`.
    fn space_vector(&self, ep: &Epoch, obs: &SpaceObserver, body: Body, corr: Corrections) -> Result<Vec3> {
        let target = body
            .naif()
            .ok_or_else(|| Error::InvalidArgument(format!("body {} has no position off Earth", body.id())))?;
        if obs.world().map(|w| w.naif()) == Some(target) {
            return Err(Error::InvalidArgument("the observer is on or at this body".into()));
        }
        let state = self.space_observer_state(ep, obs)?;
        Ok(apparent(self.pack(), ep, &state, target, corr, LightTime::Iterate)?.0)
    }

    /// Position of `body` seen from `obs` at `jd_ut`, in `zodiac`, with numerical speeds.
    pub fn space_position(
        &self,
        body: Body,
        jd_ut: f64,
        obs: &SpaceObserver,
        zodiac: Zodiac,
        corr: Corrections,
    ) -> Result<SpacePosition> {
        let world = obs.world();
        let sample = |jd: f64| -> Result<(f64, f64, f64)> {
            let ep = self.space_epoch(jd);
            let v = self.space_vector(&ep, obs, body, corr)?;
            let m = self.space_zodiac_matrix(&ep, world, zodiac)?;
            let (lon, lat, r) = to_spherical(mat_vec(&m, v));
            Ok((lon / DEG, lat / DEG, r / AU_KM))
        };
        let (lon, lat, dist) = sample(jd_ut)?;
        let (l0, b0, r0) = sample(jd_ut - SPEED_STEP)?;
        let (l1, b1, r1) = sample(jd_ut + SPEED_STEP)?;
        let dl = crate::math::diff_deg(l1, l0);
        Ok(SpacePosition {
            lon,
            lat,
            dist_au: dist,
            lon_speed: dl / (2.0 * SPEED_STEP),
            lat_speed: (b1 - b0) / (2.0 * SPEED_STEP),
            dist_speed: (r1 - r0) / (2.0 * SPEED_STEP),
        })
    }

    /// Azimuth and elevation of `body` for a surface observer (apparent, airless).
    pub fn space_horizontal(&self, body: Body, jd_ut: f64, site: &Site) -> Result<Horizontal> {
        let ep = self.space_epoch(jd_ut);
        let obs = SpaceObserver::Surface(*site);
        let v = self.space_vector(&ep, &obs, body, Corrections::APPARENT)?;
        let m = self.world_frame(&ep, site.world).to_body;
        let u = unit(mat_vec(&m, v));
        let (_, up) = site.body_fixed();
        let east = unit(cross([0.0, 0.0, 1.0], up));
        let north = cross(up, east);
        let el = dot(u, up).clamp(-1.0, 1.0).asin() / DEG;
        let az = norm_deg(dot(u, east).atan2(dot(u, north)) / DEG);
        Ok(Horizontal { azimuth: az, elevation: el })
    }

    /// House cusps and angles for a site on a world. `system` as for [`houses::houses`];
    /// Sunshine ('I', 'i') and Gauquelin ('G') are not available off Earth.
    pub fn space_houses(&self, jd_ut: f64, site: &Site, system: char, zodiac: Zodiac) -> Result<SpaceHouses> {
        if matches!(system, 'I' | 'i' | 'G' | 'g') {
            return Err(Error::InvalidArgument(format!(
                "house system '{system}' is not available off Earth"
            )));
        }
        let ep = self.space_epoch(jd_ut);
        let f = self.world_frame(&ep, site.world);
        let (_, up_bf) = site.body_fixed();
        let up = mat_t_vec(&f.to_body, up_bf);
        let z = self.space_zodiac_matrix(&ep, Some(site.world), zodiac)?;
        let a = f.spin_axis();
        let obliquity = dot(a, z[2]).clamp(-1.0, 1.0).acos() / DEG;
        // The house routines take an obliquity below 90°: past it (a world spinning against
        // the zodiac's sense), work with the opposite zodiac pole and map the cusps back as
        // directions, so their numbering still follows the sky's daily motion.
        let k = if obliquity > 90.0 { scale(z[2], -1.0) } else { z[2] };
        let node = unit(cross(a, k));
        let eps = dot(a, k).clamp(-1.0, 1.0).acos() / DEG;
        let spin_latitude = dot(up, a).clamp(-1.0, 1.0).asin() / DEG;
        let sidereal_angle = norm_deg(dot(cross(a, node), up).atan2(dot(node, up)) / DEG);
        let h = houses::houses(system, sidereal_angle, spin_latitude, eps)?;
        let other = cross(k, node);
        let to_zodiac = |x: f64| {
            let (s, c) = (x * DEG).sin_cos();
            let u = mat_vec(&z, add(scale(node, c), scale(other, s)));
            norm_deg(u[1].atan2(u[0]) / DEG)
        };
        Ok(SpaceHouses {
            houses: Houses {
                cusps: h.cusps.map(to_zodiac),
                asc: to_zodiac(h.asc),
                mc: to_zodiac(h.mc),
                armc: h.armc,
                vertex: to_zodiac(h.vertex),
                equatorial_asc: to_zodiac(h.equatorial_asc),
                co_asc_koch: to_zodiac(h.co_asc_koch),
                co_asc_munkasey: to_zodiac(h.co_asc_munkasey),
                polar_asc: to_zodiac(h.polar_asc),
            },
            sidereal_angle,
            spin_latitude,
            obliquity,
        })
    }

    /// GCRS unit vector of a zodiac longitude on the zodiac circle (latitude 0).
    pub fn space_zodiac_direction(
        &self,
        ep: &Epoch,
        world: Option<World>,
        zodiac: Zodiac,
        lon_deg: f64,
    ) -> Result<Vec3> {
        let z = self.space_zodiac_matrix(ep, world, zodiac)?;
        let (s, c) = (lon_deg * DEG).sin_cos();
        Ok(mat_t_vec(&z, [c, s, 0.0]))
    }
}

/// Mars Sol Date (sols since 1873-12-29, Mars24 / Allison & McEwen 2000) for JD TT.
pub fn mars_sol_date(jd_tt: f64) -> f64 {
    (jd_tt - 2451549.5) / 1.0274912517 + 44796.0 - 0.0009626
}

/// JD TT for a Mars Sol Date.
pub fn jd_tt_from_mars_sol_date(msd: f64) -> f64 {
    (msd - 44796.0 + 0.0009626) * 1.0274912517 + 2451549.5
}

/// Coordinated Mars Time (mean solar time at the prime meridian), hours.
pub fn mars_coordinated_time(jd_tt: f64) -> f64 {
    24.0 * mars_sol_date(jd_tt).rem_euclid(1.0)
}

/// Local mean solar time on Mars at east longitude `lon_deg`, hours.
pub fn mars_local_mean_time(jd_tt: f64, lon_deg: f64) -> f64 {
    (mars_coordinated_time(jd_tt) + lon_deg / 15.0).rem_euclid(24.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn orientation_at(w: World) -> Orientation {
        w.orientation(J2000 + 9000.3)
    }

    #[test]
    fn body_frame_is_a_rotation_whose_z_is_the_pole() {
        for w in World::ALL {
            let o = orientation_at(w);
            let m = o.icrf_to_body();
            let p = o.pole();
            let z = mat_vec(&m, p);
            assert!((z[2] - 1.0).abs() < 1e-12, "{w:?}");
            let d = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
            assert!((d - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn mars_sol_date_round_trips() {
        for jd in [2451549.5, 2460000.25, 2400000.0] {
            assert!((jd_tt_from_mars_sol_date(mars_sol_date(jd)) - jd).abs() < 1e-9);
        }
        // Mars24 epoch: MSD 44796.0 at 2000 Jan 6 00:00 TT, less the 0.0009626 sol offset.
        assert!((mars_sol_date(2451549.5) - (44796.0 - 0.0009626)).abs() < 1e-12);
        assert!(
            (mars_local_mean_time(2451549.5, 180.0) - (mars_coordinated_time(2451549.5) + 12.0) % 24.0).abs()
                < 1e-9
        );
    }
}
