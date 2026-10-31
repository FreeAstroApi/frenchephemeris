// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
// The WebAssembly build must give the Python package's numbers (same engine, same packs).
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

import { BODIES, julday, julianDayFromDate, load, revjul, signOf } from "../index.js";

const ref = JSON.parse(await readFile(new URL("./reference.json", import.meta.url), "utf8"));
const eph = await load();

// WebAssembly and native code round the last bits differently (speeds are numerical
// derivatives): compare to 1e-9 (degrees, about 4 µas) plus 1e-12 relative.
const close = (a, b, tol, what) => {
    for (let i = 0; i < b.length; i++) {
        const t = Math.max(tol, 1e-9) + 1e-12 * Math.abs(b[i]);
        assert.ok(Math.abs(a[i] - b[i]) <= t, `${what}[${i}]: ${a[i]} vs ${b[i]}`);
    }
};
const vec = (p) => [p.lon, p.lat, p.distance, p.lonSpeed, p.latSpeed, p.distanceSpeed];

test("julian day and calendar", () => {
    assert.equal(julday(1990, 5, 17, 14.5), ref.jd);
    const d = revjul(ref.jd);
    assert.deepEqual([d.year, d.month, d.day], [1990, 5, 17]);
    assert.ok(Math.abs(d.hour - 14.5) < 1e-6);
    assert.equal(julianDayFromDate(new Date(Date.UTC(2000, 0, 1, 12))), 2451545.0);
});

test("positions match the Python package", () => {
    for (const [name, want] of Object.entries(ref.bodies)) close(vec(eph.calcUt(ref.jd, name)), want, 1e-12, name);
    for (const [name, want] of Object.entries(ref.sidereal)) {
        close(vec(eph.calcUt(ref.jd, name, { sidereal: "lahiri" })), want, 1e-12, `${name} sidereal`);
    }
    close(vec(eph.calcUt(ref.jd, "mars", { equatorial: true })), ref.equatorial, 1e-12, "mars equatorial");
    close(vec(eph.calcUt(ref.jd, "moon", { topo: { lon: 2.3522, lat: 48.8566, height: 35 } })), ref.topo_moon, 1e-12, "topocentric moon");
    assert.equal(eph.calcUt(ref.jd, BODIES.sun).lon, eph.calcUt(ref.jd, "sun").lon);
});

test("houses, ayanamsa, Delta T and stars match", () => {
    const h = eph.houses(ref.jd, 48.8566, 2.3522, "P");
    close(h.cusps, ref.houses.cusps, 1e-12, "cusps");
    close([h.asc, h.mc, h.armc, h.vertex], ref.houses.ascmc.slice(0, 4), 1e-12, "angles");
    assert.ok(Math.abs(eph.ayanamsa(ref.jd, "lahiri") - ref.ayanamsa) < 1e-9);
    assert.ok(Math.abs(eph.deltaT(ref.jd) - ref.deltat_s) < 1e-9);
    const spica = eph.fixstar("Spica", ref.jd);
    // The distance rate of a star 15.8 million AU away is cancellation noise: not compared.
    close(vec(spica).slice(0, 5), ref.spica.slice(0, 5), 1e-12, "Spica");
    assert.equal(spica.name, "Spica,alf Vir");
    assert.equal(eph.fixstar("HR 5056", ref.jd).name, "Spica,alf Vir");
});

test("charts", () => {
    const c = eph.chart(ref.jd, 48.8566, 2.3522);
    assert.equal(c.bodies.sun.sign, "Taurus");
    assert.equal(c.bodies.moon.sign, "Aquarius");
    assert.equal(signOf(c.houses.asc).sign, "Libra");
    const m = eph.spaceChart(ref.jd, { world: "mars", lon: 77.45, lat: 18.44 });
    assert.ok(!("mars" in m.bodies) && "earth" in m.bodies);
    assert.ok(Math.abs(m.bodies.earth.lon - ref.mars_earth.lon) < 1e-9);
    close(m.houses.cusps, ref.mars_houses.cusps, 1e-9, "Mars cusps");
    const t = eph.marsTime(ref.jd, 77.45);
    const back = eph.jdFromMarsTime(Math.floor(t.marsSolDate + 77.45 / 360), t.localMeanSolarTime, 77.45);
    assert.ok(Math.abs(back - ref.jd) < 1e-8);
});

test("errors are thrown, not returned", () => {
    assert.throws(() => eph.calcUt(ref.jd, "vulcan"), /unknown body/);
    assert.throws(() => eph.calcUt(ref.jd, "sun", { sidereal: "nope" }), /ayanamsa/);
    assert.throws(() => eph.fixstar("Not A Star", ref.jd), /unknown fixed star/);
    assert.throws(() => eph.calcUt(julday(2600, 1, 1), "sun"));
});
