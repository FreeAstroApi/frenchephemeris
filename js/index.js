// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//
// French Ephemeris for JavaScript: the Rust engine compiled to WebAssembly, with the bundled
// mobile data set (NASA JPL DE442 refit, years 1800 to 2400, about 2 MB).
// Works in Node.js 18+, Deno, Bun and browsers (ES module; bundlers that understand
// `new URL(..., import.meta.url)`, such as Vite and webpack 5, copy the files automatically).

import initWasm, { initSync, Engine, julday as rawJulday, revjul as rawRevjul, marsSolDate, jdTtFromMarsSolDate, marsLocalMeanTime } from "./pkg/engine.js";

const WASM_URL = new URL("./pkg/engine_bg.wasm", import.meta.url);
const MAIN_PACK_URL = new URL("./data/fepk-de442-mobile-1800-2400-v1.bin", import.meta.url);
const MINOR_PACK_URL = new URL("./data/fepk-minor-mobile-1800-2400-v1.bin", import.meta.url);

/** Classic body ids, as in pyswisseph and the Python package. */
export const BODIES = Object.freeze({
    sun: 0, moon: 1, mercury: 2, venus: 3, mars: 4, jupiter: 5, saturn: 6, uranus: 7, neptune: 8, pluto: 9,
    meanNode: 10, trueNode: 11, meanLilith: 12, trueLilith: 13, earth: 14,
    chiron: 15, pholus: 16, ceres: 17, pallas: 18, juno: 19, vesta: 20,
});

const SIGNS = ["Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo", "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces"];

const FLAG = { equatorial: 1, heliocentric: 2, truePosition: 4, noNutation: 8, j2000: 16 };

const isNode = typeof process !== "undefined" && !!process.versions?.node && typeof window === "undefined";

async function readBytes(source) {
    if (source instanceof Uint8Array) return source;
    if (source instanceof ArrayBuffer) return new Uint8Array(source);
    const url = source instanceof URL ? source : new URL(String(source), import.meta.url);
    if (isNode && url.protocol === "file:") {
        const { readFile } = await import("node:fs/promises");
        return new Uint8Array(await readFile(url));
    }
    const response = await fetch(url);
    if (!response.ok) throw new Error(`French Ephemeris: could not load ${url} (${response.status})`);
    return new Uint8Array(await response.arrayBuffer());
}

let wasmReady = null;

async function ensureWasm(wasm) {
    if (!wasmReady) {
        wasmReady = (async () => {
            const bytes = await readBytes(wasm ?? WASM_URL);
            initSync({ module: bytes });
        })();
    }
    return wasmReady;
}

function bodyId(body) {
    if (typeof body === "number") return body;
    const id = BODIES[body];
    if (id === undefined) throw new Error(`French Ephemeris: unknown body ${JSON.stringify(body)}`);
    return id;
}

function flags(options = {}) {
    let f = 0;
    for (const [k, bit] of Object.entries(FLAG)) if (options[k]) f |= bit;
    return f;
}

function position(a) {
    return { lon: a[0], lat: a[1], distance: a[2], lonSpeed: a[3], latSpeed: a[4], distanceSpeed: a[5], retrograde: a[3] < 0 };
}

/** Sign name and degree within the sign of an ecliptic longitude. */
export function signOf(lon) {
    const l = ((lon % 360) + 360) % 360;
    const i = Math.floor(l / 30);
    return { sign: SIGNS[i], degree: l - 30 * i };
}

/** Julian day (UT) for a calendar date and decimal hour; Gregorian unless `gregorian` is false. */
export function julday(year, month, day, hour = 0, gregorian = true) {
    return rawJulday(year, month, day, hour, gregorian);
}

/** Julian day (UT) for a JavaScript Date. */
export function julianDayFromDate(date) {
    return date.getTime() / 86400000 + 2440587.5;
}

/** Calendar date for a Julian day: { year, month, day, hour }. */
export function revjul(jd, gregorian = true) {
    const [year, month, day, hour] = rawRevjul(jd, gregorian);
    return { year, month, day, hour };
}

/** An ephemeris engine. Create one with `load()` and reuse it. */
export class Ephemeris {
    constructor(engine) {
        this.engine = engine;
    }

    /**
     * Position of a body at a Julian day (UT).
     * options: { sidereal: "lahiri" | ..., equatorial, heliocentric, truePosition, noNutation,
     *            j2000, topo: { lon, lat, height } }
     */
    calcUt(jdUt, body, options = {}) {
        const topo = options.topo ? [options.topo.lon, options.topo.lat, options.topo.height ?? 0] : undefined;
        return position(this.engine.calcUt(jdUt, bodyId(body), flags(options), options.sidereal, topo));
    }

    /** House cusps and angles for a place (east longitude, latitude); `system` such as "P" (Placidus) or "W". */
    houses(jdUt, lat, lon, system = "P", options = {}) {
        const a = this.engine.housesUt(jdUt, lat, lon, system, options.sidereal);
        return {
            cusps: Array.from(a.slice(0, 12)),
            asc: a[12], mc: a[13], armc: a[14], vertex: a[15],
            equatorialAsc: a[16], coAscKoch: a[17], coAscMunkasey: a[18], polarAsc: a[19],
        };
    }

    /** Ayanamsa in degrees (mean, as pyswisseph's get_ayanamsa_ut, unless `trueEquinox`). */
    ayanamsa(jdUt, mode = "lahiri", trueEquinox = false) {
        return this.engine.ayanamsaUt(jdUt, mode, trueEquinox);
    }

    /** Fixed star by name or designation ("Spica", "alf Vir", "HR 5056"): position, name and V magnitude. */
    fixstar(name, jdUt, options = {}) {
        return {
            name: this.engine.fixstarName(name),
            magnitude: this.engine.fixstarMag(name),
            ...position(this.engine.fixstarUt(name, jdUt, flags(options), options.sidereal)),
        };
    }

    /** Delta T (TT minus UT) in seconds. */
    deltaT(jdUt) {
        return this.engine.deltaT(jdUt);
    }

    /**
     * A whole chart: positions of the main bodies, houses and angles.
     * options: { houseSystem: "P", sidereal, bodies: [...names] }
     */
    chart(jdUt, lat, lon, options = {}) {
        const names = options.bodies ?? ["sun", "moon", "mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto", "trueNode", "chiron"];
        const bodies = {};
        for (const name of names) {
            const p = this.calcUt(jdUt, name, { sidereal: options.sidereal });
            bodies[name] = { ...p, ...signOf(p.lon) };
        }
        return { bodies, houses: this.houses(jdUt, lat, lon, options.houseSystem ?? "P", { sidereal: options.sidereal }) };
    }

    /**
     * A chart cast from a site on another world: Mars, the Moon, or any planet out to Pluto.
     * site: { world: "mars", lon (east), lat }; options: { zodiac: "earth" | "native", houseSystem }
     */
    spaceChart(jdUt, site, options = {}) {
        const zodiac = options.zodiac ?? "earth";
        const bodies = {};
        for (const name of ["sun", "moon", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto"]) {
            if (name === site.world) continue;
            const p = position(this.engine.spacePosition(jdUt, BODIES[name], site.world, site.lon, site.lat, zodiac));
            bodies[name] = { ...p, ...signOf(p.lon) };
        }
        const h = this.engine.spaceHouses(jdUt, site.world, site.lon, site.lat, options.houseSystem ?? "P", zodiac);
        return {
            bodies,
            houses: { cusps: Array.from(h.slice(0, 12)), asc: h[12], mc: h[13], siderealAngle: h[14], obliquity: h[16] },
        };
    }

    /** Mars Sol Date and local mean solar time (hours) at east longitude `lon`, for a Julian day UT. */
    marsTime(jdUt, lon = 0) {
        const jdTt = this.engine.utToTt(jdUt);
        return { marsSolDate: marsSolDate(jdTt), localMeanSolarTime: marsLocalMeanTime(jdTt, lon) };
    }

    /** Julian day UT for a Mars sol and local mean solar time (hours) at east longitude `lon`. */
    jdFromMarsTime(sol, hours, lon = 0) {
        return this.engine.ttToUt(jdTtFromMarsSolDate(sol + (hours - lon / 15) / 24));
    }

    /** Sources and coverage of the loaded data. */
    manifest() {
        return this.engine.manifest();
    }
}

/**
 * Load the engine and its data. With no argument, uses the files shipped in this package
 * (years 1800 to 2400). Pass bytes or URLs to use other packs:
 * { wasm, mainPack, minorPack } (Uint8Array, ArrayBuffer, URL or string).
 */
export async function load(options = {}) {
    await ensureWasm(options.wasm);
    const [main, minor] = await Promise.all([
        readBytes(options.mainPack ?? MAIN_PACK_URL),
        options.minorPack === null ? Promise.resolve(undefined) : readBytes(options.minorPack ?? MINOR_PACK_URL),
    ]);
    return new Ephemeris(new Engine(main, minor));
}

export { initWasm };
