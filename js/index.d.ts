// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI

export type BodyName =
    | "sun" | "moon" | "mercury" | "venus" | "mars" | "jupiter" | "saturn" | "uranus" | "neptune" | "pluto"
    | "meanNode" | "trueNode" | "meanLilith" | "trueLilith" | "earth"
    | "chiron" | "pholus" | "ceres" | "pallas" | "juno" | "vesta";

export declare const BODIES: Readonly<Record<BodyName, number>>;

export interface Position {
    /** Longitude (or right ascension with `equatorial`), degrees. */
    lon: number;
    /** Latitude (or declination), degrees. */
    lat: number;
    /** Distance, AU. */
    distance: number;
    lonSpeed: number;
    latSpeed: number;
    distanceSpeed: number;
    retrograde: boolean;
}

export interface CalcOptions {
    /** Ayanamsa name for sidereal positions: "lahiri", "raman", "krishnamurti", "fagan_bradley", ... */
    sidereal?: string;
    equatorial?: boolean;
    heliocentric?: boolean;
    /** Geometric position (no light time, deflection or aberration). */
    truePosition?: boolean;
    noNutation?: boolean;
    j2000?: boolean;
    /** Topocentric place: east longitude, latitude (degrees), height (m). */
    topo?: { lon: number; lat: number; height?: number };
}

export interface Houses {
    cusps: number[];
    asc: number;
    mc: number;
    armc: number;
    vertex: number;
    equatorialAsc: number;
    coAscKoch: number;
    coAscMunkasey: number;
    polarAsc: number;
}

export interface SignPosition {
    sign: string;
    degree: number;
}

export interface Chart {
    bodies: Record<string, Position & SignPosition>;
    houses: Houses;
}

export interface SpaceChart {
    bodies: Record<string, Position & SignPosition>;
    houses: { cusps: number[]; asc: number; mc: number; siderealAngle: number; obliquity: number };
}

export interface Site {
    world: "mercury" | "venus" | "earth" | "moon" | "mars" | "jupiter" | "saturn" | "uranus" | "neptune" | "pluto";
    /** East longitude, degrees. */
    lon: number;
    /** Planetographic latitude, degrees. */
    lat: number;
}

export declare class Ephemeris {
    calcUt(jdUt: number, body: BodyName | number, options?: CalcOptions): Position;
    houses(jdUt: number, lat: number, lon: number, system?: string, options?: { sidereal?: string }): Houses;
    ayanamsa(jdUt: number, mode?: string, trueEquinox?: boolean): number;
    fixstar(name: string, jdUt: number, options?: CalcOptions): Position & { name: string; magnitude: number };
    deltaT(jdUt: number): number;
    chart(jdUt: number, lat: number, lon: number, options?: { houseSystem?: string; sidereal?: string; bodies?: BodyName[] }): Chart;
    spaceChart(jdUt: number, site: Site, options?: { zodiac?: "earth" | "native"; houseSystem?: string }): SpaceChart;
    marsTime(jdUt: number, lon?: number): { marsSolDate: number; localMeanSolarTime: number };
    jdFromMarsTime(sol: number, hours: number, lon?: number): number;
    manifest(): string;
}

type Bytes = Uint8Array | ArrayBuffer | URL | string;

/** Load the engine with the bundled data (years 1800 to 2400), or with the packs given. */
export declare function load(options?: { wasm?: Bytes; mainPack?: Bytes; minorPack?: Bytes | null }): Promise<Ephemeris>;

export declare function julday(year: number, month: number, day: number, hour?: number, gregorian?: boolean): number;
export declare function julianDayFromDate(date: Date): number;
export declare function revjul(jd: number, gregorian?: boolean): { year: number; month: number; day: number; hour: number };
export declare function signOf(lon: number): SignPosition;
