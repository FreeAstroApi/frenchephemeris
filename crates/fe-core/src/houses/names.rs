// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
//! House system names (`house_name`), as the classic API labels them.

/// Name of a house system code; empty for unknown or unsupported codes.
/// Case-insensitive except 'i' (Sunshine, Makransky's solution).
pub fn house_name(system: char) -> &'static str {
    let s = if system == 'i' { 'i' } else { system.to_ascii_uppercase() };
    match s {
        'A' | 'E' => "equal",
        'B' => "Alcabitius",
        'C' => "Campanus",
        'D' => "equal (MC)",
        'F' => "Carter poli-equ.",
        'G' => "Gauquelin sectors",
        'H' => "horizon/azimut",
        'I' => "Sunshine",
        'J' => "Savard-A",
        'i' => "Sunshine/alt.",
        'K' => "Koch",
        'L' => "Pullen SD",
        'M' => "Morinus",
        'N' => "equal/1=Aries",
        'O' => "Porphyry",
        'P' => "Placidus",
        'Q' => "Pullen SR",
        'R' => "Regiomontanus",
        'S' => "Sripati",
        'T' => "Polich/Page",
        'U' => "Krusinski-Pisa-Goelzer",
        'V' => "equal/Vehlow",
        'W' => "equal/ whole sign",
        'X' => "axial rotation system/Meridian houses",
        'Y' => "APC houses",
        _ => "",
    }
}
