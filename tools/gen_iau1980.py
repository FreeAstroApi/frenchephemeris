# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
# Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
"""Generate the IAU 1980 nutation table and the pinned IERS celestial-pole-offset series.

1. crates/fe-core/src/frames/iau1980/nut80_table.rs — the 106-term IAU 1980 nutation series
   (Wahr 1981; Seidelmann 1982; Explanatory Supplement 1992, Table 3.222.1), transcribed
   from ERFA's nut80.c (BSD-3, (c) NumFOCUS, derived from SOFA):
   https://raw.githubusercontent.com/liberfa/erfa/master/src/nut80.c
2. crates/fe-core/src/frames/eop_c04_iau1980.bin — daily dPsi, dEps of the IERS EOP 14 C04
   series referred to IAU 1976/1980 precession-nutation, from
   https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU1980_one_file_1962-now.txt
   Format: little-endian i32 pairs (dPsi, dEps) in micro-arcseconds, one per day from the
   first MJD (printed in frames/iau1980.rs, EOP_FIRST_MJD).

  .venv/bin/python -I tools/gen_iau1980.py <nut80.c> <eopc04_IAU1980.txt>
"""
import hashlib, os, re, struct, sys

here = os.path.dirname(os.path.abspath(__file__))
out_dir = os.path.join(here, '..', 'crates', 'fe-core', 'src', 'frames')
nut_c, eop = sys.argv[1], sys.argv[2]

src = open(nut_c).read()
rows = re.findall(r'\{\s*(-?\d+),\s*(-?\d+),\s*(-?\d+),\s*(-?\d+),\s*(-?\d+),\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+)\s*\}', src)
assert len(rows) == 106, len(rows)
lines = ['// SPDX-License-Identifier: (AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial) AND BSD-3-Clause',
         '// Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI',
         '// Data: IAU 1980 nutation series transcribed from ERFA nut80.c (BSD-3-Clause, Copyright (C) 2013-2021 NumFOCUS Foundation; derived from IAU SOFA). See NOTICE.',
         '//! IAU 1980 nutation series (Wahr 1981; Seidelmann 1982, Celest. Mech. 27, 79;',
         '//! Explanatory Supplement to the Astronomical Almanac 1992, Table 3.222.1).',
         '//! Transcribed by tools/gen_iau1980.py from ERFA nut80.c (BSD-3, (c) 2013-2023 NumFOCUS',
         '//! Foundation, derived with permission from the IAU SOFA library). Do not edit.',
         '//! Multipliers of (l, l\', F, D, Ω), then Δψ sine and Δε cosine coefficients and',
         '//! their rates, units 0.1 mas and 0.1 mas per Julian century.', '',
         'pub(crate) const NUT80: [([i8; 5], [f64; 4]); 106] = [']
for r in rows:
    m = ', '.join(r[:5]); c = ', '.join(f'{float(v):.1f}' for v in r[5:])
    lines.append(f'    ([{m}], [{c}]),')
lines.append('];')
open(os.path.join(out_dir, 'iau1980', 'nut80_table.rs'), 'w').write('\n'.join(lines) + '\n')

data = []
first = None
for line in open(eop):
    p = line.split()
    if len(p) >= 10 and p[0].isdigit() and len(p[0]) == 4:
        mjd = int(p[3])
        if first is None:
            first = mjd
        assert mjd == first + len(data), (mjd, first, len(data))
        data.append((round(float(p[8]) * 1e6), round(float(p[9]) * 1e6)))
blob = b''.join(struct.pack('<ii', a, b) for a, b in data)
open(os.path.join(out_dir, 'eop_c04_iau1980.bin'), 'wb').write(blob)
print('nut80 rows', len(rows), 'eop days', len(data), 'first MJD', first, 'last MJD', first + len(data) - 1,
      'source sha256', hashlib.sha256(open(eop, 'rb').read()).hexdigest())
