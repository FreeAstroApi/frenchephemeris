# French Ephemeris for JavaScript

The French Ephemeris engine compiled to WebAssembly: planetary positions from NASA JPL data,
24 house systems, 46 ayanamsas, 8,969 fixed stars and charts cast from Mars, in Node.js,
Deno, Bun and the browser. Same engine and same numbers as the Python package
(`pip install frenchephemeris`).

```bash
npm install frenchephemeris
```

```js
import { load, julday, signOf } from "frenchephemeris";

const eph = await load();                       // engine + bundled data, about 4.5 MB
const jd = julday(1990, 5, 17, 14.5);           // 17 May 1990, 14:30 UT

eph.calcUt(jd, "moon");                         // { lon, lat, distance, lonSpeed, ..., retrograde }
eph.calcUt(jd, "sun", { sidereal: "lahiri" });  // sidereal zodiac
eph.houses(jd, 48.8566, 2.3522, "P");           // Paris, Placidus: { cusps, asc, mc, ... }
eph.fixstar("Spica", jd);                       // { name: "Spica,alf Vir", magnitude, lon, ... }

const chart = eph.chart(jd, 48.8566, 2.3522);   // every body with its sign, and the houses
signOf(chart.houses.asc);                       // { sign: "Libra", degree: 3.56 }

eph.spaceChart(jd, { world: "mars", lon: 77.45, lat: 18.44 });   // the same moment, from Mars
```

## Data

The package ships the **mobile** data set: a compact refit of NASA JPL DE442 for years
1800 to 2400 (planets within 0.007″ of JPL, Moon within 0.4″), with Chiron, Pholus, Ceres,
Pallas, Juno and Vesta. Dates outside that range throw an error. Other packs (JPL's exact numbers,
years −13198 to 17190) can be passed to `load({ mainPack, minorPack })` as bytes or URLs; see
the [data packs](https://github.com/FreeAstroApi/frenchephemeris#data-packs).

## In the browser

`load()` fetches `engine_bg.wasm` and the data next to the module with
`new URL(..., import.meta.url)`, which Vite and webpack 5 handle. You can also host the files
yourself and pass their URLs: `load({ wasm: "/fe/engine_bg.wasm", mainPack: "/fe/main.bin", minorPack: "/fe/minor.bin" })`.

## Licence

GNU AGPL-3.0, or a commercial licence (USD 200 once per company) for closed source products:
<https://www.freeastroapi.com/frenchephemeris#licensing>. Data: NASA/JPL; see NOTICE.
