/****  Pixel-wise SUHI + UTFVI seasonal/period exports (2003–2025)
 *
 *  OUTPUT per city:
 *   1) ONE multi-band SUHI GeoTIFF
 *   2) ONE multi-band UTFVI GeoTIFF
 *
 *  Bands (for BOTH products):
 *   Annual_Full_day, Annual_Full_night, Annual_Early_day, ..., Annual_Diff_night,
 *   Summer_* (DJF), Winter_* (JJA)
 *
 *  Pixel-wise SUHI definition (urban footprint):
 *     SUHI(x) = LST(x) - mean(LST_rural over city)
 *     then masked to URBAN pixels.
 *
 *  Pixel-wise UTFVI definition (urban footprint, monthly):
 *     UTFVI(x) = (LST(x) - mean(LST_urban over city)) / mean(LST_urban over city)
 *     then masked to URBAN pixels.
 *
 *  HOMOGENEITY FIX (EE ImageCollection mean requirement):
 *   - clamp to fixed range and cast to Float for EVERY monthly image
 *     so the collection is homogeneous for .mean().
 */

///////////////////////
// SETTINGS
///////////////////////
var START_YEAR = 2003;
var END_YEAR   = 2025;

var LST_MODE = "BOTH";   // "TERRA", "AQUA", "BOTH"
var MINPX = 10;          // minimum baseline pixels required (rural for SUHI, urban for UTFVI)

// SUHI clamp range (°C)
var SUHI_CLAMP_MIN = -30;
var SUHI_CLAMP_MAX =  30;

// UTFVI clamp range (dimensionless)
var UTFVI_CLAMP_MIN = -1.0;
var UTFVI_CLAMP_MAX =  1.0;
var UTFVI_EPS = 0.1;     // guard for |urbanMean| < EPS (°C) -> invalid month

var EXPORT_FOLDER = "GEE_exports";
var SUHI_EXPORT_PREFIX  = "SUHI_PIXEL_SEASONAL_MEANS";
var UTFVI_EXPORT_PREFIX = "UTFVI_PIXEL_SEASONAL_MEANS";

///////////////////////
// CITY ASSETS
///////////////////////
var cities = [
  {name: "Adelaide",  fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Adelaide")},
  {name: "Brisbane",  fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Brisbane")},
  {name: "Melbourne", fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Melbourne")},
  {name: "Perth",     fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Perth")},
  {name: "Sydney",    fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Sydney")}
];

///////////////////////
// LAND COVER (SUE masks)
///////////////////////
var lc = ee.ImageCollection("MODIS/061/MCD12Q1").select("LC_Type1");

// Urban=13, Water=17, Rural = not(urban) and not(water)
function sueMasks(year){
  year = ee.Number(year);

  var lcSorted = lc.sort("system:time_start");

  // Find most recent LC year <= requested year (fallback to first if not found)
  var lcBefore = lcSorted
    .filter(ee.Filter.calendarRange(2001, year, "year"))
    .sort("system:time_start", false)
    .first();

  var lcFallback = ee.Image(lcSorted.first());
  var lcYear = ee.Image(ee.Algorithms.If(lcBefore, lcBefore, lcFallback));

  var urban = lcYear.eq(13);
  var water = lcYear.eq(17);
  var rural = urban.not().and(water.not());

  return ee.Image.cat([urban.rename("urbanMask"), rural.rename("ruralMask")]);
}

///////////////////////
// LST with QC bits 0–1 <= 1
///////////////////////
function maskMOD11A2(img, dayOrNight){
  img = ee.Image(img);

  var qcBand = (dayOrNight === "Day") ? "QC_Day" : "QC_Night";
  var qc = img.select(qcBand);

  // bits 0-1: mandatory QA flags
  var quality = qc.bitwiseAnd(3);
  var good = quality.lte(1);

  var lstBand = (dayOrNight === "Day") ? "LST_Day_1km" : "LST_Night_1km";
  var lstC = img.select(lstBand).multiply(0.02).subtract(273.15);

  return lstC.rename(lstBand)
    .updateMask(good)
    .copyProperties(img, img.propertyNames());
}

var terra = ee.ImageCollection("MODIS/061/MOD11A2");
var aqua  = ee.ImageCollection("MODIS/061/MYD11A2");

function monthlyLSTComposite(year, month, dayOrNight){
  year = ee.Number(year);
  month = ee.Number(month);

  var start = ee.Date.fromYMD(year, month, 1);
  var end = start.advance(1, "month");

  function monthlyMean(col){
    var masked = col.filterDate(start, end).map(function(im){
      return maskMOD11A2(im, dayOrNight);
    });
    return masked.mean().set({
      "system:time_start": start.millis(),
      "year": year,
      "month": month
    });
  }

  if (LST_MODE === "TERRA") return monthlyMean(terra);
  if (LST_MODE === "AQUA")  return monthlyMean(aqua);

  var t = monthlyMean(terra);
  var a = monthlyMean(aqua);

  return ee.ImageCollection([t, a]).mean().set({
    "system:time_start": start.millis(),
    "year": year,
    "month": month
  });
}

///////////////////////
// BUILD YEAR-MONTH LIST
///////////////////////
function ymPairs(startYear, endYear){
  var years = ee.List.sequence(startYear, endYear);
  return years.map(function(y){
    y = ee.Number(y);
    return ee.List.sequence(1, 12).map(function(m){
      return ee.Dictionary({year: y, month: ee.Number(m)});
    });
  }).flatten();
}
var pairs = ymPairs(START_YEAR, END_YEAR);

///////////////////////
// SEASON + PERIOD FILTERS
///////////////////////
function seasonMonths(season){
  season = ee.String(season);
  return ee.List(ee.Algorithms.If(
    season.compareTo("Annual").eq(0), ee.List.sequence(1, 12),
    ee.Algorithms.If(season.compareTo("Summer").eq(0), ee.List([12, 1, 2]), ee.List([6, 7, 8])) // Winter = JJA
  ));
}

function periodFilter(period){
  period = ee.String(period);
  var fFull  = ee.Filter.rangeContains("year", START_YEAR, END_YEAR);
  var fEarly = ee.Filter.rangeContains("year", 2003, 2013);
  var fLate  = ee.Filter.rangeContains("year", 2014, 2025);

  return ee.Filter(ee.Algorithms.If(
    period.compareTo("Full").eq(0), fFull,
    ee.Algorithms.If(period.compareTo("Early").eq(0), fEarly, fLate)
  ));
}

///////////////////////
// PIXEL-WISE MONTHLY SUHI FOR A CITY
//  SUHI(x) = LST(x) - mean(LST_rural_over_city)
//  masked to URBAN pixels
//
// HOMOGENEITY FIX:
//  - clamp to fixed range and cast to Float for EVERY monthly image
///////////////////////
function monthlySUHI_pixel(cityGeom, year, month, dayOrNight){
  year = ee.Number(year);
  month = ee.Number(month);

  var lst = monthlyLSTComposite(year, month, dayOrNight);
  var lstBandName = ee.String(lst.bandNames().get(0)); // LST_Day_1km or LST_Night_1km

  var masks = sueMasks(year);
  var urbanMask = masks.select("urbanMask");
  var ruralMask = masks.select("ruralMask");

  // Rural baseline stats (mean + count)
  var ruralStats = lst.updateMask(ruralMask).reduceRegion({
    reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.count(), sharedInputs: true}),
    geometry: cityGeom,
    scale: 1000,
    bestEffort: true,
    maxPixels: 1e13
  });

  var ruralMean = ee.Number(ruralStats.get(lstBandName.cat("_mean")));
  var ruralCnt  = ee.Number(ruralStats.get(lstBandName.cat("_count")));

  var ok = ruralCnt.gte(MINPX);

  // Good month
  var suhiOK = lst
    .subtract(ruralMean)
    .rename("suhi")
    .updateMask(urbanMask)
    .clamp(SUHI_CLAMP_MIN, SUHI_CLAMP_MAX)
    .toFloat();

  // Bad month: fully masked float
  var suhiBad = ee.Image.constant(0).toFloat()
    .rename("suhi")
    .updateMask(ee.Image(0));

  return ee.Image(ee.Algorithms.If(ok, suhiOK, suhiBad))
    .set({
      "system:time_start": ee.Date.fromYMD(year, month, 1).millis(),
      "year": year,
      "month": month,
      "daynight": dayOrNight
    });
}

function monthlySUHI_collection(cityGeom, dayOrNight){
  return ee.ImageCollection(pairs.map(function(d){
    d = ee.Dictionary(d);
    var y = ee.Number(d.get("year"));
    var m = ee.Number(d.get("month"));
    return monthlySUHI_pixel(cityGeom, y, m, dayOrNight);
  }));
}

// Seasonal/period pixel-wise mean SUHI
function seasonalPeriodMeanSUHI(cityGeom, season, period, dayOrNight){
  var col = monthlySUHI_collection(cityGeom, dayOrNight);
  col = col.filter(ee.Filter.inList("month", seasonMonths(season)));
  col = col.filter(periodFilter(period));
  return col.mean().rename("suhi").toFloat();
}

///////////////////////
// PIXEL-WISE MONTHLY UTFVI FOR A CITY
//  UTFVI(x) = (LST(x) - mean(LST_urban_over_city)) / mean(LST_urban_over_city)
//  masked to URBAN pixels
//
// HOMOGENEITY FIX:
//  - clamp to fixed range and cast to Float for EVERY monthly image
///////////////////////
function monthlyUTFVI_pixel(cityGeom, year, month, dayOrNight){
  year = ee.Number(year);
  month = ee.Number(month);

  var lst = monthlyLSTComposite(year, month, dayOrNight);
  var lstBandName = ee.String(lst.bandNames().get(0));

  var masks = sueMasks(year);
  var urbanMask = masks.select("urbanMask");

  // Urban baseline stats (mean + count)
  var urbanStats = lst.updateMask(urbanMask).reduceRegion({
    reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.count(), sharedInputs: true}),
    geometry: cityGeom,
    scale: 1000,
    bestEffort: true,
    maxPixels: 1e13
  });

  var urbanMean = ee.Number(urbanStats.get(lstBandName.cat("_mean")));
  var urbanCnt  = ee.Number(urbanStats.get(lstBandName.cat("_count")));

  var okCnt  = urbanCnt.gte(MINPX);
  var okMean = urbanMean.abs().gte(UTFVI_EPS);
  var ok = okCnt.and(okMean);

  var utfviOK = lst
    .subtract(urbanMean)
    .divide(urbanMean)
    .rename("utfvi")
    .updateMask(urbanMask)
    .clamp(UTFVI_CLAMP_MIN, UTFVI_CLAMP_MAX)
    .toFloat();

  var utfviBad = ee.Image.constant(0).toFloat()
    .rename("utfvi")
    .updateMask(ee.Image(0));

  return ee.Image(ee.Algorithms.If(ok, utfviOK, utfviBad))
    .set({
      "system:time_start": ee.Date.fromYMD(year, month, 1).millis(),
      "year": year,
      "month": month,
      "daynight": dayOrNight
    });
}

function monthlyUTFVI_collection(cityGeom, dayOrNight){
  return ee.ImageCollection(pairs.map(function(d){
    d = ee.Dictionary(d);
    var y = ee.Number(d.get("year"));
    var m = ee.Number(d.get("month"));
    return monthlyUTFVI_pixel(cityGeom, y, m, dayOrNight);
  }));
}

// Seasonal/period pixel-wise mean UTFVI
function seasonalPeriodMeanUTFVI(cityGeom, season, period, dayOrNight){
  var col = monthlyUTFVI_collection(cityGeom, dayOrNight);
  col = col.filter(ee.Filter.inList("month", seasonMonths(season)));
  col = col.filter(periodFilter(period));
  return col.mean().rename("utfvi").toFloat();
}

///////////////////////
// BUILD MULTI-BAND EXPORT IMAGES PER CITY
///////////////////////
function addBand(out, img, name){
  return ee.Image(out).addBands(ee.Image(img).rename(name));
}

function citySUHI_SeasonalStack(cityObj){
  var name = cityObj.name;
  var geom = cityObj.geom;

  var seasons = ["Annual", "Summer", "Winter"];
  var periods = ["Full", "Early", "Late"];

  var out = ee.Image([]);

  seasons.forEach(function(season){

    periods.forEach(function(period){
      var dayImg = seasonalPeriodMeanSUHI(geom, season, period, "Day");
      var nigImg = seasonalPeriodMeanSUHI(geom, season, period, "Night");

      out = addBand(out, dayImg, season + "_" + period + "_day");
      out = addBand(out, nigImg, season + "_" + period + "_night");
    });

    var dayLate  = seasonalPeriodMeanSUHI(geom, season, "Late",  "Day");
    var dayEarly = seasonalPeriodMeanSUHI(geom, season, "Early", "Day");
    var nigLate  = seasonalPeriodMeanSUHI(geom, season, "Late",  "Night");
    var nigEarly = seasonalPeriodMeanSUHI(geom, season, "Early", "Night");

    out = addBand(out, dayLate.subtract(dayEarly).toFloat(), season + "_Diff_day");
    out = addBand(out, nigLate.subtract(nigEarly).toFloat(), season + "_Diff_night");
  });

  return out
    .toFloat()
    .clip(geom)
    .set({
      "city": name,
      "START_YEAR": START_YEAR,
      "END_YEAR": END_YEAR,
      "early_period": "2003-2013",
      "late_period": "2014-2025",
      "SUHI_def": "LST_pixel - mean(LST_rural), masked to urban",
      "SUHI_clamp": "[" + SUHI_CLAMP_MIN + "," + SUHI_CLAMP_MAX + "]",
      "LST_MODE": LST_MODE
    });
}

function cityUTFVI_SeasonalStack(cityObj){
  var name = cityObj.name;
  var geom = cityObj.geom;

  var seasons = ["Annual", "Summer", "Winter"];
  var periods = ["Full", "Early", "Late"];

  var out = ee.Image([]);

  seasons.forEach(function(season){

    periods.forEach(function(period){
      var dayImg = seasonalPeriodMeanUTFVI(geom, season, period, "Day");
      var nigImg = seasonalPeriodMeanUTFVI(geom, season, period, "Night");

      out = addBand(out, dayImg, season + "_" + period + "_day");
      out = addBand(out, nigImg, season + "_" + period + "_night");
    });

    var dayLate  = seasonalPeriodMeanUTFVI(geom, season, "Late",  "Day");
    var dayEarly = seasonalPeriodMeanUTFVI(geom, season, "Early", "Day");
    var nigLate  = seasonalPeriodMeanUTFVI(geom, season, "Late",  "Night");
    var nigEarly = seasonalPeriodMeanUTFVI(geom, season, "Early", "Night");

    out = addBand(out, dayLate.subtract(dayEarly).toFloat(), season + "_Diff_day");
    out = addBand(out, nigLate.subtract(nigEarly).toFloat(), season + "_Diff_night");
  });

  return out
    .toFloat()
    .clip(geom)
    .set({
      "city": name,
      "START_YEAR": START_YEAR,
      "END_YEAR": END_YEAR,
      "early_period": "2003-2013",
      "late_period": "2014-2025",
      "UTFVI_def": "(LST_pixel - mean(LST_urban)) / mean(LST_urban), masked to urban",
      "UTFVI_clamp": "[" + UTFVI_CLAMP_MIN + "," + UTFVI_CLAMP_MAX + "]",
      "LST_MODE": LST_MODE
    });
}

///////////////////////
// EXPORT (SUHI + UTFVI per city)
// NOTE: We DO NOT rely on cityGeoms to avoid scope issues.
// We build {name, geom} directly from `cities`.
///////////////////////
cities.forEach(function(obj){
  var c = { name: obj.name, geom: ee.FeatureCollection(obj.fc).geometry() };

  // SUHI stack
  var suhiStack = citySUHI_SeasonalStack(c);
  Export.image.toDrive({
    image: suhiStack,
    description: SUHI_EXPORT_PREFIX + "_" + c.name + "_2003_2025",
    folder: EXPORT_FOLDER,
    fileNamePrefix: SUHI_EXPORT_PREFIX + "_" + c.name + "_2003_2025",
    region: c.geom.bounds(),
    scale: 1000,
    maxPixels: 1e13
  });

  // UTFVI stack
  var utfviStack = cityUTFVI_SeasonalStack(c);
  Export.image.toDrive({
    image: utfviStack,
    description: UTFVI_EXPORT_PREFIX + "_" + c.name + "_2003_2025",
    folder: EXPORT_FOLDER,
    fileNamePrefix: UTFVI_EXPORT_PREFIX + "_" + c.name + "_2003_2025",
    region: c.geom.bounds(),
    scale: 1000,
    maxPixels: 1e13
  });
});

print("Configured SUHI + UTFVI exports. Run them from the Tasks tab.");
