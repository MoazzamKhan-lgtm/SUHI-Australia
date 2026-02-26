/****  SUHI monthly export + seasonal/period mean rasters (2003–2025)
 *
 *  - Monthly city table export (SUHI day/night + LST means + optional EVI)
 *  - Seasonal (Annual, Summer DJF, Winter JJA) mean SUHI day/night rasters per city
 *  - Early (2003–2013), Late (2014–2025), Diff (Late−Early) + Full (2003–2025)
 *  - Also exports one raster for ALL-CITIES mean (averaging city-month entries)
 *
 *  NOTE: seasonal “rasters” are constant-value images clipped to each city geometry
 */

// ======================
// USER SETTINGS
// ======================
var START_YEAR = 2003;
var END_YEAR   = 2025;
var LST_MODE = "BOTH";   // "TERRA", "AQUA", "BOTH"
var MINPX = 10;

var ADD_EVI  = true;
var ADD_ERA5 = false;   // IMPORTANT: turn OFF to avoid OOM (unused below)

var EXPORT_DESC   = "SUHI_SUE_5cities_MONTHLY_2003_2025_noERA5";
var EXPORT_FOLDER = "GEE_exports";

// seasonal/summary raster exports
var EXPORT_SUMMARY_FOLDER = EXPORT_FOLDER; // change if you want

// ======================
// CITY ASSETS
// ======================
var cities = [
  {name: "Adelaide",  fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Adelaide")},
  {name: "Brisbane",  fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Brisbane")},
  {name: "Melbourne", fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Melbourne")},
  {name: "Perth",     fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Perth")},
  {name: "Sydney",    fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Sydney")}
];

function cityGeom(obj){
  return {name: obj.name, geom: ee.FeatureCollection(obj.fc).geometry()};
}
var cityGeoms = cities.map(cityGeom);

// ======================
// LAND COVER (SUE masks)
// ======================
var lc = ee.ImageCollection("MODIS/061/MCD12Q1").select("LC_Type1");

function sueMasks(year){
  year = ee.Number(year);
  var lcSorted = lc.sort("system:time_start");
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

// ======================
// LST MASKING (QC bits 0–1 <= 1)
// ======================
function maskMOD11A2(img, dayOrNight){
  img = ee.Image(img);
  var qcBand = (dayOrNight === "Day") ? "QC_Day" : "QC_Night";
  var qc = img.select(qcBand);
  var quality = qc.bitwiseAnd(3); // bits 0-1
  var good = quality.lte(1);
  var lstBand = (dayOrNight === "Day") ? "LST_Day_1km" : "LST_Night_1km";
  var lstC = img.select(lstBand).multiply(0.02).subtract(273.15);
  return lstC.rename(lstBand).updateMask(good).copyProperties(img, img.propertyNames());
}

var terra = ee.ImageCollection("MODIS/061/MOD11A2");
var aqua  = ee.ImageCollection("MODIS/061/MYD11A2");

function monthlyLSTComposite(year, month, dayOrNight){
  year = ee.Number(year); month = ee.Number(month);
  var start = ee.Date.fromYMD(year, month, 1);
  var end = start.advance(1, "month");

  function monthlyMean(col){
    var masked = col.filterDate(start, end).map(function(im){
      return maskMOD11A2(im, dayOrNight);
    });
    return masked.mean().set({"system:time_start": start.millis(), "year": year, "month": month});
  }

  if (LST_MODE === "TERRA") return monthlyMean(terra);
  if (LST_MODE === "AQUA")  return monthlyMean(aqua);

  var t = monthlyMean(terra);
  var a = monthlyMean(aqua);
  return ee.ImageCollection([t, a]).mean()
    .set({"system:time_start": start.millis(), "year": year, "month": month});
}

// ======================
// EVI (optional)
// ======================
var vi = ee.ImageCollection("MODIS/061/MOD13A2").select("EVI");

function monthlyEVI(year, month){
  var start = ee.Date.fromYMD(year, month, 1);
  var end = start.advance(1, "month");
  return vi.filterDate(start, end).mean().multiply(0.0001).rename("EVI")
    .set({"system:time_start": start.millis(), "year": year, "month": month});
}

// ======================
// STATS HELPERS
// ======================
function meanCountOverMask(img, geom, maskImg, scale){
  return ee.Image(img).updateMask(maskImg).reduceRegion({
    reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.count(), sharedInputs: true}),
    geometry: geom,
    scale: scale,
    bestEffort: true,
    maxPixels: 1e13
  });
}

function safeSubtract(meanA, meanB, countA, countB, minpx){
  countA = ee.Number(countA);
  countB = ee.Number(countB);
  var ok = countA.gte(minpx).and(countB.gte(minpx));
  return ee.Algorithms.If(ok, ee.Number(meanA).subtract(ee.Number(meanB)), null);
}

// ======================
// BUILD YEAR-MONTH LIST
// ======================
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

// ======================
// MONTHLY FC PER CITY
// ======================
function cityMonthlyFC(cityObj){
  var name = cityObj.name, geom = cityObj.geom;

  return ee.FeatureCollection(pairs.map(function(d){
    d = ee.Dictionary(d);
    var year = ee.Number(d.get("year"));
    var month = ee.Number(d.get("month"));
    var date = ee.Date.fromYMD(year, month, 1);

    var masks = sueMasks(year);
    var urbanMask = masks.select("urbanMask");
    var ruralMask = masks.select("ruralMask");

    var lstDay = monthlyLSTComposite(year, month, "Day");
    var lstNig = monthlyLSTComposite(year, month, "Night");

    var dayU = meanCountOverMask(lstDay, geom, urbanMask, 1000);
    var dayR = meanCountOverMask(lstDay, geom, ruralMask, 1000);
    var nigU = meanCountOverMask(lstNig, geom, urbanMask, 1000);
    var nigR = meanCountOverMask(lstNig, geom, ruralMask, 1000);

    // Means (°C)
    var uDay = dayU.get("LST_Day_1km_mean");
    var rDay = dayR.get("LST_Day_1km_mean");
    var uNig = nigU.get("LST_Night_1km_mean");
    var rNig = nigR.get("LST_Night_1km_mean");

    // Counts
    var uDayN = dayU.get("LST_Day_1km_count");
    var rDayN = dayR.get("LST_Day_1km_count");
    var uNigN = nigU.get("LST_Night_1km_count");
    var rNigN = nigR.get("LST_Night_1km_count");

    // SUHI
    var suhiDay = safeSubtract(uDay, rDay, uDayN, rDayN, MINPX);
    var suhiNig = safeSubtract(uNig, rNig, uNigN, rNigN, MINPX);

    // Export props
    var props = ee.Dictionary({
      city: name,
      year: year,
      month: month,
      date: date.format("YYYY-MM-dd"),

      // SUHI
      suhi_day_C: suhiDay,
      suhi_night_C: suhiNig,

      // LST means
      lst_urban_day_C: uDay,
      lst_rural_day_C: rDay,
      lst_urban_night_C: uNig,
      lst_rural_night_C: rNig,

      // pixel counts
      n_urban_day: uDayN,
      n_rural_day: rDayN,
      n_urban_night: uNigN,
      n_rural_night: rNigN
    });

    if (ADD_EVI){
      var eviImg = monthlyEVI(year, month);
      var eU = meanCountOverMask(eviImg, geom, urbanMask, 1000);
      var eR = meanCountOverMask(eviImg, geom, ruralMask, 1000);

      var uE  = eU.get("EVI_mean");
      var rE  = eR.get("EVI_mean");
      var uEN = eU.get("EVI_count");
      var rEN = eR.get("EVI_count");

      var dEVI = safeSubtract(uE, rE, uEN, rEN, MINPX);

      props = props.combine(ee.Dictionary({
        dEVI: dEVI,
        evi_urban: uE,
        evi_rural: rE,
        n_evi_urban: uEN,
        n_evi_rural: rEN
      }), true);
    }

    return ee.Feature(null, props).set({"system:time_start": date.millis()});
  }));
}

// Combine all cities
var monthlyAll = ee.FeatureCollection(cityGeoms.map(function(c){
  return cityMonthlyFC(c);
})).flatten();

print("Monthly preview", monthlyAll.limit(5));

// ======================
// EXPORT MONTHLY TABLE
// ======================
Export.table.toDrive({
  collection: monthlyAll,
  description: EXPORT_DESC,
  folder: EXPORT_FOLDER,
  fileNamePrefix: EXPORT_DESC,
  fileFormat: "CSV"
});

// =====================================================================
// SEASONAL + PERIOD MEAN SUHI EXPORTS AS RASTERS (CITY + ALL-CITIES)
// =====================================================================

// ----- Season filter (Annual, Summer DJF, Winter JJA) -----
function filterSeason(fc, seasonName){
  seasonName = ee.String(seasonName);

  var annual = ee.Filter.inList("month", ee.List.sequence(1, 12));
  var summer = ee.Filter.inList("month", ee.List([12, 1, 2])); // DJF
  var winter = ee.Filter.inList("month", ee.List([6, 7, 8]));  // JJA

  return ee.FeatureCollection(ee.Algorithms.If(
    seasonName.compareTo("Annual").eq(0), ee.FeatureCollection(fc).filter(annual),
    ee.Algorithms.If(seasonName.compareTo("Summer").eq(0), ee.FeatureCollection(fc).filter(summer),
      ee.FeatureCollection(fc).filter(winter)
    )
  ));
}

// ----- Period filter (Full, Early, Late) -----
function filterPeriod(fc, periodName){
  periodName = ee.String(periodName);

  var full  = ee.FeatureCollection(fc).filter(ee.Filter.rangeContains("year", START_YEAR, END_YEAR));
  var early = ee.FeatureCollection(fc).filter(ee.Filter.rangeContains("year", 2003, 2013));
  var late  = ee.FeatureCollection(fc).filter(ee.Filter.rangeContains("year", 2014, 2025));

  return ee.FeatureCollection(ee.Algorithms.If(
    periodName.compareTo("Full").eq(0), full,
    ee.Algorithms.If(periodName.compareTo("Early").eq(0), early, late)
  ));
}

// Mean of a property over a FeatureCollection (skips nulls automatically)
function meanProp(fc, propName){
  return ee.Number(
    ee.FeatureCollection(fc).reduceColumns(ee.Reducer.mean(), [propName]).get("mean")
  );
}

// Mean SUHI for a given city/season/period
function meanSUHI_forCitySeasonPeriod(cityName, seasonName, periodName){
  var fc = monthlyAll.filter(ee.Filter.eq("city", cityName));
  fc = filterSeason(fc, seasonName);
  fc = filterPeriod(fc, periodName);

  var day   = meanProp(fc, "suhi_day_C");
  var night = meanProp(fc, "suhi_night_C");

  return ee.Dictionary({day: day, night: night});
}

// Mean SUHI across ALL cities (averaging city-month entries)
function meanSUHI_allCities(seasonName, periodName){
  var fc = filterSeason(monthlyAll, seasonName);
  fc = filterPeriod(fc, periodName);

  var day   = meanProp(fc, "suhi_day_C");
  var night = meanProp(fc, "suhi_night_C");

  return ee.Dictionary({day: day, night: night});
}

// Build two bands (day/night) from dictionary of {day,night}
function buildBandsFromDict(dict, prefix){
  dict = ee.Dictionary(dict);
  var day   = ee.Number(dict.get("day"));
  var night = ee.Number(dict.get("night"));
  var b1 = ee.Image.constant(day).rename(prefix + "_day");
  var b2 = ee.Image.constant(night).rename(prefix + "_night");
  return b1.addBands(b2);
}

// Build multi-band seasonal/period SUHI image for one city geometry
function citySummaryImage(cityObj){
  var cityName = cityObj.name;
  var geom = cityObj.geom;

  var seasons = ["Annual", "Summer", "Winter"]; // DJF, JJA
  var out = ee.Image([]);

  seasons.forEach(function(seasonName){
    var full  = meanSUHI_forCitySeasonPeriod(cityName, seasonName, "Full");
    var early = meanSUHI_forCitySeasonPeriod(cityName, seasonName, "Early");
    var late  = meanSUHI_forCitySeasonPeriod(cityName, seasonName, "Late");

    var diff = ee.Dictionary({
      day:   ee.Number(late.get("day")).subtract(ee.Number(early.get("day"))),
      night: ee.Number(late.get("night")).subtract(ee.Number(early.get("night")))
    });

    out = out
      .addBands(buildBandsFromDict(full,  seasonName + "_Full"))
      .addBands(buildBandsFromDict(early, seasonName + "_Early"))
      .addBands(buildBandsFromDict(late,  seasonName + "_Late"))
      .addBands(buildBandsFromDict(diff,  seasonName + "_Diff"));
  });

  return out.toFloat().clip(geom).set({
    city: cityName,
    START_YEAR: START_YEAR,
    END_YEAR: END_YEAR,
    early_period: "2003-2013",
    late_period: "2014-2025",
    note: "Bands are constant city-mean SUHI values"
  });
}

// Build multi-band seasonal/period SUHI image for all-cities mean
function allCitiesSummaryImage(){
  var unionGeom = ee.FeatureCollection(cityGeoms.map(function(c){
    return ee.Feature(c.geom);
  })).geometry().dissolve();

  var seasons = ["Annual", "Summer", "Winter"];
  var out = ee.Image([]);

  seasons.forEach(function(seasonName){
    var full  = meanSUHI_allCities(seasonName, "Full");
    var early = meanSUHI_allCities(seasonName, "Early");
    var late  = meanSUHI_allCities(seasonName, "Late");

    var diff = ee.Dictionary({
      day:   ee.Number(late.get("day")).subtract(ee.Number(early.get("day"))),
      night: ee.Number(late.get("night")).subtract(ee.Number(early.get("night")))
    });

    out = out
      .addBands(buildBandsFromDict(full,  seasonName + "_Full"))
      .addBands(buildBandsFromDict(early, seasonName + "_Early"))
      .addBands(buildBandsFromDict(late,  seasonName + "_Late"))
      .addBands(buildBandsFromDict(diff,  seasonName + "_Diff"));
  });

  return out.toFloat().clip(unionGeom).set({
    city: "ALL_CITIES_MEAN",
    START_YEAR: START_YEAR,
    END_YEAR: END_YEAR,
    early_period: "2003-2013",
    late_period: "2014-2025",
    note: "Bands are constant mean SUHI across cities (averaging city-month entries)"
  });
}

// ======================
// EXPORT CITY SUMMARY RASTERS (5 exports)
// ======================
cityGeoms.forEach(function(c){
  var img = citySummaryImage(c);

  Export.image.toDrive({
    image: img,
    description: "SUHI_SEASONAL_MEANS_" + c.name + "_2003_2025",
    folder: EXPORT_SUMMARY_FOLDER,
    fileNamePrefix: "SUHI_SEASONAL_MEANS_" + c.name + "_2003_2025",
    region: c.geom.bounds(),
    scale: 1000,
    maxPixels: 1e13
  });
});

// ======================
// EXPORT ALL-CITIES MEAN SUMMARY RASTER (1 export)
// ======================
var imgAll = allCitiesSummaryImage();
var allRegion = ee.FeatureCollection(cityGeoms.map(function(c){
  return ee.Feature(c.geom);
})).geometry().bounds();

Export.image.toDrive({
  image: imgAll,
  description: "SUHI_SEASONAL_MEANS_ALLCITIES_2003_2025",
  folder: EXPORT_SUMMARY_FOLDER,
  fileNamePrefix: "SUHI_SEASONAL_MEANS_ALLCITIES_2003_2025",
  region: allRegion,
  scale: 1000,
  maxPixels: 1e13
});

// ======================
// OPTIONAL: EXPORT TIDY SUMMARY TABLES (city-season + all-cities)
// ======================
function summaryRow(cityName, seasonName){
  var full  = meanSUHI_forCitySeasonPeriod(cityName, seasonName, "Full");
  var early = meanSUHI_forCitySeasonPeriod(cityName, seasonName, "Early");
  var late  = meanSUHI_forCitySeasonPeriod(cityName, seasonName, "Late");

  var diffD = ee.Number(late.get("day")).subtract(ee.Number(early.get("day")));
  var diffN = ee.Number(late.get("night")).subtract(ee.Number(early.get("night")));

  return ee.Feature(null, {
    city: cityName,
    season: seasonName,
    full_day:  full.get("day"),
    full_night: full.get("night"),
    early_day: early.get("day"),
    early_night: early.get("night"),
    late_day:  late.get("day"),
    late_night: late.get("night"),
    diff_day: diffD,
    diff_night: diffN
  });
}

var seasonsList = ee.List(["Annual","Summer","Winter"]);
var cityNamesList = ee.List(cities.map(function(o){ return o.name; }));

var citySeasonRows = ee.FeatureCollection(
  cityNamesList.map(function(cityName){
    cityName = ee.String(cityName);
    return seasonsList.map(function(seasonName){
      seasonName = ee.String(seasonName);
      return summaryRow(cityName, seasonName);
    });
  }).flatten()
);

Export.table.toDrive({
  collection: citySeasonRows,
  description: "SUHI_SEASONAL_CITY_SUMMARY_2003_2025",
  folder: EXPORT_SUMMARY_FOLDER,
  fileNamePrefix: "SUHI_SEASONAL_CITY_SUMMARY_2003_2025",
  fileFormat: "CSV"
});

function summaryRowAll(seasonName){
  seasonName = ee.String(seasonName);
  var full  = meanSUHI_allCities(seasonName, "Full");
  var early = meanSUHI_allCities(seasonName, "Early");
  var late  = meanSUHI_allCities(seasonName, "Late");

  var diffD = ee.Number(late.get("day")).subtract(ee.Number(early.get("day")));
  var diffN = ee.Number(late.get("night")).subtract(ee.Number(early.get("night")));

  return ee.Feature(null, {
    city: "ALL_CITIES_MEAN",
    season: seasonName,
    full_day:  full.get("day"),
    full_night: full.get("night"),
    early_day: early.get("day"),
    early_night: early.get("night"),
    late_day:  late.get("day"),
    late_night: late.get("night"),
    diff_day: diffD,
    diff_night: diffN
  });
}

var allRows = ee.FeatureCollection(seasonsList.map(summaryRowAll));

Export.table.toDrive({
  collection: allRows,
  description: "SUHI_SEASONAL_ALLCITIES_SUMMARY_2003_2025",
  folder: EXPORT_SUMMARY_FOLDER,
  fileNamePrefix: "SUHI_SEASONAL_ALLCITIES_SUMMARY_2003_2025",
  fileFormat: "CSV"
});F