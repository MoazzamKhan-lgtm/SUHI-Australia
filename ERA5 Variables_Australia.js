var START_YEAR = 2003;
var END_YEAR   = 2025;

var EXPORT_DESC   = "ERA5_5cities_MONTHLY_2003_2025";
var EXPORT_FOLDER = "GEE_exports";

var cities = [
  {name: "Adelaide",  fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Adelaide")},
  {name: "Brisbane",  fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Brisbane")},
  {name: "Melbourne", fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Melbourne")},
  {name: "Perth",     fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Perth")},
  {name: "Sydney",    fc: ee.FeatureCollection("projects/ee-farhanmoazzam/assets/Sydney")}
];
function cityGeom(obj){ return {name: obj.name, geom: ee.FeatureCollection(obj.fc).geometry()}; }
var cityGeoms = cities.map(cityGeom);

var era5 = ee.ImageCollection("ECMWF/ERA5/HOURLY");

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

function monthlyERA5Image(year, month){
  var start = ee.Date.fromYMD(year, month, 1);
  var end = start.advance(1, "month");
  var subset = era5.filterDate(start, end);

  var t2m = subset.select("temperature_2m").mean().subtract(273.15).rename("t2m_C");
  var d2m = subset.select("dewpoint_temperature_2m").mean().subtract(273.15).rename("d2m_C");
  var u10 = subset.select("u_component_of_wind_10m").mean();
  var v10 = subset.select("v_component_of_wind_10m").mean();
  var wind = u10.pow(2).add(v10.pow(2)).sqrt().rename("wind10_ms");
  var pr = subset.select("total_precipitation").sum().multiply(1000).rename("precip_mm");
  var cloud = subset.select("total_cloud_cover").mean().rename("cloud_frac");
  var blh = subset.select("boundary_layer_height").mean().rename("blh_m");

  var rh = ee.Image(100).multiply(
    d2m.expression("exp((17.625*Td)/(243.04+Td))", {Td: d2m}).divide(
      t2m.expression("exp((17.625*T)/(243.04+T))", {T: t2m})
    )
  ).rename("rh_pct");

  return ee.Image.cat([t2m, rh, wind, pr, cloud, blh])
    .set({"system:time_start": start.millis(), "year": year, "month": month});
}

function cityERA5MonthlyFC(cityObj){
  var name = cityObj.name, geom = cityObj.geom;

  return ee.FeatureCollection(pairs.map(function(d){
    d = ee.Dictionary(d);
    var year = ee.Number(d.get("year"));
    var month = ee.Number(d.get("month"));
    var start = ee.Date.fromYMD(year, month, 1);

    var img = monthlyERA5Image(year, month);
    var stats = img.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geom,
      scale: 10000,
      bestEffort: true,
      maxPixels: 1e13
    });

    var props = ee.Dictionary({
      city: name, year: year, month: month,
      date: start.format("YYYY-MM-dd")
    }).combine(stats, true);

    return ee.Feature(null, props).set({"system:time_start": start.millis()});
  }));
}

var eraAll = ee.FeatureCollection(cityGeoms.map(function(c){ return cityERA5MonthlyFC(c); })).flatten();
print("ERA5 preview", eraAll.limit(5));

Export.table.toDrive({
  collection: eraAll,
  description: EXPORT_DESC,
  folder: EXPORT_FOLDER,
  fileFormat: "CSV"
});
