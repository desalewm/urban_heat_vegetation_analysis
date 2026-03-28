
// *******************************************************************************************************
//  Remote sensing of urban heat dynamics and the cooling effect of urban green spaces in Ethiopian cities
//  ******************************************************************************************************

// =====================================================
// 1️ STUDY AREA & GLOBAL SETTINGS
// =====================================================

// Define study cities
var cities = {
  "Adama": ee.FeatureCollection("projects/desalewmoges/assets/adama"),
  "Addis": ee.FeatureCollection("projects/desalewmoges/assets/addis"),
  "Harar": ee.FeatureCollection("projects/desalewmoges/assets/harar"),
  "Jimma": ee.FeatureCollection("projects/desalewmoges/assets/jimma")
};

// Analysis years
var years = [2021, 2022, 2023, 2024];

// Center map
Map.centerObject(cities["Addis"], 10);

// =====================================================
// 2️ PREPROCESSING FUNCTIONS
// =====================================================

// Sentinel-2 cloud masking
function maskS2(image) {
  var qa = image.select('QA60');
  var mask = qa.bitwiseAnd(1 << 10).eq(0)
    .and(qa.bitwiseAnd(1 << 11).eq(0));
  return image.updateMask(mask).divide(10000)
    .copyProperties(image, image.propertyNames());
}

// Landsat cloud masking
function maskL8(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)
    .and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask);
}

// Landsat scaling
function scaleL8(image) {
  var optical = image.select('SR_B.*').multiply(0.0000275).add(-0.2);
  var thermal = image.select('ST_B10').multiply(0.00341802).add(149.0);
  return image.addBands(optical, null, true)
              .addBands(thermal, null, true);
}

// =====================================================
// 3️ NDVI WORKFLOW (SEASONAL, MONTHLY, TIMESERIES)
// =====================================================

// NDVI calculation
function addNDVI(image, geom) {
  return image.addBands(
    image.normalizedDifference(['B8', 'B4']).rename('NDVI')
  ).clip(geom);
}

// Sentinel-2 collection (Jan–May)
function getS2(geom, start, end) {
  return ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(geom)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .filter(ee.Filter.calendarRange(1, 5, 'month'))
    .map(maskS2);
}

// ---- Seasonal NDVI ----
function seasonalNDVI(city, name, year) {
  var geom = city.geometry();
  var col = getS2(geom, year + '-01-01', year + '-12-31')
    .map(function(img){ return addNDVI(img, geom); });

  var mean = col.mean().select('NDVI');

  Map.addLayer(mean, {min:0,max:1,palette:['white','green']}, name + ' NDVI ' + year);

  Export.image.toDrive({
    image: mean,
    description: name + '_NDVI_' + year,
    region: geom,
    scale: 10,
    maxPixels: 1e13
  });
}

// ---- Monthly NDVI ----
function monthlyNDVI(city, name) {
  var geom = city.geometry();
  var col = getS2(geom, '2021-01-01', '2024-12-31')
    .map(function(img){ return addNDVI(img, geom).select('NDVI'); });

  for (var m = 1; m <= 5; m++) {
    var monthly = col.filter(ee.Filter.calendarRange(m, m, 'month')).mean();

    Export.image.toDrive({
      image: monthly,
      description: name + '_NDVI_Month_' + m,
      region: geom,
      scale: 10,
      maxPixels: 1e13
    });
  }
}

// ---- NDVI Time Series ----
function ndviTimeSeries(city, name) {
  var geom = city.geometry();

  var col = ee.ImageCollection("NASA/HLS/HLSS30/v002")
    .filterBounds(geom)
    .filterDate('2021-01-01', '2024-12-31')
    .filter(ee.Filter.calendarRange(2,5,'month'))
    .map(function(img){
      return img.normalizedDifference(['B8','B4'])
        .rename('NDVI')
        .copyProperties(img, ['system:time_start']);
    });

  var chart = ui.Chart.image.series({
    imageCollection: col,
    region: geom,
    reducer: ee.Reducer.mean(),
    scale: 100
  });

  print(name + ' NDVI Time Series', chart);
}

// =====================================================
// 4️ LST WORKFLOW (LANDSAT)
// =====================================================

// LST calculation
function computeLST(image) {
  var ndvi = image.normalizedDifference(['SR_B5','SR_B4']);
  var fv = ndvi.subtract(ndvi.reduceRegion({
    reducer: ee.Reducer.min(),
    geometry: image.geometry(),
    scale: 30
  }).values().get(0))
  .pow(2);

  var emissivity = fv.multiply(0.004).add(0.986);

  return image.expression(
    '(Tb / (1 + (0.00115 * (Tb / 1.438)) * log(em))) - 273.15',
    {Tb: image.select('ST_B10'), em: emissivity}
  ).rename('LST');
}

// Landsat collection
function getL8(geom, year) {
  return ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
    .filterBounds(geom)
    .filterDate(year+'-01-01', year+'-12-31')
    .filter(ee.Filter.calendarRange(1,5,'month'))
    .filter(ee.Filter.lt('CLOUD_COVER', 30))
    .map(maskL8)
    .map(scaleL8);
}

// Seasonal LST
function seasonalLST(city, name, year) {
  var geom = city.geometry();

  var col = getL8(geom, year).map(computeLST);
  var mean = col.mean();

  Map.addLayer(mean, {min:20,max:40,palette:['blue','red']}, name + ' LST ' + year);

  Export.image.toDrive({
    image: mean,
    description: name + '_LST_' + year,
    region: geom,
    scale: 30,
    maxPixels: 1e13
  });
}

// =====================================================
// 5️ MODIS LST (DAY/NIGHT TIMESERIES & MAPS)
// =====================================================

function getMODIS() {
  return ee.ImageCollection("MODIS/061/MOD11A1")
    .filterDate('2021-01-01','2024-12-31')
    .filter(ee.Filter.calendarRange(2,5,'month'))
    .map(function(img){
      return ee.Image.cat(
        img.select('LST_Day_1km').multiply(0.02).subtract(273.15),
        img.select('LST_Night_1km').multiply(0.02).subtract(273.15)
      ).rename(['Day','Night'])
      .copyProperties(img, ['system:time_start']);
    });
}

// Export MODIS time series
function modisExport(city, name) {
  var geom = city.geometry();
  var col = getMODIS();

  var table = col.map(function(img){
    return ee.Feature(null, {
      date: ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'),
      day: img.select('Day').reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geom,
        scale: 1000
      }).get('Day')
    });
  });

  Export.table.toDrive({
    collection: table,
    description: name + '_MODIS_LST',
    fileFormat: 'CSV'
  });
}

// =====================================================
// 6️ EXECUTION PIPELINE
// =====================================================

Object.keys(cities).forEach(function(name){
  var city = cities[name];

  years.forEach(function(year){
    seasonalNDVI(city, name, year);
    seasonalLST(city, name, year);
  });

  monthlyNDVI(city, name);
  ndviTimeSeries(city, name);
  modisExport(city, name);
});