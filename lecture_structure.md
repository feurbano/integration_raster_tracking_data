# Concept of the lecture

Introduce the EURODEER project as a study case where a spatial database is the core element and the engine of scientific activities in wildlife tracking. Illustrate technical and not that technical lessons learnt from the development and management of this large wildlife tracking database. Show 2 cases of integration of tracking data (vector) and environmental information (raster) as an additional demonstration of the potentiality of the PostgreSQL/PostGIS platform.


## Introduction to Eurodeer (20 minutes)

## DEMONSTRATION 1: Analyzing movement data with a (raster) environmental layer (20 minutes, see the code in a separate document)
* General introductions
* Visualize all animals
* Load a background
* Create trajectory
* Visualize trajectory
* Create a Convex Hull
* Visualize Convex Hull
* Limit to a specific animal
* Import land cover
* Visualize land cover
* Intersect points with land cover for a single animal
* Show stats for points
* Intersect convex hull with land cover for a single animal
* Show stats for convex hull
* Compare stats points/convex hull
* Intersect monthly convex hulls of an animal with land cover
* Intersect points for female/male with land cover

## DEMONSTRATION 2: Analyzing location data with a time series of environmental layers (10 minutes, see the code in a separate document)
* General introductions
* Load and timestamp a raster (NDVI) time series
* Extraction of a NDVI value for a point/time
* Extraction of a NDVI time series of values of a given fix
* Extraction of the NDVI value for a fix as temporal interpolation of the 2 closest images
* Extraction of the NDVI values for a set of fixes as temporal interpolation of the 2 closest images for animal 782
* Plot raster time series stored in PostgreSQL/PostGIS from R
* Future directions: Google Earth Engine

## Space for discussion and questions - (10 minutes)
