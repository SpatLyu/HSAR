# Dataset of properties in the municipality of Athens

A dataset of apartments in the municipality of Athens for 2017. Point
location of the properties is given together with their main
characteristics and the distance to the closest metro/train station.

## Usage

``` r
properties
```

## Format

An object of class `sf` (inherits from `data.frame`) with 1000 rows and
7 columns.

## Details

An sf object of 1000 points with the following 6 variables.

- id:

  An unique identifier for each property.

- size:

  The size of the property (unit: square meters)

- price:

  The asking price (unit: euros)

- prpsqm:

  The asking price per squre meter (unit: euroes/square meter).

- age:

  Age of property in 2017 (unit: years).

- dist_metro:

  The distance to closest train/metro station (unit: meters).
