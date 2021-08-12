[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5183824.svg)](https://doi.org/10.5281/zenodo.5183824)

# junco
files for analysis and visualization of data on junco bone microstructure

requires a single spreadsheet of specimen and bone morphology data
stored on Dryad; DOI to follow

exports 8 spreadsheets and 5 figures to the directory indicated by the user
these files match the tables and figures in the paper; link to follow

Junco hyemalis bone microstructure raw data table
Written 2021-06-22 by Leeann Louis

## General notes:

  * All specimens are included twice in the spreadsheet: once for data on the humerus, and once for data on the femur.
  * There are some missing values for elevation of collection (11 out of 87 specimens).

## Full definitions for all data fields:

* Museum: Institution from which the specimen was loaned. Here it is always the American Museum of Natural History (AMNH) in New York, NY, US
* SpecNum: Specimen number of the skeleteon in the museum
* Name: Full species name of the specimen (Genus species subspecies)
* Subspp: Subspecies name
* Sex: Specimen sex
* Location: Details on where the specimen was collected
* Elev: Altitude (in m) at which the specimen was collected, if available.
* Date: Month, day, and year at which the specimen was collected
* Day: Julian date: the integer day of the year where 1 corresponds to January 1st and 365 corresponds to December 31st (in a non-leap year)
* Mass: Weight at collection in grams
* Side: Right or left wing / leg from which the bone was taken
* Bone: Femur or humerus
* CtDOI: DOI number for cortical computed tomography data available at [MorphoSource](https://www.morphosource.org/)
* TbDOI: DOI number for trabecular computed tomography data available at [MorphoSource](https://www.morphosource.org/)
* Length: Bone length in mm, taken during the scout-view
* TbBV: Trabecular bone volume (mm^3): Volume of mineralized trabecular bone tissue within the total volume analyzed
* BVTV: Trabecular bone volume fraction (%): Ratio of trabecular bone divided by the total volume of tissue analyzed
* ConnDens: Connectivity density (mm^-3): The degree of connectivity between trabecular, normalized by the total volume assessed
* SMI: Structure model index (SMI): Trabecular structure assessed from 0 to 3, where 0 is plate-like and 3 is rod-like
* TbN: Trabecular number (mm^-1): Average number of trabeculae per unit length. Calculated as the inverse of the mean distance between the mid-axes of the trabeculae
* TbTh: Trabecular thickness (mm): Average thickness of a trabecula. Determined by fitting a sphere to the trabecula and averaging the diameter of the sphere for all trabeculae
* TbSp: Trabecular separation (mm): Average distance between trabeculae. Found using the same process as trabecular thickness but spheres are fit to the empty space rather than the trabeculae
* CtAr: Cortical bone area (mm^2): Average cross-sectional area within the tissue that is mineralized bone
* TAr: Tissue area (mm^2): Cross-sectional area within the outer edge of the diaphysis (e.g. the whole bone cross section), including both mineralized bone and interior bone marrow. Averaged across all cross-sections in the analyzed volume
* MaAr: Marrow area (mm^2): Average cross-sectional area of marrow tissue, equal to T.Ar – Bn.Ar
* CtArTAr: Cortical bone area fraction (%): Fraction of total tissue area occupied by mineralized bone
* CtTh: Cortical thickness (mm): Thickness of the cortical bone averaged across all cross-sections
* J: Polar area moment of inertia (mm^4) Geometric resistance of the cortical bone to twisting. For a hollow cylinder, this is , J = π(D4 – d4)/32, where D is the outer diameter of the circle, and d is the inner diameter
* CtTMD: Cortical tissue mineral density (mg cm^-3): Average amount of mineralized tissue within a volume of cortical bone, calculated from the grayscale value and converted to a physical density using a calibration phantom. This includes small (5 μm) pores containing cells and their dendrites, but not larger pores containing blood vessels
