## Function to simulate a tree using the HiGeoSSE process. Depends on 'tree.classe' and 'diversitree:::default.argnames.classe" functions of the package diversitree.

## Entries will be a matrix with exactly 7 rows in the same order as the arguments of the GeoSSE model. The number of columns will inform the number of hidden states.
## Also a squared matrix with dimensions equal to a multiplier of 3 (ordered A, B, AB). Multiples of 3 will inform the number of hidden states.
## Function need to privide a option to just return a list with empty parameters for the own function in the correct format. This will help a lot when using this.
## Function will return the same vectors that the SimulateHisse does. However, returning the $results in the same format will be hard. Will just return a phylogeny with tip states maybe.

SimulateHiGeoSSE <- function(){

}
