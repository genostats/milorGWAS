# this is just a test for reading a .dose file (or dosages in a VCF)
read.dosage.file <- function(filename) read_dose_file(path.expand(filename))

# cette fonction vérifie le format de tout le fichier
dim.dosage.file <- function(filename) dose_file_dim(path.expand(filename))

# celle-ci ne lit que la première ligne
nb.inds.dosage.file <- function(filename) nb_inds_dose_file(path.expand(filename))

# cette fonction lit les samples names (if any)
samples.dosage.file <- function(filename) samples_dose_file(path.expand(filename))

