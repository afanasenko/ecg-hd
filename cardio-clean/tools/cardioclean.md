

Package cardio-clean contains single command **cardioclean** to process ECG files

*Example:*

$ cardioclean [-h] [-c \<CONFIG>] -i \<SOURCE> -o \<DEST>

- CONFIG

Yaml config file (Use sample config.yaml found in the repository as a start point)

default: config.yaml in current directory.

- SOURCE

Path to source ECG record in MIT format. Relative or full path to .hea or .dat file should be specified

- DEST

Base name for output record. Cardioclean generates a (.dat + .hea) file pair in current directory with specified DEST base name.

