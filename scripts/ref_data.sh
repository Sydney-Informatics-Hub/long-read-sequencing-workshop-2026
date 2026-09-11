cd ~
mkdir -p data/ref
cd data/ref

mkdir amrfinderplus_db
cd amrfinderplus_db
amrfinder_update -d amrfinderplus_db

cd ../
plassembler-exec  plassembler download -d plasmid_db_plassembler

mkdir -p busco
cd busco
ln -s  /cvmfs/data.galaxyproject.org/byhand/busco/all-odb12-2026-03-20-145944/lineages/bacteria_odb12 bacteria_odb12.2

cd ../
ln -s  /cvmfs/data.galaxyproject.org/managed/kraken2_databases/kalamari kalamari

