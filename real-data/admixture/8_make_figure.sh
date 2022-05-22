# Created 'pong' environment with python=2.7, activated and ran 'pip install pong'
d="/project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/"
r=`echo ${d}admixture_results/`

# Create file map
cd $r
rm -f pong.filemap
for K in `seq 2 4`;do
  for rep in `seq 1 30`;do
    echo -e K${K}_rep${rep}"\t"${K}"\t"K${K}_rep${rep}.Q >> pong.filemap
  done
done

# Create 'ind2pop' file
awk '{print $1}' ${d}1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG_LD-LD_200kb-1kbwin-0.3.fam > pong.ind2pop

# Create population order file
rm -f pong.poporder
for P in YRI MSL LWK GWD ESN ACB ASW GENOA_AA CLM MXL PEL PUR GENOA_EA GBR CEU IBS TSI;do
  echo $P >> pong.poporder
done

# Create colour choice file
rm -f pong.colours
for col in "#A5170E" "#F7F056" "#437DBF" "#EE8026" "#4EB265" "#882E72";do
  echo $col >> pong.colours
done

# Open tunnel on LOCAL MACHINE
ssh -N -L 4000:localhost:4000 sgopalan@discovery.usc.edu
# (session will hang - open new tab and log in to discovery, then:)
cd /project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/admixture_results/
# (created 'pong' environment with python=2.7, activated and ran 'pip install pong')
conda activate pong
pong -m pong.filemap -i pong.ind2pop -n pong.poporder -l pong.colours --port 4000
reated 'pong' environment with python=2.7, activated and ran 'pip install pong'
d="/project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/"
r=`echo ${d}admixture_results/`

# Create file map
cd $r
K=3
rm -f GENOA_pong.filemap
for rep in `seq 1 30`;do
    echo -e K${K}_rep${rep}"\t"${K}"\t"GENOA-only_K${K}_rep${rep}.Q >> GENOA_pong.filemap
done

# Create population order file
rm -f GENOA_pong.poporder
for P in GENOA_AA GENOA_EA;do
  echo $P >> GENOA_pong.poporder
done

# Create colour choice file
rm -f GENOA_pong.colours
for col in "#A5170E" "#F7F056" "#437DBF";do
  echo $col >> GENOA_pong.colours
done

# Open tunnel on LOCAL MACHINE
ssh -N -L 4000:localhost:4000 sgopalan@discovery.usc.edu
# (enter password and do 2-factor - session will hang - open new tab and log in to discovery, then:)
cd /project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/admixture_results/
# (created 'pong' environment with python=2.7, activated and ran 'pip install pong')
conda activate pong
pong -m GENOA_pong.filemap -i GENOA_pong.ind2pop -n GENOA_pong.poporder -l GENOA_pong.colours --port 4000


