step 1:
fit constitutive behaviors of materials.
"fitMaterialProperties.py" 

step 2:
generate geometry and material settings.
"python3 setGeomAndMatr.py"

step 3:
generate cae and inp files.
"abaqus cae noGUI=genCAEInABAQUS.py -- .\configuration.json" (win)
"abaqus cae noGUI=genCAEInABAQUS.py -- ./configuration.json" (linux)


step 4:
run abaqus simulation using inp file. Go to 'tension' folder.
"abaqus job=tension inp=tension.inp cpus=4 interactive" (only applicable in case of rather small model)
"sbatch launcherOftensionOfRope.sh" (submit the job to UniLu HPC)

step 5:
run postprocessing in abaqus
"abaqus cae noGUI=postProcessInAbq.py -- ./pickleSetup ./pickleWireAggregate" (linux)

step 6:
download the work directory from unilu hpc to local destination of required and run the following command
python3 postProcessOutAbq.py ./tension/results.npz "/home/lichen/Dropbox/ropemodels/scripts/ropeModel/ropeExperimentalResults/Sample Nr.1/03. Rope/Rope_8.305mm_OZ_LFx6.4_Strain_Force_Stress.xlsx" ./FU.png


APPENDIX

commands to run on UniLu HPC:


available packgaes on UniLu HPC using 'module keyword packageName':
lang/Python: lang/Python/2.7.18-GCCcore-10.2.0, lang/Python/3.8.6-GCCcore-10.2.0
lang/SciPy-bundle: lang/SciPy-bundle/2020.11-foss-2020b-Python-2.7.18, lang/SciPy-bundle/2020.11-foss-2020b, ...
cae/ABAQUS: cae/ABAQUS/2021-hotfix-2207

module load cae/ABAQUS/2021-hotfix-2207
module load lang/Python
module load lang/SciPy-bundle

------------------------------------------------------------
Third party Python packages:
python3 -m pip install --user shapely
python3 -m pip install --user sympy
python3 -m pip install --user matplotlib
