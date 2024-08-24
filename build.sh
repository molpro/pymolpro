# This script will install pymolpro

conda install -y -c conda-forge --file requirements.txt
python -m pip install --no-deps --force-reinstall --verbose . 

