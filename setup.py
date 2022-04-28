import subprocess
import os
import setuptools

# Installation of MULTI

# Define pymulti_path. This path will store the code of nicole and any other temporary files in the calculation.
pymulti_path = '{}/.pymulti/'.format(os.environ['HOME'])
multi_path = pymulti_path + 'files/multi/mul23/'

# Create the folder according to pymulti_path
if not(os.path.isdir(pymulti_path)):
    os.mkdir(pymulti_path)
    
# Copy the files folder into working directory
if not(os.path.isdir(pymulti_path + 'files')):
    os.mkdir(pymulti_path + 'files')
rm_status = subprocess.run(['rm', '-r', pymulti_path + 'files/'], stdout=subprocess.PIPE)
cp_status = subprocess.run(['cp', '-r', 'pymulti/files', pymulti_path + 'files'], stdout=subprocess.PIPE)

# Install MULTI
#   Make the double precision version of MULTI
mke_dbl = subprocess.run(['./make_mul23_dbl.sh', '-f'],  cwd=multi_path + 'source/', stdout=subprocess.PIPE)
#   Install dp version of MULTI
install = subprocess.run(['make', 'mul23gus.x'],  cwd=multi_path + 'source_dp/', stdout=subprocess.PIPE)
install = subprocess.run(['make', 'mul23lus.x'],  cwd=multi_path + 'source_dp/', stdout=subprocess.PIPE)

# Check if NICOLE is in the folder
if not(os.path.isfile(multi_path+'source_dp/mul23gus.x')) or not(os.path.isfile(multi_path+'source_dp/mul23lus.x')):
    raise ValueError("MULTI is not installed correctly!")
else:
    print('Successfully installed MULTI!')
        
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
      name='pymulti',
      version='0.0.1',
      description='The python wrapper to run NLTE spectra synthesis code MULTI.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/MingjieJian/pymulti',
      author='Mingjie Jian',
      author_email='ssaajianmingjie@gmail.com',
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Framework :: IPython",
        "Operating System :: OS Independent",
        "Development Status :: 1 - Planning",
        "Topic :: Scientific/Engineering :: Astronomy"
      ],
      python_requires=">=3.5",
      packages=setuptools.find_packages(),
      install_requires=[
          'numpy >= 1.18.0',
          'pandas >= 1.0.0',
          'matplotlib >= 3.1.0',
          'mendeleev >= 0.6.0',
          'scipy >= 1.4.0',
          'astropy >= 4.0',
          'spectres',
          'tqdm',
          'pymoog',
          'rulerwd',
      ],
      include_package_data=True,  
    #   package_data={'': ['moog_nosm/moog_nosm_FEB2017/']},
      zip_safe=False)