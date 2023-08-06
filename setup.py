import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
  name = 'HaploDynamics',
  packages=setuptools.find_packages(),
  version = 'v0.1-beta.0',
  license='gpl-3.0',
  author = 'Remy Tuyeras',
  author_email = 'rtuyeras@gmail.com',
  description = 'A python library to develop genomic data simulators',
  long_description = long_description,
  long_description_content_type = "text/markdown",
  url = 'https://github.com/remytuyeras/HaploDynamics',
  download_url = 'https://github.com/remytuyeras/HaploDynamics/archive/refs/tags/v0.1-beta.0.tar.gz',
  keywords = ['Simulator', 'Genomics', 'Genomic', 'Microarray','SNP chip','VCF', 'Linkage disequilibrium', 'Hardy-Weinberg equilibrium'],
  install_requires=['scipy'],
  classifiers=[
    'Development Status :: 4 - Beta',  #"3 - Alpha", "4 - Beta" or "5 - Production/Stable"
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.9',
  ],
)
