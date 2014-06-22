# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os

readme_filename = os.path.join(os.path.dirname(__file__), 'README.md')
with open(readme_filename, 'r') as f:
  readme = f.read()

try:
  import pypandoc
  readme = pypandoc.convert(readme, to='rst', format='md')
except:
  print "Conversion of long_description from markdown to reStructuredText failed, skipping..."
  pass


from setuptools import setup

if __name__ == '__main__':
    setup(
        name='immuno',
        version="0.2",
        description="Evaluate immunotherapy from given set of variants",
        url="https://github.com/hammerlab/immuno",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            'numpy>=1.7',
            'pandas>=0.13.1',
            'scikit-learn>=0.14.1',
            'biopython',
            'dsltools',
            'datacache',
            'epitopes',
        ],
        long_description=readme,
        packages=['immuno'],
        package_data = {'immuno' : ['data/*txt']},
        include_package_data = True
    )
