# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from immuno import load_vcf

def test_vcf_to_dataframe():
    vcf_file = 'example.vcf'
    df = load_vcf.vcf_to_dataframe(vcf_file)
    assert df is not None
    assert(len(df) == 3)
    assert(len(df.columns) == 8)

if __name__ == '__main__':
  from dsltools import testing_helpers
  testing_helpers.run_local_tests()


