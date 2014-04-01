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

import appdirs
import logging

logging.basicConfig(
    format="[%(levelname)s %(filename)s:%(lineno)d %(funcName)s] %(message)s",
    level=logging.DEBUG)

def peptide_substrings(full_peptide, window_length):
    n = len(full_peptide)
    return [full_peptide[i:i+window_length]
            for i in xrange(n + 1 - window_length)]
