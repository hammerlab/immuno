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

class ImmunoPipeline():
  def __init__(self):
    self._generators = []
    self._scorers = []

  def add_scorer(self, element):
    self._scorers.append(element)

  def add_generator(self, element):
    self._scorers.append(element)

  def score(self, data):
    for element in self._scorers:
      data = element.apply(data)

    return data

  def generate(self):
    for element in self._generators:
      data = self._generators(self._epitopes, self._alleles)

    return data


class PipelineElement(object):
  def __init__(self, name):
    self.name = name
    pass

  def apply(self, data):
    data[self.name] = self._apply(data)
    return data

  def _apply(self, data):
    pass

  def verify(self):
    pass