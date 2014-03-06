"""

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""

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