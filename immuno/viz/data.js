var seed = 1;
Math.random = function() {
    var x = Math.sin(seed++) * 10000;
    return x - Math.floor(x);
}

function randomSample(lst, n) {
  return _.map(_.range(n), function(n) { return _.sample(lst, n, true); });
}

  ///////////////////////////
 //  Generate mock data.  //
///////////////////////////
var NUM_PEPTIDES = 50;
var ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
var peptides = _.range(NUM_PEPTIDES).map(function() {
  var length = _.random(55, 65),
      sequence = randomSample(ALPHABET, length).join(''),
      gene = randomSample(ALPHABET + "1234567890", _.random(3, 6)).join(''),
      mutStart = _.random(length-2),
      mutation = "" + randomSample(ALPHABET, 1) + mutStart + sequence[mutStart];
  return {sequence: sequence, gene: gene, length: length,
          mutStart: mutStart, mutEnd: mutStart+1, mutation: mutation};
});

_.flatMap = _.compose(_.flatten, _.map);

// Add epitopes & associated info to peptides.
window.DATA = _.map(peptides, function(peptide) {
  peptide.epitopes = _.flatMap(_.range(9, 16), function(epitopeLength) {
    return _.map(_.range(peptide.length-epitopeLength), function(startPosition) {
      var scores = {'HLA-A*34:01': {percentile: _.random(100), bindingScore: _.random(0, 1000)},
                    'HLA-A*02:01': {percentile: _.random(100), bindingScore: _.random(0, 1000)},
                    'HLA-B*12:09': {percentile: _.random(100), bindingScore: _.random(0, 1000)},
                    'HLA-B*07:21': {percentile: _.random(100), bindingScore: _.random(0, 1000)},
                    'HLA-C*19:13': {percentile: _.random(100), bindingScore: _.random(0, 1000)},
                    'HLA-C*11:11': {percentile: _.random(100), bindingScore: _.random(0, 1000)}};
      return {start: startPosition, length: epitopeLength, scores: scores, sequence: peptide.sequence.slice(startPosition, startPosition+epitopeLength+1)};
    });
  });
  return peptide;
});
