(function() {
'use strict';

// TODO:
//     1. Show all peptides & epitopes (don't cut off at static height).
//     2. Support sorting with inspector view.
//     3. Support different sorting functions.
//     4. Display of immungenicity (like binding score bar chart).
//     5. Hover over score bars for more info.

var PEPTIDE_HEIGHT = 20,
    EPITOPE_LEN = 6,
    EPITOPE_INFO_HEIGHT = 125,
    MAX_PEPTIDE_LENGTH = 66,
    WIDTH = 1000,
    HEIGHT = 12000,
    GENE_WIDTH = 100;


var aaX = d3.scale.ordinal()
        .domain(d3.range(0, MAX_PEPTIDE_LENGTH))
        .rangeBands([0, WIDTH-GENE_WIDTH]),
    svg = d3.select('#peptides')
      .append('svg')
        .attr('width', WIDTH)
        .attr('height', HEIGHT)
      .append('g')
         .attr('transform', 'translate(0,'+PEPTIDE_HEIGHT+')');


function main(data) {
  renderPeptides(data);
  initializeThresholdSlider(data);
}

function renderPeptides(data) {
  var peptides = svg.selectAll('.peptide')
    .data(data, function(d) { return d.sequence; });

  peptides
    .enter().append('g')
      .attr('class', 'peptide')
      .attr('transform', function(d, i) {
        return 'translate(' + GENE_WIDTH + ',' +
          (i * (PEPTIDE_HEIGHT + 5)) + ')';
      });

  peptides
    .transition().duration(600)
      .attr('transform', function(d, i) {
        return 'translate(' + GENE_WIDTH + ',' +
          (i * (PEPTIDE_HEIGHT + 5)) + ')';
      });

  peptides.exit().remove();

  renderGenes(peptides);
  renderHighlightEpitopes(peptides);
  renderPeptideSequences(peptides);
  initializePeptideClickBoxes(peptides);
}

function renderGenes(peptides) {
  peptides.selectAll('.gene')
      .data(function(d,i) { return [d.gene]; })
    .enter().append('text')
      .attr('class', 'gene')
      .attr('dx', -GENE_WIDTH)
      .text(function(d) { return d; });
}

function renderHighlightEpitopes(peptides) {
  var epitopes = peptides.selectAll('.epitope')
      .data(function(d, i) {
        // TODO(ihodes): Is this a good way to handle this? This shows fewer
        //               than you'd normally see, so you just see the strongest
        //               in the overview/peptide list; once the list is
        //               exploded, all the epitopes with one MHC allelse above
        //               the score are shown.
        return epitopesAbove(3, d.epitopes, getThreshold());
      }, function(d, i) { return [d.start, d.sequence]; });

  epitopes
    .enter().append('g')
      .attr('transform', function(d, i) {
        return 'translate(' + aaX(d.start) + ',0)';
      })
      .attr('class', 'epitope')
    .append('rect')
      .attr('class', 'passing-epitope')
      .attr('width', function(d, i) {
        return aaX(d.start + d.length + 1) - aaX(d.start);
      })
      .attr('height', 16)
      .attr('x', -2)
      .attr('y', -12);

  epitopes
    .exit()
      .remove();

  // TODO(ihodes): can we not do these, or do them less?
  peptides.selectAll('.sequence').each(function() {
    this.parentNode.appendChild(this);
  });
  peptides.selectAll('.clickbox').each(function() {
    this.parentNode.appendChild(this);
  });
}

function renderPeptideSequences(peptides) {
  var sequence = peptides.selectAll('.sequence')
      .data(function(d) { return d.sequence.split(''); });

  sequence
    .enter().append('text')
      .text(function(d, i) { return d;})
      .attr('class', function(d, i) {
        if (i == this.parentNode.__data__.mutStart
            && i < this.parentNode.__data__.mutEnd) {
          return 'sequence mutation';
        } else {
          return 'sequence';
        }
      })
      .attr('width', 12)
      .attr('height', 12)
      .attr('y', 0)
      .attr('x', function(d, i) {
        return aaX(i);
      });
}

// This is a box we can click, to allow easy selecting of peptides.
function initializePeptideClickBoxes(peptides) {
  peptides.append('rect')
      .attr('class', 'clickbox')
      .attr('width', WIDTH)
      .attr('height', PEPTIDE_HEIGHT)
      .attr('y', -PEPTIDE_HEIGHT)
      .attr('opacity', 0)
      .on('click', function(d, i) {
        explodeEpitopes(peptides, this.parentNode, this, d);
      });
}

function explodeEpitopes(peptides, peptideEl, peptideClickBox, peptideData) {
  d3.select('#close-inspector')
      .style('display', 'inline')
      .on('click', function() {
        d3.select(this).style('display', 'none');
        collapseEpitopes(peptides, peptideEl, peptideClickBox);
      });

  d3.select(peptideClickBox)
      .on('click', function(d, i) {
        var xPos = d3.event.x - GENE_WIDTH - 6,
            acidIdx = aaX.domain()[d3.bisect(aaX.range(), xPos) - 1],
            epitopes = peptideData.epitopes;

        d3.selectAll('.epitope-info-window').remove();
        highlightAcid(acidIdx, peptideEl, peptideClickBox, peptideData);

        epitopes = epitopesOverlapping(acidIdx, epitopes);
        epitopes = sortEpitopes(epitopes, getThreshold());
        renderEpitopesWithScores(epitopes, peptideEl);
      });

  peptides.attr('class', 'unselected peptide');
  d3.select(peptideEl).attr('class', 'selected peptide');
  d3.selectAll('.unselected')
    .transition().duration(800).ease(d3.ease('quad'))
      .attr('transform', function(d, i) {
        return 'translate(' + ((2 * WIDTH) * (i%2 == 0 ? 1 : -1)) + ', '+
          (i * PEPTIDE_HEIGHT) + ')';
      });

  d3.select(peptideEl)
    .transition().duration(700)
      .attr('transform', function() {
        return 'translate(' + GENE_WIDTH + ', 5)'
      });

  renderEpitopesWithScores(peptideData.epitopes, peptideEl);

  // TODO(ihodes): Not sure this is what we want.
  // epitopes
  //     .on('mouseover', function(d, i) {
  //       showToolTip(d3.event.x, d3.event.y);
  //     })
  //     .on('mouseout', function(d, i) {
  //       d3.select('#tooltip').style('display', 'none');
  //     });
}

function showToolTip(x, y) {
  d3.select('#tooltip')
      .style('display', 'block')
      .style('-webkit-transform', function() {
        return 'translate3d( ' + (x-100) + 'px, ' + y + 'px, 0px)';
      });
}


function collapseEpitopes(peptides, peptideEl, peptideClickBox) {
  peptides
      .attr('class', 'peptide')
    .transition().duration(800)
      .attr('transform', function(d, i) {
        return 'translate(' + GENE_WIDTH + ',' + (i*PEPTIDE_HEIGHT) + ')';
      });

  d3.selectAll('.epitope-info-window').remove();

  var scores = d3.selectAll('.ep-score-chart').remove()

  var epitopes = d3.select(peptideEl).selectAll('.epitope');
  epitopes.transition().duration(500)
      .attr('transform', function(d, i) {
        return 'translate(' + aaX(d.start) + ', 0)';
      });
  epitopes.selectAll('.ep-sequence')
      .remove();

  renderHighlightEpitopes(peptides);

  epitopes.selectAll('rect')
      .on('mouseover', null)
      .on('mouseout', null)
      .on('click', null);

  d3.selectAll('.highlight').remove();
  d3.select(peptideClickBox)
      .on('click', function(d, i) {
        explodeEpitopes(peptides, this.parentNode, this, d);
      })
      .on('mousemove', null);

  peptideEl.appendChild(peptideClickBox);
}

function renderEpitopesWithScores(epitopesList, peptideEl) {
  var epitopes = d3.select(peptideEl).selectAll('.epitope')
      .data(sortEpitopes(epitopesList, getThreshold()), function(d, i) {
        return [d.start, d.sequence];
      });

  epitopes
    .enter().append('g')
      .attr('class', 'epitope')
      .attr('transform', function(d, i) {
        return 'translate(' + aaX(d.start) + ',0)';
      });

  // Draw epitope sequences.
  epitopes.selectAll('.ep-sequence')
      .data(function(d) { return d.sequence.split(''); })
    .enter().append('text')
      .text(function(d, i) { return d; })
      .style('font-size', '1.1em')
      .attr('class', 'ep-sequence')
      .attr('width', 12)
      .attr('height', 12)
      .attr('y', 0)
      .attr('x', function(d, i) {
        return aaX(i);
      });

  epitopes
      .filter(function(d) {
        // !... is to ensure we don't keep adding highlight rects to the same
        // epitopes.
        return passesEpitopeThreshold(d, getThreshold()) &&
          !d3.select(this).select('.passing-epitope').node();
      })
    .append('rect')
      .attr('class', 'passing-epitope')
      .attr('width', function(d, i) {
        return aaX(d.start + d.length + 1) - aaX(d.start);
      })
      .attr('height', 16)
      .attr('x', -2)
      .attr('y', -12);

  epitopes
    .append('rect')
      .attr('class', 'epitope-highlight')
      .attr('x', function(d, i) {
        return -aaX(d.start)-GENE_WIDTH;
      })
      .attr('y', -PEPTIDE_HEIGHT+8)
      .attr('width', WIDTH+GENE_WIDTH)
      .attr('height', PEPTIDE_HEIGHT-4)
      .attr('opacity', 0)
      .attr('fill', '#ddd')
      .on('mouseover', function(d, i) {
        d3.select(this).transition().attr('opacity', 0.25);
      })
      .on('mouseout', function(d, i) {
        d3.select(this).transition().attr('opacity', 0.0);
      })
      .on('click', function(d, i) {
        // TODO(ihodes): Not just remove the info window, but push back up
        //               all the epitopes below it.
        // TODO(ihodes): Chart/display binding and immugenicity, and whatever
        //               other metainformation would be cool.
        d3.selectAll('.epitope-info-window').remove();

        var epitopeEl = this.parentNode,
            peptide = d3.select(epitopeEl.parentNode),
            epSequence = d.sequence,
            yOffset = epitopeEl.getCTM().f; // f is the y translate
                                            // element in the element's
                                            // transformation matrix.

        // Moves down, by the size of the info window, all epitopes which
        // are rendered at a lower Y position than the epitope clicked on.
        peptide
          .selectAll('.epitope')
            .filter(function(d, i) {
              return this.getCTM().f > yOffset;
            })
          .transition()
            .attr('transform', function(d, i) {
              return 'translate(' + aaX(d.start) + ','
                + (this.getCTM().f + EPITOPE_INFO_HEIGHT - 20) + ')';
            });

        d3.select(epitopeEl)
          .append('rect')
            .attr('class', 'epitope-info-window')
            .attr('y', 0)
            .attr('x', function(d, i) {
              return -aaX(d.start) - GENE_WIDTH;
            })
            .attr('width', WIDTH+GENE_WIDTH)
            .attr('height', EPITOPE_INFO_HEIGHT);
      });


  epitopes
    .transition().duration(800)
      .attr('transform', function(d, i) {
        return 'translate(' + aaX(d.start) + ',' + (i*15 + PEPTIDE_HEIGHT) + ')';
      });

  epitopes
    .exit()
      .remove();

  var scores = epitopes
    .append('g')
      .attr('class', 'ep-score-chart')
      .attr('transform', function(d, i) {
        return 'translate(' + (-aaX(d.start)-GENE_WIDTH) + ',0)';
      })

  scores.selectAll('scoreBar')
      .data(function(d) { return _.reduce(d.scores, function(acc, v, k) {
        v.mch = k;
        acc.push(v);
        return acc;
      }, []); })
    .enter().append('rect')
      .attr('x', function(d, i) { return i*10; })
      .attr('y', function(d, i) { return -d.percentile/10; })
      .attr('height', function(d, i) { return d.percentile/10; })
      .attr('width', 10);

  // NB maybe too noisy?
  // scores
  //   .append('line')
  //     .attr('class', 'ep-score-threshold')
  //     .attr('x1', 0)
  //     .attr('x2', 60)
  //     .attr('y1', -getThreshold()/10)
  //     .attr('y2', -getThreshold()/10);
}

function highlightAcid(acidIdx, peptideEl, peptideClickBox, peptideData) {
  d3.selectAll('.highlight')
    .transition()
      .attr('opacity', 0).remove();

  d3.select(peptideClickBox.parentNode)
    .append('rect')
      .attr('x', aaX(acidIdx)-3)
      .attr('y', -12)
      .attr('height', HEIGHT)
      .attr('width', 15)
      .attr('class', 'highlight')
      .attr('fill', '#adadad')
      .attr('opacity', 0.0)
      .on('click', function() {
        d3.select(this).remove();
        var epitopes = sortEpitopes(peptideData.epitopes, getThreshold());
        renderEpitopesWithScores(epitopes, peptideEl);
      })
    .transition()
      .attr('opacity', 0.25);

  var epitopes = d3.select(peptideEl).selectAll('.epitope')
      .each(function(){
        this.parentNode.appendChild(this);
      });
}

function initializeThresholdSlider(data) {
  d3.select('#threshold')
      .on('input', function() {
        d3.select('#tval').text(this.value);
        d3.select('#controls .suffix').text(ordinalSuffix(this.value));
      })
      .on('change', function() {
        var thresh = parseInt(this.value);
        d3.select('#tval').text(this.value)
        d3.select('#controls .suffix').text(ordinalSuffix(this.value));

        data = data.sort(function(a,b){
          var aScore = peptideBindingScore(a, thresh),
          bScore = peptideBindingScore(b, thresh);
          return aScore > bScore ? -1 : 1;
        });

        renderPeptides(data);
      });
}

// Returns list of epitopes which overlap idx
function epitopesOverlapping(idx, epitopes) {
  return _.filter(epitopes, function(epitope) {
    return epitope.start <= idx && (epitope.start + epitope.length) >= idx;
  });
}

function getThreshold() {
  return parseInt(d3.select('#threshold').node().value);
}

// Returns number of passing epitopes (HLA with percentile >= threshhold).
  function peptideBindingScore(peptide, threshold) {
  return _.reduce(peptide.epitopes, function(acc, epitope) {
    return acc + epitopeBindingScore(epitope, threshold)
  }, 0);
}

// Return the number of binding scores above a given percentile threshold.
function epitopeBindingScore(epitope, threshold) {
  return _.reduce(epitope.scores, function(acc, score, hla) {
    if (score.percentile >= threshold)
      return acc + 1;
    else
      return acc;
  }, 0);
}

function sortEpitopes(epitopes, threshold) {
  return epitopes.sort(function(a, b) {
    var aScore = epitopeBindingScore(a, threshold),
        bScore = epitopeBindingScore(b, threshold);
    return aScore > bScore ? -1 : 1;
  });
}

// Returns all epitopes with percentile scores above a threshold.
function epitopesAbove(numAbove, epitopes, threshold) {
  return _.reduce(epitopes, function(acc, epitope) {
    var num = epitopeBindingScore(epitope, threshold);
    if (num >= numAbove) acc.push(epitope);
    return acc;
  }, []);
}

function passesEpitopeThreshold(epitope, threshold) {
  var num = epitopeBindingScore(epitope, threshold);
  return num > 0;
}

function ordinalSuffix(n) {
  var sn = String(n);
  if (sn.slice(-2)[0] == '1')
    return 'th';
  switch (sn[sn.length-1]) {
    case '1':
      return 'st';
    case '2':
      return 'nd';
    case '3':
      return 'rd';
    default:
      return 'th';
  }
}

main(window.DATA);

})();
