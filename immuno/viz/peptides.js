(function() {
'use strict';

// TODO:
//     1. (?) Display of immungenicity (like binding score bar chart).
//     2. (?) Hover over score bars for more info.
//     3. (?) Dynamic SVG height (pay attention to perceptual diffs then).

var WIDTH = 1000,
    HEIGHT = 12000,
    PEPTIDE_HEIGHT = 20,
    ACID_DIM = 12,
    EPITOPE_INFO_HEIGHT = 125,  // global
    SLIDER_TYPE = 'percentile', // global
    SLIDER_BINDING_SCORE = 500, // global
    SLIDER_PERCENTILE = 8,      // global
    GENE_LETTER_WIDTH = 14,
    GENE_WIDTH,
    EPITOPE_INFO_HEIGHT = 125,
    acidX = d3.scale.ordinal();


function main(data) {
  GENE_WIDTH = d3.max(data, function(d) { return d.gene.length; }) * GENE_LETTER_WIDTH;

  acidX
    .rangeBands([0, WIDTH - GENE_WIDTH])
    .domain(d3.range(0, d3.max(data, function(d) {
    return d.sequence.length;
  })));



  d3.select('#peptides')
      .append('svg')
        .attr('width', WIDTH)
        .attr('height', HEIGHT)
      .append('g')
         .attr('transform', 'translate(0,' + PEPTIDE_HEIGHT + ')')
         .attr('id', 'svg');

  renderPeptides(data);
  initializeSliderHandler(data);
}

function renderPeptides(data) {
  var peptides = d3.select('#svg')
    .selectAll('.peptide')
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

  peptides
      .call(renderGenes)
      .call(renderHighlightEpitopes)
      .call(renderPeptideSequences)
      .call(initializePeptideHandlers);
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
        return epitopesPassing(getSliderAttr(), 3, getSliderValue(), d.epitopes);
      }, _epitope_key_fn);

  epitopes
    .enter().append('g')
      .attr('transform', function(d, i) {
        return 'translate(' + acidX(d.start) + ',0)';
      })
      .attr('class', 'epitope')
    .call(renderPassingEpitopeHighlights);

  epitopes
    .exit()
      .remove();

  // Ensure epitopes are below the letters of the sequence and the click box.
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
        if (i >= this.parentNode.__data__.mutStart
            && i < this.parentNode.__data__.mutEnd) {
          return 'sequence mutation';
        } else {
          return 'sequence';
        }
      })
      .attr('width', ACID_DIM)
      .attr('height', ACID_DIM)
      .attr('y', 0)
      .attr('x', function(d, i) {
        return acidX(i);
      });
}

function initializePeptideHandlers(peptides) {
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
  initializeEpitopeSliderHandler(peptideData.epitopes, peptideEl);

  d3.select('#close-inspector')
      .style('display', 'inline')
      .on('click', function() {
        d3.select(this).style('display', 'none');
        collapseEpitopes(peptides, peptideEl, peptideClickBox);
      });

  d3.select(peptideClickBox)
      .on('click', function(d, i) {
        var xPos = d3.event.x - GENE_WIDTH - 6,
            acidIdx = acidX.domain()[d3.bisect(acidX.range(), xPos) - 1],
            epitopes = peptideData.epitopes;

        highlightAcid(acidIdx, peptideEl, peptideClickBox, peptideData);

        epitopes = epitopesOverlapping(acidIdx, epitopes);
        epitopes = sortEpitopes(getSliderAttr(), getSliderValue(), epitopes);
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
}

function collapseEpitopes(peptides, peptideEl, peptideClickBox) {
  resetEpitopeInfoWindow();
  initializeSliderHandler(peptides.data());

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
        return 'translate(' + acidX(d.start) + ', 0)';
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
  resetEpitopeInfoWindow();

  var epitopes = d3.select(peptideEl).selectAll('.epitope')
      .data(sortEpitopes(getSliderAttr(), getSliderValue(), epitopesList), function(d, i) {
        return [d.start, d.sequence];
      });

  epitopes
    .enter().append('g')
      .attr('class', 'epitope')
      .attr('transform', function(d, i) {
        return 'translate(' + acidX(d.start) + ',0)';
      });

  epitopes
      .call(renderEpitopeSequence)
      .call(renderEpitopeScores)
      .call(renderPassingEpitopeHighlights)
      .call(initializeEpitopeHandlers);

  epitopes
    .transition().duration(800)
      .attr('transform', function(d, i) {
        return 'translate(' + acidX(d.start) + ',' + (i*15 + PEPTIDE_HEIGHT) + ')';
      });

  epitopes
    .exit()
      .remove();
}

function renderEpitopeScores(epitopes) {
  var scores = epitopes
    .append('g')
      .attr('class', 'ep-score-chart')
      .attr('transform', function(d, i) {
        return 'translate(' + (-acidX(d.start)-GENE_WIDTH) + ',0)';
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
}

function renderEpitopeSequence(epitopes) {
  epitopes.selectAll('.ep-sequence')
      .data(function(d) { return d.sequence.split(''); })
    .enter().append('text')
      .text(function(d, i) { return d; })
      .style('font-size', '1.1em')
      .attr('class', 'ep-sequence')
      .attr('width', ACID_DIM)
      .attr('height', ACID_DIM)
      .attr('y', 0)
      .attr('x', function(d, i) {
        return acidX(i);
      });

}

function renderPassingEpitopeHighlights(epitopes) {
  epitopes
      .selectAll('.passing-epitope')
    .remove();

  epitopes
    .filter(function(d) {
        // !... is to ensure we don't keep adding highlight rects to the same
        // epitopes.
      var alleles = allelesBelow(getSliderAttr(), getSliderValue(), d);
      return (alleles > 0) && !d3.select(this).select('.passing-epitope').node();
    })
    .append('rect')
      .attr('class', 'passing-epitope')
      .attr('width', function(d, i) {
        return acidX(d.start + d.length) - acidX(d.start) + ACID_DIM;
      })
      .attr('height', 16)
      .attr('x', -2)
      .attr('y', -ACID_DIM);
}

function initializeEpitopeHandlers(epitopes) {
  epitopes
    .append('rect')
      .attr('class', 'epitope-highlight')
      .attr('x', function(d, i) {
        return -acidX(d.start)-GENE_WIDTH;
      })
      .attr('y', -PEPTIDE_HEIGHT+8)
      .attr('width', WIDTH+GENE_WIDTH)
      .attr('height', PEPTIDE_HEIGHT-4)
      .attr('opacity', 0)
      .attr('fill', '#ddd')
      .on('click', function(d, i) {
        resetEpitopeInfoWindow();
        renderEpitopeInfoWindow(this.parentNode, d);
      });
}

function renderEpitopeInfoWindow(epitopeEl, epitopeData) {
  var peptide = d3.select('.selected.peptide'),
      epitopeShiftGroup = peptide.append('g').attr('id', 'epitope-shift-group'),
      yOffset = epitopeEl.getCTM().f;

  d3.select(epitopeEl).select('.epitope-highlight').transition().attr('opacity', .25);

  peptide
    .selectAll('.epitope')
    .filter(function(d, i) {
      return this.getCTM().f > yOffset;
    })
      .each(function(d, i) {
        epitopeShiftGroup.node().appendChild(this);
      });

  epitopeShiftGroup
      .attr('transform', function(d, i) {
        return 'translate(0,' + (PEPTIDE_HEIGHT + EPITOPE_INFO_HEIGHT) + ')';
      });

  d3.select(epitopeEl)
    .append('g')
      .attr('id', 'epitope-info-window')
      .attr('transform', function(d, i) {
        return 'translate(' + (-acidX(d.start) - GENE_WIDTH) + ', '
          + PEPTIDE_HEIGHT + ')';
      })
      .call(renderEpitopeInfoWindowChart);
}

function resetEpitopeInfoWindow(peptide) {
  var peptide = d3.select('.selected.peptide'),
      epitopeShiftGroup = d3.select('#epitope-shift-group');
  d3.selectAll('.epitope-highlight').attr('opacity', 0)
  d3.select('#epitope-info-window').remove();
  epitopeShiftGroup
      .attr('transform', function() {
        return 'translate(0, 0)';
      })
    .selectAll('.epitope')
      .each(function(d, i) {
        peptide.node().appendChild(this);
      });
  epitopeShiftGroup.remove();
}

function renderEpitopeInfoWindowChart(epitopeInfoWindow) {
  var epitopeData = this.node().__data__,
      scores = _.reduce(epitopeData.scores, function(acc, v, k) {
        v.name = k;
        acc.push(v);
        return acc;
      }, []);

  var scoreBody = epitopeInfoWindow
    .append('foreignObject')
      .attr('width', WIDTH+GENE_WIDTH)
      .attr('height', EPITOPE_INFO_HEIGHT)
    .append('xhtml:body')
      .attr('class', 'epitope-body'),
      tbl = scoreBody.append('table'),
      header = tbl.append('thead')
    .append('tr');

  header.append('td').html('HLA Allele');
  header.append('td').html('Binding Score');
  header.append('td').html('Percentile');

  var rows = tbl.selectAll('tr')
      .data(scores)
    .enter().append('tr');

  rows.append('td')
      .html(function(d) { return d.name; });
  rows.append('td')
      .html(function(d) { return d.bindingScore; });
  rows.append('td')
      .html(function(d) { return d.percentile; });
}

function highlightAcid(acidIdx, peptideEl, peptideClickBox, peptideData) {
  d3.selectAll('.highlight')
    .transition()
      .attr('opacity', 0)
      .remove();

  d3.select(peptideClickBox.parentNode)
    .append('rect')
      .attr('x', acidX(acidIdx)-3)
      .attr('y', -ACID_DIM)
      .attr('height', HEIGHT)
      .attr('width', 15)
      .attr('class', 'highlight')
      .attr('fill', '#adadad')
      .attr('opacity', 0.0)
      .on('click', function() {
        d3.select(this).remove();
        var epitopes = sortEpitopes(getSliderAttr(), getSliderValue(), peptideData.epitopes);
        renderEpitopesWithScores(epitopes, peptideEl);
      })
    .transition()
      .attr('opacity', 0.25);

  var epitopes = d3.select(peptideEl).selectAll('.epitope')
      .each(function() {
        this.parentNode.appendChild(this);
      });
}


// Helpers for xxxSliderHandlers
function setSliderText() {
  if (SLIDER_TYPE === 'percentile') {
    d3.select('#tval').text(SLIDER_PERCENTILE);
    d3.select('#suffix').text(ordinalSuffix(SLIDER_PERCENTILE));
  } else if (SLIDER_TYPE === 'ic50') {
    d3.select('#tval').text(SLIDER_BINDING_SCORE);
    d3.select('#suffix').text('mMol');
  }
}

function sortPeptides(peptides, attr, threshold) {
  return peptides.sort(function(a,b){
    var passing = _.partial(peptidesPassingEpitopes, attr, threshold);
    return passing(a) > passing(b) ? -1 : 1;
  });
}

function initializeEpitopeSliderHandler(epitopes, peptideEl) {
  d3.select('input[value=percentile]')
    .on('click', function() {
      SLIDER_TYPE = 'percentile';
      d3.select('#slider')
        .attr('min', 1)
        .attr('max', 99)
      .node().value = 100 - SLIDER_PERCENTILE; // attr('value', ..) wasn't working.
      d3.select('span#slider-type').text('Binding score percentile');
      setSliderText();

      epitopes = sortEpitopes('percentile', SLIDER_PERCENTILE, epitopes);
      renderEpitopesWithScores(epitopes, peptideEl);
    });

  d3.select('input[value=ic50]')
    .on('click', function() {
      SLIDER_TYPE = 'ic50';
      d3.select('#slider')
        .attr('min', 1)
        .attr('max', 2500)
      .node().value = SLIDER_BINDING_SCORE; // attr('value', ..) wasn't working.
      d3.select('span#slider-type').text('Binding score IC50');
      setSliderText();

      epitopes = sortEpitopes('bindingScore', SLIDER_BINDING_SCORE, epitopes);
      renderEpitopesWithScores(epitopes, peptideEl);
    });

  d3.select('#slider')
      .on('input', function() {
        if (SLIDER_TYPE === 'percentile') {
          SLIDER_PERCENTILE = 100 - parseInt(this.value);
        } else if (SLIDER_TYPE === 'ic50') {
          SLIDER_BINDING_SCORE = parseInt(this.value);
        }

        setSliderText();
      })
      .on('change', function() {
        if (SLIDER_TYPE === 'percentile') {
          SLIDER_PERCENTILE = 100 - parseInt(this.value);
          epitopes = sortEpitopes('percentile', SLIDER_PERCENTILE, epitopes);
        } else if (SLIDER_TYPE === 'ic50') {
          SLIDER_BINDING_SCORE = parseInt(this.value);
          epitopes = sortEpitopes('bindingScore', SLIDER_BINDING_SCORE, epitopes);
        }

        renderEpitopesWithScores(epitopes, peptideEl);
        setSliderText();
      });
}

function initializeSliderHandler(peptides) {
  d3.select('input[value=percentile]')
    .on('click', function() {
      SLIDER_TYPE = 'percentile';
      d3.select('#slider')
        .attr('min', 1)
        .attr('max', 99)
      .node().value = 100 - SLIDER_PERCENTILE; // attr('value', ..) wasn't working.
      d3.select('span#slider-type').text('Binding score percentile');
      setSliderText();

      peptides = sortPeptides(peptides, 'percentile', SLIDER_PERCENTILE);
      renderPeptides(peptides);
    });

  d3.select('input[value=ic50]')
    .on('click', function() {
      SLIDER_TYPE = 'ic50';
      d3.select('#slider')
        .attr('min', 1)
        .attr('max', 2500)
      .node().value = SLIDER_BINDING_SCORE; // attr('value', ..) wasn't working.
      d3.select('span#slider-type').text('Binding score IC50');
      setSliderText();

      peptides = sortPeptides(peptides, 'bindingScore', SLIDER_BINDING_SCORE);
      renderPeptides(peptides);
    });

  d3.select('#slider')
      .on('input', function() {
        if (SLIDER_TYPE === 'percentile') {
          SLIDER_PERCENTILE = 100 - parseInt(this.value);
        } else if (SLIDER_TYPE === 'ic50') {
          SLIDER_BINDING_SCORE = parseInt(this.value);
        }

        setSliderText();
      })
      .on('change', function() {
        if (SLIDER_TYPE === 'percentile') {
          SLIDER_PERCENTILE = 100 - parseInt(this.value);
          peptides = sortPeptides(peptides, 'percentile', SLIDER_PERCENTILE);
        } else if (SLIDER_TYPE === 'ic50') {
          SLIDER_BINDING_SCORE = parseInt(this.value);

          peptides = sortPeptides(peptides, 'bindingScore', SLIDER_BINDING_SCORE);
        }

        renderPeptides(peptides);
        setSliderText();
      });
}




/************************************************
 ** Utilities not specific to d3 visualiztion. **
 ************************************************
 */


function getSliderValue() {
  return SLIDER_TYPE == 'percentile' ? SLIDER_PERCENTILE : SLIDER_BINDING_SCORE;
}

function getSliderAttr() {
  if (SLIDER_TYPE === 'percentile') return 'percentile';
  else if (SLIDER_TYPE === 'ic50') return 'bindingScore';
}

// Returns list of epitopes which overlap idx
function epitopesOverlapping(idx, epitopes) {
  return _.filter(epitopes, function(epitope) {
    return epitope.start <= idx && (epitope.start + epitope.length) >= idx;
  });
}

// Returns number of passing epitopes
function peptidesPassingEpitopes(peptide, attr, threshold) {
  return _.reduce(peptide.epitopes, function(acc, epitope) {
    return acc + allelesBelow(epitope, attr, threshold)
  }, 0);
}

// Return the number of alleles epitope has with attr below threshold
function allelesBelow(attr, threshold, epitope) {
  return _.reduce(epitope.scores, function(acc, score) {
    if (score[attr] < threshold)
      return acc + 1;
    else
      return acc;
  }, 0);
}

function sortEpitopes(attr, threshold, epitopes) {
  var qual = _.partial(allelesBelow, attr, threshold);
  return epitopes.sort(function(a, b) {
    return qual(a) > qual(b) ? -1 : 1;
  });
}

// Returns all epitopes with number of alleles below a binding score percentile
// threshold
function epitopesPassing(attr, numAbove, threshold, epitopes) {
  var alleles = _.partial(allelesBelow, attr, threshold);
  return _.filter(epitopes, function(epitope) {
    return alleles(epitope) >= numAbove;
  });
}

function ordinalSuffix(n) {
  var sn = String(n);
  if (sn.length > 1 && sn.slice(-2)[0] == '1')
    return 'th';

  switch (sn.slice(-1)) {
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

function _epitope_key_fn(d, i) {
  return [d.start, d.sequence];
}


  ///////////////
 //  Export:  //
///////////////

window.peptides = window.peptides || {};
window.peptides.run = main;

})();
