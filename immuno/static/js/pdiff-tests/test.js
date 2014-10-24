var pdifftest = require('./lib/pdifftest');

// 1 = number of tests in this suite.
casper.test.begin('Peptide list compare', 1, function(test) {
  // This sequence just configures the test, it doesn't run it.
  pdifftest.startTest(casper, 'index.html', function() {
    // ... click things ...

    var screenshot = pdifftest.takeScreenshot(casper, {"top":0,"left":0,"width":1000,"height":2000});
    pdifftest.checkAgainstGolden(casper, screenshot, 'pdiff-tests/golden/peptides-list.png');
  });

  // Until you call "test.done()", the test will hang!
  casper.run(function() {
    test.done();
  });
});
