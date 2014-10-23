var pdifftest = require('./lib/pdifftest');

casper.test.begin('Clicked peptide exploded epitopes compare', 1, function(test) {
  // This sequence just configures the test, it doesn't run it.
  pdifftest.startTest(casper, 'index.html', function() {
    casper.click('.clickbox');
    casper.wait(2000, function() {
      var screenshot = pdifftest.takeScreenshot(casper, {"top":0,"left":0,"width":1000,"height":2000});
      pdifftest.checkAgainstGolden(casper, screenshot, 'pdiff-tests/golden/exploded-epitopes.png');
    });
  });

  // Until you call "test.done()", the test will hang!
  casper.run(function() {
    test.done();
  });
});
