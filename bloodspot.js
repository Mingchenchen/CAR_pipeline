
// create casper instance
var casper = require("casper").create({

	loadImages:  false,
	loadPlugins: false,
	pageSettings: {
		webSecurityEnabled: false
	},
	verbose: true,
	loglevel: 'debug'

});

// utility imports
require("utils").dump(casper.cli.options);

// variables
var gene = casper.cli.get("gene");
var fs = require('fs');
var fname = gene + '.csv';
var save = fs.pathJoin(fs.workingDirectory, 'bloodspot', fname);

// init casper
casper.start();

// user agent assignment
casper.userAgent('Mozilla/5.0 (Macintosh; Intel Mac OS X)');

// collection function
casper.thenOpen('http://servers.binf.ku.dk/bloodspot/?gene=' + gene + '&dataset=MERGED_AML', function() {

	this.echo("[ page: " + this.getTitle() + " ]");

	this.wait(1000, function() {
		var csv = this.page.evaluate(function() {
			return window.getDataInFormatCSV();
		});
		//console.log(csv);
		fs.write(save, csv, 'w');
	});

});

// run
casper.run(function() {

	this.echo('[ finished ]');
	this.exit();

});
