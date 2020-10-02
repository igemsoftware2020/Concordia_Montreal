//===========================================================================
//						search.js
//				  Developed by Maher Hassanain
//						 June 2020
//===========================================================================
var express = require('express');
var router = express.Router();
const mongo = require('mongodb');
var assert = require('assert');
var url = 'mongodb://localhost:27017';
// var url = 'mongodb://mongo.ccb.lab:27017';
const collection = "geneResults";
const secondCollection = "metaData";
/* GET users listing. */
router.get('/', function(req, res, next) {
    res.send('respond with a resource');
});



module.exports = router;
