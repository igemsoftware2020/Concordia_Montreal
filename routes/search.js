//===========================================================================
// © Copyright 2020 iGEM Concordia, Maher Hassanain, Benjamin Clark, Hajar Mouddene, Grecia Orlano 
// This file is part of AstroBio.
//
//     AstroBio is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     AstroBio is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with AstroBio.  If not, see <https://www.gnu.org/licenses/>.
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
