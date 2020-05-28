var express = require('express');
var router = express.Router();
const mongo = require('mongodb');
var assert = require('assert');
var url = 'mongodb://localhost:27017';
const collection = "yeastGenes";

/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('index', { title: 'Home Page' });
});

// testInsertToDB();

function testInsertToDB() {
  var test = {
    Name: "Maher",
    Major : "IT"
  };
  mongo.connect(url, function(err,client) {
    assert.equal(null,err);
    var db = client.db('igemConcordia2020');
    db.collection(collection).insertOne(test, function(err,res) {
      assert.equal(null,err);
      console.log("Inserted to DB successfully");
      client.close();
    });
  });
}



module.exports = router;
