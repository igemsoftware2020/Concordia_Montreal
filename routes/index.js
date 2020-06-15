var express = require('express');
var router = express.Router();
const mongo = require('mongodb');
var assert = require('assert');
var url = 'mongodb://localhost:27017';
const collection = "yeastGenes";

/* GET home page. */
router.get('/', function(req, res, next) {
  // insert_GSE4136_Yeast_Data();
  // insertNewColumn_GSE4136();
  addAccessionGSE4136();
  mongo.connect(url,function(err, client) {
    assert.equal(null,err);
    var db = client.db('igemConcordia2020');
    db.collection(collection).find({}).toArray((err, documents) => {
      if(err) {
        console.log(err);
      } else {
        console.log(documents.length);
        // console.log("test");
        res.render('index', { title: 'Home Page' , test: documents });
      }
    });
  })
});
router.post('/', function(req,res) {
  var temp =req.body.txt;
  JSON.stringify(temp);
  // console.log(temp);
  if(temp !== null && temp !== '' && temp !== undefined) {
    mongo.connect(url, function(err,client) {
      assert.equal(null,err);
      var db = client.db('igemConcordia2020');
      db.collection(collection).find({ "Platform_ORF" : temp }).toArray((err,documents) => {
        if(err){
          console.log(err);
          client.close();
        } else {
          if(documents === undefined || documents.length === 0 || documents.length === null) {
            console.log("No data for this ORF");
            console.log("Now render to new page and display no data available");
            client.close();
            var message = "Nothing was found";
            res.render('search', {test: message});
          } else {
            console.log("Data available!");
            console.log("Now render to new page and display data");
            // console.log(documents);
            client.close();
            res.render('search', {test: documents});
          }
        }
      })
    });
  }
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
function insertNewColumn_GSE4136(){
  mongo.connect(url,function(err,client) {
    assert.equal(null,err);
    var db = client.db('igemConcordia2020');
    db.collection(collection).update({},
        {$set : {"AccessionNumber":"GSE4136" }},
        {upsert:false,
          multi:true});
    console.log("added new column")
    // Client.close();
  })
}
function addAccessionGSE4136(){
  mongo.connect(url, function(err, client) {
    assert.equal(null,err);
    var db = client.db('igemConcordia2020');
    db.collection(collection).count((err, count) => {
      if(err) {
        console.log("Error in fetching from DB");
        client.close();
      } else {
        if (count !== 0 && !err) {
          db.collection(collection).find({ "AccessionNumber" : "GSE4136" }).toArray((err,documents) => {
            if(err){
                console.log(err);
                client.close();
            } else {
              if(documents === undefined || documents.length == null || documents.length === 0) {
                console.log("AccessionNumber is not in DB! add column!");
                insertNewColumn_GSE4136();
                client.close();
              } else {
                console.log("column already exists for this study");
                client.close();
              }
            }
          })
        } else {
          console.log("GSE4136 column is already added");
          client.close();
        }
      }
    });
  });
}



module.exports = router;
