//===========================================================================
//						index.js
//				  Developed by Maher Hassanain
//						 June 2020
//===========================================================================
var express = require('express');
var router = express.Router();
const mongo = require('mongodb');
var assert = require('assert');
var url = 'mongodb://localhost:27017';
const collection = "yeastGenes";
const secondCollection = "metaData";
var fs = require('fs');
var fasta = require('bionode-fasta');
// var blast = require('blastjs');
//
// var type = 'nucl';
// var fileIn = './FASTA//sample.fasta';
// var outPath = './';
// var name = 'example';


// blast.makeDB(type, fileIn, outPath, name, function(err){
//   if(err){
//     console.error(err);
//   } else {
//     console.log('database created at', outPath);
//   }
// });


// fs.createReadStream('./FASTA/sample.fasta')
//     .pipe(fasta())
//     .pipe(process.stdout);
//
// fasta.obj('./FASTA/sample.fasta').on('data', console.log);


// fasta.obj('./input.fasta').on('data', console.log)


/* GET home page. */
router.get('/', function(req, res, next) {
  // insert_GSE4136_Yeast_Data();
  // insertNewColumn_GSE4136();

  // addAccessionGSE4136();

  // displayDB();
  // UN-COMMENT BEFORE LAUNCH

  convertFloatlogFC();

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
  var temp = req.body.txt;
  var logfc = req.body.sel1;
  var pval = req.body.adjpval;
  var orf = req.body.platORF;
  JSON.stringify(temp);
  JSON.stringify(logfc);
  JSON.stringify(pval);
  JSON.stringify(orf);
  // console.log(temp);
  // console.log(logfc);
  // console.log(pval);

  if(pval !== null && pval !== '' && pval !== undefined && temp !== null && temp !== '' && temp !== undefined) {
    // search for all
    console.log("All data filled up");

    if(logfc === "Up") {
      mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        // let flt = pval;
        db.collection(collection).find({
          // "Gene.symbol" : temp
          $and : [{"Gene.symbol" : temp} , {"logFC" : {$gt : 0.0}} , {"adj.P.Val" : {$lt : parseFloat(pval)}}]
        }).sort({ "logFC" : -1 }).toArray((err,documents) => {
          if(err){
            console.log(err);
            client.close();
          } else {
            if(documents === undefined || documents.length === 0 || documents.length === null) {
              console.log("No data found for this gene");
              console.log("Now render to new page and display no data available");
              client.close();
              var message = "Nothing was found";
              res.render('search', {data: message});
            } else {
              console.log("Data available!");
              console.log("Now render to new page and display data");
              let temporary = documents[0].meta_data;
              console.log(temporary);
              db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
                if(err) {
                  console.log("error in fetching meta data");
                } else {
                  console.log("Meta-data fetched successfully");
                  // console.log(doc[0]);
                  client.close();
                  res.render('search', {data: documents , metaData : doc});
                }
              });
              // console.log(documents);
            }
          }
        })
      });
    } else if(logfc === "Down") {
      mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find({
          // "Gene.symbol" : temp
          $and : [{"Gene.symbol" : temp} , {"logFC" : {$lt : 0.0}} , {"adj.P.Val" : {$lt : parseFloat(pval)}}]
        }).sort({"logFC" : 1}).toArray((err,documents) => {
          if(err){
            console.log(err);
            client.close();
          } else {
            if(documents === undefined || documents.length === 0 || documents.length === null) {
              console.log("No data found for this gene");
              console.log("Now render to new page and display no data available");
              client.close();
              var message = "Nothing was found";
              res.render('search', {data: message});
            } else {
              console.log("Data available!");
              console.log("Now render to new page and display data");
              let temporary = documents[0].meta_data;
              console.log(temporary);
              db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
                if(err) {
                  console.log("error in fetching meta data");
                } else {
                  console.log("Meta-data fetched successfully");
                  // console.log(doc[0]);
                  client.close();
                  res.render('search', {data: documents , metaData : doc});
                }
              });
              // console.log(documents);
            }
          }
        })
      });
    } else {
      // do nothing
    }
  } else if (temp !== null && temp !== '' && temp !== undefined && !pval) {
    // search based on gene name
    if(logfc === "Up") {
      //find positive
      // console.log("Uppppp");
      mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find({
          // "Gene.symbol" : temp
          $and : [{"Gene.symbol" : temp} , {"logFC" : {$gt : 0.0}}]
        }).sort({"logFC" : -1}).toArray((err,documents) => {
          if(err){
            console.log(err);
            client.close();
          } else {
            if(documents === undefined || documents.length === 0 || documents.length === null) {
              console.log("No data found for this gene");
              console.log("Now render to new page and display no data available");
              client.close();
              var message = "Nothing was found";
              res.render('search', {data: message});
            } else {
              console.log("Data available!");
              console.log("Now render to new page and display data");
              let temporary = documents[0].meta_data;
              console.log(temporary);
              db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
                if(err) {
                  console.log("error in fetching meta data");
                } else {
                  console.log("Meta-data fetched successfully");
                  // console.log(doc[0]);
                  client.close();
                  res.render('search', {data: documents , metaData : doc});
                }
              });
              // console.log(documents);
            }
          }
        })
      });
    } else if (logfc === "Down") {
      // find negative
      // console.log("Dowwwwwwwn");
      mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find({
          // "Gene.symbol" : temp
          $and : [{"Gene.symbol" : temp} , {"logFC" : {$lt : 0.0}}]
        }).sort({"logFC" : 1}).toArray((err,documents) => {
          if(err){
            console.log(err);
            client.close();
          } else {
            if(documents === undefined || documents.length === 0 || documents.length === null) {
              console.log("No data found for this gene");
              console.log("Now render to new page and display no data available");
              client.close();
              var message = "Nothing was found";
              res.render('search', {data: message});
            } else {
              console.log("Data available!");
              console.log("Now render to new page and display data");
              let temporary = documents[0].meta_data;
              console.log(temporary);
              db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
                if(err) {
                  console.log("error in fetching meta data");
                } else {
                  console.log("Meta-data fetched successfully");
                  // console.log(doc[0]);
                  client.close();
                  res.render('search', {data: documents , metaData : doc});
                }
              });
              // console.log(documents);
            }
          }
        })
      });
    } else {
      // do Nothing, if someone wants to inject their own input
    }
  } else if(pval !== null && pval !== '' && pval !== undefined && !temp) {
    // search based on pval limit
    console.log("no gene name --- search for PVAL ONLY");
    if(logfc === "Up") {
      mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find({
          // "Gene.symbol" : temp
          $and : [ {"adj.P.Val" : {$lt : parseFloat(pval)}} , {"logFC" : {$gt : 0.0}}]
        }).sort({"logFC" : -1}).limit(100).toArray((err,documents) => {
          if(err){
            console.log(err);
            client.close();
          } else {
            if(documents === undefined || documents.length === 0 || documents.length === null) {
              console.log("No data found for this gene");
              console.log("Now render to new page and display no data available");
              client.close();
              var message = "Nothing was found";
              res.render('search', {data: message});
            } else {
              console.log("Data available!");
              console.log("Now render to new page and display data");
              let temporary = documents[0].meta_data;
              console.log(temporary);
              db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
                if(err) {
                  console.log("error in fetching meta data");
                } else {
                  console.log("Meta-data fetched successfully");
                  // console.log(doc[0]);
                  client.close();
                  res.render('search', {data: documents , metaData : doc});
                }
              });
              // console.log(documents);
            }
          }
        })
      });
    } else if (logfc === "Down") {
      mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find({
          // "Gene.symbol" : temp
          $and : [ {"adj.P.Val" : {$lt : parseFloat(pval)}} , {"logFC" : {$lt : 0.0}}]
        }).sort({"logFC" : 1}).limit(100).toArray((err,documents) => {
          if(err){
            console.log(err);
            client.close();
          } else {
            if(documents === undefined || documents.length === 0 || documents.length === null) {
              console.log("No data found for this gene");
              console.log("Now render to new page and display no data available");
              client.close();
              var message = "Nothing was found";
              res.render('search', {data: message});
            } else {
              console.log("Data available!");
              console.log("Now render to new page and display data");
              let temporary = documents[0].meta_data;
              console.log(temporary);
              db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
                if(err) {
                  console.log("error in fetching meta data");
                } else {
                  console.log("Meta-data fetched successfully");
                  // console.log(doc[0]);
                  client.close();
                  res.render('search', {data: documents , metaData : doc});
                }
              });
              // console.log(documents);
            }
          }
        })
      });
    }  else {
      // console.log("else last")
    }
  } else {
    // No input specified, use the logFC dropdown as mean to view data
    // console.log("last else");
    if(logfc === "Up") {
      mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find({
          // "Gene.symbol" : temp
          "logFC" : {$gt : 0.0}
        }).sort({"logFC" : -1}).limit(100).toArray((err,documents) => {
          if(err){
            console.log(err);
            client.close();
          } else {
            if(documents === undefined || documents.length === 0 || documents.length === null) {
              console.log("No data found for this gene");
              console.log("Now render to new page and display no data available");
              client.close();
              var message = "Nothing was found";
              res.render('search', {data: message});
            } else {
              console.log("Data available!");
              console.log("Now render to new page and display data");
              let temporary = documents[0].meta_data;
              console.log(temporary);
              db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
                if(err) {
                  console.log("error in fetching meta data");
                } else {
                  console.log("Meta-data fetched successfully");
                  // console.log(doc[0]);
                  client.close();
                  res.render('search', {data: documents , metaData : doc});
                }
              });
              // console.log(documents);
            }
          }
        })
      });
    } else if (logfc === "Down") {
      mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find({
          // "Gene.symbol" : temp
          "logFC" : {$lt : 0.0}
        }).sort({"logFC" : 1}).limit(100).toArray((err,documents) => {
          if(err){
            console.log(err);
            client.close();
          } else {
            if(documents === undefined || documents.length === 0 || documents.length === null) {
              console.log("No data found for this gene");
              console.log("Now render to new page and display no data available");
              client.close();
              var message = "Nothing was found";
              res.render('search', {data: message});
            } else {
              console.log("Data available!");
              console.log("Now render to new page and display data");
              let temporary = documents[0].meta_data;
              console.log(temporary);
              db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
                if(err) {
                  console.log("error in fetching meta data");
                } else {
                  console.log("Meta-data fetched successfully");
                  // console.log(doc[0]);
                  client.close();
                  res.render('search', {data: documents , metaData : doc});
                }
              });
              // console.log(documents);
            }
          }
        })
      });
    } else {
      //nothing end of logfc if and else
    }
  }
  // console.log(temp);

  // if(temp !== null && temp !== '' && temp !== undefined) {
  //   mongo.connect(url, function(err,client) {
  //     assert.equal(null,err);
  //     var db = client.db('igemConcordia2020');
  //     db.collection(collection).find({ "Gene.symbol" : temp }).toArray((err,documents) => {
  //       if(err){
  //         console.log(err);
  //         client.close();
  //       } else {
  //         if(documents === undefined || documents.length === 0 || documents.length === null) {
  //           console.log("No data found for this gene");
  //           console.log("Now render to new page and display no data available");
  //           client.close();
  //           var message = "Nothing was found";
  //           res.render('search', {data: message});
  //         } else {
  //           console.log("Data available!");
  //           console.log("Now render to new page and display data");
  //           let temporary = documents[0].meta_data;
  //           console.log(temporary);
  //           db.collection(secondCollection).find({"Link" : temporary}).toArray((err,doc) => {
  //             if(err) {
  //               console.log("error in fetching meta data");
  //             } else {
  //               console.log("Meta-data fetched successfully");
  //               // console.log(doc[0]);
  //               client.close();
  //               res.render('search', {data: documents , metaData : doc});
  //             }
  //           });
  //           // console.log(documents);
  //         }
  //       }
  //     })
  //   });
  // }

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
                // insertNewColumn_GSE4136();
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
function convertFloatlogFC() {
    // need to make a function to check if data is already changed to float or not, if already float then skip
    mongo.connect(url, function(err, client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find().forEach(function(data) {
            db.collection(collection).update({
                "_id" :data._id,
                "moop" : data.moop
            }, {
                "$set" : {
                    "logFC" : parseFloat(data.logFC),
                    "adj.P.Val" : parseFloat(data.adj.P.Val),
                    "P.Value" : parseFloat(data.P.Value),
                    "t" : parseFloat(data.t),
                    "B" : parseFloat(data.B)
                }
            });
        })
    });
}
function displayDB() {
    mongo.connect(url, function(err,client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).find({}).sort({"adj.P.Val" : 1}).toArray((err,document) => {
            if(err) {
                console.log("error: " + err);
            } else {
                console.log("Data fetched and ready to display");
                for(var i = 0; i <document.length; i++) {
                    console.log(document[i].Gene.symbol);
                    console.log(document[i].adj.P.Val);
                    console.log(document[i].logFC);
                }
            }
        })
    })
}

// let NCBIAPI = 'https://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=mouse[orgn]';





module.exports = router;
