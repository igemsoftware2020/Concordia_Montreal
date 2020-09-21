//===========================================================================
//						GSE40468_meta.js
//				  Developed by Maher Hassanain
//						 July 2020
//===========================================================================
const mongo = require('mongodb');
var assert = require('assert');
var url = 'mongodb://localhost:27017';
const collection = "metaData";
const DSE4136DATA = require('../JSON/GSE40648_meta_Data');


function insert_GSE40648_Bacteria_Data_Meta(){
    mongo.connect(url, function(err, client) {
        assert.equal(null,err);
        var db = client.db('igemConcordia2020');
        db.collection(collection).count((err, count) => {
            if(err) {
                console.log("Error in fetching from DB");
            } else {
                if (count === 0 && !err) {
                    console.log("Collection for GSE4136 meta-data is empty, insert to DB");
                    db.collection(collection).insertMany(DSE4136DATA, function(err, res) {
                        if (err) {
                            console.log("Error in inserting GSE4136 Data to DB");
                            console.log(err);
                        } else {
                            console.log("GSE4136 Data inserted successfully to DB");
                            client.close();
                        }
                    });
                } else {
                    console.log("GSE4136 Data is already filled, keep as it is");
                    client.close();
                }
            }
        });
    });
}

module.exports = insert_GSE40648_Bacteria_Data_Meta();