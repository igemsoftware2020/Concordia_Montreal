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
const mongo = require('mongodb');
var assert = require('assert');
var url = 'mongodb://localhost:27017';
// var url = 'mongodb://mongo.ccb.lab:27017';
const collection = "metaData";
const DSE4136DATA = require('../JSON/GSE105058_meta_Data');


function insert_GSE105058_Plant_Data_Meta(){
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

module.exports = insert_GSE105058_Plant_Data_Meta();
