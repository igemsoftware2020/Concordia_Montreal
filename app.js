//===========================================================================
// © Copyright 2020 iGEM Concordia, Maher Hassanain, Benjamin Clark, Hajar Mouddene, Grecia Olano
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
var createError = require('http-errors');
var express = require('express');
var path = require('path');
var cookieParser = require('cookie-parser');
var logger = require('morgan');
var fs = require('fs');
const CSVToJSON = require('csvtojson');
var indexRouter = require('./routes/index');
// Routes
var literatureRouter = require('./routes/literature');
var communityRouter = require('./routes/community');
var researchRouter = require('./routes/research');
// const sample = require('./routes/sample');
// GSE - section
const gse4136_5th = require('./routes/GSE4136_5th');
const gse4136_25th = require('./routes/GSE4136_25th');
const gse40468 = require('./routes/GSE40468');
const gse95388 = require('./routes/GSE95388');
const gse50881 = require('./routes/GSE50881');
const gse105058 = require('./routes/GSE105058');
const gse90166 = require('./routes/GSE90166');

// META DATA section
const gse4136_5th_meta = require('./routes/GSE4136_meta_5th');
const gse4136_25th_meta = require('./routes/GSE4136_meta_25th');
const gse40468_meta = require('./routes/GSE40468_meta');
const gse95388_meta = require('./routes/GSE95388_meta');
const gse50881_meta = require('./routes/GSE50881_meta');
const gse105058_meta = require('./routes/GSE105058_meta');
const gse90166_meta = require('./routes/GSE90166_meta');


// console.log(sample);
var app = express();

// view engine setup
app.set('views', path.join(__dirname, 'views'));
app.set('view engine', 'ejs');

app.use(logger('dev'));
app.use(express.json());
app.use(express.urlencoded({ extended: false }));
app.use(cookieParser());
app.use(express.static(path.join(__dirname, 'public')));

app.use('/', indexRouter);
app.use('/literature', literatureRouter);
app.use('/community', communityRouter);
app.use('/research', researchRouter);


// catch 404 and forward to error handler
app.use(function(req, res, next) {
  next(createError(404));
});

// error handler
app.use(function(err, req, res, next) {
  // set locals, only providing error in development
  res.locals.message = err.message;
  res.locals.error = req.app.get('env') === 'development' ? err : {};

  // render the error page
  res.status(err.status || 500);
  res.render('error');
});

// SAMPLE DATA CONVERT TO JSON

// CSVToJSON().fromFile("./de_expression/GSE4136_YEAST_EXAMPLE_DE.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)"
//
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE4136_YEAST.json', data);
// });


// UP-TO-DATE DATA FROM GSE 4136 5TH AND 25TH GENERATIONS

// CSVToJSON().fromFile("./de_expression/GSE4136_5thGen.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)",
//     "EGEOD" : "GSE4136",
//     "Organism" : "Yeast",
//     "Species" : "S.Cervisiae",
//     "Strain" : "BY4743",
//     "StudyType" : "HARV",
//     "AssayType" : "Microarray",
//     "Gen" : "5",
//     "meta_data" : "0"
//
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE4136_5thGen_Data.json', data);
// });

// CSVToJSON().fromFile("./de_expression/GSE4136_25thGen.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)",
//     "EGEOD" : "GSE4136",
//     "Organism" : "Yeast",
//     "Species" : "S.Cervisiae",
//     "Strain" : "BY4743",
//     "StudyType" : "HARV",
//     "AssayType" : "Microarray",
//     "Gen" : "5",
//     "meta_data" : "0"
//
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE4136_25thGen_Data.json', data);
// });

// CSVToJSON().fromFile("./de_expression/GSE95388.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)",
//     "EGEOD" : "GSE4136",
//     "Organism" : "Yeast",
//     "Species" : "S.Cervisiae",
//     "Strain" : "BY4743",
//     "StudyType" : "HARV",
//     "AssayType" : "Microarray",
//     "Gen" : "5",
//     "meta_data" : "0"
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE95388_Data.json', data);
// });


// CSVToJSON().fromFile("./de_expression/GSE40648_Ecoli.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)",
//     "EGEOD" : "GSE4136",
//     "Organism" : "Yeast",
//     "Species" : "S.Cervisiae",
//     "Strain" : "BY4743",
//     "StudyType" : "HARV",
//     "AssayType" : "Microarray",
//     "Gen" : "5",
//     "meta_data" : "0"
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE40648_Data.json', data);
// });

// CSVToJSON().fromFile("./de_expression/GSE64468.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)",
//     "EGEOD" : "GSE4136",
//     "Organism" : "Yeast",
//     "Species" : "S.Cervisiae",
//     "Strain" : "BY4743",
//     "StudyType" : "HARV",
//     "AssayType" : "Microarray",
//     "Gen" : "5",
//     "meta_data" : "0"
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE64468_Data.json', data);
// });

// CSVToJSON().fromFile("./de_expression/GSE50881.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)",
//     "EGEOD" : "GSE4136",
//     "Organism" : "Yeast",
//     "Species" : "S.Cervisiae",
//     "Strain" : "BY4743",
//     "StudyType" : "HARV",
//     "AssayType" : "Microarray",
//     "Gen" : "5",
//     "meta_data" : "0"
//
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE50881_Data.json', data);
// });


// CSVToJSON().fromFile("./de_expression/GSE4136_YEAST_EXAMPLE_DE.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)"
//
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE4136_YEAST.json', data);
// });


// UP-TO-DATE DATA FROM GSE 4136 5TH AND 25TH GENERATIONS

// CSVToJSON().fromFile("./de_expression/GSE105058.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)",
//     "EGEOD" : "GSE4136",
//     "Organism" : "Yeast",
//     "Species" : "S.Cervisiae",
//     "Strain" : "BY4743",
//     "StudyType" : "HARV",
//     "AssayType" : "Microarray",
//     "Gen" : "5",
//     "meta_data" : "0"
//
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE105058_Data.json', data);
// });



// UP-TO-DATE DATA FROM GSE 4136 5TH AND 25TH GENERATIONS
//
// CSVToJSON().fromFile("./de_expression/GSE90166.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "adj.P.Val" : "0.000137538508078003",
//     "P.Value" : "2.29827749126127e-08",
//     "t" : "-31.3039813536524",
//     "B" : "9.60606491713567",
//     "logFC" : "-3.64977447963713",
//     "Gene_symbol" : "AAH1",
//     "Platform_ORF" : "YNL141W",
//     "Gene_title" : "adenine deaminase",
//     "GO_Function" : "adenine deaminase activity///adenine deaminase activity///adenine deaminase activity///deaminase activity///hydrolase activity///metal ion binding///zinc ion binding",
//     "GO_Process" : "adenine catabolic process///adenine catabolic process///adenine catabolic process///hypoxanthine salvage///hypoxanthine salvage///hypoxanthine salvage///nucleotide metabolic process///purine ribonucleoside monophosphate biosynthetic process///purine-containing compound salvage///purine-containing compound salvage",
//     "GO_Component" : "cytoplasm///cytoplasm///nucleus///nucleus",
//     "Chromosome_annotation" : "Chromosome XIV, NC_001146.8 (359596..360639)",
//     "EGEOD" : "GSE4136",
//     "Organism" : "Yeast",
//     "Species" : "S.Cervisiae",
//     "Strain" : "BY4743",
//     "StudyType" : "HARV",
//     "AssayType" : "Microarray",
//     "Gen" : "5",
//     "meta_data" : "0"
//
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE90166_Data.json', data);
// });

// META-DATA CONVERSION TO JSON

// CSVToJSON().fromFile("./de_expression/GSE4136_meta_25thGen.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "Count": "1",
//     "Accesssions": "GSM94606",
//     "Treatment": "Control sample 5th generation_rep1",
//     "Description": "This sample is a control for comparison against samples grown in low-shear modeled microgravity (LSMMG).",
//     "Link": "1",
//     "Experimenter": "maher,,hassanain",
//     "Contact": "kmcinnerney@montana.edu",
//     "Title": "Yeast Genomic Expression Patterns in Response to Low-Shear Modeled Microgravity",
//     "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4136",
//     "PMIDs": "17201921",
//     "Institute": "MConcordia University",
//     "Design": "Four conditions are compared with three replicates each:  yeast grown in low-shear modeled microgravity (HARV  etc..",
//     "PlatformID": "GPL2529",
//     "Type": "Expression profiling by array",
//     "Summary": "The goal of this study was to assess whether low shear-modeled microgravity (LSMMG) etc.."
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE4136_meta_25thGen_Data.json', data);
// });

// CSVToJSON().fromFile("./de_expression/GSE4136_meta_5thGen.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "Count" : "1",
//     "Accesssions" : "GSM94606",
//     "Treatment" : "Control sample 5th generation_rep1",
//     "Description" : "This sample is a control for comparison against samples grown in low-shear modeled microgravity (LSMMG).",
//     "Link" : "0",
//     "Experimenter": "maher,,hassanain",
//     "Contact": "kmcinnerney@montana.edu",
//     "Title": "Yeast Genomic Expression Patterns in Response to Low-Shear Modeled Microgravity",
//     "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4136",
//     "PMIDs": "17201921",
//     "Institute": "MConcordia University",
//     "Design": "Four conditions are compared with three replicates each:  yeast grown in low-shear modeled microgravity (HARV  etc..",
//     "PlatformID": "GPL2529",
//     "Type": "Expression profiling by array",
//     "Summary": "The goal of this study was to assess whether low shear-modeled microgravity (LSMMG) etc.."
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE4136_meta_5thGen_Data.json', data);
// });


// CSVToJSON().fromFile("./de_expression/GSE95388_meta_.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "Count" : "1",
//     "Accesssions" : "GSM94606",
//     "Treatment" : "Control sample 5th generation_rep1",
//     "Description" : "This sample is a control for comparison against samples grown in low-shear modeled microgravity (LSMMG).",
//     "Link" : "0",
//     "Experimenter": "maher,,hassanain",
//     "Contact": "kmcinnerney@montana.edu",
//     "Title": "Yeast Genomic Expression Patterns in Response to Low-Shear Modeled Microgravity",
//     "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4136",
//     "PMIDs": "17201921",
//     "Institute": "MConcordia University",
//     "Design": "Four conditions are compared with three replicates each:  yeast grown in low-shear modeled microgravity (HARV  etc..",
//     "PlatformID": "GPL2529",
//     "Type": "Expression profiling by array",
//     "Summary": "The goal of this study was to assess whether low shear-modeled microgravity (LSMMG) etc.."
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE95388_meta_Data.json', data);
// });

// CSVToJSON().fromFile("./de_expression/GSE40648_meta_ecoli.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "Count" : "1",
//     "Accesssions" : "GSM94606",
//     "Treatment" : "Control sample 5th generation_rep1",
//     "Description" : "This sample is a control for comparison against samples grown in low-shear modeled microgravity (LSMMG).",
//     "Link" : "0",
//     "Experimenter": "maher,,hassanain",
//     "Contact": "kmcinnerney@montana.edu",
//     "Title": "Yeast Genomic Expression Patterns in Response to Low-Shear Modeled Microgravity",
//     "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4136",
//     "PMIDs": "17201921",
//     "Institute": "MConcordia University",
//     "Design": "Four conditions are compared with three replicates each:  yeast grown in low-shear modeled microgravity (HARV  etc..",
//     "PlatformID": "GPL2529",
//     "Type": "Expression profiling by array",
//     "Summary": "The goal of this study was to assess whether low shear-modeled microgravity (LSMMG) etc.."
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE40648_meta_Data.json', data);
// });

// CSVToJSON().fromFile("./de_expression/GSE50881_meta_y.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "Count" : "1",
//     "Accesssions" : "GSM94606",
//     "Treatment" : "Control sample 5th generation_rep1",
//     "Description" : "This sample is a control for comparison against samples grown in low-shear modeled microgravity (LSMMG).",
//     "Link" : "0",
//     "Experimenter": "maher,,hassanain",
//     "Contact": "kmcinnerney@montana.edu",
//     "Title": "Yeast Genomic Expression Patterns in Response to Low-Shear Modeled Microgravity",
//     "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4136",
//     "PMIDs": "17201921",
//     "Institute": "MConcordia University",
//     "Design": "Four conditions are compared with three replicates each:  yeast grown in low-shear modeled microgravity (HARV  etc..",
//     "PlatformID": "GPL2529",
//     "Type": "Expression profiling by array",
//     "Summary": "The goal of this study was to assess whether low shear-modeled microgravity (LSMMG) etc.."
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE50881_meta_Data.json', data);
// });


// CSVToJSON().fromFile("./de_expression/GSE64468_meta_.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "Count" : "1",
//     "Accesssions" : "GSM94606",
//     "Treatment" : "Control sample 5th generation_rep1",
//     "Description" : "This sample is a control for comparison against samples grown in low-shear modeled microgravity (LSMMG).",
//     "Link" : "0",
//     "Experimenter": "maher,,hassanain",
//     "Contact": "kmcinnerney@montana.edu",
//     "Title": "Yeast Genomic Expression Patterns in Response to Low-Shear Modeled Microgravity",
//     "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4136",
//     "PMIDs": "17201921",
//     "Institute": "MConcordia University",
//     "Design": "Four conditions are compared with three replicates each:  yeast grown in low-shear modeled microgravity (HARV  etc..",
//     "PlatformID": "GPL2529",
//     "Type": "Expression profiling by array",
//     "Summary": "The goal of this study was to assess whether low shear-modeled microgravity (LSMMG) etc.."
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE64468_meta_Data.json', data);
// });

// CSVToJSON().fromFile("./de_expression/GSE105058_meta_Atha.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "Count" : "1",
//     "Accesssions" : "GSM94606",
//     "Treatment" : "Control sample 5th generation_rep1",
//     "Description" : "This sample is a control for comparison against samples grown in low-shear modeled microgravity (LSMMG).",
//     "Link" : "0",
//     "Experimenter": "maher,,hassanain",
//     "Contact": "kmcinnerney@montana.edu",
//     "Title": "Yeast Genomic Expression Patterns in Response to Low-Shear Modeled Microgravity",
//     "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4136",
//     "PMIDs": "17201921",
//     "Institute": "MConcordia University",
//     "Design": "Four conditions are compared with three replicates each:  yeast grown in low-shear modeled microgravity (HARV  etc..",
//     "PlatformID": "GPL2529",
//     "Type": "Expression profiling by array",
//     "Summary": "The goal of this study was to assess whether low shear-modeled microgravity (LSMMG) etc.."
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE105058_meta_Data.json', data);
// });
//
// CSVToJSON().fromFile("./de_expression/GSE90166_meta_Smut.csv").then(source => {
//   // console.log(source);
//   source.push({
//     "Count" : "1",
//     "Accesssions" : "GSM94606",
//     "Treatment" : "Control sample 5th generation_rep1",
//     "Description" : "This sample is a control for comparison against samples grown in low-shear modeled microgravity (LSMMG).",
//     "Link" : "0",
//     "Experimenter": "maher,,hassanain",
//     "Contact": "kmcinnerney@montana.edu",
//     "Title": "Yeast Genomic Expression Patterns in Response to Low-Shear Modeled Microgravity",
//     "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4136",
//     "PMIDs": "17201921",
//     "Institute": "MConcordia University",
//     "Design": "Four conditions are compared with three replicates each:  yeast grown in low-shear modeled microgravity (HARV  etc..",
//     "PlatformID": "GPL2529",
//     "Type": "Expression profiling by array",
//     "Summary": "The goal of this study was to assess whether low shear-modeled microgravity (LSMMG) etc.."
//   });
//   var data = JSON.stringify(source);
//   fs.writeFileSync('./JSON/GSE90166_meta_Data.json', data);
// });

module.exports = app;
