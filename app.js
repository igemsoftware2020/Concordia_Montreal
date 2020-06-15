var createError = require('http-errors');
var express = require('express');
var path = require('path');
var cookieParser = require('cookie-parser');
var logger = require('morgan');
var fs = require('fs');
const CSVToJSON = require('csvtojson');
var indexRouter = require('./routes/index');
const sample = require('./routes/sample');
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



module.exports = app;
